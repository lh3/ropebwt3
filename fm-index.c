#include <string.h>
#include <assert.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "rle.h"
#include "kthread.h"

/***********************
 * Encoding conversion *
 ***********************/

rld_t *rb3_enc_plain2rld(int64_t len, const uint8_t *bwt, int cbits)
{
	int64_t i, i0;
	rld_t *e;
	rlditr_t ei;
	if (cbits <= 0) cbits = 3;
	e = rld_init(RB3_ASIZE, cbits); // TODO: check alphabet
	rld_itr_init(e, &ei, 0);
	for (i0 = 0, i = 1; i <= len; ++i) {
		if (i == len || bwt[i0] != bwt[i]) {
			rld_enc(e, &ei, i - i0, bwt[i0]);
			i0 = i;
		}
	}
	rld_enc_finish(e, &ei);
	return e;
}

rld_t *rb3_enc_fmr2fmd(mrope_t *r, int cbits, int is_free)
{
	rld_t *e;
	rlditr_t ei;
	mritr_t ri;
	const uint8_t *block;

	if (cbits <= 0) cbits = 3;
	e = rld_init(RB3_ASIZE, cbits);
	mr_itr_first(r, &ri, 1);
	rld_itr_init(e, &ei, 0);
	while ((block = mr_itr_next_block(&ri)) != 0) {
		const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
		while (q < end) {
			int c = 0;
			int64_t l;
			rle_dec1(q, c, l);
			rld_enc(e, &ei, l, c);
		}
	}
	if (is_free) mr_destroy(r);
	rld_enc_finish(e, &ei);
	return e;
}

mrope_t *rb3_enc_plain2fmr(int64_t len, const uint8_t *bwt, int max_nodes, int block_len)
{
	int64_t i, i0, off, cnt[256], c6;
	int32_t c;
	mrope_t *r;
	rpcache_t cache;

	if (max_nodes <= 0) max_nodes = ROPE_DEF_MAX_NODES;
	if (block_len <= 0) block_len = ROPE_DEF_BLOCK_LEN;
	memset(cnt, 0, 256 * sizeof(int64_t));
	for (i = 0; i < len; ++i) ++cnt[bwt[i]];
	for (c = 6, c6 = 0; c < 256; ++c) c6 += cnt[c];
	assert(c6 == 0);

	memset(&cache, 0, sizeof(rpcache_t));
	r = mr_init(max_nodes, block_len, MR_SO_IO);
	for (c = 0, off = 0; c < 6; ++c) {
		for (i0 = 0, i = 1; i <= cnt[c]; ++i) {
			if (i == cnt[c] || bwt[off+i0] != bwt[off+i]) {
				rope_insert_run(r->r[c], i0, bwt[off+i0], i - i0, &cache);
				i0 = i;
			}
		}
		off += cnt[c];
	}
	return r;
}

mrope_t *rb3_enc_fmd2fmr(rld_t *e, int max_nodes, int block_len, int is_free)
{
	mrope_t *r;
	rlditr_t itr;
	int64_t l, off;
	int c, a;
	rpcache_t cache;

	if (max_nodes <= 0) max_nodes = ROPE_DEF_MAX_NODES;
	if (block_len <= 0) block_len = ROPE_DEF_BLOCK_LEN;
	r = mr_init(max_nodes, block_len, MR_SO_IO);
	memset(&cache, 0, sizeof(rpcache_t));

	rld_itr_init(e, &itr, 0);
	a = 0, off = 0;
	while ((l = rld_dec(e, &itr, &c, is_free)) > 0) {
		if (off + l <= e->cnt[a+1]) { // <= is important
			rope_insert_run(r->r[a], off - e->cnt[a], c, l, &cache);
			off += l;
		} else {
			while (off < e->cnt[a+1]) {
				int64_t t = e->cnt[a+1] - off;
				rope_insert_run(r->r[a], off - e->cnt[a], c, t, &cache);
				off += t, l -= off, ++a;
			}
		}
	}
	if (is_free) rld_destroy(e);
	return r;
}

/*********
 * Merge *
 *********/

void rb3_mg_rank1(const rb3_fmi_t *fa, const rb3_fmi_t *fb, int64_t *rb, int64_t p)
{
	int64_t aca[RB3_ASIZE+1], acb[RB3_ASIZE+1], ka, kb;
	int c, last_c = 0;
	rb3_fmi_get_acc(fa, aca);
	rb3_fmi_get_acc(fb, acb);
	ka = aca[1], kb = p;
	while (1) {
		int64_t oa[RB3_ASIZE], ob[RB3_ASIZE];
		c = rb3_fmi_rank1a(fb, kb, ob);
		rb[kb] = (ka + kb) << 6 | c << 3 | last_c;
		last_c = c;
		if (c == 0) break;
		kb = acb[c] + ob[c];
		rb3_fmi_rank1a(fa, ka, oa);
		ka = aca[c] + oa[c];
	}
}

typedef struct {
	const rb3_fmi_t *fa, *fb;
	int64_t *rb;
} mgaux_t;

static void worker(void *data, long k, int tid)
{
	mgaux_t *a = (mgaux_t*)data;
	rb3_mg_rank1(a->fa, a->fb, a->rb, k);
}

void rb3_mg_rank(const rb3_fmi_t *fa, const rb3_fmi_t *fb, int64_t *rb, int n_threads)
{
	int64_t k, acb[RB3_ASIZE+1];
	rb3_fmi_get_acc(fb, acb);
	if (n_threads > 1) {
		mgaux_t a;
		a.fa = fa, a.fb = fb, a.rb = rb;
		kt_for(n_threads, worker, &a, acb[1]);
	} else {
		for (k = 0; k < acb[1]; ++k)
			rb3_mg_rank1(fa, fb, rb, k);
	}
}

void rb3_fmi_merge(mrope_t *r, rb3_fmi_t *fb, int n_threads, int free_fb)
{
	rb3_fmi_t fa;
	rpcache_t cache;
	int64_t *rb, i, aca[RB3_ASIZE+1], acb[RB3_ASIZE+1];
	int c;

	rb3_fmi_init(&fa, 0, r);
	rb3_fmi_get_acc(&fa, aca);
	rb3_fmi_get_acc(fb, acb);
	rb = RB3_MALLOC(int64_t, acb[RB3_ASIZE]);
	rb3_mg_rank(&fa, fb, rb, n_threads);
	if (free_fb) rb3_fmi_destroy(fb);

	memset(&cache, 0, sizeof(rpcache_t));
	for (i = 0, c = 0; i < acb[RB3_ASIZE]; ++i) {
		int64_t k = rb[i]>>6;
		int b = rb[i]&7, c = rb[i]>>3&7;
		rope_insert_run(r->r[b], k - (aca[b] + acb[b]), c, 1, &cache);
	}
	free(rb);
}
