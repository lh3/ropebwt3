#include <string.h>
#include <assert.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "rle.h"
#include "kthread.h"
#include "kalloc.h"

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
		while (l > 0 && off + l > e->cnt[a+1] && off <= e->cnt[a+1]) {
			int64_t t = e->cnt[a+1] - off;
			if (t > 0) rope_insert_run(r->r[a], off - e->cnt[a], c, t, &cache);
			off += t, l -= t, ++a;
		}
		if (off + l <= e->cnt[a+1]) { // <= is important
			rope_insert_run(r->r[a], off - e->cnt[a], c, l, &cache);
			off += l;
			while (a < RB3_ASIZE && off == e->cnt[a+1]) ++a;
		}
	}
	if (is_free) rld_destroy(e);
	return r;
}

/**********************
 * Parallel plain2fmr *
 **********************/

typedef struct {
	const uint8_t *bwt;
	const int64_t *cnt;
	mrope_t *r;
} p2fmr_aux_t;

static void worker_p2fmr(void *data, long c, int tid)
{
	p2fmr_aux_t *a = (p2fmr_aux_t*)data;
	int64_t i, i0, off;
	rope_t *r = a->r->r[c];
	rpcache_t cache;
	memset(&cache, 0, sizeof(rpcache_t));
	for (i = 0, off = 0; i < c; ++i)
		off += a->cnt[i];
	for (i0 = 0, i = 1; i <= a->cnt[c]; ++i) {
		if (i == a->cnt[c] || a->bwt[off+i0] != a->bwt[off+i]) {
			rope_insert_run(r, i0, a->bwt[off+i0], i - i0, &cache);
			i0 = i;
		}
	}
}

mrope_t *rb3_enc_plain2fmr(int64_t len, const uint8_t *bwt, int max_nodes, int block_len, int32_t n_threads)
{
	p2fmr_aux_t aux;
	int64_t i, cnt[256], c6;
	int32_t c;

	if (max_nodes <= 0) max_nodes = ROPE_DEF_MAX_NODES;
	if (block_len <= 0) block_len = ROPE_DEF_BLOCK_LEN;
	memset(cnt, 0, 256 * sizeof(int64_t));
	for (i = 0; i < len; ++i) ++cnt[bwt[i]];
	for (c = 6, c6 = 0; c < 256; ++c) c6 += cnt[c];
	assert(c6 == 0);

	aux.bwt = bwt, aux.cnt = cnt;
	aux.r = mr_init(max_nodes, block_len, MR_SO_IO);
	if (n_threads > 6) n_threads = 6; // up to 6 threads
	if (n_threads > 1) {
		kt_for(n_threads, worker_p2fmr, &aux, 6);
	} else {
		for (c = 0; c < 6; ++c)
			worker_p2fmr(&aux, c, c);
	}
	return aux.r;
}

/******************************
 * Calculate rank for merging *
 ******************************/

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
} mgrank_aux_t;

static void worker_cal_rank(void *data, long k, int tid)
{
	mgrank_aux_t *a = (mgrank_aux_t*)data;
	rb3_mg_rank1(a->fa, a->fb, a->rb, k);
}

void rb3_mg_rank(const rb3_fmi_t *fa, const rb3_fmi_t *fb, int64_t *rb, int n_threads)
{
	int64_t k, acb[RB3_ASIZE+1];
	rb3_fmi_get_acc(fb, acb);
	if (n_threads > 1) {
		mgrank_aux_t a;
		a.fa = fa, a.fb = fb, a.rb = rb;
		kt_for(n_threads, worker_cal_rank, &a, acb[1]);
	} else {
		for (k = 0; k < acb[1]; ++k)
			rb3_mg_rank1(fa, fb, rb, k);
	}
}

/*********
 * Merge *
 *********/

typedef struct {
	const int64_t *rank;
	const int64_t *aca, *acb;
	mrope_t *r;
} mgins_aux_t;

static void worker_mgins(void *data, long c, int tid)
{
	mgins_aux_t *a = (mgins_aux_t*)data;
	int64_t i;
	rope_t *r = a->r->r[c];
	rpcache_t cache;
	memset(&cache, 0, sizeof(rpcache_t));
	for (i = a->acb[c]; i < a->acb[c+1]; ++i) {
		int64_t x = a->rank[i];
		assert((x&7) == c);
		rope_insert_run(r, (x>>6) - (a->aca[c] + a->acb[c]), x>>3&7, 1, &cache);
	}
}

void rb3_fmi_merge(mrope_t *r, rb3_fmi_t *fb, int n_threads, int free_fb)
{
	rb3_fmi_t fa;
	int64_t *rb, aca[RB3_ASIZE+1], acb[RB3_ASIZE+1];
	int c;
	mgins_aux_t aux;

	rb3_fmi_init(&fa, 0, r);
	rb3_fmi_get_acc(&fa, aca);
	rb3_fmi_get_acc(fb, acb);
	rb = RB3_MALLOC(int64_t, acb[RB3_ASIZE]);
	rb3_mg_rank(&fa, fb, rb, n_threads);
	if (free_fb) rb3_fmi_destroy(fb);
	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] caculated ranks for %ld symbols\n", __func__, rb3_realtime(), rb3_percent_cpu(), (long)acb[RB3_ASIZE]);

	aux.rank = rb, aux.aca = aca, aux.acb = acb, aux.r = r;
	if (n_threads > 1) {
		kt_for(n_threads < 6? n_threads : 6, worker_mgins, &aux, 6);
	} else {
		for (c = 0; c < 6; ++c)
			worker_mgins(&aux, c, c);
	}
	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] inserted %ld symbols\n", __func__, rb3_realtime(), rb3_percent_cpu(), (long)acb[RB3_ASIZE]);
	free(rb);
}

/***************
 * Exact match *
 ***************/

void rb3_fmd_extend(const rb3_fmi_t *f, const rb3_sai_t *ik, rb3_sai_t ok[RB3_ASIZE], int is_back)
{
	int64_t tk[RB3_ASIZE], tl[RB3_ASIZE];
	int c;
	is_back = !!is_back; // 0 or 1
	rb3_fmi_rank2a(f, ik->x[!is_back], ik->x[!is_back] + ik->size, tk, tl);
	for (c = 0; c < RB3_ASIZE; ++c) {
		ok[c].x[!is_back] = f->acc[c] + tk[c];
		ok[c].size = (tl[c] -= tk[c]);
	}
	ok[0].x[is_back] = ik->x[is_back];
	ok[4].x[is_back] = ok[0].x[is_back] + tl[0];
	ok[3].x[is_back] = ok[4].x[is_back] + tl[4];
	ok[2].x[is_back] = ok[3].x[is_back] + tl[3];
	ok[1].x[is_back] = ok[2].x[is_back] + tl[2];
	ok[5].x[is_back] = ok[1].x[is_back] + tl[1];
}

static void rb3_sai_reverse(rb3_sai_t *a, int64_t l)
{
	int64_t i;
	rb3_sai_t t;
	for (i = 0; i < l>>1; ++i)
		t = a[i], a[i] = a[l - 1 - i], a[l - 1 - i] = t;
}

int64_t rb3_fmd_smem1(void *km, const rb3_fmi_t *f, int64_t min_occ, int64_t min_len, int64_t len, const uint8_t *q, int64_t x, rb3_sai_v *mem, rb3_sai_v *curr, rb3_sai_v *prev)
{
	int64_t i, j, ret;
	rb3_sai_t ik, ok[6];
	rb3_sai_v *swap;
	size_t oldn = mem->n;

	assert(len <= INT32_MAX); // this can be relaxed if we define a new struct for mem
	rb3_fmd_set_intv(f, q[x], &ik);
	ik.info = x + 1;
	if (ik.size == 0) return x + 1;
	for (i = x + 1, curr->n = 0; i < len; ++i) { // forward extension
		int c = rb3_comp(q[i]);
		rb3_fmd_extend(f, &ik, ok, 0);
		if (ok[c].size != ik.size) {
			Kgrow(km, rb3_sai_t, curr->a, curr->n, curr->m);
			curr->a[curr->n++] = ik;
			if (ok[c].size < min_occ) break;
		}
		ik = ok[c]; ik.info = i + 1;
	}
	if (i == len) {
		Kgrow(km, rb3_sai_t, curr->a, curr->n, curr->m);
		curr->a[curr->n++] = ik;
	}
	rb3_sai_reverse(curr->a, curr->n);
	ret = curr->a[0].info;
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= -1; --i) { // backward extension
		int c = i < 0? 0 : q[i];
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			rb3_sai_t *p = &prev->a[j];
			rb3_fmd_extend(f, p, ok, 1);
			if (c == 0 || ok[c].size < min_occ) {
				if (curr->n == 0 && (int32_t)p->info - i - 1 >= min_len && (mem->n == oldn || i + 1 < mem->a[mem->n-1].info>>32)) {
					rb3_sai_t *q;
					Kgrow(km, rb3_sai_t, mem->a, mem->n, mem->m);
					q = &mem->a[mem->n++];
					*q = *p; q->info |= (int64_t)(i + 1)<<32;
				}
			} else if (curr->n == 0 || ok[c].size != curr->a[curr->n-1].size) {
				ok[c].info = p->info;
				Kgrow(km, rb3_sai_t, curr->a, curr->n, curr->m);
				curr->a[curr->n++] = ok[c];
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}

	rb3_sai_reverse(&mem->a[oldn], mem->n - oldn);
	return ret;
}

int64_t rb3_fmd_smem(void *km, const rb3_fmi_t *f, int64_t len, const uint8_t *q, rb3_sai_v *mem, int64_t min_occ, int64_t min_len)
{
	int64_t x = 0;
	rb3_sai_v curr = {0,0,0}, prev = {0,0,0};
	mem->n = 0;
	do {
		x = rb3_fmd_smem1(km, f, min_occ, min_len, len, q, x, mem, &curr, &prev);
	} while (x < len);
	kfree(km, curr.a);
	kfree(km, prev.a);
	return mem->n;
}

/*******************
 * Other utilities *
 *******************/

int64_t rb3_fmi_retrieve(const rb3_fmi_t *f, int64_t k, kstring_t *s)
{
	int64_t i, ok[RB3_ASIZE], acc[RB3_ASIZE+1];
	int c;
	s->l = 0;
	rb3_fmi_get_acc(f, acc);
	if (k < 0 || k >= acc[RB3_ASIZE]) return -1;
	while ((c = rb3_fmi_rank1a(f, k, ok)) > 0) {
		RB3_GROW(char, s->s, s->l + 1, s->m);
		s->s[s->l++] = "$ACGTN"[c];
		k = acc[c] + ok[c];
	}
	s->s[s->l] = 0;
	for (i = 0; i < s->l>>1; ++i) // reverse
		c = s->s[i], s->s[i] = s->s[s->l - 1 - i], s->s[s->l - 1 - i] = c;
	return k;
}
