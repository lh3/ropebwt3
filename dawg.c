#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "dawg.h"
#include "kalloc.h"
#include "libsais.h" // for libsais()
#include "io.h" // for rb3_nt6_table[]
#define kh_packed
#include "khashl-km.h"

/*******************
 * Lightweight BWT *
 *******************/

static uint32_t bwtl_cnt_table[256];

void rb3_bwtl_init(void)
{
	int i, c;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (c = 0; c != 4; ++c)
			x |= (((i&3) == c) + ((i>>2&3) == c) + ((i>>4&3) == c) + (i>>6 == c)) << (c<<3);
		bwtl_cnt_table[i] = x;
	}
}

rb3_bwtl_t *rb3_bwtl_gen(void *km, int len, const uint8_t *seq)
{
	rb3_bwtl_t *b;
	int32_t i;
	b = Kcalloc(km, rb3_bwtl_t, 1); // allocate long-term memory first to reduce memory fragmentation
	b->km = km;
	b->seq_len = len;
	b->bwt_size = (len + 15) / 16;
	b->bwt = Kcalloc(km, uint32_t, b->bwt_size);
	b->n_occ = (len + 16) / 16 * 4;
	b->occ = Kcalloc(km, int32_t, b->n_occ);

	{ // calculate b->bwt
		uint8_t *s, *s8;
		s8 = Kcalloc(km, uint8_t, len);
		for (i = 0; i < len; ++i) {
			s8[i] = rb3_nt6_table[seq[i]];
			if (s8[i] == 5) s8[i] = 1; // NB: convert ambiguous bases to A
		}
		b->sa = Kcalloc(km, int32_t, len + 1);
		libsais(s8, &b->sa[1], len, 0, 0);
		b->sa[0] = len;
		s = Kcalloc(km, uint8_t, len + 1);
		for (i = 0; i <= len; ++i) {
			if (b->sa[i] == 0) b->primary = i;
			else s[i] = s8[b->sa[i] - 1] - 1;
		}
		kfree(km, s8);
		for (i = b->primary; i < len; ++i) s[i] = s[i + 1];
		for (i = 0; i < len; ++i)
			b->bwt[i>>4] |= s[i] << ((15 - (i&15)) << 1);
		kfree(km, s);
	}
	{ // calculate b->occ
		int32_t c[4];
		memset(c, 0, 16);
		for (i = 0; i < len; ++i) {
			if (i % 16 == 0)
				memcpy(b->occ + (i/16) * 4, c, 16);
			++c[rb3_bwtl_B0(b, i)];
		}
		if (i % 16 == 0)
			memcpy(b->occ + (i/16) * 4, c, 16);
		memcpy(&b->acc[1], c, 16);
		b->acc[0] = 1;
		for (i = 1; i < 5; ++i) b->acc[i] += b->acc[i-1];
	}
	return b;
}

void rb3_bwtl_rank1a(const rb3_bwtl_t *bwt, int32_t k, int32_t cnt[4])
{
	uint32_t x, b;
	if (k > bwt->primary) --k; // because $ is not in bwt
	memcpy(cnt, bwt->occ + (k>>4<<2), 16);
	if (k % 16 == 0) return;
	--k;
	b = bwt->bwt[k>>4] & ~((1U<<((~k&15)<<1)) - 1);
	x = bwtl_cnt_table[b&0xff] + bwtl_cnt_table[b>>8&0xff] + bwtl_cnt_table[b>>16&0xff] + bwtl_cnt_table[b>>24];
	x -= 15 - (k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

void rb3_bwtl_rank2a(const rb3_bwtl_t *bwt, int32_t k, int32_t l, int32_t cntk[4], int32_t cntl[4])
{
	rb3_bwtl_rank1a(bwt, k, cntk);
	rb3_bwtl_rank1a(bwt, l, cntl);
}

void rb3_bwtl_destroy(rb3_bwtl_t *bwt)
{
	kfree(bwt->km, bwt->occ);
	kfree(bwt->km, bwt->bwt);
	kfree(bwt->km, bwt->sa);
	kfree(bwt->km, bwt);
}

/*********************
 * DAWG Construction *
 *********************/

typedef struct {
	int32_t deg, cnt, id;
} deg_cell_t;

KHASHL_MAP_INIT(KH_LOCAL, sw_deg_t, sw_deg, uint64_t, deg_cell_t, kh_hash_uint64, kh_eq_generic)

static sw_deg_t *sw_cal_deg(void *km, const rb3_bwtl_t *bwt) // calculate the in-degree of each node in DAWG
{
	sw_deg_t *h;
	int32_t n = 0, m = 16, c, absent;
	uint64_t *a;
	khint_t itr;

	h = sw_deg_init2(km);
	sw_deg_resize(h, bwt->seq_len + 1); // preallocate for efficiency
	itr = sw_deg_put(h, bwt->seq_len + 1, &absent); // put the root
	kh_val(h, itr).deg = 0;

	a = Kmalloc(km, uint64_t, m);
	a[n++] = bwt->seq_len + 1; // the root interval is [0,bwt->seq_len + 1)
	while (n > 0) { // count the in-degree of each node in DAWG
		uint64_t x = a[--n]; // pop
		int32_t rlo[4], rhi[4];
		rb3_bwtl_rank2a(bwt, x>>32, (int32_t)x, rlo, rhi);
		for (c = 3; c >= 0; --c) { // traverse children
			uint64_t key;
			int32_t lo = bwt->acc[c] + rlo[c];
			int32_t hi = bwt->acc[c] + rhi[c];
			if (lo == hi) continue;
			key = (uint64_t)lo << 32 | hi;
			itr = sw_deg_put(h, key, &absent);
			if (absent) {
				kh_val(h, itr).deg = 0;
				Kgrow(km, uint64_t, a, n, m);
				a[n++] = key;
			}
			kh_val(h, itr).deg++;
		}
	}
	kfree(km, a);
	return h;
}

rb3_dawg_t *rb3_dawg_gen(void *km, const rb3_bwtl_t *q) // generate DAWG
{
	khint_t itr;
	sw_deg_t *h;
	rb3_dawg_t *g;
	int32_t i, off_pre = 0, n_a = 0, id = 0;
	uint64_t *a;
	rb3_dawg_node_t *p;

	h = sw_cal_deg(km, q);
	g = Kcalloc(km, rb3_dawg_t, 1); // allocate the DAWG upfront
	g->bwt = q;
	g->n_node = kh_size(h);
	g->node = Kcalloc(km, rb3_dawg_node_t, g->n_node);
	kh_foreach(h, itr) {
		g->n_pre += kh_val(h, itr).deg; // n_pre is sum of in-degrees across all nodes
		kh_val(h, itr).cnt = kh_val(h, itr).id = 0;
	}
	g->pre = Kcalloc(km, int32_t, g->n_pre);

	a = Kmalloc(km, uint64_t, g->n_node); // this is way over-allocated, but that is fine
	p = &g->node[id++];
	p->lo = 0, p->hi = q->seq_len + 1, p->n_pre = 0, p->pre = g->pre;
	a[n_a++] = q->seq_len + 1; // the root interval
	while (n_a > 0) { // topological sorting; this while loop is different from the one in sw_cal_deg()
		uint64_t x = a[--n_a]; // pop
		int32_t rlo[4], rhi[4], c;
		rb3_bwtl_rank2a(q, x>>32, (int32_t)x, rlo, rhi);
		for (c = 3; c >= 0; --c) { // traverse children
			uint64_t key;
			int32_t lo = q->acc[c] + rlo[c];
			int32_t hi = q->acc[c] + rhi[c];
			if (lo == hi) continue;
			key = (uint64_t)lo << 32 | hi;
			itr = sw_deg_get(h, key);
			assert(itr != kh_end(h));
			kh_val(h, itr).cnt++;
			if (kh_val(h, itr).cnt == kh_val(h, itr).deg) { // the predecessors being visited
				kh_val(h, itr).id = id;
				p = &g->node[id++];
				p->lo = lo, p->hi = hi, p->c = c + 1, p->n_pre = 0, p->pre = &g->pre[off_pre]; // c+1 for the nt6 encoding
				off_pre += kh_val(h, itr).deg;
				a[n_a++] = key;
			}
		}
	}
	assert(id == g->n_node && off_pre == g->n_pre);
	kfree(km, a);

	for (i = 0; i < g->n_node; ++i) { // populate predecessors
		int32_t rlo[4], rhi[4], c;
		rb3_bwtl_rank2a(q, g->node[i].lo, g->node[i].hi, rlo, rhi);
		for (c = 0; c < 4; ++c) { // traverse i's children
			int32_t lo = q->acc[c] + rlo[c];
			int32_t hi = q->acc[c] + rhi[c];
			if (lo == hi) continue;
			itr = sw_deg_get(h, (uint64_t)lo << 32 | hi);
			p = &g->node[kh_val(h, itr).id];
			p->pre[p->n_pre++] = i;
		}
	}
	sw_deg_destroy(h);

	if (rb3_dbg_flag & RB3_DBG_DAWG) { // for debugging
		for (i = 0; i < g->n_node; ++i) {
			rb3_dawg_node_t *p = &g->node[i];
			int j;
			fprintf(stderr, "DG\t%d\t[%d,%d)\t", i, p->lo, p->hi);
			for (j = 0; j < p->n_pre; ++j) {
				if (j) fputc(',', stderr);
				fprintf(stderr, "%d", p->pre[j]);
			}
			fputc('\n', stderr);
		}
	}
	return g;
}

rb3_dawg_t *rb3_dawg_gen_linear(void *km, int32_t len, const uint8_t *seq)
{
	int32_t i;
	rb3_dawg_t *g;
	g = Kcalloc(km, rb3_dawg_t, 1);
	g->n_node = len + 1;
	g->node = Kcalloc(km, rb3_dawg_node_t, g->n_node);
	g->n_pre = len;
	g->pre = Kcalloc(km, int32_t, g->n_pre);
	g->node[0].n_pre = 0, g->node[0].c = -1, g->node[0].pre = g->pre;
	g->node[0].lo = len, g->node[0].hi = -1;
	for (i = 0; i < len; ++i) {
		rb3_dawg_node_t *p = &g->node[i + 1];
		p->lo = len - 1 - i, p->hi = -1;
		p->c = rb3_nt6_table[seq[p->lo]];
		p->n_pre = 1;
		p->pre = &g->pre[i];
		p->pre[0] = i;
	}
	return g;
}

void rb3_dawg_destroy(void *km, rb3_dawg_t *g)
{
	kfree(km, g->pre); kfree(km, g->node); kfree(km, g);
}
