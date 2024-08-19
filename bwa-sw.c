#include <string.h>
#include <stdio.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "align.h"
#include "kalloc.h"
#include "dawg.h"
#define kh_packed
#include "khashl-km.h"
#include "ksort.h" // for binary heap

#define reverse_lt(a, b) ((a) > (b)) // sorting will be in the descending order; this is intentional
KSORT_INIT(rb3_64, uint64_t, reverse_lt)
KSORT_INIT(rb3_32, int32_t, reverse_lt)

void rb3_swopt_init(rb3_swopt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->n_best = 25;
	opt->min_sc = 30;
	opt->match = 1, opt->mis = 3;
	opt->gap_open = 3, opt->gap_ext = 1;
	opt->end_len = 11;
	opt->min_mem_len = 0;
	opt->r2cache_size = 0x10000;
}

// NB: don't change these values!
#define SW_FROM_H    0
#define SW_FROM_E    1
#define SW_FROM_F    2
#define SW_FROM_OPEN 0
#define SW_FROM_EXT  1

#define SW_F_UNSET (0xfffffff) // 28 bits

typedef struct { // 48 bytes
	int32_t H, E, F;
	uint32_t H_from:2, E_from:1, F_from:1, F_from_off:28;
	uint32_t H_from_pos, E_from_pos;
	int32_t rlen, qlen;
	int64_t lo, hi;
} sw_cell_t;

typedef struct {
	int32_t n;
	sw_cell_t *a;
} sw_row_t;

#define sw_cell_hash(x) (kh_hash_uint64((x).lo) + kh_hash_uint64((x).hi))
#define sw_cell_eq(x, y) ((x).lo == (y).lo && (x).hi == (y).hi)
KHASHL_SET_INIT(KH_LOCAL, sw_candset_t, sw_candset, sw_cell_t, sw_cell_hash, sw_cell_eq)

static void sw_push_state(int32_t last_op, int32_t op, int c, rb3_swrst_t *rst, int32_t len_only)
{
	if (!len_only) {
		rst->rseq[rst->rlen] = c;
		if (last_op == op)
			rst->cigar[rst->n_cigar - 1] += 1U<<4;
		else
			rst->cigar[rst->n_cigar++] = 1U<<4 | op;
	} else {
		rst->n_cigar += last_op == op? 0 : 1;
	}
	if (op == 7 || op == 8) rst->qlen++, rst->rlen++;
	else if (op == 1) rst->qlen++;
	else if (op == 2) rst->rlen++;
}

static void sw_backtrack_core(const rb3_swopt_t *opt, const rb3_fmi_t *f, const rb3_dawg_t *g, const sw_row_t *row, uint32_t pos, rb3_swrst_t *rst, int32_t len_only)
{ // this is adapted from ns_backtrack() in miniprot
	int32_t n_col = opt->n_best, last = 0, last_op = -1;
	rst->n_cigar = rst->rlen = rst->qlen = 0;
	while (pos > 0) {
		int32_t r = pos / n_col, c;
		const sw_cell_t *p = &row[r].a[pos%n_col];
		int32_t x = (int32_t)p->H_from | (int32_t)p->E_from<<2 | (int32_t)p->F_from<<3;
		int32_t state = last == 0? x&0x3 : last;
		int32_t ext = state == 1 || state == 2? x>>(state+1)&1 : 0; // gap extension or not
		int32_t op = state; // cigar operator
		if (rb3_dbg_flag & RB3_DBG_BT)
			fprintf(stderr, "BT\t%d\t%d\t%d\n", r, pos%n_col, p->H);
		for (c = 1; c < 7; ++c)
			if (f->acc[c] > p->lo) break;
		--c; // this is the reference base
		if (state == SW_FROM_H) {
			op = c == g->node[r].c? 7 : 8; // 7 for "=" and 8 for "X"
			pos = p->H_from_pos;
		} else if (state == SW_FROM_E) {
			assert(p->E > 0 && p->E_from_pos != UINT32_MAX);
			pos = p->E_from_pos;
		} else if (state == SW_FROM_F) {
			assert(p->F > 0 && p->F_from_off != SW_F_UNSET);
			pos = r * n_col + p->F_from_off;
		}
		sw_push_state(last_op, op, c, rst, len_only);
		last_op = op;
		last = (state == 1 || state == 2) && ext? state : 0;
	}
}

static void sw_backtrack(const rb3_swopt_t *opt, const rb3_fmi_t *f, const rb3_dawg_t *g, const sw_row_t *row, uint32_t pos, rb3_swrst_t *rst)
{
	int32_t k;
	const rb3_dawg_node_t *p;
	const sw_cell_t *q;

	// get CIGAR
	if (rst->score == 0) return;
	sw_backtrack_core(opt, f, g, row, pos, rst, 1);
	rst->rseq = RB3_CALLOC(uint8_t, rst->rlen);
	rst->cigar = RB3_CALLOC(uint32_t, rst->n_cigar);
	sw_backtrack_core(opt, f, g, row, pos, rst, 0);

	// calculate block length and matching length
	rst->mlen = rst->blen = 0;
	for (k = 0; k < rst->n_cigar; ++k) {
		int32_t op = rst->cigar[k]&0xf, len = rst->cigar[k]>>4;
		rst->blen += len;
		if (op == 7) rst->mlen += len;
	}

	// get positions
	p = &g->node[pos / opt->n_best];
	q = &row[pos / opt->n_best].a[pos % opt->n_best];
	rst->lo = q->lo, rst->hi = q->hi;
	rst->n_qoff = p->hi - p->lo;
	rst->qoff = RB3_CALLOC(int32_t, p->hi - p->lo);
	for (k = p->lo; k < p->hi; ++k)
		rst->qoff[k - p->lo] = g->bwt->sa[k];
}

static sw_cell_t *sw_update_candset(sw_candset_t *h, const sw_cell_t *p)
{
	khint_t k;
	int absent;
	k = sw_candset_put(h, *p, &absent);
	if (!absent) {
		sw_cell_t *q = &kh_key(h, k);
		if (q->E < p->E) q->E = p->E, q->E_from = p->E_from, q->E_from_pos = p->E_from_pos;
		if (q->F < p->F) q->F = p->F, q->F_from = p->F_from; // NB: F_from_off is populated differently
		if (q->H < p->H) {
			q->H = p->H, q->H_from = p->H_from;
			q->rlen = p->rlen, q->qlen = p->qlen;
			if (p->H_from == SW_FROM_H)
				q->H_from_pos = p->H_from_pos; // TODO: is this correct?
		}
	}
	return &kh_key(h, k);
}

static inline int32_t sw_heap_insert1(uint64_t *heap, int32_t max, int32_t *sz, uint32_t score, uint32_t id)
{
	uint64_t x = (uint64_t)score<<32 | id;
	if (*sz < max) {
		heap[(*sz)++] = x;
		ks_heapup_rb3_64(*sz, heap);
		return 1;
	} else if (x > heap[0]) {
		heap[0] = x;
		ks_heapdown_rb3_64(0, *sz, heap);
		return 1;
	}
	return 0;
}

static void sw_track_F(void *km, const rb3_fmi_t *f, void *rc, sw_candset_t *h, sw_row_t *row)
{ // compute F_from_off at row
	int32_t j, n_F = 0;
	for (j = 0; j < row->n; ++j)
		if (row->a[j].F > 0)
			++n_F;
	if (n_F == 0) return; // no F is calculated
	sw_candset_clear(h);
	for (j = 0; j < row->n; ++j) { // collect cells where F_from_off needs to be calculated
		int absent;
		if (row->a[j].F == 0) continue;
		sw_cell_t r = row->a[j];
		r.H = j; // reuse "H" for index
		sw_candset_put(h, r, &absent);
	}
	for (j = 0; j < row->n - 1; ++j) {
		sw_cell_t r = row->a[j];
		int64_t clo[6], chi[6];
		int c;
		khint_t k;
		rb3_fmi_rank2a_cached(f, rc, r.lo, r.hi, clo, chi);
		for (c = 1; c < 6; ++c) {
			r.lo = f->acc[c] + clo[c];
			r.hi = f->acc[c] + chi[c];
			if (r.lo == r.hi) continue;
			k = sw_candset_get(h, r);
			if (k != kh_end(h)) {
				int32_t i = kh_key(h, k).H;
				if (row->a[i].F_from_off == SW_F_UNSET)
					row->a[i].F_from_off = j;
			}
		}
	}
	for (j = 0; j < row->n - 1; ++j) {
		if (row->a[j].F_from_off == SW_F_UNSET) { // this may happen if the parent cell is not in row; this happens!
			assert(row->a[j].H_from != SW_FROM_F); // but this shouldn't happen
			row->a[j].F = 0; // prevent backtrack to F
		}
	}
}

static void sw_core(void *km, const rb3_swopt_t *opt, const rb3_fmi_t *f, const rb3_dawg_t *g, const rb3_bwtl_t *bwt, rb3_swrst_t *rst)
{
	uint32_t best_pos = 0;
	int32_t i, c, n_col = opt->n_best, m_fstack = opt->n_best * 3;
	int32_t *ks_a, ks_m;
	sw_cell_t *cell, *fstack, *p;
	sw_row_t *row;
	int64_t clo[6], chi[6];
	uint64_t *heap;
	sw_candset_t *h;
	void *rc = 0;

	rc = rb3_r2cache_init(km, opt->r2cache_size);
	cell = Kcalloc(km, sw_cell_t, g->n_node * n_col); // this is the backtracking matrix
	row = Kcalloc(km, sw_row_t, g->n_node);
	for (i = 0; i < g->n_node; ++i)
		row[i].a = &cell[i * n_col];
	p = &row[0].a[row[0].n++]; // point to the first cell
	p->lo = 0, p->hi = f->acc[6];
	p->H_from = SW_FROM_H;
	rst->score = 0;

	fstack = Kcalloc(km, sw_cell_t, m_fstack);
	heap = Kcalloc(km, uint64_t, opt->n_best);
	h = sw_candset_init2(km);
	sw_candset_resize(h, opt->n_best * 4);
	ks_m = opt->n_best * 3;
	ks_a = Kmalloc(km, int32_t, ks_m);
	for (i = 1; i < g->n_node; ++i) { // traverse all nodes in the DAWG in the topological order
		const rb3_dawg_node_t *t = &g->node[i];
		sw_row_t *ri = &row[i];
		int32_t j, k, heap_sz, max_min_sc = 0;
		khint_t itr;
		sw_candset_clear(h);

		// calculate max_min_sc; ignore a cell if its score can't reach max_min_sc
		if (t->n_pre > 1) { // only relevant if there are multiple predecessors
			int32_t n_cell = 0;
			for (j = 0; j < t->n_pre; ++j)
				n_cell += row[t->pre[j]].n;
			if (n_cell > opt->n_best) { // only relevant if there are enough cells
				int32_t l = 0;
				Kgrow(km, int32_t, ks_a, n_cell, ks_m);
				for (j = 0, max_min_sc = 0; j < t->n_pre; ++j) {
					int32_t pid = t->pre[j];
					for (k = 0; k < row[pid].n; ++k)
						ks_a[l++] = row[pid].a[k].H;
				}
				max_min_sc = ks_ksmall_rb3_32(n_cell, ks_a, opt->n_best);
			}
			max_min_sc -= opt->gap_open + opt->gap_ext > opt->mis? opt->gap_open + opt->gap_ext : opt->mis;
			if (max_min_sc < 0) max_min_sc = 0;
		}

		// compute E and H
		for (j = 0; j < t->n_pre; ++j) { // traverse all the predecessors
			int32_t pid = t->pre[j]; // parent/predecessor ID
			if (row[pid].n == 0) continue;
			for (k = 0; k < row[pid].n; ++k) {
				sw_cell_t r;
				p = &row[pid].a[k];
				if (p->H + opt->match < max_min_sc) continue; // this node can't reach opt->n_best
				memset(&r, 0, sizeof(sw_cell_t));
				r.F_from_off = SW_F_UNSET;
				// calculate E
				if (p->H - opt->gap_open > p->E)
					r.E_from = SW_FROM_OPEN, r.E = p->H - opt->gap_open;
				else
					r.E_from = SW_FROM_EXT,  r.E = p->E;
				r.E -= opt->gap_ext;
				if (r.E > 0 && r.E >= max_min_sc && p->qlen >= opt->end_len) { // add to row
					r.lo = p->lo, r.hi = p->hi;
					r.H = r.E;
					r.H_from = SW_FROM_E;
					r.E_from_pos = pid * n_col + k, r.H_from_pos = UINT32_MAX;
					r.rlen = p->rlen, r.qlen = p->qlen + 1;
					sw_update_candset(h, &r);
				}
				// calculate H
				rb3_fmi_rank2a_cached(f, rc, p->lo, p->hi, clo, chi);
				r.H_from = SW_FROM_H, r.H_from_pos = pid * n_col + k;
				r.E = 0, r.E_from_pos = UINT32_MAX;
				for (c = 1; c < 6; ++c) {
					int32_t sc = c == t->c? opt->match : -opt->mis;
					if (p->H + sc <= 0 || p->H + sc < max_min_sc) continue;
					if (c != t->c && p->qlen < opt->end_len) continue;
					r.lo = f->acc[c] + clo[c];
					r.hi = f->acc[c] + chi[c];
					if (r.lo == r.hi) continue;
					r.H = p->H + sc;
					r.rlen = p->rlen + 1, r.qlen = p->qlen + 1;
					sw_update_candset(h, &r);
				}
			}
		}
		ri->n = 0;
		if (kh_size(h) == 0) continue;

		// find top-n hits
		heap_sz = 0;
		kh_foreach(h, itr) {
			sw_heap_insert1(heap, opt->n_best, &heap_sz, kh_key(h, itr).H, itr);
		}
		ks_heapsort_rb3_64(heap_sz, heap);
		ri->n = heap_sz;
		for (j = 0; j < ri->n; ++j)
			ri->a[j] = kh_key(h, (uint32_t)heap[j]);
		for (j = 0; j < heap_sz>>1; ++j) { // reverse heap[] such that it remains a heap
			uint64_t tmp = heap[j];
			heap[j] = heap[heap_sz - j - 1];
			heap[heap_sz - j - 1] = tmp;
		}

		if (p->qlen >= opt->end_len) { // update F
			int32_t n_fstack = 0;
			for (j = ri->n - 1; j >= 0; --j)
				if (ri->a[j].H > opt->gap_open + opt->gap_ext)
					fstack[n_fstack++] = ri->a[j];
			while (n_fstack > 0) {
				sw_cell_t r, z = fstack[--n_fstack];
				int32_t min = heap_sz < opt->n_best? 0 : heap[0]>>32;
				memset(&r, 0, sizeof(sw_cell_t));
				r.H_from_pos = r.E_from_pos = UINT32_MAX, r.F_from_off = SW_F_UNSET;
				if (z.H - opt->gap_open > z.F)
					r.F_from = SW_FROM_OPEN, r.F = z.H - opt->gap_open;
				else
					r.F_from = SW_FROM_EXT,  r.F = z.F;
				r.F -= opt->gap_ext;
				r.H = r.F, r.H_from = SW_FROM_F;
				r.rlen = p->rlen + 1, r.qlen = p->qlen;
				if (r.H <= min) continue;
				rb3_fmi_rank2a_cached(f, rc, z.lo, z.hi, clo, chi);
				for (c = 1; c < 6; ++c) {
					sw_cell_t *q;
					r.lo = f->acc[c] + clo[c];
					r.hi = f->acc[c] + chi[c];
					if (r.lo == r.hi) continue;
					q = sw_update_candset(h, &r);
					sw_heap_insert1(heap, opt->n_best, &heap_sz, r.H, UINT32_MAX);
					if (r.H - opt->gap_ext > min) {
						Kgrow(km, sw_cell_t, fstack, n_fstack, m_fstack);
						fstack[n_fstack++] = *q;
					}
				}
			}
		}

		heap_sz = 0;
		kh_foreach(h, itr) {
			sw_heap_insert1(heap, opt->n_best, &heap_sz, kh_key(h, itr).H, itr);
		}
		ks_heapsort_rb3_64(heap_sz, heap);
		assert(heap_sz > 0);
		ri->n = heap_sz;
		for (j = 0; j < ri->n; ++j)
			ri->a[j] = kh_key(h, (uint32_t)heap[j]);
		sw_track_F(km, f, rc, h, &row[i]);
		if (ri->a[0].H > rst->score)
			rst->score = ri->a->H, best_pos = i * n_col;

		// for debugging
		if (rb3_dbg_flag & RB3_DBG_SW) { // NB: single-threaded only
			fprintf(stderr, "SW\t%d\t[%d,%d)\t%d\t", i, t->lo, t->hi, ri->n);
			for (j = 0; j < t->n_pre; ++j) {
				if (j) fputc(',', stderr);
				fprintf(stderr, "%d", t->pre[j]);
			}
			fputc('\t', stderr);
			for (j = 0; j < ri->n; ++j) {
				if (j) fputc(',', stderr);
				fprintf(stderr, "%d(%d)", ri->a[j].H, ri->a[j].qlen - ri->a[j].rlen);
			}
			fputc('\n', stderr);
		}
	}
	if (rst->score >= opt->min_sc)
		sw_backtrack(opt, f, g, row, best_pos, rst);
	else rst->cigar = 0, rst->rseq = 0, rst->qoff = 0;

	kfree(km, ks_a);
	kfree(km, fstack);
	sw_candset_destroy(h);
	kfree(km, heap);
	kfree(km, row);
	kfree(km, cell);
	rb3_r2cache_destroy(rc);
}

void rb3_sw(void *km, const rb3_swopt_t *opt, const rb3_fmi_t *f, int len, const uint8_t *seq, rb3_swrst_t *rst)
{
	rb3_bwtl_t *q;
	rb3_dawg_t *g;
	rst->score = rst->n_qoff = rst->n_cigar = rst->rlen = rst->qlen = rst->blen = rst->mlen = 0;
	rst->lo_pos = rst->lo_sid = -1;
	if (opt->min_mem_len > 0 && opt->min_mem_len > opt->end_len) {
		if (!rb3_fmd_smem_present(f, len, seq, opt->min_mem_len))
			return;
	}
	q = rb3_bwtl_gen(km, len, seq);
	g = rb3_dawg_gen(km, q);
	sw_core(km, opt, f, g, q, rst);
	rb3_dawg_destroy(km, g); // this doesn't deallocate q
	rb3_bwtl_destroy(q);
	if (f->ssa)
		rst->lo_pos = rb3_ssa(f, f->ssa, rst->lo, &rst->lo_sid);
}

void rb3_swrst_free(rb3_swrst_t *rst)
{
	free(rst->rseq); free(rst->cigar); free(rst->qoff);
}
