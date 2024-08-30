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
	opt->gap_open = 5, opt->gap_ext = 2;
	opt->end_len = 11;
	opt->e2e_drop = -1; // disabled by default
	opt->min_mem_len = 0;
	opt->r2cache_size = 0x10000;
}

#define SW_FROM_H    0
#define SW_FROM_E    1
#define SW_FROM_F    2
#define SW_FROM_OPEN 0
#define SW_FROM_EXT  1

#define SW_F_UNSET (0x3ffffff) // 26 bits

typedef struct {
	int32_t H, E, F;
	uint32_t flt:1, H_from:2, E_from:1, F_from:1, F_from_off:26, F_off_set:1;
	uint32_t H_from_pos, E_from_pos;
	int32_t rlen, qlen;
	int64_t lo, hi, lo_rc; // (lo, lo_rc, hi-lo) is the SA bi-interval
} sw_cell_t;

typedef struct {
	int32_t n;
	sw_cell_t *a;
} sw_row_t;

#define sw_cell_hash(x) (kh_hash_uint64((x).lo) + kh_hash_uint64((x).hi))
#define sw_cell_eq(x, y) ((x).lo == (y).lo && (x).hi == (y).hi)
KHASHL_SET_INIT(KH_LOCAL, sw_candset_t, sw_candset, sw_cell_t, sw_cell_hash, sw_cell_eq)

/*************
 * Backtrack *
 *************/

static void sw_push_state(int32_t last_op, int32_t op, int c, rb3_swhit_t *hit, int32_t len_only)
{
	if (!len_only) { // generating the CIGAR
		hit->rseq[hit->rlen] = c;
		if (last_op == op)
			hit->cigar[hit->n_cigar - 1] += 1U<<4;
		else
			hit->cigar[hit->n_cigar++] = 1U<<4 | op;
	} else { // only update the CIGAR length
		hit->n_cigar += last_op == op? 0 : 1;
	}
	if (op == 7 || op == 8) hit->qlen++, hit->rlen++;
	else if (op == 1) hit->qlen++;
	else if (op == 2) hit->rlen++;
}

static int32_t sw_backtrack1_core(const rb3_swopt_t *opt, const rb3_fmi_t *f, const rb3_dawg_t *g, const sw_row_t *row, uint32_t pos, rb3_swhit_t *hit, int32_t len_only)
{ // this is adapted from ns_backtrack() in miniprot
	int32_t n_col = opt->n_best, last = 0, last_op = -1, ed = 0;
	hit->score = row[pos / n_col].a[pos % n_col].H;
	hit->n_cigar = hit->rlen = hit->qlen = 0;
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
			pos = p->H_from_pos, ed += (op == 8);
		} else if (state == SW_FROM_E) {
			assert(p->E > 0 && p->E_from_pos != UINT32_MAX);
			pos = p->E_from_pos, ++ed;
		} else if (state == SW_FROM_F) {
			assert(p->F > 0 && p->F_off_set);
			pos = r * n_col + p->F_from_off, ++ed;
		}
		sw_push_state(last_op, op, c, hit, len_only);
		last_op = op;
		last = (state == 1 || state == 2) && ext? state : 0;
	}
	return ed;
}

static void sw_backtrack1(const rb3_swopt_t *opt, const rb3_fmi_t *f, const rb3_dawg_t *g, const sw_row_t *row, uint32_t pos, rb3_swhit_t *hit)
{
	int32_t k;
	const rb3_dawg_node_t *p;
	const sw_cell_t *q;

	// get CIGAR
	sw_backtrack1_core(opt, f, g, row, pos, hit, 1); // compute length without allocation
	hit->rseq = RB3_CALLOC(uint8_t, hit->rlen);
	hit->cigar = RB3_CALLOC(uint32_t, hit->n_cigar);
	sw_backtrack1_core(opt, f, g, row, pos, hit, 0);

	// calculate block length and matching length
	hit->mlen = hit->blen = 0;
	for (k = 0; k < hit->n_cigar; ++k) {
		int32_t op = hit->cigar[k]&0xf, len = hit->cigar[k]>>4;
		hit->blen += len;
		if (op == 7) hit->mlen += len;
	}

	// get query positions
	p = &g->node[pos / opt->n_best];
	q = &row[pos / opt->n_best].a[pos % opt->n_best];
	hit->lo = q->lo, hit->hi = q->hi;
	if (p->hi >= 0) { // [p->lo, p->hi) is a SA interval on the query
		hit->n_qoff = p->hi - p->lo;
		hit->qoff = RB3_CALLOC(int32_t, p->hi - p->lo);
		for (k = p->lo; k < p->hi; ++k)
			hit->qoff[k - p->lo] = g->bwt->sa[k];
	} else { // p->lo is the actual position on the query
		hit->n_qoff = 1;
		hit->qoff = RB3_CALLOC(int32_t, 1);
		hit->qoff[0] = p->lo;
	}

	// get reference position for the first hit in the SA interval
	hit->lo_pos = hit->lo_sid = -1;
	if (f->ssa) hit->lo_pos = rb3_ssa(f, f->ssa, hit->lo, &hit->lo_sid);
}

static void sw_cell_dedup(void *km, sw_row_t *row)
{ // mark a cell to be filtered if [lo_rc,lo_rc+(hi-lo)) is contained in a cell at a higher score
	int32_t i, j, k, *a;
	if (row->n <= 1) return; // no need
	a = Kmalloc(km, int32_t, row->n);
	a[0] = 0;
	for (i = k = 1; i < row->n; ++i) {
		sw_cell_t *p = &row->a[i];
		for (j = 0; j < k; ++j) {
			sw_cell_t *q = &row->a[a[j]];
			if (q->lo_rc <= p->lo_rc && q->lo_rc + (q->hi - q->lo) >= p->lo_rc + (p->hi - p->lo))
				break;
		}
		if (j == k) a[k++] = i;
		else p->flt = 1;
	}
	kfree(km, a);
}

static void sw_backtrack(const rb3_swopt_t *opt, const rb3_fmi_t *f, const rb3_dawg_t *g, const sw_row_t *row, int32_t qlen, uint32_t best_pos, rb3_swrst_t *r, rb3_hapdiv_t *a)
{
	int32_t i, n_col = opt->n_best;
	if (opt->flag & (RB3_SWF_E2E|RB3_SWF_HAPDIV)) { // end-to-end mode
		const sw_row_t *p = &row[g->n_node - 1]; // last row, i.e. the end of the query sequence
		int32_t n = 0, H0;
		rb3_swhit_t tmp;
		if (p->n == 0) return; // do nothing if the alignment doesn't reach the end
		memset(&tmp, 0, sizeof(rb3_swhit_t));
		H0 = p->a[0].H;
		for (i = 0, n = 0; i < p->n; ++i) { // count hits
			const sw_cell_t *q = &p->a[i];
			if (!q->flt && q->H_from == SW_FROM_H && q->H >= opt->min_sc && (opt->e2e_drop < 0 || H0 - q->H <= opt->e2e_drop))
				++n;
		}
		if (n == 0) return;
		if (r) r->n = n, r->a = RB3_CALLOC(rb3_swhit_t, n);
		if (a) {
			memset(a, 0, sizeof(*a));
			a->n_al = n;
		}
		for (i = 0, n = 0; i < p->n; ++i) { // backtrack
			const sw_cell_t *q = &p->a[i];
			uint32_t pos = (g->n_node - 1) * n_col + i;
			if (!q->flt && q->H_from == SW_FROM_H && q->H >= opt->min_sc && (opt->e2e_drop < 0 || H0 - q->H <= opt->e2e_drop)) {
				if (r) { // get full alignment
					sw_backtrack1(opt, f, g, row, pos, &r->a[n++]);
				} else if (a) { // get summary information
					int32_t ed;
					ed = sw_backtrack1_core(opt, f, g, row, pos, &tmp, 1);
					a->max_ed = a->max_ed > ed? a->max_ed : ed;
					ed = ed < RB2_SW_MAX_ED? ed : RB2_SW_MAX_ED;
					a->n_hap[ed] += q->hi - q->lo;
				}
			}
		}
	} else { // local mode; TODO: support split alignment
		r->n = 1;
		r->a = RB3_CALLOC(rb3_swhit_t, r->n);
		sw_backtrack1(opt, f, g, row, best_pos, &r->a[0]);
	}
}

/************************
 * Filling the "matrix" *
 ************************/

static sw_cell_t *sw_update_candset(sw_candset_t *h, const sw_cell_t *p)
{
	khint_t k;
	int absent;
	k = sw_candset_put(h, *p, &absent);
	if (!absent) {
		sw_cell_t *q = &kh_key(h, k);
		q->rlen = q->rlen > p->rlen? q->rlen : p->rlen; // furthest extension
		q->qlen = q->qlen > p->qlen? q->qlen : p->qlen;
		if (q->E < p->E) q->E = p->E, q->E_from = p->E_from, q->E_from_pos = p->E_from_pos;
		if (q->F < p->F) q->F = p->F, q->F_from = p->F_from; // NB: F_from_off is populated differently
		if (q->H < p->H) {
			q->H = p->H, q->H_from = p->H_from;
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

static void sw_track_F(void *km, const rb3_fmi_t *f, void *rc, sw_candset_t *h, rb3_u128_t *fpar, sw_row_t *row)
{ // compute F_from_off at row
	int32_t j;
	sw_candset_clear(h);
	for (j = 0; j < row->n; ++j) {
		int absent;
		sw_cell_t r = row->a[j];
		r.H = j; // reuse "H" for index
		sw_candset_put(h, r, &absent);
	}
	for (j = 0; j < row->n; ++j) {
		sw_cell_t r, *p = &row->a[j];
		khint_t k;
		if (p->F == 0 || p->F_from_off == SW_F_UNSET) continue;
		r.lo = fpar[p->F_from_off].x, r.hi = fpar[p->F_from_off].y;
		k = sw_candset_get(h, r);
		if (k != kh_end(h)) {
			p->F_from_off = kh_key(h, k).H, p->F_off_set = 1;
		} else { // this may happen if the parent cell is not in row; this happens!
			assert(p->H_from != SW_FROM_F); // but this shouldn't happen
			p->F_from_off = SW_F_UNSET; // prevent backtrack to F
		}
	}
}

#define sw_cell2sai(cell, sai) ((sai)->x[0] = (cell)->lo, (sai)->x[1] = (cell)->lo_rc, (sai)->size = (cell)->hi - (cell)->lo)
#define sw_sai2cell(sai, cell) ((cell)->lo = (sai)->x[0], (cell)->hi = (sai)->x[0] + (sai)->size, (cell)->lo_rc = (sai)->x[1])

static void sw_core(void *km, const rb3_swopt_t *opt, const rb3_fmi_t *f, const rb3_dawg_t *g, int32_t qlen, rb3_swrst_t *rst, rb3_hapdiv_t *anno)
{
	uint32_t best_pos = 0;
	int32_t i, c, n_col = opt->n_best, m_fstack, m_fpar, best_score;
	int32_t *ks_a, ks_m;
	sw_cell_t *cell, *fstack, *p;
	sw_row_t *row;
	uint64_t *heap;
	sw_candset_t *h;
	rb3_u128_t *fpar;
	void *rc = 0;

	if (rst) rst->n = 0, rst->a = 0;
	rc = rb3_r2cache_init(km, opt->r2cache_size);
	cell = Kcalloc(km, sw_cell_t, g->n_node * n_col); // this is the backtracking matrix
	row = Kcalloc(km, sw_row_t, g->n_node);
	for (i = 0; i < g->n_node; ++i)
		row[i].a = &cell[i * n_col];
	p = &row[0].a[row[0].n++]; // point to the first cell
	p->lo = 0, p->hi = f->acc[6], p->lo_rc = 0; // the SA bi-interval of an empty string, the root
	p->H_from = SW_FROM_H;
	best_score = 0;

	m_fstack = m_fpar = opt->n_best * 3; // fstack and fpar are temporary arrays for computing the keeping track of the F state
	fstack = Kcalloc(km, sw_cell_t, m_fstack);
	fpar = Kcalloc(km, rb3_u128_t, m_fpar);
	heap = Kcalloc(km, uint64_t, opt->n_best);
	h = sw_candset_init2(km);
	sw_candset_resize(h, opt->n_best * 4);
	ks_m = opt->n_best * 3; // ks_a is used for computing k-small
	ks_a = Kmalloc(km, int32_t, ks_m);
	for (i = 1; i < g->n_node; ++i) { // traverse all nodes in the DAWG in the topological order
		const rb3_dawg_node_t *t = &g->node[i];
		sw_row_t *ri = &row[i];
		int32_t j, k, heap_sz, max_min_sc = 0, n_fpar = 0;
		rb3_sai_t ik, ok[RB3_ASIZE];
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

		// compute H and E
		for (j = 0; j < t->n_pre; ++j) { // traverse all the predecessors
			int32_t pid = t->pre[j]; // parent/predecessor ID
			if (row[pid].n == 0) continue;
			for (k = 0; k < row[pid].n; ++k) {
				sw_cell_t r;
				p = &row[pid].a[k];
				if (p->H + opt->match < max_min_sc) continue; // this cell can't reach opt->n_best
				memset(&r, 0, sizeof(sw_cell_t));
				r.F_from_off = SW_F_UNSET;
				// calculate H
				r.H_from = SW_FROM_H, r.H_from_pos = pid * n_col + k, r.E_from_pos = UINT32_MAX;
				sw_cell2sai(p, &ik);
				rb3_fmd_extend_cached(f, rc, &ik, ok, 1);
				for (c = 1; c < 6; ++c) {
					int32_t sc = c == t->c? opt->match : -opt->mis;
					if (ok[c].size == 0) continue;
					if (p->H + sc <= 0 || p->H + sc < max_min_sc) continue;
					if (c != t->c && p->qlen < opt->end_len) continue;
					sw_sai2cell(&ok[c], &r);
					r.H = p->H + sc;
					r.rlen = p->rlen + 1, r.qlen = p->qlen + 1;
					sw_update_candset(h, &r);
				}
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

		if (p->qlen >= opt->end_len) { // update F; TODO: this algorithm is not good and even is not really correct
			int32_t n_fstack = 0;
			n_fpar = 0;
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
				r.rlen = z.rlen + 1, r.qlen = z.qlen;
				if (r.H <= min) continue;
				sw_cell2sai(&z, &ik);
				rb3_fmd_extend_cached(f, rc, &ik, ok, 1);
				for (c = 1; c < 6; ++c) {
					sw_cell_t *q;
					if (ok[c].size == 0) continue;
					sw_sai2cell(&ok[c], &r);
					q = sw_update_candset(h, &r);
					if (q->F == r.F) { // if q->F > r.F, this deletion is not really added
						sw_heap_insert1(heap, opt->n_best, &heap_sz, r.H, UINT32_MAX);
						Kgrow(km, rb3_u128_t, fpar, n_fpar, m_fpar);
						fpar[n_fpar].x = z.lo, fpar[n_fpar].y = z.hi;
						q->F_from_off = n_fpar++;
						if (r.H - opt->gap_ext > min) {
							Kgrow(km, sw_cell_t, fstack, n_fstack, m_fstack);
							fstack[n_fstack++] = *q;
						}
					}
				}
			}
		}

		heap_sz = 0;
		kh_foreach(h, itr) { // rebuild the heap
			sw_heap_insert1(heap, opt->n_best, &heap_sz, kh_key(h, itr).H, itr);
		}
		ks_heapsort_rb3_64(heap_sz, heap);
		assert(heap_sz > 0);
		ri->n = heap_sz;
		for (j = 0; j < ri->n; ++j)
			ri->a[j] = kh_key(h, (uint32_t)heap[j]);
		if (n_fpar > 0) sw_track_F(km, f, rc, h, fpar, ri); // compute F_from_off for backtrack
		if (ri->a[0].H > best_score)
			best_score = ri->a->H, best_pos = i * n_col;
		if (i == g->n_node - 1) sw_cell_dedup(km, ri); // dedup the last cell

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
	kfree(km, ks_a);
	kfree(km, fpar);
	kfree(km, fstack);
	sw_candset_destroy(h);
	kfree(km, heap);
	rb3_r2cache_destroy(rc);

	if (best_score >= opt->min_sc)
		sw_backtrack(opt, f, g, row, qlen, best_pos, rst, anno);

	kfree(km, row);
	kfree(km, cell);
}

/*****************
 * External APIs *
 *****************/

void rb3_sw(void *km, const rb3_swopt_t *opt, const rb3_fmi_t *f, int len, const uint8_t *seq, rb3_swrst_t *rst)
{
	rb3_bwtl_t *q = 0;
	rb3_dawg_t *g;
	if (opt->min_mem_len > 0 && opt->min_mem_len > opt->end_len) {
		if (!rb3_fmd_smem_present(f, len, seq, opt->min_mem_len))
			return;
	}
	if (opt->flag & RB3_SWF_E2E) {
		g = rb3_dawg_gen_linear(km, len, seq);
	} else {
		q = rb3_bwtl_gen(km, len, seq);
		g = rb3_dawg_gen(km, q);
	}
	sw_core(km, opt, f, g, len, rst, 0);
	rb3_dawg_destroy(km, g); // this doesn't deallocate q
	if (q) rb3_bwtl_destroy(q);
}

void rb3_hapdiv(void *km, const rb3_swopt_t *opt, const rb3_fmi_t *f, int len, const uint8_t *seq, rb3_hapdiv_t *hd)
{
	rb3_dawg_t *g;
	g = rb3_dawg_gen_linear(km, len, seq);
	sw_core(km, opt, f, g, len, 0, hd);
	rb3_dawg_destroy(km, g);
}

void rb3_swrst_free(rb3_swrst_t *rst)
{
	int32_t i;
	for (i = 0; i < rst->n; ++i) {
		free(rst->a[i].rseq); free(rst->a[i].cigar); free(rst->a[i].qoff);
	}
	free(rst->a);
}
