#ifndef RB3_FM_INDEX_H
#define RB3_FM_INDEX_H

#include <stdint.h>
#include <string.h>
#include "rb3priv.h"
#include "rld0.h"
#include "mrope.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RB3_ASIZE 6

typedef enum { RB3_PLAIN, RB3_FMD, RB3_FMR, RB3_TREE } rb3_fmt_t;

typedef struct {
	int32_t is_fmd;
	rld_t *e;
	mrope_t *r;
	int64_t acc[RB3_ASIZE+1];
} rb3_fmi_t;

typedef struct {
	int64_t x[2]; // 0: start of the interval, backward; 1: forward
	int64_t size;
	int64_t info;
} rb3_sai_t;

typedef struct {
	int32_t n_best;
	int32_t min_sc;
	int32_t match, mis, ambi;
	int32_t gap_open, gap_ext;
} rb3_swopt_t;

typedef struct { size_t n, m; rb3_sai_t *a; } rb3_sai_v;

rld_t *rb3_enc_plain2rld(int64_t len, const uint8_t *bwt, int cbits);
rld_t *rb3_enc_fmr2fmd(mrope_t *r, int cbits, int is_free);
mrope_t *rb3_enc_plain2fmr(int64_t len, const uint8_t *bwt, int max_nodes, int block_len, int32_t n_threads);
mrope_t *rb3_enc_fmd2fmr(rld_t *e, int max_nodes, int block_len, int is_free);

void *rb3_r2cache_init(void *km, int32_t max);
void rb3_r2cache_destroy(void *rc_);
void rb3_fmi_rank2a_cached(const rb3_fmi_t *fmi, void *rc_, int64_t k, int64_t l, int64_t ok[6], int64_t ol[6]);

void rb3_mg_rank(const rb3_fmi_t *fa, const rb3_fmi_t *fb, int64_t *rb, int n_threads);
void rb3_fmi_merge(mrope_t *r, rb3_fmi_t *fb, int n_threads, int free_fb);

int64_t rb3_fmi_retrieve(const rb3_fmi_t *f, int64_t k, kstring_t *s);
void rb3_fmd_extend(const rb3_fmi_t *f, const rb3_sai_t *ik, rb3_sai_t ok[RB3_ASIZE], int is_back);
int64_t rb3_fmd_smem(void *km, const rb3_fmi_t *f, int64_t len, const uint8_t *q, rb3_sai_v *mem, int64_t min_occ, int64_t min_len);
int64_t rb3_fmd_gmem(void *km, const rb3_fmi_t *f, int64_t len, const uint8_t *q, rb3_sai_v *mem, int64_t min_occ, int64_t min_len);

void rb3_swopt_init(rb3_swopt_t *opt);
void rb3_sw(void *km, const rb3_swopt_t *opt, const rb3_fmi_t *f, int len, const uint8_t *seq);

static inline int rb3_comp(int c)
{
	return c >= 1 && c <= 4? 5 - c : c;
}

static inline void rb3_fmd_set_intv(const rb3_fmi_t *f, int c, rb3_sai_t *ik)
{
	ik->x[0] = f->acc[c], ik->size = f->acc[c+1] - f->acc[c], ik->x[1] = f->acc[rb3_comp(c)], ik->info = 0;
}

static inline int64_t rb3_fmi_get_acc(const rb3_fmi_t *fmi, int64_t acc[RB3_ASIZE+1])
{
	if (fmi->is_fmd) {
		memcpy(acc, fmi->e->cnt, (RB3_ASIZE+1) * sizeof(int64_t));
		return fmi->e->cnt[RB3_ASIZE];
	} else return mr_get_ac(fmi->r, acc);
}

static inline void rb3_fmi_init(rb3_fmi_t *f, rld_t *e, mrope_t *r)
{
	if (e) f->is_fmd = 1, f->e = e, f->r = 0;
	else f->is_fmd = 0, f->e = 0, f->r = r;
	rb3_fmi_get_acc(f, f->acc);
}

static inline void rb3_fmi_rank2a(const rb3_fmi_t *fmi, int64_t k, int64_t l, int64_t *ok, int64_t *ol)
{
	if (fmi->is_fmd) rld_rank2a(fmi->e, k, l, (uint64_t*)ok, (uint64_t*)ol);
	else mr_rank2a(fmi->r, k, l, ok, ol);
}

static inline int rb3_fmi_rank1a(const rb3_fmi_t *fmi, int64_t k, int64_t *ok)
{
	return fmi->is_fmd? rld_rank1a(fmi->e, k, (uint64_t*)ok) : mr_rank1a(fmi->r, k, ok);
}

static inline void rb3_fmi_destroy(rb3_fmi_t *fmi)
{
	if (fmi->is_fmd) rld_destroy(fmi->e);
	else mr_destroy(fmi->r);
}

static inline void rb3_fmi_restore(rb3_fmi_t *fmi, const char *fn, int use_mmap)
{
	fmi->r = 0, fmi->e = 0;
	fmi->e = use_mmap? rld_restore_mmap(fn) : rld_restore(fn);
	if (fmi->e == 0) {
		fmi->r = mr_restore_file(fn);
		fmi->is_fmd = 0;
	} else fmi->is_fmd = 1;
	rb3_fmi_get_acc(fmi, fmi->acc);
}

static inline int64_t rb3_fmi_extend1(const rb3_fmi_t *f, int64_t *k, int64_t *l, int c)
{
	int64_t tk[RB3_ASIZE], tl[RB3_ASIZE];
	rb3_fmi_rank2a(f, *k, *l, tk, tl);
	*k = f->acc[c] + tk[c];
	*l = f->acc[c] + tl[c];
	return *l - *k;
}

#ifdef __cplusplus
}
#endif

#endif
