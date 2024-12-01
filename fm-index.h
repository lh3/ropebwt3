#ifndef RB3_FM_INDEX_H
#define RB3_FM_INDEX_H

#include <stdint.h>
#include <string.h>
#include "rb3priv.h"
#include "rld0.h"
#include "mrope.h"
#include "io.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RB3_ASIZE 6

typedef enum { RB3_PLAIN, RB3_FMD, RB3_FMR, RB3_TREE, RB3_BRE } rb3_fmt_t;

#define RB3_LOAD_MMAP  0x1
#define RB3_LOAD_SSA   0x2
#define RB3_LOAD_SID   0x4
#define RB3_LOAD_ALL   (RB3_LOAD_SSA|RB3_LOAD_SID)

typedef struct {
	int64_t x[2]; // 0: start of the interval, backward; 1: forward
	int64_t size;
	int64_t info;
} rb3_sai_t;

typedef struct {
	int32_t ms;
	int32_t ss;
	int64_t m, n_ssa;
	uint64_t *r2i; // rank -> index
	uint64_t *ssa; // sampled suffix array
} rb3_ssa_t;

typedef struct {
	int64_t sid, pos;
} rb3_pos_t;

typedef struct {
	int32_t is_fmd;
	rld_t *e;
	mrope_t *r;
	rb3_ssa_t *ssa;
	rb3_sid_t *sid;
	int64_t acc[RB3_ASIZE+1];
} rb3_fmi_t;

typedef struct { size_t n, m; rb3_sai_t *a; } rb3_sai_v;

rld_t *rb3_enc_plain2rld(int64_t len, const uint8_t *bwt, int cbits);
rld_t *rb3_enc_fmr2fmd(mrope_t *r, int cbits, int is_free);
mrope_t *rb3_enc_plain2fmr(int64_t len, const uint8_t *bwt, int max_nodes, int block_len, int32_t n_threads);
mrope_t *rb3_enc_fmd2fmr(rld_t *e, int max_nodes, int block_len, int is_free);

void *rb3_r2cache_init(void *km, int32_t max);
void rb3_r2cache_destroy(void *rc_);
void rb3_fmi_rank2a_cached(const rb3_fmi_t *fmi, void *rc_, int64_t k, int64_t l, int64_t ok[6], int64_t ol[6]);
void rb3_fmd_extend_cached(const rb3_fmi_t *f, void *rc, const rb3_sai_t *ik, rb3_sai_t ok[RB3_ASIZE], int is_back);

void rb3_mg_rank(const rb3_fmi_t *fa, const rb3_fmi_t *fb, int64_t *rb, int n_threads);
void rb3_mg_rank_plain(const rb3_fmi_t *fa, int64_t len, const uint8_t *seq, int64_t *rb, int64_t acc[RB3_ASIZE+1], int n_threads);
void rb3_fmi_merge(mrope_t *r, rb3_fmi_t *fb, int n_threads, int free_fb);
void rb3_fmi_merge_plain(mrope_t *r, int64_t len, const uint8_t *seq, int n_threads);

int64_t rb3_fmi_get_r(const rb3_fmi_t *f);
int64_t rb3_fmi_get_acc(const rb3_fmi_t *fmi, int64_t acc[RB3_ASIZE+1]);
int64_t rb3_fmi_retrieve(const rb3_fmi_t *f, int64_t k, kstring_t *s);
void rb3_fmd_extend(const rb3_fmi_t *f, const rb3_sai_t *ik, rb3_sai_t ok[RB3_ASIZE], int is_back);
int64_t rb3_fmd_smem(void *km, const rb3_fmi_t *f, int64_t len, const uint8_t *q, rb3_sai_v *mem, int64_t min_occ, int64_t min_len);
int64_t rb3_fmd_smem_TG(void *km, const rb3_fmi_t *f, int64_t len, const uint8_t *q, rb3_sai_v *mem, int64_t min_occ, int64_t min_len);
int32_t rb3_fmd_smem_present(const rb3_fmi_t *f, int64_t len, const uint8_t *q, int64_t min_len);

int64_t rb3_ssa(const rb3_fmi_t *f, const rb3_ssa_t *sa, int64_t k, int64_t *si);
int64_t rb3_ssa_multi(void *km, const rb3_fmi_t *f, const rb3_ssa_t *ssa, int64_t lo, int64_t hi, int64_t max_sa, rb3_pos_t *sa);
void rb3_ssa_destroy(rb3_ssa_t *sa);
int rb3_ssa_dump(const rb3_ssa_t *sa, const char *fn);
rb3_ssa_t *rb3_ssa_restore(const char *fn);
rb3_ssa_t *rb3_ssa_gen(const rb3_fmi_t *f, int ssa_shift, int n_threads);

int rb3_fmi_load_all(rb3_fmi_t *f, const char *fn, int32_t load_flag);

static inline int rb3_comp(int c)
{
	return c >= 1 && c <= 4? 5 - c : c;
}

static inline void rb3_fmd_set_intv(const rb3_fmi_t *f, int c, rb3_sai_t *ik)
{
	ik->x[0] = f->acc[c], ik->size = f->acc[c+1] - f->acc[c], ik->x[1] = f->acc[rb3_comp(c)], ik->info = 0;
}

static inline void rb3_fmi_init(rb3_fmi_t *f, rld_t *e, mrope_t *r)
{
	if (e) f->is_fmd = 1, f->e = e, f->r = 0;
	else f->is_fmd = 0, f->e = 0, f->r = r;
	f->ssa = 0;
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

static inline void rb3_fmi_free(rb3_fmi_t *fmi)
{
	if (fmi->is_fmd) rld_destroy(fmi->e);
	else mr_destroy(fmi->r);
	if (fmi->ssa) rb3_ssa_destroy(fmi->ssa);
	if (fmi->sid) rb3_sid_destroy(fmi->sid);
	fmi->e = 0, fmi->r = 0, fmi->ssa = 0;
}

static inline void rb3_fmi_restore(rb3_fmi_t *fmi, const char *fn, int use_mmap)
{
	fmi->r = 0, fmi->e = 0, fmi->ssa = 0, fmi->sid = 0;
	fmi->e = use_mmap? rld_restore_mmap(fn) : rld_restore(fn);
	if (fmi->e == 0) {
		fmi->r = mr_restore_file(fn);
		fmi->is_fmd = 0;
	} else fmi->is_fmd = 1;
	if (fmi->e == 0 && fmi->r == 0) return;
	rb3_fmi_get_acc(fmi, fmi->acc);
}

static inline int rb3_fmi_is_symmetric(const rb3_fmi_t *f)
{
	return ((f->acc[1]&1) == 0 && f->acc[2] - f->acc[1] == f->acc[5] - f->acc[4] && f->acc[3] - f->acc[2] == f->acc[4] - f->acc[3]);
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
