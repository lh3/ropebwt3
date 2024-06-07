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
} rb3_fmi_t;

rld_t *rb3_enc_plain2rld(int64_t len, const uint8_t *bwt, int cbits);
rld_t *rb3_enc_fmr2fmd(mrope_t *r, int cbits, int is_free);
mrope_t *rb3_enc_plain2fmr(int64_t len, const uint8_t *bwt, int max_nodes, int block_len, int32_t n_threads);
mrope_t *rb3_enc_fmd2fmr(rld_t *e, int max_nodes, int block_len, int is_free);

void rb3_mg_rank(const rb3_fmi_t *fa, const rb3_fmi_t *fb, int64_t *rb, int n_threads);
void rb3_fmi_merge(mrope_t *r, rb3_fmi_t *fb, int n_threads, int free_fb);

int64_t rb3_fmi_retrieve(const rb3_fmi_t *f, int64_t k, kstring_t *s);

static inline void rb3_fmi_init(rb3_fmi_t *f, rld_t *e, mrope_t *r)
{
	if (e) f->is_fmd = 1, f->e = e, f->r = 0;
	else f->is_fmd = 0, f->e = 0, f->r = r;
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

static inline int64_t rb3_fmi_get_acc(const rb3_fmi_t *fmi, int64_t acc[RB3_ASIZE+1])
{
	if (fmi->is_fmd) {
		memcpy(acc, fmi->e->cnt, (RB3_ASIZE+1) * sizeof(int64_t));
		return fmi->e->cnt[RB3_ASIZE];
	} else return mr_get_ac(fmi->r, acc);
}

static inline void rb3_fmi_destroy(rb3_fmi_t *fmi)
{
	if (fmi->is_fmd) rld_destroy(fmi->e);
	else mr_destroy(fmi->r);
}

static inline void rb3_fmi_restore(rb3_fmi_t *fmi, const char *fn)
{
	fmi->r = mr_restore_file(fn);
	if (fmi->r == 0) {
		fmi->e = rld_restore(fn);
		fmi->is_fmd = 1;
	} else fmi->is_fmd = 0;
}

#ifdef __cplusplus
}
#endif

#endif
