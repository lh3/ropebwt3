#ifndef RB3_FM_INDEX_H
#define RB3_FM_INDEX_H

#include <stdint.h>
#include <string.h>
#include "rld0.h"
#include "mrope.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RB3_ASIZE 6

typedef enum { RB3_PLAIN, RB3_FMD, RB3_FMR } rb3_fmt_t;

typedef struct {
	int32_t is_fmd;
	mrope_t *r;
	rld_t *e;
} rb3_fmi_t;

rld_t *rb3_enc_plain2rld(int64_t len, const uint8_t *bwt);
mrope_t *rb3_enc_plain2fmr(int64_t len, const uint8_t *bwt, int max_nodes, int block_len);

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

#ifdef __cplusplus
}
#endif

#endif
