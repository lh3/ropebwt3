#ifndef RB3_ALIGN_H
#define RB3_ALIGN_H

#include "fm-index.h"

typedef struct {
	int32_t n_best;
	int32_t min_sc, end_len;
	int32_t match, mis;
	int32_t gap_open, gap_ext;
	int32_t r2cache_size;
} rb3_swopt_t;

typedef struct {
	int32_t score;
	int32_t qlen, rlen;
	int32_t n_cigar;
	int32_t n_qoff, blen, mlen;
	int64_t lo, hi;
	uint8_t *rseq;
	uint32_t *cigar;
	int32_t *qoff;
} rb3_swrst_t;

void rb3_swopt_init(rb3_swopt_t *opt);
void rb3_sw(void *km, const rb3_swopt_t *opt, const rb3_fmi_t *f, int len, const uint8_t *seq, rb3_swrst_t *rst);
void rb3_swrst_free(rb3_swrst_t *rst);

#endif
