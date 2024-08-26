#ifndef RB3_DAWG_H
#define RB3_DAWG_H

#include <stdint.h>

typedef struct {
	void *km;
	int32_t seq_len, bwt_size, n_occ;
	int32_t primary;
	int32_t *occ, *sa, acc[5];
	uint32_t *bwt;
} rb3_bwtl_t;

#define rb3_bwtl_B0(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

typedef struct {
	int32_t n_pre, c;
	int32_t lo, hi;
	int32_t *pre;
} rb3_dawg_node_t;

typedef struct {
	int32_t n_node, n_pre;
	rb3_dawg_node_t *node;
	int32_t *pre;
	const rb3_bwtl_t *bwt;
} rb3_dawg_t;

void rb3_bwtl_init(void);
rb3_bwtl_t *rb3_bwtl_gen(void *km, int len, const uint8_t *seq);
void rb3_bwtl_rank1a(const rb3_bwtl_t *bwt, int32_t k, int32_t cnt[4]);
void rb3_bwtl_rank2a(const rb3_bwtl_t *bwt, int32_t k, int32_t l, int32_t cntk[4], int32_t cntl[4]);
void rb3_bwtl_destroy(rb3_bwtl_t *bwt);

rb3_dawg_t *rb3_dawg_gen(void *km, const rb3_bwtl_t *q);
rb3_dawg_t *rb3_dawg_gen_linear(void *km, int32_t len, const uint8_t *seq);
void rb3_dawg_destroy(void *km, rb3_dawg_t *g);

#endif
