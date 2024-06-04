#ifndef RB3_ENCODE_H
#define RB3_ENCODE_H

#include "rld0.h"
#include "mrope.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RB3_ASIZE 6

typedef enum { RB3_PLAIN, RB3_FMD, RB3_FMR } rb3_fmt_t;

rld_t *rb3_enc_plain2rld(int64_t len, const uint8_t *bwt);
mrope_t *rb3_enc_plain2fmr(int64_t len, const uint8_t *bwt, int max_nodes, int block_len);

#ifdef __cplusplus
}
#endif

#endif
