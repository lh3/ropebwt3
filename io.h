#ifndef RB3_IO_H
#define RB3_IO_H

#include <stdint.h>
#include "rb3priv.h" // for kstring_t

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int64_t tot_len;
	int64_t n_seq;
	char **name;
	int32_t *len;
} rb3_seqlist_t;

struct rb3_seqio_s;
typedef struct rb3_seqio_s rb3_seqio_t;

extern const uint8_t rb3_nt6_table[128];

rb3_seqio_t *rb3_seq_open(const char *fn, int is_line);
void rb3_seq_close(rb3_seqio_t *fp);
int64_t rb3_seq_read(rb3_seqio_t *fp, kstring_t *seq, int64_t max_len, int is_for, int is_rev);
char *rb3_seq_read1(rb3_seqio_t *fp, int64_t *len, const char **name);

void rb3_char2nt6(int64_t l, uint8_t *s);
void rb3_revcomp6(int64_t l, uint8_t *s);
void rb3_reverse_all(int64_t len, uint8_t *seq);

void rb3_sprintf_lite(kstring_t *s, const char *fmt, ...);

rb3_seqlist_t *rb3_sl_read(const char *fn);
void rb3_sl_destroy(rb3_seqlist_t *sl);

#ifdef __cplusplus
}
#endif

#endif
