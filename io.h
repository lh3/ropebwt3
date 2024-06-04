#ifndef RB3_IO_H
#define RB3_IO_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct rb3_seqio_s;
typedef struct rb3_seqio_s rb3_seqio_t;

rb3_seqio_t *rb3_seq_open(const char *fn, int is_line);
void rb3_seq_close(rb3_seqio_t *fp);
int64_t rb3_seq_read(rb3_seqio_t *fp, kstring_t *seq, int64_t max_len, int is_for, int is_rev);

#ifdef __cplusplus
}
#endif

#endif
