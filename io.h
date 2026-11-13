#ifndef RB3_IO_H
#define RB3_IO_H

#include <stdint.h>
#include <stdio.h>
#include "rb3priv.h" // for kstring_t

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int64_t tot_len;
	int64_t n_seq;
	char **name;    // NULL when the table is mapped from a .sid file
	int32_t *len;   // points into the mapping when the table is mapped
	// set only when the table is mapped from a .sid file (rb3_sid_map())
	const char *names;    // the names, each ending in a NUL byte
	const uint64_t *off;  // name i starts at names + off[i]; n_seq + 1 entries
	void *map;            // the whole mapping
	size_t map_len;
} rb3_sid_t;

static inline const char *rb3_sid_name(const rb3_sid_t *sl, int64_t i)
{
	return sl->map? sl->names + sl->off[i] : sl->name[i];
}

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

int64_t rb3_sprintf_lite(kstring_t *s, const char *fmt, ...);
void rb3_str_append(kstring_t *s, const char *st, const char *en);

rb3_sid_t *rb3_sid_read(const char *fn);
rb3_sid_t *rb3_sid_map(const char *fn);
int64_t rb3_sid_dump(const char *fn, FILE *fp);
void rb3_sid_destroy(rb3_sid_t *sl);

#ifdef __cplusplus
}
#endif

#endif
