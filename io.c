#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "rb3priv.h"
#include "io.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static const uint8_t rb3_nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

static inline void rb3_char2nt6(int64_t l, uint8_t *s)
{
	int64_t i;
	for (i = 0; i < l; ++i)
		s[i] = s[i] < 128? rb3_nt6_table[s[i]] : 5;
}

static inline void rb3_revcomp6(int64_t l, uint8_t *s)
{
	int64_t i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
		s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
		s[i] = tmp;
	}
	if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}

struct rb3_seqio_s {
	int32_t is_line;
	kseq_t *fx;
	kstream_t *fl;
	gzFile fp;
};

rb3_seqio_t *rb3_seq_open(const char *fn, int is_line)
{
	rb3_seqio_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (f == 0) return 0;
	fp = RB3_CALLOC(rb3_seqio_t, 1);
	fp->fp = f;
	fp->is_line = !!is_line;
	if (is_line) fp->fl = ks_init(f);
	else fp->fx = kseq_init(f);
	return fp;
}

void rb3_seq_close(rb3_seqio_t *fp)
{
	if (fp == 0) return;
	if (fp->is_line) ks_destroy(fp->fl);
	else kseq_destroy(fp->fx);
	gzclose(fp->fp);
	free(fp);
}

int64_t rb3_seq_read(rb3_seqio_t *fp, kstring_t *seq, int64_t max_len, int is_for, int is_rev)
{
	int64_t n_seq = 0;
	assert(is_for || is_rev);
	seq->l = 0;
	if (fp->is_line) { // TODO: implement this
		abort();
	} else {
		int32_t ret;
		while ((ret = kseq_read(fp->fx)) >= 0) {
			rb3_char2nt6(fp->fx->seq.l, (uint8_t*)fp->fx->seq.s);
			if (is_for) {
				RB3_GROW(char, seq->s, seq->l + fp->fx->seq.l + 1, seq->m);
				memcpy(&seq->s[seq->l], fp->fx->seq.s, fp->fx->seq.l + 1);
				seq->l += fp->fx->seq.l + 1;
				++n_seq;
			}
			if (is_rev) {
				rb3_revcomp6(fp->fx->seq.l, (uint8_t*)fp->fx->seq.s);
				RB3_GROW(char, seq->s, seq->l + fp->fx->seq.l + 1, seq->m);
				memcpy(&seq->s[seq->l], fp->fx->seq.s, fp->fx->seq.l + 1);
				seq->l += fp->fx->seq.l + 1;
				++n_seq;
			}
			if (max_len > 0 && seq->l > max_len) break;
		}
		if (ret < -1 && rb3_verbose >= 1)
			fprintf(stderr, "ERROR: FASTX parsing error (code %d)\n", ret);
	}
	return n_seq;
}
