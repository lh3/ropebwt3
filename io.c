#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <zlib.h>
#include "rb3priv.h"
#include "io.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

const uint8_t rb3_nt6_table[128] = {
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

static inline void rb3_reverse(int64_t l, uint8_t *s)
{
	int64_t i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		s[l-1-i] = s[i];
		s[i] = tmp;
	}
}

struct rb3_seqio_s {
	int32_t is_line;
	kseq_t *fx;
	kstream_t *fl;
	gzFile fp;
	kstring_t line_buf;
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
	free(fp->line_buf.s);
	if (fp->is_line) ks_destroy(fp->fl);
	else kseq_destroy(fp->fx);
	gzclose(fp->fp);
	free(fp);
}

static int64_t rb3_seq_add(kstring_t *seq, int is_for, int is_rev, int64_t l, char *s)
{
	int64_t n_added = 0;
	rb3_char2nt6(l, (uint8_t*)s);
	if (is_for) {
		RB3_GROW(char, seq->s, seq->l + l + 1, seq->m);
		memcpy(&seq->s[seq->l], s, l + 1); // this includes the trailing NULL
		seq->l += l + 1;
		++n_added;
	}
	if (is_rev) {
		rb3_revcomp6(l, (uint8_t*)s);
		RB3_GROW(char, seq->s, seq->l + l + 1, seq->m);
		memcpy(&seq->s[seq->l], s, l + 1);
		seq->l += l + 1;
		++n_added;
	}
	return n_added;
}

int64_t rb3_seq_read(rb3_seqio_t *fp, kstring_t *seq, int64_t max_len, int is_for, int is_rev)
{
	int64_t n_seq = 0;
	int32_t ret;
	assert(is_for || is_rev);
	seq->l = 0;
	if (fp->is_line) {
		int dret;
		while ((ret = ks_getuntil(fp->fl, KS_SEP_LINE, &fp->line_buf, &dret)) >= 0) {
			n_seq += rb3_seq_add(seq, is_for, is_rev, fp->line_buf.l, fp->line_buf.s);
			if (max_len > 0 && seq->l > max_len) break;
		}
	} else {
		while ((ret = kseq_read(fp->fx)) >= 0) {
			n_seq += rb3_seq_add(seq, is_for, is_rev, fp->fx->seq.l, fp->fx->seq.s);
			if (max_len > 0 && seq->l > max_len) break;
		}
		if (ret < -1 && rb3_verbose >= 1)
			fprintf(stderr, "ERROR: FASTX parsing error (code %d)\n", ret);
	}
	return n_seq;
}

const char *rb3_seq_read1(rb3_seqio_t *fp, int64_t *len, const char **name)
{
	int ret, dret;
	const char *s = 0;
	*len = 0;
	if (fp->is_line) {
		ret = ks_getuntil(fp->fl, KS_SEP_LINE, &fp->line_buf, &dret);
		*len = fp->line_buf.l;
		s = fp->line_buf.s;
		if (name) *name = 0;
	} else {
		ret = kseq_read(fp->fx);
		*len = fp->fx->seq.l;
		s = fp->fx->seq.s;
		if (name) *name = fp->fx->name.s;
	}
	return ret < 0? 0 : s;
}

void rb3_reverse_all(int64_t len, uint8_t *seq)
{
	int64_t i, j;
	for (i = j = 0; j < len; ++j) {
		if (seq[j] == 0) {
			if (j > i) rb3_reverse(j - i, &seq[i]);
			i = j + 1;
		}
	}
}

/**********************
 * Simplified sprintf *
 **********************/

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = RB3_REALLOC(char, s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

void rb3_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[32]; // for integer to string conversion
	const char *p, *q;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) str_copy(s, q, p);
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 'l' && *(p+1) == 'd') {
				int c, i, l = 0;
				unsigned long x;
				c = va_arg(ap, long);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				++p;
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				str_copy(s, r, r + strlen(r));
			} else if (*p == 'c') {
				str_enlarge(s, 1);
				s->s[s->l++] = va_arg(ap, int);
			} else {
				fprintf(stderr, "ERROR: unrecognized type '%%%c'\n", *p);
				abort();
			}
			q = p + 1;
		}
	}
	if (p > q) str_copy(s, q, p);
	va_end(ap);
	s->s[s->l] = 0;
}
