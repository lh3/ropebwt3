#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <zlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include "rb3priv.h"
#include "io.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

const uint8_t rb3_nt6_table[128] = {
    0, 1, 2, 3,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

void rb3_char2nt6(int64_t l, uint8_t *s)
{
	int64_t i;
	for (i = 0; i < l; ++i)
		s[i] = s[i] < 128? rb3_nt6_table[s[i]] : 5;
}

void rb3_revcomp6(int64_t l, uint8_t *s)
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

char *rb3_seq_read1(rb3_seqio_t *fp, int64_t *len, const char **name)
{
	int ret, dret;
	char *s = 0;
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

/***********
 * seqlist *
 ***********/

// One line of a .len file holds a name and a length, separated by a space or a
// tab; anything after them is ignored. Returns the length, or -1 for a line
// that holds no sequence. The reader and the .sid writer both use it, so they
// keep the same sequences in the same order.
static int64_t sid_parse(char *s, char **name)
{
	int32_t i;
	char *p, *q;
	int64_t len = -1;
	*name = 0;
	for (p = q = s, i = 0;; ++p) {
		if (*p == ' ' || *p == '\t' || *p == 0) {
			int32_t c = *p;
			*p = 0;
			if (i == 0) {
				*name = q;
			} else if (i == 1) {
				len = atol(q);
			}
			++i, q = p + 1;
			if (c == 0 || i == 2) break;
		}
	}
	return i == 2 && len > 0? len : -1;
}

rb3_sid_t *rb3_sid_read(const char *fn)
{
	rb3_sid_t *sl;
	gzFile fp;
	kstream_t *ks;
	int32_t l, dret;
	int64_t m_seq = 0;
	kstring_t str = {0,0,0};

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	sl = RB3_CALLOC(rb3_sid_t, 1);
	while ((l = ks_getuntil(ks, KS_SEP_LINE, &str, &dret)) >= 0) {
		char *name;
		int64_t len = sid_parse(str.s, &name);
		if (len > 0) {
			assert(len <= INT32_MAX);
			RB3_GROW(char*, sl->name, sl->n_seq, m_seq);
			sl->len = RB3_REALLOC(int32_t, sl->len, m_seq);
			sl->name[sl->n_seq] = rb3_strdup(name);
			sl->len[sl->n_seq++] = len;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return sl;
}

/* rb3_sid_dump() writes the names and lengths of a .len file (gzip'd or not)
 * as a .sid file, which rb3_sid_map() maps instead of parsing. Little-endian,
 * like the .ssa:
 *   0    "SID\1" and 4 bytes of 0
 *   8    n_seq
 *   16   len[n_seq], 32-bit, then 4 bytes of 0 if n_seq is odd
 *        off[n_seq+1], 64-bit: name i starts at byte off[i] of the names
 *        the names, each ending in a NUL byte; off[n_seq] is their total size
 * It returns the number of sequences, or -1 on a read, write or memory error. */
int64_t rb3_sid_dump(const char *fn, FILE *fo)
{
	static const char magic[8] = { 'S', 'I', 'D', 1, 0, 0, 0, 0 };
	static const char zero[4] = { 0, 0, 0, 0 };
	gzFile fp;
	kstream_t *ks;
	int32_t l, dret;
	int zerr;
	int64_t n = 0, m = 1<<20, ret = -1;
	uint64_t l_names = 0, m_names = 1<<24;
	int32_t *len;
	uint64_t *off;
	char *names;
	kstring_t str = {0,0,0};

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) return -1;
	ks = ks_init(fp);
	len = RB3_MALLOC(int32_t, m);
	off = RB3_MALLOC(uint64_t, m + 1);
	names = RB3_MALLOC(char, m_names);
	if (len == 0 || off == 0 || names == 0) goto end;
	while ((l = ks_getuntil(ks, KS_SEP_LINE, &str, &dret)) >= 0) {
		char *name;
		int64_t ln = sid_parse(str.s, &name);
		uint64_t nl;
		if (ln < 0) continue;
		assert(ln <= INT32_MAX);
		nl = strlen(name) + 1;
		if (n == m) {
			int32_t *len2;
			uint64_t *off2;
			m += m>>1;
			if ((len2 = RB3_REALLOC(int32_t, len, m)) == 0) goto end;
			len = len2;
			if ((off2 = RB3_REALLOC(uint64_t, off, m + 1)) == 0) goto end;
			off = off2;
		}
		if (l_names + nl > m_names) {
			char *names2;
			m_names = l_names + nl;
			m_names += m_names>>1;
			if ((names2 = RB3_REALLOC(char, names, m_names)) == 0) goto end;
			names = names2;
		}
		len[n] = ln, off[n++] = l_names;
		memcpy(names + l_names, name, nl);
		l_names += nl;
	}
	gzerror(fp, &zerr);
	if (zerr != Z_OK) goto end; // a truncated or damaged gzip file
	off[n] = l_names;
	if (fwrite(magic, 1, 8, fo) != 8 || fwrite(&n, 8, 1, fo) != 1) goto end;
	if (n > 0 && fwrite(len, 4, n, fo) != (size_t)n) goto end;
	if ((n & 1) && fwrite(zero, 1, 4, fo) != 4) goto end;
	if (fwrite(off, 8, n + 1, fo) != (size_t)(n + 1)) goto end;
	if (l_names > 0 && fwrite(names, 1, l_names, fo) != l_names) goto end;
	if (fflush(fo) != 0) goto end;
	ret = n;
end:
	free(len); free(off); free(names); free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return ret;
}

/* rb3_sid_map() maps a .sid file read-only. It returns NULL when the file is
 * missing, too short, or its sizes disagree with its header, as a truncated
 * copy would. */
rb3_sid_t *rb3_sid_map(const char *fn)
{
	int fd;
	struct stat st;
	uint8_t *a;
	uint64_t n, o_off, o_names, l_names, size;
	rb3_sid_t *sl;

	if ((fd = open(fn, O_RDONLY)) < 0) return 0;
	if (fstat(fd, &st) < 0 || st.st_size < 24) {
		close(fd);
		return 0;
	}
	size = st.st_size;
	a = (uint8_t*)mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	close(fd);
	if (a == (uint8_t*)MAP_FAILED) return 0;
	memcpy(&n, a + 8, 8);
	if (memcmp(a, "SID\1\0\0\0\0", 8) != 0 || n > (size - 24) / 12) goto fail; // a sequence takes 12 bytes or more
	o_off = 16 + (n * 4 + 7) / 8 * 8;
	o_names = o_off + (n + 1) * 8;
	if (o_names > size) goto fail;
	memcpy(&l_names, a + o_off + n * 8, 8);
	if (l_names != size - o_names || (l_names > 0 && a[size - 1] != 0)) goto fail;
	sl = RB3_CALLOC(rb3_sid_t, 1);
	sl->n_seq = n;
	sl->len = (int32_t*)(a + 16);
	sl->off = (const uint64_t*)(a + o_off);
	sl->names = (const char*)(a + o_names);
	sl->map = a, sl->map_len = size;
	return sl;
fail:
	munmap(a, size);
	return 0;
}

void rb3_sid_destroy(rb3_sid_t *sl)
{
	int64_t i;
	if (sl == 0) return;
	if (sl->map) {
		munmap(sl->map, sl->map_len);
	} else {
		for (i = 0; i < sl->n_seq; ++i)
			free(sl->name[i]);
		free(sl->name); free(sl->len);
	}
	free(sl);
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

void rb3_str_append(kstring_t *s, const char *st, const char *en)
{
	str_copy(s, st, en);
}

int64_t rb3_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[32]; // for integer to string conversion
	const char *p, *q;
	int64_t len = 0;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) {
				len += p - q;
				if (s) str_copy(s, q, p);
			}
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				len += l;
				if (s) {
					str_enlarge(s, l);
					for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				}
			} else if (*p == 'l' && *(p+1) == 'd') {
				int i, l = 0;
				long int c;
				unsigned long x;
				c = va_arg(ap, long);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				len += l;
				if (s) {
					str_enlarge(s, l);
					for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				}
				++p;
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				len += l;
				if (s) {
					str_enlarge(s, l);
					for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				}
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				int l;
				l = strlen(r);
				len += l;
				if (s) str_copy(s, r, r + l);
			} else if (*p == 'c') {
				++len;
				if (s) {
					str_enlarge(s, 1);
					s->s[s->l++] = va_arg(ap, int);
				}
			} else {
				fprintf(stderr, "ERROR: unrecognized type '%%%c'\n", *p);
				abort();
			}
			q = p + 1;
		}
	}
	if (p > q) {
		len += p - q;
		if (s) str_copy(s, q, p);
	}
	va_end(ap);
	if (s) s->s[s->l] = 0;
	return len;
}
