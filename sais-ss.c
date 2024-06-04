#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rb3priv.h"
#include "io.h"
#include "libsais.h"
#include "libsais64.h"
#include "ketopt.h"

static const int32_t sais_extra_len = 10000;

void rb3_build_sais32(int64_t n_seq, int64_t len, char *seq, int n_threads)
{
	int32_t i, k = 0, *tmp, *SA;
	assert(len + sais_extra_len <= INT32_MAX);
	tmp = RB3_MALLOC(int32_t, len);
	for (i = 0; i < len; ++i)
		tmp[i] = seq[i]? n_seq + seq[i] : ++k;
	SA = RB3_MALLOC(int32_t, len + sais_extra_len);
#ifdef LIBSAIS_OPENMP
	if (n_threads > 1)
		libsais_int_omp(tmp, SA, len, n_seq + 6, sais_extra_len, n_threads);
	else
		libsais_int(tmp, SA, len, n_seq + 6, sais_extra_len);
#else
	libsais_int(tmp, SA, len, n_seq + 6, sais_extra_len);
#endif
	for (i = 0; i < len; ++i) {
		int32_t k = SA[i] == 0? len - 1 : SA[i] - 1;
		tmp[i] = seq[k];
	}
	for (i = 0; i < len; ++i)
		seq[i] = tmp[i];
	free(SA); free(tmp);
}

void rb3_build_sais64(int64_t n_seq, int64_t len, char *seq, int n_threads)
{
	int64_t i, k = 0, *tmp, *SA;
	tmp = RB3_MALLOC(int64_t, len);
	for (i = 0; i < len; ++i)
		tmp[i] = seq[i]? n_seq + seq[i] : ++k;
	SA = RB3_MALLOC(int64_t, len + sais_extra_len);
#ifdef LIBSAIS_OPENMP
	if (n_threads > 1)
		libsais64_long_omp(tmp, SA, len, n_seq + 6, sais_extra_len, n_threads);
	else
		libsais64_long(tmp, SA, len, n_seq + 6, sais_extra_len);
#else
	libsais64_long(tmp, SA, len, n_seq + 6, sais_extra_len);
#endif
	for (i = 0; i < len; ++i) {
		int64_t k = SA[i] == 0? len - 1 : SA[i] - 1;
		tmp[i] = seq[k];
	}
	for (i = 0; i < len; ++i)
		seq[i] = tmp[i];
	free(SA); free(tmp);
}

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: ropebwt3 sais [options] <arguments>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -F         no forward strand\n");
	fprintf(fp, "  -R         no reverse strand\n");
	fprintf(fp, "  -6         force to use 64-bit integers\n");
#ifdef LIBSAIS_OPENMP
	fprintf(fp, "  -t INT     number of threads\n");
#endif
	return fp == stdout? 0 : 1;
}

int main_sais(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t c, n_threads = 1, is_for = 1, is_rev = 1, is_line = 0, use64 = 0;
	rb3_fmt_t fmt = RB3_PLAIN;
	rb3_seqio_t *fp;
	kstring_t seq = {0,0,0};
	int64_t n_seq = 0;

	while ((c = ketopt(&o, argc, argv, 1, "t:FRL6", 0)) >= 0) {
		if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'F') is_for = 0;
		else if (c == 'R') is_rev = 0;
		else if (c == 'L') is_line = 1;
		else if (c == '6') use64 = 1;
	}
	if (argc == o.ind) return usage(stdout);

	fp = rb3_seq_open(argv[o.ind], is_line);
	n_seq = rb3_seq_read(fp, &seq, 0, is_for, is_rev);
	if (use64 == 0 && seq.l + sais_extra_len >= INT32_MAX) use64 = 1;
	if (use64) rb3_build_sais64(n_seq, seq.l, seq.s, n_threads);
	else rb3_build_sais32(n_seq, seq.l, seq.s, n_threads);
	if (fmt == RB3_PLAIN) {
		int64_t i;
		for (i = 0; i < seq.l; ++i)
			seq.s[i] = "$ACGTN"[(uint8_t)seq.s[i]];
		puts(seq.s);
	}
	free(seq.s);
	rb3_seq_close(fp);
	return 0;
}
