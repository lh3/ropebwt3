#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rb3priv.h"
#include "libsais16.h"
#include "libsais16x64.h"

static const int32_t sais_extra_len = 10000;

void rb3_build_sais32x32(int64_t n_seq, int64_t len, char *seq, int n_threads)
{
	int32_t i, k = 0, *tmp, *SA;
	assert(len + sais_extra_len <= INT32_MAX);
	tmp = RB3_MALLOC(int32_t, len);
	for (i = 0; i < len; ++i)
		tmp[i] = seq[i]? n_seq + seq[i] : ++k;
	SA = RB3_MALLOC(int32_t, len + sais_extra_len);
#ifdef LIBSAIS_OPENMP
	if (n_threads > 1)
		libsais16_int_omp(tmp, SA, len, n_seq + 6, sais_extra_len, n_threads);
	else
		libsais16_int(tmp, SA, len, n_seq + 6, sais_extra_len);
#else
	libsais16_int(tmp, SA, len, n_seq + 6, sais_extra_len);
#endif
	for (i = 0; i < len; ++i) {
		int32_t k = SA[i] == 0? len - 1 : SA[i] - 1;
		tmp[i] = seq[k];
	}
	for (i = 0; i < len; ++i)
		seq[i] = tmp[i];
	free(SA); free(tmp);
}

void rb3_build_sais64x64(int64_t n_seq, int64_t len, char *seq, int n_threads)
{
	int64_t i, k = 0, *tmp, *SA;
	tmp = RB3_MALLOC(int64_t, len);
	for (i = 0; i < len; ++i)
		tmp[i] = seq[i]? n_seq + seq[i] : ++k;
	SA = RB3_MALLOC(int64_t, len + sais_extra_len);
#ifdef LIBSAIS_OPENMP
	if (n_threads > 1)
		libsais16x64_long_omp(tmp, SA, len, n_seq + 6, sais_extra_len, n_threads);
	else
		libsais16x64_long(tmp, SA, len, n_seq + 6, sais_extra_len);
#else
	libsais16x64_long(tmp, SA, len, n_seq + 6, sais_extra_len);
#endif
	for (i = 0; i < len; ++i) {
		int64_t k = SA[i] == 0? len - 1 : SA[i] - 1;
		tmp[i] = seq[k];
	}
	for (i = 0; i < len; ++i)
		seq[i] = tmp[i];
	free(SA); free(tmp);
}

void rb3_build_sais16x64(int64_t n_seq, int64_t len, char *seq, int n_threads)
{
	int64_t i, k = 0, *SA;
	uint16_t *tmp;
	assert(n_seq + 5 < UINT16_MAX);
	tmp = RB3_MALLOC(uint16_t, len);
	for (i = 0; i < len; ++i)
		tmp[i] = seq[i]? n_seq + seq[i] : ++k;
	SA = RB3_MALLOC(int64_t, len + sais_extra_len);
#ifdef LIBSAIS_OPENMP
	if (n_threads > 1)
		libsais16x64_omp(tmp, SA, len, sais_extra_len, 0, n_threads);
	else
		libsais16x64(tmp, SA, len, sais_extra_len, 0);
#else
	libsais16x64(tmp, SA, len, sais_extra_len, 0);
#endif
	for (i = 0; i < len; ++i) {
		int64_t k = SA[i] == 0? len - 1 : SA[i] - 1;
		tmp[i] = seq[k];
	}
	for (i = 0; i < len; ++i)
		seq[i] = tmp[i];
	free(SA); free(tmp);
}

void rb3_build_sais(int64_t n_seq, int64_t len, char *seq, int n_threads)
{
	if (n_seq + 5 < UINT16_MAX)
		rb3_build_sais16x64(n_seq, len, seq, n_threads);
	else if (len + sais_extra_len >= INT32_MAX)
		rb3_build_sais64x64(n_seq, len, seq, n_threads);
	else
		rb3_build_sais32x32(n_seq, len, seq, n_threads);
}
