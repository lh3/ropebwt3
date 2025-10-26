#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rb3priv.h"
#include "libsais.h"
#include "libsais64.h"

static const int32_t sais_extra_len = 10000;

void rb3_build_sais32(int64_t n_seq, int64_t len, char *seq, int n_threads)
{
	int32_t i, *SA;
	uint8_t *T = (uint8_t*)seq;
	SA = RB3_MALLOC(int32_t, len + sais_extra_len);
#ifdef LIBSAIS_OPENMP
	if (n_threads > 1)
		libsais_gsa_omp(T, SA, len, sais_extra_len, 0, n_threads);
	else
		libsais_gsa(T, SA, len, sais_extra_len, 0);
#else
	libsais_gsa(T, SA, len, sais_extra_len, 0);
#endif
	for (i = 0; i < len; ++i)
		SA[i] = T[SA[i] == 0? len - 1 : SA[i] - 1];
	for (i = 0; i < len; ++i)
		T[i] = SA[i];
	free(SA);
}

void rb3_build_sais64(int64_t n_seq, int64_t len, char *seq, int n_threads)
{
	int64_t i, *SA;
	uint8_t *T = (uint8_t*)seq;
	SA = RB3_MALLOC(int64_t, len + sais_extra_len);
#ifdef LIBSAIS_OPENMP
	if (n_threads > 1)
		libsais64_gsa_omp(T, SA, len, sais_extra_len, 0, n_threads);
	else
		libsais64_gsa(T, SA, len, sais_extra_len, 0);
#else
	libsais64_gsa(T, SA, len, sais_extra_len, 0);
#endif
	for (i = 0; i < len; ++i)
		SA[i] = T[SA[i] == 0? len - 1 : SA[i] - 1];
	for (i = 0; i < len; ++i)
		T[i] = SA[i];
	free(SA);
}

void rb3_build_sais(int64_t n_seq, int64_t len, char *seq, int n_threads)
{
	if (len + sais_extra_len >= INT32_MAX)
		rb3_build_sais64(n_seq, len, seq, n_threads);
	else
		rb3_build_sais32(n_seq, len, seq, n_threads);
}
