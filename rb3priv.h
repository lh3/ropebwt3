#ifndef RB3PRIV_H
#define RB3PRIV_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef KSTRING_T // same as kstring_t in kseq.h
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#define RB3_MALLOC(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define RB3_CALLOC(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define RB3_REALLOC(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

#define RB3_GROW(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = RB3_REALLOC(type, ptr, (__m)); \
		} \
	} while (0)

#define RB3_GROW0(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			size_t old_m = (__m); \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = RB3_REALLOC(type, ptr, (__m)); \
			memset((ptr) + old_m, 0, ((__m) - old_m) * sizeof(type)); \
		} \
	} while (0)

extern int rb3_verbose;

// in sys.c
double rb3_cputime(void);
double rb3_realtime(void);
double rb3_percent_cpu(void);
long rb3_peakrss(void);
int64_t rb3_parse_num(const char *str);

// in sais-ss.c
void rb3_build_sais(int64_t n_seq, int64_t len, char *seq, int n_threads, int force64);

#ifdef __cplusplus
}
#endif

#endif
