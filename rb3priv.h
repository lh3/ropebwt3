#ifndef RB3PRIV_H
#define RB3PRIV_H

#include <stddef.h>
#include <stdint.h>

#define RB3_DBG_DAWG   0x1
#define RB3_DBG_SW     0x2
#define RB3_DBG_QNAME  0x4
#define RB3_DBG_BT     0x8

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

extern int rb3_verbose, rb3_dbg_flag;

// in misc.c
void rb3_init(void);

double rb3_cputime(void);
double rb3_realtime(void);
double rb3_percent_cpu(void);
long rb3_peakrss(void);

int64_t rb3_parse_num(const char *str);
char *rb3_strdup(const char *src);

// in sais-ss.c
void rb3_build_sais(int64_t n_seq, int64_t len, char *seq, int n_threads);

#ifdef __cplusplus
}
#endif

#endif
