#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "kalloc.h"
#include "kthread.h"
#include "ketopt.h"
#include "ksort.h"

/********************
 * ssa construction *
 ********************/

typedef struct { size_t n, m; uint64_t *a; } uint64_v;

static void ssa_gen1(void *km, const rb3_fmi_t *f, rb3_ssa_t *sa, int64_t k, uint64_v *buf)
{
	int32_t c, mask = (1<<sa->ss) - 1;
	int64_t ok[RB3_ASIZE], k0 = k, l = 0;
	size_t i;
	buf->n = 0;
	do {
		++l;
		c = rb3_fmi_rank1a(f, k, ok);
		k = f->acc[c] + ok[c];
		if (c) {
			if (((k - f->acc[1]) & mask) == 0) {
				int64_t x = (k - f->acc[1]) >> sa->ss;
				assert(x < sa->n_ssa);
				sa->ssa[x] = l;
				Kgrow(km, uint64_t, buf->a, buf->n, buf->m);
				buf->a[buf->n++] = x;
			}
		} else sa->r2i[k] = k0;
	} while (c);
	for (i = 0; i < buf->n; ++i)
		sa->ssa[buf->a[i]] = (l - 1 - sa->ssa[buf->a[i]]) << sa->ms | k0;
}

typedef struct {
	void **km;
	const rb3_fmi_t *f;
	rb3_ssa_t *sa;
	uint64_v *buf;
} worker_t;

static void worker(void *data, long i, int tid)
{
	worker_t *w = (worker_t*)data;
	ssa_gen1(w->km[tid], w->f, w->sa, i, &w->buf[tid]);
}

rb3_ssa_t *rb3_ssa_gen(const rb3_fmi_t *f, int ssa_shift, int n_threads)
{
	rb3_ssa_t *sa;
	worker_t *w;
	int i;

	sa = RB3_CALLOC(rb3_ssa_t, 1);
	sa->ss = ssa_shift;
	sa->m = f->acc[1];
	for (sa->ms = 1; 1LL<<sa->ms < sa->m; ++sa->ms) {}
	sa->n_ssa = (f->acc[RB3_ASIZE] - f->acc[1] + (1LL<<sa->ss) - 1LL) >> sa->ss;
	sa->r2i = RB3_CALLOC(uint64_t, sa->m);
	sa->ssa = RB3_CALLOC(uint64_t, sa->n_ssa);

	w = RB3_CALLOC(worker_t, 1);
	w->km = RB3_CALLOC(void*, n_threads);
	for (i = 0; i < n_threads; ++i)
		w->km[i] = km_init();
	w->buf = RB3_CALLOC(uint64_v, n_threads);
	w->sa = sa, w->f = f;

	kt_for(n_threads, worker, w, sa->m);

	for (i = 0; i < n_threads; ++i)
		km_destroy(w->km[i]);
	free(w->km); free(w->buf); free(w);
	return sa;
}

void rb3_ssa_destroy(rb3_ssa_t *sa)
{
	if (sa == 0) return;
	if (sa->map) munmap(sa->map, sa->map_len);
	else free(sa->r2i), free(sa->ssa);
	free(sa);
}

/**************
 * Compute SA *
 **************/

int64_t rb3_ssa(const rb3_fmi_t *f, const rb3_ssa_t *sa, int64_t k, int64_t *si)
{
	int32_t c, mask = (1<<sa->ss) - 1;
	int64_t x = 0;
	int64_t ok[RB3_ASIZE];
	uint64_t v;
	*si = -1;
	if (k >= f->acc[6]) return -1;
	while (k < f->acc[1] || ((k - f->acc[1]) & mask)) {
		++x;
		c = rb3_fmi_rank1a(f, k, ok);
		k = f->acc[c] + ok[c];
		if (c == 0) {
			*si = rb3_ssa_r2i(sa, k);
			return x - 1;
		}
	}
	k = (k - f->acc[1]) >> sa->ss;
	v = rb3_ssa_val(sa, k);
	*si = v & ((1ULL<<sa->ms) - 1);
	return x + (v >> sa->ms);
}

typedef struct {
	int64_t off;
	int64_t lo, hi;
} ssa_intv_t;

#define intv_lt(x, y) ((x).hi - (x).lo < (y).hi - (y).lo)
KSORT_INIT(ssa_intv, ssa_intv_t, intv_lt)

typedef struct {
	int64_t n_sa, max_sa, n0;
	int32_t n_a, m_a;
	ssa_intv_t *a;
	rb3_pos_t *sa;
	void *km;
} ssa_aux_t;

static inline void ssa_add_intv1(ssa_aux_t *aux, int64_t lo, int64_t hi, int64_t off)
{
	Kgrow(aux->km, ssa_intv_t, aux->a, aux->n_a, aux->m_a);
	aux->a[aux->n_a].off = off, aux->a[aux->n_a].lo = lo, aux->a[aux->n_a].hi = hi;
	aux->n_a++;
	ks_heapup_ssa_intv(aux->n_a, aux->a);
}

static int32_t ssa_add_intv(const rb3_ssa_t *ssa, ssa_aux_t *aux, int64_t lo, int64_t hi, int64_t off)
{
	int64_t m = aux->n0;
	int64_t k = ((lo - m) >> ssa->ss << ssa->ss) + m;
	if (aux->n_sa == aux->max_sa) return -1;
	for (; k < hi; k += 1LL << ssa->ss) {
		int64_t l = (k - m) >> ssa->ss;
		uint64_t v;
		if (k < lo) continue;
		assert(l < ssa->n_ssa && aux->n_sa < aux->max_sa);
		v = rb3_ssa_val(ssa, l);
		aux->sa[aux->n_sa].sid = v & ((1LL << ssa->ms) - 1);
		aux->sa[aux->n_sa].pos = off + (v >> ssa->ms);
		aux->n_sa++;
		if (aux->n_sa == aux->max_sa) return -1;
		if (lo < k) ssa_add_intv1(aux, lo, k, off);
		lo = k + 1;
	}
	ssa_add_intv1(aux, lo, hi, off);
	return 0;
}

int64_t rb3_ssa_multi(void *km, const rb3_fmi_t *f, const rb3_ssa_t *ssa, int64_t lo, int64_t hi, int64_t max_sa, rb3_pos_t *sa)
{
	ssa_aux_t aux;
	int64_t ok[RB3_ASIZE], ol[RB3_ASIZE];
	if (max_sa == 0 || lo >= hi) return 0;
	memset(&aux, 0, sizeof(aux));
	aux.max_sa = max_sa < hi - lo? max_sa : hi - lo;
	aux.m_a = 256, aux.n_a = 0;
	aux.a = Kmalloc(km, ssa_intv_t, aux.m_a);
	aux.km = km, aux.sa = sa, aux.n0 = f->acc[1];
	ssa_add_intv(ssa, &aux, lo, hi, 0);
	while (aux.n_a > 0 && aux.n_sa < aux.max_sa) {
		int64_t l;
		int32_t c;
		ssa_intv_t x = aux.a[0];
		--aux.n_a;
		if (aux.n_a > 0) { // maintain heap
			aux.a[0] = aux.a[aux.n_a];
			ks_heapdown_ssa_intv(0, aux.n_a, aux.a);
		}
		rb3_fmi_rank2a(f, x.lo, x.hi, ok, ol);
		for (l = ok[0]; l < ol[0]; ++l) { // reaching sentinels
			aux.sa[aux.n_sa].sid = rb3_ssa_r2i(ssa, l);
			aux.sa[aux.n_sa].pos = x.off;
			aux.n_sa++;
			if (aux.n_sa == aux.max_sa) goto end_ssa_multi;
		}
		for (c = 1; c < 6; ++c)
			if (ok[c] < ol[c])
				ssa_add_intv(ssa, &aux, f->acc[c] + ok[c], f->acc[c] + ol[c], x.off + 1);
	}
end_ssa_multi:
	kfree(km, aux.a);
	return aux.n_sa;
}

/***********
 * ssa I/O *
 ***********/

int rb3_ssa_dump(const rb3_ssa_t *sa, const char *fn)
{
	uint32_t y;
	FILE *fp;
	fp = fn && strcmp(fn, "-")? fopen(fn, "wb") : fdopen(1, "wb");
	if (fp == 0) return -1;
	fwrite("SSA\1", 1, 4, fp);
	y = sa->ss; fwrite(&y, 4, 1, fp);
	y = sa->ms; fwrite(&y, 4, 1, fp);
	fwrite(&sa->m, 8, 1, fp);
	fwrite(&sa->n_ssa, 8, 1, fp);
	if (sa->map) { // the arrays as they are in the mapped file
		fwrite(sa->r2i_m, 8, sa->m, fp);
		fwrite(sa->ssa_m, 8, sa->n_ssa, fp);
	} else {
		fwrite(sa->r2i, 8, sa->m, fp);
		fwrite(sa->ssa, 8, sa->n_ssa, fp);
	}
	fclose(fp);
	return 0;
}

rb3_ssa_t *rb3_ssa_restore(const char *fn)
{
	FILE *fp;
	uint32_t y;
	char magic[4];
	rb3_ssa_t *sa;

	fp = fn && strcmp(fn, "-")? fopen(fn, "rb") : fdopen(0, "rb");
	if (fp == 0) return 0;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, "SSA\1", 4) != 0) return 0; // wrong magic
	sa = RB3_CALLOC(rb3_ssa_t, 1);
	fread(&y, 4, 1, fp); sa->ss = y;
	fread(&y, 4, 1, fp); sa->ms = y;
	fread(&sa->m, 8, 1, fp);
	fread(&sa->n_ssa, 8, 1, fp);
	sa->r2i = RB3_CALLOC(uint64_t, sa->m);
	sa->ssa = RB3_CALLOC(uint64_t, sa->n_ssa);
	if (sa->ssa == 0 || sa->r2i == 0) {
		free(sa->r2i); free(sa->ssa); free(sa);
		return 0;
	}
	fread(sa->r2i, 8, sa->m, fp);
	fread(sa->ssa, 8, sa->n_ssa, fp);
	fclose(fp);
	return sa;
}

/* rb3_ssa_map() maps a .ssa file read-only, as rb3_ssa_dump() wrote it, instead
 * of reading it into memory. On a RAM disk the pages are the RAM disk's own, so
 * every process shares one copy and nothing is read at load. It returns NULL
 * when the file is missing, or its size isn't 28 + 8 * (m + n_ssa) bytes, as a
 * truncated copy's would not be. */
rb3_ssa_t *rb3_ssa_map(const char *fn)
{
	int fd;
	struct stat st;
	uint8_t *a;
	uint32_t ss, ms;
	int64_t m, n_ssa;
	uint64_t n;
	rb3_ssa_t *sa;

	if ((fd = open(fn, O_RDONLY)) < 0) return 0;
	if (fstat(fd, &st) < 0 || st.st_size < 28) {
		close(fd);
		return 0;
	}
	a = (uint8_t*)mmap(0, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
	close(fd);
	if (a == (uint8_t*)MAP_FAILED) return 0;
	memcpy(&ss, a + 4, 4);
	memcpy(&ms, a + 8, 4);
	memcpy(&m, a + 12, 8);
	memcpy(&n_ssa, a + 20, 8);
	n = (st.st_size - 28) / 8; // the number of 64-bit entries the file holds
	// rb3_ssa() masks with (1<<ss) - 1 in an int, and ssa_add_intv() with
	// (1LL<<ms) - 1, so larger shifts would overflow
	if (memcmp(a, "SSA\1", 4) != 0 || ss > 30 || ms < 1 || ms > 62 || (st.st_size - 28) % 8 != 0
		|| m < 0 || n_ssa < 0 || (uint64_t)m > n || (uint64_t)n_ssa != n - m) {
		munmap(a, st.st_size);
		return 0;
	}
	sa = RB3_CALLOC(rb3_ssa_t, 1);
	sa->ss = ss, sa->ms = ms, sa->m = m, sa->n_ssa = n_ssa;
	sa->r2i_m = a + 28;
	sa->ssa_m = a + 28 + 8 * m;
	sa->map = a, sa->map_len = st.st_size;
	return sa;
}

/*******************
 * main() function *
 *******************/

int main_ssa(int argc, char *argv[])
{
	int c, n_threads = 4, ssa_shift = 8;
	rb3_ssa_t *sa;
	rb3_fmi_t f;
	char *fn = 0;
	ketopt_t o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, "t:s:o:", 0)) >= 0) {
		if (c == 't') n_threads = atoi(o.arg);
		else if (c == 's') ssa_shift = atoi(o.arg);
		else if (c == 'o') fn = o.arg;
	}
	if (argc == o.ind) {
		fprintf(stderr, "Usage: ropebwt3 ssa [options] <in.fmd>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t INT     number of threads [%d]\n", n_threads);
		fprintf(stderr, "  -s INT     sample rate one SA per 2**INT bases [%d]\n", ssa_shift);
		fprintf(stderr, "  -o FILE    output to file [stdout]\n");
		return 1;
	}
	rb3_fmi_restore(&f, argv[o.ind], 0);
	if (f.e == 0 && f.r == 0) {
		fprintf(stderr, "[E::%s] failed to load the FM-index\n", __func__);
		return 1;
	}
	sa = rb3_ssa_gen(&f, ssa_shift, n_threads);

	rb3_ssa_dump(sa, fn);
	rb3_fmi_free(&f);
	rb3_ssa_destroy(sa);
	return 0;
}
