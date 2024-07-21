#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "kalloc.h"
#include "kthread.h"
#include "ketopt.h"

typedef struct { size_t n, m; uint64_t  *a; } uint64_v;

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
	free(sa->r2i); free(sa->ssa); free(sa);
}

int64_t rb3_ssa(const rb3_fmi_t *f, const rb3_ssa_t *sa, int64_t k, int64_t *si)
{
	int32_t c, mask = (1<<sa->ss) - 1;
	int64_t x = 0;
	int64_t ok[RB3_ASIZE];
	*si = -1;
	if (k >= f->acc[6]) return -1;
	while (k < f->acc[1] || ((k - f->acc[1]) & mask)) {
		++x;
		c = rb3_fmi_rank1a(f, k, ok);
		k = f->acc[c] + ok[c];
		if (c == 0) {
			*si = sa->r2i[k];
			return x - 1;
		}
	}
	k = (k - f->acc[1]) >> sa->ss;
	*si = sa->ssa[k] & ((1ULL<<sa->ms) - 1);
	return x + (sa->ssa[k] >> sa->ms);
}

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
	fwrite(sa->r2i, 8, sa->m, fp);
	fwrite(sa->ssa, 8, sa->n_ssa, fp);
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

int main_ssa(int argc, char *argv[])
{
	int c, n_threads = 4, ssa_shift = 6;
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
		fprintf(stderr, "  -s INT     sample rate one SA per 1<<INT bases [%d]\n", ssa_shift);
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
	rb3_fmi_destroy(&f);
	rb3_ssa_destroy(sa);
	return 0;
}
