#include <stdio.h>
#include <string.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "ketopt.h"

#define RB3_VERSION "3.0pre-r38"

int main_build(int argc, char *argv[]);
int main_merge(int argc, char *argv[]);
int main_get(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: ropebwt3 <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  build      construct a BWT\n");
	fprintf(fp, "  merge      merge BWTs\n");
	fprintf(fp, "  get        retrieve i-th sequence from BWT\n");
	fprintf(fp, "  version    print the version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	rb3_realtime();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "build") == 0) ret = main_build(argc-1, argv+1);
	else if (strcmp(argv[1], "merge") == 0) ret = main_merge(argc-1, argv+1);
	else if (strcmp(argv[1], "get") == 0) ret = main_get(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		printf("%s\n", RB3_VERSION);
		return 0;
	} else {
		fprintf(stderr, "ERROR: unknown command '%s'\n", argv[1]);
		return 1;
	}

	if (rb3_verbose >= 3 && argc > 2 && ret == 0) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, RB3_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, rb3_realtime(), rb3_cputime(), rb3_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}

int main_merge(int argc, char *argv[])
{
	int32_t c, i, n_threads = 1;
	ketopt_t o = KETOPT_INIT;
	mrope_t *r;
	char *fn_tmp = 0;

	while ((c = ketopt(&o, argc, argv, 1, "t:o:S:", 0)) >= 0) {
		if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'o') freopen(o.arg, "wb", stdout);
		else if (c == 'S') fn_tmp = o.arg;
	}
	if (argc - o.ind < 2) {
		fprintf(stdout, "Usage: ropebwt3 merge [options] <base.fmr> <other1.fmr> [...]\n");
		fprintf(stdout, "Options:\n");
		fprintf(stdout, "  -t INT     number of threads [%d]\n", n_threads);
		fprintf(stdout, "  -o FILE    output FMR to FILE [stdout]\n");
		fprintf(stderr, "  -S FILE    save the current index to FILE after each input file []\n");
		return 1;
	}

	r = mr_restore_file(argv[o.ind]);
	if (r == 0) {
		if (rb3_verbose >= 1)
			fprintf(stderr, "ERROR: failed to load FMR file '%s'\n", argv[o.ind]);
		return 1;
	}
	for (i = o.ind + 1; i < argc; ++i) {
		rb3_fmi_t fb;
		rb3_fmi_restore(&fb, argv[i]);
		if (fb.e == 0 && fb.r == 0) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: failed to load FMR/FMD file '%s'\n", argv[i]);
			break;
		}
		rb3_fmi_merge(r, &fb, n_threads, 1);
		if (fn_tmp) {
			FILE *fp;
			fp = fopen(fn_tmp, "w");
			if (fp != 0) mr_dump(r, fp);
			fclose(fp);
			if (rb3_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] saved the current index to '%s'\n", __func__, rb3_realtime(), rb3_percent_cpu(), fn_tmp);
		}
	}
	mr_dump(r, stdout);
	mr_destroy(r);
	return 0;
}

int main_get(int argc, char *argv[])
{
	int32_t c, i;
	ketopt_t o = KETOPT_INIT;
	rb3_fmi_t fmi;
	kstring_t s = {0,0,0};

	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) { }
	if (argc - o.ind < 2) {
		fprintf(stdout, "Usage: ropebwt3 get <idx.fmr> <int> [...]\n");
		return 0;
	}
	rb3_fmi_restore(&fmi, argv[o.ind]);
	if (fmi.e == 0 && fmi.r == 0) {
		if (rb3_verbose >= 1)
			fprintf(stderr, "ERROR: failed to load index file '%s'\n", argv[o.ind]);
		return 1;
	}
	for (i = o.ind + 1; i < argc; ++i) {
		int64_t k, r;
		k = atol(argv[i]);
		r = rb3_fmi_retrieve(&fmi, k, &s);
		if (r >= 0) {
			printf(">%ld %ld\n", (long)k, (long)r);
			puts(s.s);
		}
	}
	free(s.s);
	rb3_fmi_destroy(&fmi);
	return 0;
}
