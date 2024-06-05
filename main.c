#include <stdio.h>
#include <string.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "ketopt.h"

#define RB3_VERSION "3.0pre-r9"

int main_sais(int argc, char *argv[]);
int main_merge_fmi(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: ropebwt3 <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  sais       construct multi-string SA with libsais\n");
	fprintf(fp, "  merge-fmi  merge fm-indices\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	rb3_realtime();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "sais") == 0) main_sais(argc-1, argv+1);
	else if (strcmp(argv[1], "merge-fmi") == 0) main_merge_fmi(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		printf("%s\n", RB3_VERSION);
		return 0;
	} else {
		fprintf(stderr, "ERROR: unknown command '%s'\n", argv[1]);
		return 1;
	}

	if (rb3_verbose >= 3 && argc > 2) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, RB3_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, rb3_realtime(), rb3_cputime(), rb3_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}

int main_merge_fmi(int argc, char *argv[])
{
	int32_t c, i, n_threads = 1;
	ketopt_t o = KETOPT_INIT;
	mrope_t *r;

	while ((c = ketopt(&o, argc, argv, 1, "t:o:", 0)) >= 0) {
		if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'o') freopen(o.arg, "wb", stdout);
	}
	if (argc - o.ind < 2) {
		fprintf(stdout, "Usage: ropebwt3 merge-fmi [options] <base.fmr> <other1.fmd> [...]\n");
		fprintf(stdout, "Options:\n");
		fprintf(stdout, "  -t INT     number of threads [%d]\n", n_threads);
		fprintf(stdout, "  -o FILE    output FMR to FILE [stdout]\n");
		return 0;
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
	}
	mr_dump(r, stdout);
	mr_destroy(r);
	return 0;
}
