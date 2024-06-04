#include <stdio.h>
#include <string.h>
#include "rb3priv.h"

#define RB3_VERSION "3.0pre-r7"

int main_sais(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: ropebwt3 <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  sais     construct multi-string SA with libsais\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	rb3_realtime();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "sais") == 0) main_sais(argc-1, argv+1);
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
