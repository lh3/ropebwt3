#include <stdio.h>
#include <string.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "io.h"
#include "ketopt.h"

#define RB3_VERSION "3.0pre-r44"

int main_build(int argc, char *argv[]);
int main_merge(int argc, char *argv[]);
int main_get(int argc, char *argv[]);
int main_suffix(int argc, char *argv[]);
int main_match(int argc, char *argv[]);
int main_fa2line(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: ropebwt3 <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  build      construct a BWT\n");
	fprintf(fp, "  merge      merge BWTs\n");
	fprintf(fp, "  match      find supermaximal exact matches (requring both strands)\n");
	fprintf(fp, "  suffix     find the longest matching suffix (aka backward search)\n");
	fprintf(fp, "  get        retrieve the i-th sequence from BWT\n");
	fprintf(fp, "  fa2line    convert FASTX to lines\n");
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
	else if (strcmp(argv[1], "match") == 0) ret = main_match(argc-1, argv+1);
	else if (strcmp(argv[1], "suffix") == 0) ret = main_suffix(argc-1, argv+1);
	else if (strcmp(argv[1], "get") == 0) ret = main_get(argc-1, argv+1);
	else if (strcmp(argv[1], "fa2line") == 0) ret = main_fa2line(argc-1, argv+1);
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

int main_suffix(int argc, char *argv[])
{
	int32_t c, j, is_line = 0;
	ketopt_t o = KETOPT_INIT;
	rb3_fmi_t fmi;
	kstring_t out = {0,0,0};
	int64_t rec_num = 0;

	while ((c = ketopt(&o, argc, argv, 1, "L", 0)) >= 0) {
		if (c == 'L') is_line = 1;
	}
	if (argc - o.ind < 2) {
		fprintf(stdout, "Usage: ropebwt3 suffix [options] <idx.fmr> <seq.fa> [...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -L        one sequence per line in the input\n");
		return 0;
	}
	rb3_fmi_restore(&fmi, argv[o.ind]);
	if (fmi.e == 0 && fmi.r == 0) {
		if (rb3_verbose >= 1)
			fprintf(stderr, "ERROR: failed to load index file '%s'\n", argv[o.ind]);
		return 1;
	}
	for (j = o.ind + 1; j < argc; ++j) {
		const char *s, *name;
		int64_t i, len;
		rb3_seqio_t *fp;
		fp = rb3_seq_open(argv[j], is_line);
		while ((s = rb3_seq_read1(fp, &len, &name)) != 0) {
			int64_t k = 0, l = fmi.acc[RB3_ASIZE], last_size = 0;
			out.l = 0;
			for (i = len - 1; i >= 0; --i) {
				int c = s[i];
				int64_t size;
				c = c < 128 && c >= 0? rb3_nt6_table[c] : 5;
				size = rb3_fmi_extend1(&fmi, &k, &l, c);
				if (size == 0) break;
				last_size = size;
			}
			if (name) rb3_sprintf_lite(&out, "%s");
			else rb3_sprintf_lite(&out, "%ld", rec_num);
			rb3_sprintf_lite(&out, "\t%ld\t%ld\t%ld\n", i + 1, len, last_size);
			fputs(out.s, stdout);
		}
		rb3_seq_close(fp);
	}
	free(out.s);
	rb3_fmi_destroy(&fmi);
	return 0;
}

int main_match(int argc, char *argv[])
{
	int32_t c, j, is_line = 0;
	ketopt_t o = KETOPT_INIT;
	rb3_fmi_t fmi;
	kstring_t out = {0,0,0};
	int64_t rec_num = 0, min_occ = 1, min_len = 1;

	while ((c = ketopt(&o, argc, argv, 1, "Ll:c:", 0)) >= 0) {
		if (c == 'L') is_line = 1;
		else if (c == 'l') min_len = atol(o.arg);
		else if (c == 'c') min_occ = atol(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stdout, "Usage: ropebwt3 match [options] <idx.fmr> <seq.fa> [...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT    min SMEM length [%ld]\n", (long)min_len);
		fprintf(stderr, "  -s INT    min interval size [%ld]\n", (long)min_occ);
		fprintf(stderr, "  -L        one sequence per line in the input\n");
		return 0;
	}
	rb3_fmi_restore(&fmi, argv[o.ind]);
	if (fmi.e == 0 && fmi.r == 0) {
		if (rb3_verbose >= 1)
			fprintf(stderr, "ERROR: failed to load index file '%s'\n", argv[o.ind]);
		return 1;
	}
	if (fmi.acc[2] - fmi.acc[1] != fmi.acc[5] - fmi.acc[4] || fmi.acc[3] - fmi.acc[2] != fmi.acc[4] - fmi.acc[3]) {
		if (rb3_verbose >= 1)
			fprintf(stderr, "ERROR: #A != #T or #C != $G\n");
		return 1;
	}
	for (j = o.ind + 1; j < argc; ++j) {
		const char *name;
		char *seq;
		int64_t len;
		rb3_sai_v mem = {0,0,0};
		rb3_seqio_t *fp;
		fp = rb3_seq_open(argv[j], is_line);
		while ((seq = rb3_seq_read1(fp, &len, &name)) != 0) {
			int64_t i;
			rb3_char2nt6(len, (uint8_t*)seq);
			rb3_fmd_smem(0, &fmi, len, (uint8_t*)seq, &mem, min_occ, min_len);
			for (i = 0; i < mem.n; ++i) {
				rb3_sai_t *p = &mem.a[i];
				int32_t st = p->info>>32, en = (int32_t)p->info;
				out.l = 0;
				if (name) rb3_sprintf_lite(&out, "%s");
				else rb3_sprintf_lite(&out, "%ld", rec_num);
				rb3_sprintf_lite(&out, "\t%ld\t%ld\t%ld\n", st, en, (long)p->size);
				fputs(out.s, stdout);
			}
		}
		free(mem.a);
		rb3_seq_close(fp);
	}
	free(out.s);
	rb3_fmi_destroy(&fmi);
	return 0;
}

int main_fa2line(int argc, char *argv[])
{
	int32_t c, j, no_rev = 0;
	ketopt_t o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, "R", 0)) >= 0) {
		if (c == 'R') no_rev = 1;
	}
	if (argc - o.ind < 1) {
		fprintf(stdout, "Usage: ropebwt3 fa2line [options] <seq.fa> [...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -R        no reverse strand\n");
		return 0;
	}
	for (j = o.ind; j < argc; ++j) {
		char *seq;
		const char *name;
		int64_t len;
		rb3_seqio_t *fp;
		fp = rb3_seq_open(argv[j], 0);
		while ((seq = rb3_seq_read1(fp, &len, &name)) != 0) {
			int64_t i;
			rb3_char2nt6(len, (uint8_t*)seq);
			for (i = 0; i < len; ++i) seq[i] = "\nACGTX"[(uint8_t)seq[i]];
			puts(seq);
			if (!no_rev) {
				rb3_char2nt6(len, (uint8_t*)seq);
				rb3_revcomp6(len, (uint8_t*)seq);
				for (i = 0; i < len; ++i) seq[i] = "\nACGTX"[(uint8_t)seq[i]];
				puts(seq);
			}
		}
		rb3_seq_close(fp);
	}
	return 0;
}
