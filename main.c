#include <stdio.h>
#include <string.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "io.h"
#include "ketopt.h"

#define RB3_VERSION "3.9-r280-dirty"

int main_build(int argc, char *argv[]);
int main_merge(int argc, char *argv[]);
int main_get(int argc, char *argv[]);
int main_ssa(int argc, char *argv[]);
int main_suffix(int argc, char *argv[]);
int main_search(int argc, char *argv[]);
int main_kount(int argc, char *argv[]);
int main_fa2line(int argc, char *argv[]);
int main_fa2kmer(int argc, char *argv[]);
int main_plain2fmd(int argc, char *argv[]);
int main_stat(int argc, char *argv[]);

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: ropebwt3 <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  Search:\n");
	fprintf(fp, "    sw         find local alignment\n");
	fprintf(fp, "    mem        find maximal exact matches\n");
	fprintf(fp, "    hapdiv     haplotype diversity with sliding k-mers\n");
	fprintf(fp, "    suffix     find the longest matching suffix\n");
	fprintf(fp, "  Construction:\n");
	fprintf(fp, "    build      construct a BWT\n");
	fprintf(fp, "    merge      merge BWTs\n");
	fprintf(fp, "    plain2fmd  convert BWT in plain text to FMD\n");
	fprintf(fp, "    ssa        generate sampled suffix array\n");
	fprintf(fp, "  Miscellaneous:\n");
	fprintf(fp, "    get        retrieve the i-th sequence from BWT\n");
	fprintf(fp, "    stat       basic statistics of BWT\n");
	fprintf(fp, "    kount      count (high-occurrence) k-mers\n");
	fprintf(fp, "    fa2line    convert FASTX to lines\n");
	fprintf(fp, "    fa2kmer    extract k-mers from FASTX\n");
	fprintf(fp, "    version    print the version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	int ret = 0;
	rb3_init();
	if (argc == 1) return usage(stdout);
	else if (strcmp(argv[1], "search") == 0) ret = main_search(argc-1, argv+1);
	else if (strcmp(argv[1], "sw") == 0) ret = main_search(argc-1, argv+1);
	else if (strcmp(argv[1], "mem") == 0) ret = main_search(argc-1, argv+1);
	else if (strcmp(argv[1], "hapdiv") == 0) ret = main_search(argc-1, argv+1);
	else if (strcmp(argv[1], "build") == 0) ret = main_build(argc-1, argv+1);
	else if (strcmp(argv[1], "merge") == 0) ret = main_merge(argc-1, argv+1);
	else if (strcmp(argv[1], "ssa") == 0) ret = main_ssa(argc-1, argv+1);
	else if (strcmp(argv[1], "stat") == 0) ret = main_stat(argc-1, argv+1);
	else if (strcmp(argv[1], "suffix") == 0) ret = main_suffix(argc-1, argv+1);
	else if (strcmp(argv[1], "get") == 0) ret = main_get(argc-1, argv+1);
	else if (strcmp(argv[1], "kount") == 0) ret = main_kount(argc-1, argv+1);
	else if (strcmp(argv[1], "fa2line") == 0) ret = main_fa2line(argc-1, argv+1);
	else if (strcmp(argv[1], "fa2kmer") == 0) ret = main_fa2kmer(argc-1, argv+1);
	else if (strcmp(argv[1], "plain2fmd") == 0) ret = main_plain2fmd(argc-1, argv+1);
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
	rb3_fmi_t fmi;
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

	rb3_fmi_restore(&fmi, argv[o.ind], 0);
	r = fmi.is_fmd? rb3_enc_fmd2fmr(fmi.e, 0, 0, 1) : fmi.r;
	if (r == 0) {
		if (rb3_verbose >= 1)
			fprintf(stderr, "ERROR: failed to load FMR file '%s'\n", argv[o.ind]);
		return 1;
	}
	for (i = o.ind + 1; i < argc; ++i) {
		rb3_fmi_t fb;
		rb3_fmi_restore(&fb, argv[i], 0);
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
	rb3_fmi_restore(&fmi, argv[o.ind], 0);
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
	rb3_fmi_free(&fmi);
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
	rb3_fmi_restore(&fmi, argv[o.ind], 0);
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
			++rec_num;
			out.l = 0;
			for (i = len - 1; i >= 0; --i) {
				int c = s[i];
				int64_t size;
				c = c < 128 && c >= 0? rb3_nt6_table[c] : 5;
				size = rb3_fmi_extend1(&fmi, &k, &l, c);
				if (size == 0) break;
				last_size = size;
			}
			if (name) rb3_sprintf_lite(&out, "%s", name);
			else rb3_sprintf_lite(&out, "seq%ld", rec_num);
			rb3_sprintf_lite(&out, "\t%ld\t%ld\t%ld\n", i + 1, len, last_size);
			fputs(out.s, stdout);
		}
		rb3_seq_close(fp);
	}
	free(out.s);
	rb3_fmi_free(&fmi);
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

int main_fa2kmer(int argc, char *argv[])
{
	int32_t c, j, kmer = 151, step = 50;
	ketopt_t o = KETOPT_INIT;
	kstring_t out = {0,0,0};
	while ((c = ketopt(&o, argc, argv, 1, "k:w:", 0)) >= 0) {
		if (c == 'k') kmer = atoi(o.arg);
		else if (c == 'w') step = atoi(o.arg);
		else {
			fprintf(stderr, "ERROR: unknown option\n");
			return 1;
		}
	}
	if (argc - o.ind < 1) {
		fprintf(stdout, "Usage: ropebwt3 fa2kmer [options] <seq.fa> [...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT      k-mer size [%d]\n", kmer);
		fprintf(stderr, "  -w INT      step size [%d]\n", step);
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
			for (i = 0; i < len; i += step) {
				int64_t en = i + step + kmer > len? len : i + kmer; // special treatment for the last k-mer
				out.l = 0;
				rb3_sprintf_lite(&out, ">%s:%ld-%ld\n", name, (long)(i + 1), (long)en);
				fwrite(out.s, 1, out.l, stdout);
				fwrite(&seq[i], 1, en - i, stdout);
				fputc('\n', stdout);
				if (en == len) break;
			}
		}
		rb3_seq_close(fp);
	}
	free(out.s);
	return 0;
}

int main_plain2fmd(int argc, char *argv[])
{
	int32_t c, j;
	uint8_t buf[0x10000];
	ketopt_t o = KETOPT_INIT;
	rld_t *e;
	rlditr_t ei;
	while ((c = ketopt(&o, argc, argv, 1, "o:", 0)) >= 0) {
		if (c == 'o') freopen(o.arg, "wb", stdout);
	}
	if (argc - o.ind < 1) {
		fprintf(stdout, "Usage: ropebwt3 plain2fmd [-o output.fmd] <in.txt>\n");
		return 0;
	}
	e = rld_init(RB3_ASIZE, 3);
	rld_itr_init(e, &ei, 0);
	for (j = o.ind; j < argc; ++j) {
		FILE *fp;
		int32_t i, len;
		fp = strcmp(argv[j], "-") == 0? stdin : fopen(argv[j], "r");
		if (fp == 0) continue; // TODO: change this to an error
		while ((len = fread(buf, 1, 0x10000, fp)) > 0) {
			for (i = 0; i < len; ++i) {
				int c = buf[i] == '\n' || buf[i] == '$'? 0 : buf[i] >= 128? 5 : rb3_nt6_table[buf[i]];
				rld_enc(e, &ei, 1, c);
			}
		}
		fclose(fp);
	}
	rld_enc_finish(e, &ei);
	rld_dump(e, "-");
	return 0;
}

typedef struct {
	int64_t k, l;
	int d, c;
} count_pair64_t;

typedef struct {
	rb3_fmi_t fmi;
	uint64_t s_top;
	count_pair64_t *stack;
	int64_t ok[6], ol[6];
	count_pair64_t top;
} count_aux_t;

int main_kount(int argc, char *argv[])
{
	int c, min_occ = 100, depth = 51, i, n = 0;
	count_pair64_t *p;
	count_aux_t *aux;
	char *str;
	kstring_t out = {0,0,0};
	ketopt_t o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, "k:m:", 0)) >= 0) {
		if (c == 'k') depth = atol(o.arg);
		else if (c == 'm') min_occ = atol(o.arg);
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: ropebwt3 kount [options] <in1.fmd> [in2.fmd [...]]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT       k-mer length [%d]\n", depth);
		fprintf(stderr, "  -m INT       min k-mer occurrence [%d]\n", min_occ);
		return 1;
	}

	n = argc - o.ind;
	aux = RB3_CALLOC(count_aux_t, n);
	for (i = 0; i < n; ++i) {
		count_aux_t *q = &aux[i];
		rb3_fmi_restore(&q->fmi, argv[o.ind + i], 0);
		if (q->fmi.e == 0 && q->fmi.r == 0) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: failed to load index '%s'\n", argv[o.ind + i]);
			return 1; // FIXME: potential memory leak
		}
		q->stack = RB3_CALLOC(count_pair64_t, (depth + 1) * 4);
		p = &q->stack[q->s_top++];
		p->k = 0, p->l = q->fmi.acc[RB3_ASIZE], p->d = p->c = 0;
	}
	str = RB3_CALLOC(char, depth + 1);
	str[depth] = 0;
	while (1) {
		int a;
		for (i = 0; i < n; ++i) {
			if (aux[i].s_top == 0) break;
			aux[i].top = aux[i].stack[--aux[i].s_top];
		}
		if (i < n) break;
		if (aux->top.d > 0) str[depth - aux->top.d] = "$ACGTN"[aux->top.c];
		for (i = 0; i < n; ++i)
			rb3_fmi_rank2a(&aux[i].fmi, aux[i].top.k, aux[i].top.l, aux[i].ok, aux[i].ol);
		for (a = 1; a <= 4; ++a) {
			for (i = 0; i < n; ++i)
				if (aux[i].ol[a] - aux[i].ok[a] >= min_occ) break;
			if (i == n) continue;
			str[depth - aux->top.d - 1] = "$ACGTN"[a];
			if (aux->top.d != depth - 1) {
				for (i = 0; i < n; ++i) {
					count_aux_t *q = &aux[i];
					p = &q->stack[q->s_top++];
					p->k = q->fmi.acc[a] + q->ok[a];
					p->l = q->fmi.acc[a] + q->ol[a];
					p->d = aux->top.d + 1;
					p->c = a;
				}
			} else {
				out.l = 0;
				rb3_sprintf_lite(&out, "%s", str);
				for (i = 0; i < n; ++i)
					rb3_sprintf_lite(&out, "\t%ld", (long)(aux[i].ol[a] - aux[i].ok[a]));
				puts(out.s);
			}
		}
	}
	free(str);
	for (i = 0; i < n; ++i) {
		free(aux[i].stack);
		rb3_fmi_free(&aux[i].fmi);
	}
	free(aux);
	return 0;
}

int main_stat(int argc, char *argv[])
{
	int32_t c, use_mmap = 0;
	ketopt_t o = KETOPT_INIT;
	rb3_fmi_t fmi;

	while ((c = ketopt(&o, argc, argv, 1, "M", 0)) >= 0) { }
	if (argc - o.ind == 0) {
		fprintf(stdout, "Usage: ropebwt3 stat [-M] <idx.fmd>\n");
		return 0;
	}
	rb3_fmi_restore(&fmi, argv[o.ind], use_mmap);
	if (fmi.e == 0 && fmi.r == 0) {
		if (rb3_verbose >= 1)
			fprintf(stderr, "ERROR: failed to load index file '%s'\n", argv[o.ind]);
		return 1;
	}
	printf("%ld sequences\n", (long)fmi.acc[1]);
	printf("%ld symbols\n", (long)fmi.acc[6]);
	printf("%ld runs\n", (long)rb3_fmi_get_r(&fmi));
	printf("%ld A\n", (long)(fmi.acc[2] - fmi.acc[1]));
	printf("%ld C\n", (long)(fmi.acc[3] - fmi.acc[2]));
	printf("%ld G\n", (long)(fmi.acc[4] - fmi.acc[3]));
	printf("%ld T\n", (long)(fmi.acc[5] - fmi.acc[4]));
	printf("%ld N\n", (long)(fmi.acc[6] - fmi.acc[5]));
	rb3_fmi_free(&fmi);
	return 0;
}
