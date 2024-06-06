#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "io.h"
#include "rld0.h"
#include "ketopt.h"

#define RB3_BF_NO_FOR  0x1
#define RB3_BF_NO_REV  0x2
#define RB3_BF_SAIS64  0x4
#define RB3_BF_LINE    0x8
#define RB3_BF_USE_RB2 0x10

typedef struct {
	int64_t flag;
	rb3_fmt_t fmt;
	int32_t n_threads;
	int32_t block_len;
	int32_t max_nodes;
	int32_t sort_order;
	int64_t batch_size;
} rb3_bopt_t;

void rb3_bopt_init(rb3_bopt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->n_threads = 1;
	opt->fmt = RB3_PLAIN;
	opt->block_len = ROPE_DEF_BLOCK_LEN;
	opt->max_nodes = ROPE_DEF_MAX_NODES;
	opt->batch_size = 7000000000LL;
	opt->sort_order = MR_SO_IO;
}

static int usage_build(FILE *fp, const rb3_bopt_t *opt)
{
	fprintf(fp, "Usage: ropebwt3 build [options] <in.fa> [...]\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Algorithm:\n");
	fprintf(fp, "    -m NUM      batch size [7G]\n");
	fprintf(fp, "    -t INT      number of threads [%d]\n", opt->n_threads);
	fprintf(fp, "    -l INT      leaf block size in B+-tree [%d]\n", opt->block_len);
	fprintf(fp, "    -n INT      max number children per internal node [%d]\n", opt->max_nodes);
	fprintf(fp, "    -6          force to use 64-bit integers for libsais\n");
	fprintf(fp, "    -2          use the ropebwt2 algorithm (libsais by default)\n");
	fprintf(fp, "    -s          build BWT in the reverse lexicographical order (RLO; force -2)\n");
	fprintf(fp, "    -r          build BWT in RCLO (force -2)\n");
	fprintf(fp, "  Input:\n");
	fprintf(fp, "    -i FILE     read existing index from FILE []\n");
	fprintf(fp, "    -L          one sequence per line in the input\n");
	fprintf(fp, "    -F          no forward strand\n");
	fprintf(fp, "    -R          no reverse strand\n");
	fprintf(fp, "  Output:\n");
	fprintf(fp, "    -o FILE     output to FILE [stdout]\n");
	fprintf(fp, "    -d          dump in the fermi-delta format (FMD)\n");
	fprintf(fp, "    -b          dump in the ropebwt format (FMR)\n");
	fprintf(fp, "    -T          output the index in the Newick format (for debugging)\n");
	return fp == stdout? 0 : 1;
}

int main_build(int argc, char *argv[])
{
	rb3_bopt_t opt;
	kstring_t seq = {0,0,0};
	int32_t c, i;
	ketopt_t o = KETOPT_INIT;
	mrope_t *r = 0;
	char *fn_in = 0;

	rb3_bopt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "l:n:m:t:62sri:LFRo:dbT", 0)) >= 0) {
		// algorithm
		if (c == 'm') opt.batch_size = rb3_parse_num(o.arg);
		else if (c == 't') opt.n_threads = atoi(o.arg);
		else if (c == 'l') opt.block_len = atoi(o.arg);
		else if (c == 'n') opt.max_nodes = atoi(o.arg);
		else if (c == '6') opt.flag |= RB3_BF_SAIS64;
		else if (c == '2') opt.flag |= RB3_BF_USE_RB2;
		else if (c == 's') opt.flag |= RB3_BF_USE_RB2, opt.sort_order = MR_SO_RLO;
		else if (c == 'r') opt.flag |= RB3_BF_USE_RB2, opt.sort_order = MR_SO_RCLO;
		// input
		else if (c == 'i') fn_in = o.arg;
		else if (c == 'L') opt.flag |= RB3_BF_LINE;
		else if (c == 'F') opt.flag |= RB3_BF_NO_FOR;
		else if (c == 'R') opt.flag |= RB3_BF_NO_REV;
		// output
		else if (c == 'o') freopen(o.arg, "wb", stdout);
		else if (c == 'd') opt.fmt = RB3_FMD;
		else if (c == 'b') opt.fmt = RB3_FMR;
		else if (c == 'T') opt.fmt = RB3_TREE;
	}
	if (argc == o.ind && fn_in == 0) return usage_build(stdout, &opt);

	if (fn_in) {
		rb3_fmi_t fmi;
		rb3_fmi_restore(&fmi, fn_in);
		if (fmi.e == 0 && fmi.r == 0) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open index file '%s'\n", fn_in);
			return 1;
		} else if (fmi.is_fmd) {
			r = rb3_enc_fmd2fmr(fmi.e, opt.max_nodes, opt.block_len, 1);
		} else r = fmi.r;
	}

	for (i = o.ind; i < argc; ++i) {
		rb3_seqio_t *fp;
		int64_t n_seq = 0;
		fp = rb3_seq_open(argv[o.ind], !!(opt.flag&RB3_BF_LINE));
		if (fp == 0) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s'\n", argv[i]);
			continue;
		}
		while ((n_seq = rb3_seq_read(fp, &seq, opt.batch_size, !(opt.flag&RB3_BF_NO_FOR), !(opt.flag&RB3_BF_NO_REV))) > 0) {
			if (opt.flag & RB3_BF_USE_RB2) { // use the ropebwt2 algorithm
				if (r == 0) r = mr_init(opt.max_nodes, opt.block_len, opt.sort_order);
				rb3_reverse_all(seq.l, (uint8_t*)seq.s);
				mr_insert_multi(r, seq.l, (uint8_t*)seq.s, (opt.n_threads > 1));
			} else { // use libsais
				rb3_build_sais(n_seq, seq.l, seq.s, opt.n_threads, !!(opt.flag&RB3_BF_SAIS64));
				if (r == 0) {
					r = rb3_enc_plain2fmr(seq.l, (uint8_t*)seq.s, opt.max_nodes, opt.block_len);
				} else {
					mrope_t *p;
					rb3_fmi_t fmi;
					p = rb3_enc_plain2fmr(seq.l, (uint8_t*)seq.s, opt.max_nodes, opt.block_len);
					rb3_fmi_init(&fmi, 0, p);
					rb3_fmi_merge(r, &fmi, opt.n_threads, 1);
				}
			}
		}
		rb3_seq_close(fp);
	}
	free(seq.s);
	if (r == 0) return 1;

	if (opt.fmt == RB3_FMR) {
		mr_dump(r, stdout);
	} else if (opt.fmt == RB3_FMD) {
		rld_t *e;
		e = rb3_enc_fmr2fmd(r, 0, 1); // most of r is deallocated here
		rld_dump(e, "-");
		rld_destroy(e);
		r = 0;
	} else if (opt.fmt == RB3_PLAIN) {
		mr_print_bwt(r, stdout);
	} else if (opt.fmt == RB3_TREE) {
		mr_print_tree(r, stdout);
	}
	if (r) mr_destroy(r);
	return 0;
}
