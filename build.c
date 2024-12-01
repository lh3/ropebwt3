#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "io.h"
#include "rld0.h"
#include "rle.h"
#include "bre.h"
#include "ketopt.h"
#include "kthread.h"

#define RB3_BF_NO_FOR     0x1
#define RB3_BF_NO_REV     0x2
#define RB3_BF_LINE       0x4
#define RB3_BF_USE_RB2    0x8

typedef struct {
	int64_t flag;
	rb3_fmt_t fmt;
	int32_t n_threads;
	int32_t sais_threads;
	int32_t block_len;
	int32_t max_nodes;
	int32_t sort_order;
	int64_t batch_size;
} rb3_bopt_t;

void rb3_bopt_init(rb3_bopt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->n_threads = 4;
	opt->sais_threads = 0;
	opt->fmt = RB3_PLAIN;
	opt->block_len = ROPE_DEF_BLOCK_LEN;
	opt->max_nodes = ROPE_DEF_MAX_NODES;
	opt->batch_size = 7000000000LL;
	opt->sort_order = MR_SO_IO;
}

typedef struct {
	int64_t len;
	uint8_t *bwt;
} step_t;

typedef struct {
	const rb3_bopt_t *opt;
	int64_t id;
	rb3_seqio_t *fp;
	mrope_t *r;
} pipeline_t;

static void *worker_pipeline(void *shared, int step, void *in)
{
	pipeline_t *p = (pipeline_t*)shared;
	step_t *t = (step_t*)in;
	if (step == 0) {
		kstring_t seq = {0,0,0};
		int64_t n_seq;
		seq.m = 0x100000;
		seq.s = RB3_MALLOC(char, seq.m + 1);
		n_seq = rb3_seq_read(p->fp, &seq, p->opt->batch_size, !(p->opt->flag&RB3_BF_NO_FOR), !(p->opt->flag&RB3_BF_NO_REV));
		if (n_seq > 0) {
			int32_t n_threads = p->id == 0? p->opt->n_threads : p->opt->sais_threads;
			if (rb3_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] read %ld symbols\n", __func__, rb3_realtime(), rb3_percent_cpu(), (long)seq.l);
			t = RB3_CALLOC(step_t, 1);
			rb3_build_sais(n_seq, seq.l, seq.s, n_threads);
			if (rb3_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] constructed partial BWT for %ld symbols\n", __func__, rb3_realtime(), rb3_percent_cpu(), (long)seq.l);
			t->len = seq.l, t->bwt = (uint8_t*)seq.s;
			p->id++;
			return t;
		} else free(seq.s);
	} else if (step == 1) {
		int32_t n_threads = p->opt->n_threads - p->opt->sais_threads;
		if (p->r == 0) p->r = rb3_enc_plain2fmr(t->len, t->bwt, p->opt->max_nodes, p->opt->block_len, n_threads);
		else rb3_fmi_merge_plain(p->r, t->len, t->bwt, n_threads);
		if (rb3_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] encoded/merged the partial BWT for %ld symbols\n", __func__, rb3_realtime(), rb3_percent_cpu(), (long)t->len);
		free(t->bwt); free(t);
	}
	return 0;
}

static void mr_print_bre(mrope_t *r, const char *fn)
{
	mritr_t ri;
	const uint8_t *block;
	bre_file_t *f;
	bre_hdr_t h;

	bre_hdr_init(&h, BRE_AT_DNA6, 2);
	f = bre_open_write(fn, &h);
	mr_itr_first(r, &ri, 1);
	while ((block = mr_itr_next_block(&ri)) != 0) {
		const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
		while (q < end) {
			int c = 0;
			int64_t l;
			rle_dec1(q, c, l);
			bre_write(f, c, l);
		}
	}
	bre_close(f);
	mr_destroy(r);
}

static int usage_build(FILE *fp, const rb3_bopt_t *opt)
{
	fprintf(fp, "Usage: ropebwt3 build [options] <in.fa> [...]\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Algorithm:\n");
	fprintf(fp, "    -m NUM      batch size [7G]\n");
	fprintf(fp, "    -t INT      total number of threads [%d]\n", opt->n_threads);
	fprintf(fp, "    -p INT      #threads for sais and run sais and merge together (more RAM) [%d]\n", opt->sais_threads);
	fprintf(fp, "    -l INT      leaf block size in B+-tree [%d]\n", opt->block_len);
	fprintf(fp, "    -n INT      max number children per internal node [%d]\n", opt->max_nodes);
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
	fprintf(fp, "    -e          dump in the BRE format\n");
	fprintf(fp, "    -T          output the index in the Newick format (for debugging)\n");
	fprintf(fp, "    -S FILE     save the current index to FILE after each input file []\n");
	return fp == stdout? 0 : 1;
}

int main_build(int argc, char *argv[])
{
	rb3_bopt_t opt;
	kstring_t seq = {0,0,0};
	int32_t c, i;
	ketopt_t o = KETOPT_INIT;
	mrope_t *r = 0;
	char *fn_in = 0, *fn_tmp = 0;

	rb3_bopt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "l:n:m:t:2sri:LFRo:dbTS:p:e", 0)) >= 0) {
		// algorithm
		if (c == 'm') opt.batch_size = rb3_parse_num(o.arg);
		else if (c == 't') opt.n_threads = atoi(o.arg);
		else if (c == 'p') opt.sais_threads = atoi(o.arg);
		else if (c == 'l') opt.block_len = atoi(o.arg);
		else if (c == 'n') opt.max_nodes = atoi(o.arg);
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
		else if (c == 'e') opt.fmt = RB3_BRE;
		else if (c == 'S') fn_tmp = o.arg;
	}
	if (argc == o.ind && fn_in == 0)
		return usage_build(stderr, &opt);

	if (fn_in) {
		rb3_fmi_t fmi;
		rb3_fmi_restore(&fmi, fn_in, 0);
		if (fmi.e == 0 && fmi.r == 0) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open index file '%s'\n", fn_in);
			return 1;
		} else if (fmi.is_fmd) {
			r = rb3_enc_fmd2fmr(fmi.e, opt.max_nodes, opt.block_len, 1);
		} else r = fmi.r;
		if (rb3_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded the index from file '%s'\n", __func__, rb3_realtime(), rb3_percent_cpu(), fn_in);
	}

	if (argc - o.ind == 1 && opt.sais_threads > 0 && opt.n_threads - opt.sais_threads > 0) {
		rb3_seqio_t *fp;
		pipeline_t p;
		fp = rb3_seq_open(argv[o.ind], !!(opt.flag&RB3_BF_LINE));
		if (fp == 0) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s'\n", argv[o.ind]);
			goto end_build;
		}
		memset(&p, 0, sizeof(p));
		p.opt = &opt, p.fp = fp, p.r = r;
		kt_pipeline(2, worker_pipeline, &p, 2);
		r = p.r;
		rb3_seq_close(fp);
		goto end_build;
	}

	for (i = o.ind; i < argc; ++i) {
		rb3_seqio_t *fp;
		int64_t n_seq = 0;
		fp = rb3_seq_open(argv[i], !!(opt.flag&RB3_BF_LINE));
		if (fp == 0) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s'\n", argv[i]);
			continue;
		}
		while ((n_seq = rb3_seq_read(fp, &seq, opt.batch_size, !(opt.flag&RB3_BF_NO_FOR), !(opt.flag&RB3_BF_NO_REV))) > 0) {
			if (rb3_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] read %ld symbols from file '%s'\n", __func__, rb3_realtime(), rb3_percent_cpu(), (long)seq.l, argv[i]);
			if (opt.flag & RB3_BF_USE_RB2) { // use the ropebwt2 algorithm
				if (r == 0) r = mr_init(opt.max_nodes, opt.block_len, opt.sort_order);
				rb3_reverse_all(seq.l, (uint8_t*)seq.s);
				mr_insert_multi(r, seq.l, (uint8_t*)seq.s, (opt.n_threads > 1));
				if (rb3_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] inserted %ld symbols\n", __func__, rb3_realtime(), rb3_percent_cpu(), (long)seq.l);
			} else { // use libsais
				rb3_build_sais(n_seq, seq.l, seq.s, opt.n_threads);
				if (rb3_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] constructed partial BWT for %ld symbols\n", __func__, rb3_realtime(), rb3_percent_cpu(), (long)seq.l);
				if (r == 0) {
					r = rb3_enc_plain2fmr(seq.l, (uint8_t*)seq.s, opt.max_nodes, opt.block_len, opt.n_threads);
					if (rb3_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] encoded the partial BWT for %ld symbols\n", __func__, rb3_realtime(), rb3_percent_cpu(), (long)seq.l);
				} else {
					rb3_fmi_merge_plain(r, seq.l, (uint8_t*)seq.s, opt.n_threads);
					if (rb3_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] merged the partial BWT for %ld symbols\n", __func__, rb3_realtime(), rb3_percent_cpu(), (long)seq.l);
				}
			}
		}
		rb3_seq_close(fp);
		if (fn_tmp) {
			FILE *fp;
			fp = fopen(fn_tmp, "w");
			if (fp != 0) mr_dump(r, fp);
			fclose(fp);
			if (rb3_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] saved the current index to '%s'\n", __func__, rb3_realtime(), rb3_percent_cpu(), fn_tmp);
		}
	}
	free(seq.s);

end_build:
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
	} else if (opt.fmt == RB3_BRE) {
		mr_print_bre(r, "-");
		r = 0;
	}
	if (r) mr_destroy(r);
	return 0;
}
