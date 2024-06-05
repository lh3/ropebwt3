#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "io.h"
#include "rld0.h"
#include "libsais.h"
#include "libsais64.h"
#include "ketopt.h"

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: ropebwt3 build [options] <in.fa> [...]\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -t INT     number of threads\n");
	fprintf(fp, "  -F         no forward strand\n");
	fprintf(fp, "  -R         no reverse strand\n");
	fprintf(fp, "  -6         force to use 64-bit integers\n");
	fprintf(fp, "  -o FILE    output to FILE [stdout]\n");
	fprintf(fp, "  -d         dump in the fermi-delta format (FMD)\n");
	fprintf(fp, "  -b         dump in the ropebwt format (FMR)\n");
	return fp == stdout? 0 : 1;
}

int main_build(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t c, n_threads = 1, is_for = 1, is_rev = 1, is_line = 0, use64 = 0;
	rb3_fmt_t fmt = RB3_PLAIN;
	rb3_seqio_t *fp;
	kstring_t seq = {0,0,0};
	int64_t n_seq = 0;

	while ((c = ketopt(&o, argc, argv, 1, "t:FRL6o:db", 0)) >= 0) {
		if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'o') freopen(o.arg, "wb", stdout);
		else if (c == 'F') is_for = 0;
		else if (c == 'R') is_rev = 0;
		else if (c == 'L') is_line = 1;
		else if (c == 'd') fmt = RB3_FMD;
		else if (c == 'b') fmt = RB3_FMR;
		else if (c == '6') use64 = 1;
	}
	if (argc == o.ind) return usage(stdout);

	fp = rb3_seq_open(argv[o.ind], is_line);
	n_seq = rb3_seq_read(fp, &seq, 0, is_for, is_rev);
	rb3_build_sais(n_seq, seq.l, seq.s, n_threads, use64);
	if (fmt == RB3_PLAIN) {
		int64_t i;
		for (i = 0; i < seq.l; ++i)
			seq.s[i] = "$ACGTN"[(uint8_t)seq.s[i]];
		puts(seq.s);
	} else if (fmt == RB3_FMD) {
		rld_t *e;
		e = rb3_enc_plain2rld(seq.l, (uint8_t*)seq.s);
		rld_dump(e, "-");
		rld_destroy(e);
	} else if (fmt == RB3_FMR) {
		mrope_t *r;
		r = rb3_enc_plain2fmr(seq.l, (uint8_t*)seq.s, 0, 0);
		mr_dump(r, stdout);
		mr_destroy(r);
	}
	free(seq.s);
	rb3_seq_close(fp);
	return 0;
}
