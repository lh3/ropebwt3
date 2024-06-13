#include "fm-index.h"
#include "rb3priv.h"
#include "io.h"
#include "ketopt.h"

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
			++rec_num;
			rb3_char2nt6(len, (uint8_t*)seq);
			rb3_fmd_smem(0, &fmi, len, (uint8_t*)seq, &mem, min_occ, min_len);
			for (i = 0; i < mem.n; ++i) {
				rb3_sai_t *p = &mem.a[i];
				int32_t st = p->info>>32, en = (int32_t)p->info;
				out.l = 0;
				if (name) rb3_sprintf_lite(&out, "%s", name);
				else rb3_sprintf_lite(&out, "seq%ld", rec_num);
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

