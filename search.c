#include "fm-index.h"
#include "align.h"
#include "rb3priv.h"
#include "io.h"
#include "ketopt.h"
#include "kthread.h"
#include "kalloc.h"

typedef enum { RB3_SA_MEM_TG, RB3_SA_MEM_ORI, RB3_SA_GREEDY, RB3_SA_SW } rb3_search_algo_t;

typedef struct {
	int32_t n_threads, no_kalloc, write_rs;
	rb3_search_algo_t algo;
	int64_t min_occ, min_len;
	int64_t batch_size;
	rb3_swopt_t swo;
} rb3_mopt_t;

void rb3_mopt_init(rb3_mopt_t *opt)
{
	memset(opt, 0, sizeof(rb3_mopt_t));
	opt->n_threads = 4;
	opt->min_occ = 1;
	opt->min_len = 19;
	opt->batch_size = 100000000;
	opt->algo = RB3_SA_MEM_TG;
	rb3_swopt_init(&opt->swo);
}

typedef struct mp_tbuf_s {
	void *km;
	rb3_sai_v mem; // this is allocated from km
} m_tbuf_t;

typedef struct {
	char *name;
	uint8_t *seq;
	int64_t id;
	int32_t len, n_mem;
	rb3_sai_t *mem;
	rb3_swrst_t rst;
} m_seq_t;

typedef struct {
	const rb3_mopt_t *opt;
	int64_t id;
	rb3_fmi_t fmi;
	rb3_seqio_t *fp;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
	int32_t n_seq;
	m_seq_t *seq;
	m_tbuf_t *buf;
} step_t;

static void worker_for(void *data, long i, int tid)
{
	step_t *t = (step_t*)data;
	const pipeline_t *p = t->p;
	m_seq_t *s = &t->seq[i];
	m_tbuf_t *b = &t->buf[tid];
	if (rb3_dbg_flag & RB3_DBG_QNAME)
		fprintf(stderr, "Q\t%s\t%d\n", s->name, tid);
	if (p->opt->algo == RB3_SA_SW) { // BWA-SW
		rb3_sw(b->km, &p->opt->swo, &p->fmi, s->len, s->seq, &s->rst);
	} else { // MEM algorithms
		rb3_char2nt6(s->len, s->seq);
		b->mem.n = 0;
		if (p->opt->algo == RB3_SA_MEM_TG)
			rb3_fmd_smem_TG(b->km, &p->fmi, s->len, s->seq, &b->mem, p->opt->min_occ, p->opt->min_len);
		else if (p->opt->algo == RB3_SA_GREEDY)
			rb3_fmd_gmem(b->km, &p->fmi, s->len, s->seq, &b->mem, p->opt->min_occ, p->opt->min_len);
		else if (p->opt->algo == RB3_SA_MEM_ORI)
			rb3_fmd_smem(b->km, &p->fmi, s->len, s->seq, &b->mem, p->opt->min_occ, p->opt->min_len);
		s->n_mem = b->mem.n;
		s->mem = RB3_MALLOC(rb3_sai_t, s->n_mem);
		memcpy(s->mem, b->mem.a, s->n_mem * sizeof(rb3_sai_t));
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	pipeline_t *p = (pipeline_t*)shared;
	step_t *t = (step_t*)in;
	int32_t i, j;
	if (step == 0) {
		const char *name;
		char *ss;
		int64_t len, tot = 0;
		int32_t n_seq = 0, m_seq = 0;
		m_seq_t *seq = 0;
		while ((ss = rb3_seq_read1(p->fp, &len, &name)) != 0) {
			m_seq_t *s;
			RB3_GROW0(m_seq_t, seq, n_seq, m_seq);
			s = &seq[n_seq++];
			s->name = name? rb3_strdup(name) : 0;
			s->seq = (uint8_t*)rb3_strdup(ss);
			s->len = len;
			s->id = p->id++;
			s->mem = 0, s->n_mem = 0;
			tot += len;
			if (tot >= p->opt->batch_size)
				break;
		}
		if (n_seq > 0) {
			t = RB3_CALLOC(step_t, 1);
			t->p = p;
			t->seq = seq;
			t->n_seq = n_seq;
			t->buf = RB3_CALLOC(m_tbuf_t, p->opt->n_threads);
			for (i = 0; i < p->opt->n_threads; ++i)
				t->buf[i].km = p->opt->no_kalloc? 0 : km_init();
			return t;
		}
	} else if (step == 1) {
		kt_for(p->opt->n_threads, worker_for, in, t->n_seq);
		return in;
	} else if (step == 2) {
		kstring_t out = {0,0,0};
		for (i = 0; i < p->opt->n_threads; ++i)
			km_destroy(t->buf[i].km);
		free(t->buf);
		for (j = 0; j < t->n_seq; ++j) {
			m_seq_t *s = &t->seq[j];
			free(s->seq);
			if (p->opt->algo == RB3_SA_SW) { // BWA-SW
				rb3_swrst_t *r = &s->rst;
				int32_t k;
				out.l = 0;
				if (r->score > 0 && r->n_qoff > 0) {
					if (s->name) rb3_sprintf_lite(&out, "%s", s->name);
					else rb3_sprintf_lite(&out, "seq%ld", s->id + 1);
					rb3_sprintf_lite(&out, "\t%d\t%d\t%d", s->len, r->qoff[0], r->qoff[0] + r->qlen);
					if (p->fmi.ssa) {
						int64_t pos, sid;
						pos = rb3_ssa(&p->fmi, p->fmi.ssa, r->lo, &sid);
						rb3_sprintf_lite(&out, "\t%c", "+-"[sid&1]);
						if (p->fmi.sid)
							rb3_sprintf_lite(&out, "\t%s\t%d", p->fmi.sid->name[sid>>1], p->fmi.sid->len[sid>>1]);
						else
							rb3_sprintf_lite(&out, "\t%ld\t%d", sid, r->rlen);
						rb3_sprintf_lite(&out, "\t%ld\t%ld", pos, pos + r->rlen);
					} else {
						rb3_sprintf_lite(&out, "\t*\t*\t%d\t*\t*", r->rlen);
					}
					rb3_sprintf_lite(&out, "\t%d\t%d\t0", r->mlen, r->blen);
					rb3_sprintf_lite(&out, "\tAS:i:%d\tqh:i:%d\trh:i:%ld\tcg:Z:", r->score, r->n_qoff, (long)(r->hi - r->lo));
					for (k = 0; k < r->n_cigar; ++k)
						rb3_sprintf_lite(&out, "%d%c", r->cigar[k]>>4, "MIDNSHP=X"[r->cigar[k]&0xf]);
					if (p->opt->write_rs) {
						rb3_sprintf_lite(&out, "\trs:Z:");
						for (k = 0; k < r->rlen; ++k)
							rb3_sprintf_lite(&out, "%c", "$ACGTN"[r->rseq[k]]);
					}
					rb3_sprintf_lite(&out, "\n");
					fputs(out.s, stdout);
				}
				rb3_swrst_free(r);
			} else {
				for (i = 0; i < s->n_mem; ++i) {
					rb3_sai_t *q = &s->mem[i];
					int32_t st = q->info>>32, en = (int32_t)q->info;
					out.l = 0;
					if (s->name) rb3_sprintf_lite(&out, "%s", s->name);
					else rb3_sprintf_lite(&out, "seq%ld", s->id + 1);
					rb3_sprintf_lite(&out, "\t%ld\t%ld\t%ld\n", st, en, (long)q->size);
					fputs(out.s, stdout);
				}
			}
			free(s->name);
			free(s->mem);
		}
		free(out.s);
		free(t->seq);
		if (rb3_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences\n", __func__, rb3_realtime(), rb3_percent_cpu(), t->n_seq);
		free(t);
	}
	return 0;
}

static ko_longopt_t long_options[] = {
	{ "no-ssa",          ko_no_argument,       301 },
	{ "seq",             ko_no_argument,       302 },
	{ "no-kalloc",       ko_no_argument,       501 },
	{ "dbg-dawg",        ko_no_argument,       502 },
	{ "dbg-sw",          ko_no_argument,       503 },
	{ "dbg-qname",       ko_no_argument,       504 },
	{ 0, 0, 0 }
};

int main_search(int argc, char *argv[]) // "sw" and "mem" share the same CLI
{
	int32_t c, j, is_line = 0, ret, load_flag = 0, no_ssa = 0;
	rb3_mopt_t opt;
	pipeline_t p;
	ketopt_t o = KETOPT_INIT;

	rb3_mopt_init(&opt);
	p.opt = &opt, p.id = 0;
	while ((c = ketopt(&o, argc, argv, 1, "Ll:c:t:K:MgdwN:A:B:O:E:C:m:k:", long_options)) >= 0) {
		if (c == 'L') is_line = 1;
		else if (c == 'g') opt.algo = RB3_SA_GREEDY;
		else if (c == 'w') opt.algo = RB3_SA_MEM_ORI;
		else if (c == 'd') opt.algo = RB3_SA_SW, load_flag |= RB3_LOAD_ALL;
		else if (c == 'l') opt.min_len = atol(o.arg);
		else if (c == 'c') opt.min_occ = atol(o.arg);
		else if (c == 't') opt.n_threads = atoi(o.arg);
		else if (c == 'K') opt.batch_size = rb3_parse_num(o.arg);
		else if (c == 'N') opt.swo.n_best = atoi(o.arg);
		else if (c == 'M') load_flag |= RB3_LOAD_MMAP;
		else if (c == 'A') opt.swo.match = atoi(o.arg);
		else if (c == 'B') opt.swo.mis = atoi(o.arg);
		else if (c == 'O') opt.swo.gap_open = atoi(o.arg);
		else if (c == 'E') opt.swo.gap_ext = atoi(o.arg);
		else if (c == 'C') opt.swo.r2cache_size = rb3_parse_num(o.arg);
		else if (c == 'm') opt.swo.min_sc = atoi(o.arg);
		else if (c == 'k') opt.swo.end_len = atoi(o.arg);
		else if (c == 301) no_ssa = 1;
		else if (c == 302) opt.write_rs = 1;
		else if (c == 501) opt.no_kalloc = 1;
		else if (c == 502) rb3_dbg_flag |= RB3_DBG_DAWG;
		else if (c == 503) rb3_dbg_flag |= RB3_DBG_SW;
		else if (c == 504) rb3_dbg_flag |= RB3_DBG_QNAME;
		else {
			fprintf(stderr, "ERROR: unknown option\n");
			return 1;
		}
	}
	if (strcmp(argv[0], "sw") == 0) {
		opt.algo = RB3_SA_SW;
		if (!no_ssa) load_flag |= RB3_LOAD_ALL;
	}
	if (argc - o.ind < 2) {
		fprintf(stdout, "Usage: ropebwt3 %s [options] <idx.fmr> <seq.fa> [...]\n", argv[0]);
		fprintf(stderr, "Options:\n");
		if (strcmp(argv[0], "mem") == 0 || strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -l INT      min MEM length [%ld]\n", (long)opt.min_len);
			fprintf(stderr, "  -s INT      min interval size [%ld]\n", (long)opt.min_occ);
			fprintf(stderr, "  -g          find greedy MEMs (faster but not always SMEMs)\n");
			fprintf(stderr, "  -w          use the original MEM algorithm (for testing)\n");
		}
		if (strcmp(argv[0], "search") == 0)
			fprintf(stderr, "  -d          use BWA-SW for local alignment\n");
		if (strcmp(argv[0], "sw") == 0 || strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -N INT      keep up to INT hits per DAWG node [%d]\n", opt.swo.n_best);
			fprintf(stderr, "  -k INT      initiate alignment with INT-mer matches [%d]\n", opt.swo.end_len);
			fprintf(stderr, "  -m INT      min alignment score [%d]\n", opt.swo.min_sc);
			fprintf(stderr, "  -A INT      match score [%d]\n", opt.swo.match);
			fprintf(stderr, "  -B INT      mismatch penalty [%d]\n", opt.swo.mis);
			fprintf(stderr, "  -O INT      gap open penalty [%d]\n", opt.swo.gap_open);
			fprintf(stderr, "  -E INT      gap extension penalty; a k-long gap costs O+k*E [%d]\n", opt.swo.gap_ext);
			fprintf(stderr, "  -C NUM      size of the ranking cache [%d]\n", opt.swo.r2cache_size);
			fprintf(stderr, "  --seq       write reference sequence to the rs tag\n");
			fprintf(stderr, "  --no-ssa    ignore the sampled suffix array\n");
		}
		fprintf(stderr, "  -t INT      number of threads [%d]\n", opt.n_threads);
		fprintf(stderr, "  -L          one sequence per line in the input\n");
		fprintf(stderr, "  -K NUM      query batch size [100m]\n");
		fprintf(stderr, "  -M          use mmap to load FMD\n");
		return 0;
	}
	ret = rb3_fmi_load_all(&p.fmi, argv[o.ind], load_flag);
	if (ret < 0) return 1;
	if (opt.algo != RB3_SA_SW) {
		if (!rb3_fmi_is_symmetric(&p.fmi)) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: BWT doesn't contain both strands, which is required for MEM\n");
			return 1;
		}
	}
	for (j = o.ind + 1; j < argc; ++j) {
		p.fp = rb3_seq_open(argv[j], is_line);
		kt_pipeline(2, worker_pipeline, &p, 3);
		rb3_seq_close(p.fp);
	}
	rb3_fmi_free(&p.fmi);
	return 0;
}
