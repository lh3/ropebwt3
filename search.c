#include "fm-index.h"
#include "align.h"
#include "rb3priv.h"
#include "io.h"
#include "ketopt.h"
#include "kthread.h"
#include "kalloc.h"

typedef enum { RB3_SA_MEM_TG, RB3_SA_MEM_ORI, RB3_SA_SW, RB3_SA_HAPDIV } rb3_search_algo_t;

#define RB3_MF_NO_KALLOC   0x1
#define RB3_MF_WRITE_RS    0x2
#define RB3_MF_WRITE_UNMAP 0x4
#define RB3_MF_WRITE_COV   0x8

typedef struct {
	uint32_t flag;
	int32_t n_threads, min_gap_len, hapdiv_k, hapdiv_w;
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
	opt->hapdiv_k = 101;
	opt->hapdiv_w = 50;
	opt->batch_size = 100000000;
	opt->algo = RB3_SA_MEM_TG;
	rb3_swopt_init(&opt->swo);
}

typedef struct mp_tbuf_s {
	void *km;
	int32_t n_gap, m_gap;
	uint64_t *gap;
	rb3_sai_v mem; // this is allocated from km
} m_tbuf_t;

typedef struct {
	char *name;
	uint8_t *seq;
	int64_t id;
	int32_t len, n_mem, n_gap;
	uint64_t *gap;
	rb3_sai_t *mem;
} m_seq_t;

typedef struct {
	const rb3_mopt_t *opt;
	int64_t id;
	rb3_fmi_t fmi;
	rb3_seqio_t *fp;
} pipeline_t;

typedef struct {
	int32_t id, offset;
	rb3_hapdiv_t r;
} m_hapdiv_t;

typedef struct {
	const pipeline_t *p;
	int32_t n_seq, n_hapdiv;
	m_seq_t *seq;
	rb3_swrst_t *rst;
	m_hapdiv_t *hapdiv;
	m_tbuf_t *buf;
} step_t;

static void worker_for_seq(void *data, long i, int tid)
{
	step_t *t = (step_t*)data;
	const pipeline_t *p = t->p;
	m_seq_t *s = &t->seq[i];
	m_tbuf_t *b = &t->buf[tid];
	if (rb3_dbg_flag & RB3_DBG_QNAME)
		fprintf(stderr, "Q\t%s\t%d\n", s->name, tid);
	rb3_char2nt6(s->len, s->seq);
	if (p->opt->algo == RB3_SA_SW) { // BWA-SW
		rb3_sw(b->km, &p->opt->swo, &p->fmi, s->len, s->seq, &t->rst[i]);
	} else { // MEM algorithms
		b->mem.n = 0;
		if (p->opt->algo == RB3_SA_MEM_TG)
			rb3_fmd_smem_TG(b->km, &p->fmi, s->len, s->seq, &b->mem, p->opt->min_occ, p->opt->min_len);
		else if (p->opt->algo == RB3_SA_MEM_ORI)
			rb3_fmd_smem(b->km, &p->fmi, s->len, s->seq, &b->mem, p->opt->min_occ, p->opt->min_len);
		s->n_mem = b->mem.n;
		s->mem = RB3_MALLOC(rb3_sai_t, s->n_mem);
		memcpy(s->mem, b->mem.a, s->n_mem * sizeof(rb3_sai_t));
		if (p->opt->min_gap_len > 0) {
			int32_t i, last = 0;
			b->n_gap = 0;
			Kgrow(b->km, uint64_t, b->gap, b->mem.n + 1, b->m_gap);
			for (i = 0; i < b->mem.n; ++i) {
				int32_t st = b->mem.a[i].info>>32, en = (int32_t)b->mem.a[i].info;
				if (st > last) {
					if (st - last >= p->opt->min_gap_len)
						b->gap[b->n_gap++] = (uint64_t)last<<32 | st;
					last = en;
				} else last = last > en? last : en;
			}
			if (s->len - last >= p->opt->min_gap_len)
				b->gap[b->n_gap++] = (uint64_t)last<<32 | s->len;
			s->n_gap = b->n_gap;
			s->gap = RB3_MALLOC(uint64_t, s->n_gap);
			memcpy(s->gap, b->gap, s->n_gap * 8);
		}
	}
}

static void worker_for_hapdiv(void *data, long i, int tid)
{
	step_t *t = (step_t*)data;
	const pipeline_t *p = t->p;
	m_hapdiv_t *a = &t->hapdiv[i];
	rb3_hapdiv(t->buf[tid].km, &p->opt->swo, &p->fmi, p->opt->hapdiv_k, &t->seq[a->id].seq[a->offset], &a->r);
}

static inline void write_name(kstring_t *out, const m_seq_t *s)
{
	if (s->name) rb3_sprintf_lite(out, "%s", s->name);
	else rb3_sprintf_lite(out, "seq%ld", s->id + 1);
}

static void write_paf(kstring_t *out, const rb3_fmi_t *f, const rb3_swhit_t *h, const m_seq_t *s, int32_t write_rs)
{
	int32_t k;
	write_name(out, s);
	rb3_sprintf_lite(out, "\t%d\t%d\t%d", s->len, h->qoff[0], h->qoff[0] + h->qlen);
	if (h->lo_pos >= 0 && h->lo_sid >= 0) {
		if (f->sid) { // sequence names and lengths are available
			int64_t rlen = f->sid->len[h->lo_sid>>1];
			rb3_sprintf_lite(out, "\t%c\t%s\t%ld", "+-"[h->lo_sid&1], f->sid->name[h->lo_sid>>1], (long)rlen);
			if ((h->lo_sid&1) == 0) // forward strand
				rb3_sprintf_lite(out, "\t%ld\t%ld", h->lo_pos, h->lo_pos + h->rlen);
			else // reverse strand
				rb3_sprintf_lite(out, "\t%ld\t%ld", rlen - (h->lo_pos + h->rlen), rlen - h->lo_pos);
		} else { // sequence names and lengths are not available
			rb3_sprintf_lite(out, "\t+\t%ld\t*\t%ld\t%ld", h->lo_sid, h->lo_pos, h->lo_pos + h->rlen); // always on the forward strand
		}
	} else {
		rb3_sprintf_lite(out, "\t*\t*\t%d\t*\t*", h->rlen);
	}
	rb3_sprintf_lite(out, "\t%d\t%d\t0", h->mlen, h->blen);
	rb3_sprintf_lite(out, "\tAS:i:%d\tqh:i:%d\trh:i:%ld\tcg:Z:", h->score, h->n_qoff, (long)(h->hi - h->lo));
	for (k = 0; k < h->n_cigar; ++k)
		rb3_sprintf_lite(out, "%d%c", h->cigar[k]>>4, "MIDNSHP=X"[h->cigar[k]&0xf]);
	if (write_rs) {
		rb3_sprintf_lite(out, "\trs:Z:");
		for (k = 0; k < h->rlen; ++k)
			rb3_sprintf_lite(out, "%c", "$ACGTN"[h->rseq[k]]);
	}
	rb3_sprintf_lite(out, "\n");
}

static void write_per_seq(step_t *t)
{
	const pipeline_t *p = t->p;
	int32_t i, j;
	kstring_t out = {0,0,0};
	for (j = 0; j < t->n_seq; ++j) {
		m_seq_t *s = &t->seq[j];
		free(s->seq);
		out.l = 0;
		if (p->opt->algo == RB3_SA_SW) { // BWA-SW
			rb3_swrst_t *r = &t->rst[j];
			if (r->n > 0) { // mapped
				for (i = 0; i < r->n; ++i) {
					out.l = 0;
					write_paf(&out, &p->fmi, &r->a[i], s, p->opt->flag & RB3_MF_WRITE_RS);
					fputs(out.s, stdout);
				}
			} else if (p->opt->flag & RB3_MF_WRITE_UNMAP) { // unmapped
				write_name(&out, s);
				rb3_sprintf_lite(&out, "\t%d\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\n", s->len);
				fputs(out.s, stdout);
			}
			rb3_swrst_free(r);
		} else if (p->opt->min_gap_len > 0) { // output regions not covered by long MEMs
			for (i = 0; i < s->n_gap; ++i) {
				int32_t st = s->gap[i]>>32, en = (int32_t)s->gap[i];
				out.l = 0;
				write_name(&out, s);
				rb3_sprintf_lite(&out, "\t%d\t%d\t%d\n", st, en, s->len);
				fputs(out.s, stdout);
			}
		} else if (p->opt->flag & RB3_MF_WRITE_COV) { // output breadth of coverage
			int32_t st0 = 0, en0 = 0, cov = 0;
			for (i = 0; i < s->n_mem; ++i) {
				rb3_sai_t *q = &s->mem[i];
				int32_t st = q->info>>32, en = (int32_t)q->info;
				if (st > en0) {
					cov += en0 - st0;
					st0 = st, en0 = en;
				} else en0 = en0 > en? en0 : en;
			}
			cov += en0 - st0;
			if (cov > 0) {
				out.l = 0;
				write_name(&out, s);
				rb3_sprintf_lite(&out, "\t%d\t%d\n", s->len, cov);
				fputs(out.s, stdout);
			}
		} else { // output long MEMs
			for (i = 0; i < s->n_mem; ++i) {
				rb3_sai_t *q = &s->mem[i];
				int32_t st = q->info>>32, en = (int32_t)q->info;
				out.l = 0;
				write_name(&out, s);
				rb3_sprintf_lite(&out, "\t%d\t%d\t%ld\n", st, en, (long)q->size);
				fputs(out.s, stdout);
			}
		}
		free(s->name); free(s->mem); free(s->gap);
	}
	free(out.s);
	free(t->rst);
}

static void write_hapdiv(step_t *t)
{
	int32_t j, ed;
	const m_hapdiv_t *p;
	kstring_t out = {0,0,0};
	for (j = 0; j < t->n_seq; ++j)
		free(t->seq[j].seq);
	if (t->n_hapdiv == 0) return;
	p = &t->hapdiv[0];
	for (j = 1; j <= t->n_hapdiv; ++j) {
		const m_hapdiv_t *q = t->hapdiv + j;
		if (j == t->n_hapdiv || p->id != q->id || memcmp(&p->r, &q->r, sizeof(p->r)) != 0) {
			m_seq_t *s = &t->seq[p->id];
			out.l = 0;
			write_name(&out, s);
			rb3_sprintf_lite(&out, "\t%d\t%d\t%d\t%d", p->offset, t->hapdiv[j-1].offset + t->p->opt->hapdiv_k, p->r.n_al, p->r.max_ed);
			for (ed = 0; ed <= RB2_SW_MAX_ED; ++ed)
				rb3_sprintf_lite(&out, "\t%d", p->r.n_hap[ed]);
			puts(out.s);
			p = q;
		}
	}
	for (j = 0; j < t->n_seq; ++j)
		free(t->seq[j].name);
	free(out.s);
	free(t->hapdiv);
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	pipeline_t *p = (pipeline_t*)shared;
	step_t *t = (step_t*)in;
	int32_t i;
	if (step == 0) {
		const char *name;
		char *ss;
		int64_t len, tot = 0;
		int32_t n_seq = 0, m_seq = 0;
		m_seq_t *seq = 0;
		while ((ss = rb3_seq_read1(p->fp, &len, &name)) != 0) { // read sequences
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
		if (n_seq > 0) { // construct a step_t object
			t = RB3_CALLOC(step_t, 1);
			t->p = p;
			t->seq = seq;
			t->n_seq = n_seq;
			if (p->opt->algo == RB3_SA_HAPDIV) { // the hapdiv mode
				int32_t j, n_hapdiv = 0;
				for (i = 0; i < n_seq; ++i)
					n_hapdiv += seq[i].len < p->opt->hapdiv_k? 0 : (seq[i].len - p->opt->hapdiv_k) / p->opt->hapdiv_w + 1;
				t->n_hapdiv = n_hapdiv;
				t->hapdiv = RB3_CALLOC(m_hapdiv_t, n_hapdiv);
				for (i = 0, n_hapdiv = 0; i < n_seq; ++i)
					for (j = 0; j + p->opt->hapdiv_k <= seq[i].len; j += p->opt->hapdiv_w)
						t->hapdiv[n_hapdiv].id = i, t->hapdiv[n_hapdiv++].offset = j;
				assert(n_hapdiv == t->n_hapdiv);
			} else { // per-sequence mode (sw, mem, gap and coverage)
				t->rst = RB3_CALLOC(rb3_swrst_t, n_seq);
			}
			t->buf = RB3_CALLOC(m_tbuf_t, p->opt->n_threads);
			for (i = 0; i < p->opt->n_threads; ++i)
				t->buf[i].km = p->opt->flag & RB3_MF_NO_KALLOC? 0 : km_init();
			return t;
		}
	} else if (step == 1) {
		if (p->opt->algo == RB3_SA_HAPDIV)
			kt_for(p->opt->n_threads, worker_for_hapdiv, in, t->n_hapdiv);
		else
			kt_for(p->opt->n_threads, worker_for_seq, in, t->n_seq);
		return in;
	} else if (step == 2) {
		for (i = 0; i < p->opt->n_threads; ++i)
			km_destroy(t->buf[i].km);
		free(t->buf);
		if (p->opt->algo == RB3_SA_HAPDIV)
			write_hapdiv(t);
		else
			write_per_seq(t);
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
	{ "gap",             ko_required_argument, 303 },
	{ "cov",             ko_no_argument,       304 },
	{ "old-mem",         ko_no_argument,       305 },
	{ "no-kalloc",       ko_no_argument,       501 },
	{ "dbg-dawg",        ko_no_argument,       502 },
	{ "dbg-sw",          ko_no_argument,       503 },
	{ "dbg-qname",       ko_no_argument,       504 },
	{ "dbg-bt",          ko_no_argument,       505 },
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
	while ((c = ketopt(&o, argc, argv, 1, "Ll:c:t:K:MdN:A:B:O:E:C:m:k:uj:ey:a:w:", long_options)) >= 0) {
		if (c == 'L') is_line = 1;
		else if (c == 'a') opt.algo = RB3_SA_HAPDIV, opt.hapdiv_k = atoi(o.arg);
		else if (c == 'w') opt.algo = RB3_SA_HAPDIV, opt.hapdiv_w = atoi(o.arg);
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
		else if (c == 'j') opt.swo.min_mem_len = atoi(o.arg);
		else if (c == 'e') opt.swo.flag |= RB3_SWF_E2E, opt.swo.end_len = 1;
		else if (c == 'y') opt.swo.e2e_drop = atoi(o.arg);
		else if (c == 'u') opt.flag |= RB3_MF_WRITE_UNMAP;
		else if (c == 301) no_ssa = 1;
		else if (c == 302) opt.flag |= RB3_MF_WRITE_RS;
		else if (c == 303) opt.min_gap_len = rb3_parse_num(o.arg);
		else if (c == 304) opt.flag |= RB3_MF_WRITE_COV;
		else if (c == 305) opt.algo = RB3_SA_MEM_ORI;
		else if (c == 501) opt.flag |= RB3_MF_NO_KALLOC;
		else if (c == 502) rb3_dbg_flag |= RB3_DBG_DAWG;
		else if (c == 503) rb3_dbg_flag |= RB3_DBG_SW;
		else if (c == 504) rb3_dbg_flag |= RB3_DBG_QNAME;
		else if (c == 505) rb3_dbg_flag |= RB3_DBG_BT;
		else {
			fprintf(stderr, "ERROR: unknown option\n");
			return 1;
		}
	}
	if (strcmp(argv[0], "sw") == 0) {
		opt.algo = RB3_SA_SW;
		if (!no_ssa) load_flag |= RB3_LOAD_ALL;
	} else if (strcmp(argv[0], "hapdiv") == 0) {
		opt.algo = RB3_SA_HAPDIV, opt.swo.end_len = 1;
	}
	if (opt.algo == RB3_SA_HAPDIV)
		opt.swo.flag |= RB3_SWF_E2E | RB3_SWF_HAPDIV;
	if (argc - o.ind < 2) {
		fprintf(stdout, "Usage: ropebwt3 %s [options] <idx.fmr> <seq.fa> [...]\n", argv[0]);
		fprintf(stderr, "Options:\n");
		if (strcmp(argv[0], "mem") == 0 || strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -l INT      min MEM length [%ld]\n", (long)opt.min_len);
			fprintf(stderr, "  -c INT      min interval size [%ld]\n", (long)opt.min_occ);
			fprintf(stderr, "  --old-mem   use the original MEM algorithm (for testing)\n");
			fprintf(stderr, "  --gap=NUM   output regions >=NUM that are not covered by MEMs [%d]\n", opt.min_gap_len);
			fprintf(stderr, "  --cov       output breadth of coverage\n");
		}
		if (strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -d          use BWA-SW for local alignment\n");
		}
		if (strcmp(argv[0], "hapdiv") == 0 || strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -a INT      annotate sliding INT-mers [%d]\n", opt.hapdiv_k);
			fprintf(stderr, "  -w INT      k-mer step size for annotation [%d]\n", opt.hapdiv_w);
		}
		if (strcmp(argv[0], "sw") == 0 || strcmp(argv[0], "hapdiv") == 0 || strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -N INT      keep up to INT hits per DAWG node [%d]\n", opt.swo.n_best);
			fprintf(stderr, "  -m INT      min alignment score [%d]\n", opt.swo.min_sc);
			fprintf(stderr, "  -A INT      match score [%d]\n", opt.swo.match);
			fprintf(stderr, "  -B INT      mismatch penalty [%d]\n", opt.swo.mis);
			fprintf(stderr, "  -O INT      gap open penalty [%d]\n", opt.swo.gap_open);
			fprintf(stderr, "  -E INT      gap extension penalty; a k-long gap costs O+k*E [%d]\n", opt.swo.gap_ext);
			fprintf(stderr, "  -C NUM      size of the ranking cache [%d]\n", opt.swo.r2cache_size);
			fprintf(stderr, "  -y INT      ignore secondary hits scored INT lower than the best [%d]\n", opt.swo.e2e_drop);
		}
		if (strcmp(argv[0], "sw") == 0 || strcmp(argv[0], "search") == 0) {
			fprintf(stderr, "  -e          end-to-end mode (forcing -k to 1)\n");
			fprintf(stderr, "  -j INT      min MEM length to initiate alignment [%d]\n", opt.swo.min_mem_len);
			fprintf(stderr, "  -k INT      require INT-mer match at the end of alignment [%d]\n", opt.swo.end_len);
			fprintf(stderr, "  -u          write unmapped queries to PAF\n");
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
	if (!rb3_fmi_is_symmetric(&p.fmi)) {
		if (rb3_verbose >= 1)
			fprintf(stderr, "ERROR: BWT doesn't contain both strands\n");
		return 1;
	}
	for (j = o.ind + 1; j < argc; ++j) {
		p.fp = rb3_seq_open(argv[j], is_line);
		if (p.fp == 0) {
			if (rb3_verbose >= 1)
				fprintf(stderr, "ERROR: failed to load the sequence file '%s'\n", argv[j]);
			break;
		}
		kt_pipeline(2, worker_pipeline, &p, 3);
		rb3_seq_close(p.fp);
	}
	rb3_fmi_free(&p.fmi);
	return 0;
}
