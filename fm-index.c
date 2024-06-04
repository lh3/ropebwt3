#include <string.h>
#include <assert.h>
#include "rb3priv.h"
#include "fm-index.h"

rld_t *rb3_enc_plain2rld(int64_t len, const uint8_t *bwt)
{
	int64_t i, i0;
	rld_t *e;
	rlditr_t ei;
	e = rld_init(RB3_ASIZE, 3); // TODO: check alphabet
	rld_itr_init(e, &ei, 0);
	for (i0 = 0, i = 1; i <= len; ++i) {
		if (i == len || bwt[i0] != bwt[i]) {
			rld_enc(e, &ei, i - i0, bwt[i0]);
			i0 = i;
		}
	}
	rld_enc_finish(e, &ei);
	return e;
}

mrope_t *rb3_enc_plain2fmr(int64_t len, const uint8_t *bwt, int max_nodes, int block_len)
{
	int64_t i, i0, off, cnt[256], c6;
	int32_t c;
	mrope_t *r;
	rpcache_t cache;

	if (max_nodes <= 0) max_nodes = ROPE_DEF_MAX_NODES;
	if (block_len <= 0) block_len = ROPE_DEF_BLOCK_LEN;

	memset(cnt, 0, 256 * sizeof(int64_t));
	for (i = 0; i < len; ++i) ++cnt[bwt[i]];
	for (c = 6, c6 = 0; c < 256; ++c) c6 += cnt[c];
	assert(c6 == 0);

	memset(&cache, 0, sizeof(rpcache_t));
	r = mr_init(max_nodes, block_len, MR_SO_IO);
	for (c = 0, off = 0; c < 6; ++c) {
		for (i0 = 0, i = 1; i <= cnt[c]; ++i) {
			if (i == cnt[c] || bwt[off+i0] != bwt[off+i]) {
				rope_insert_run(r->r[c], i0, bwt[off+i0], i - i0, &cache);
				i0 = i;
			}
		}
		off += cnt[c];
	}
	return r;
}
