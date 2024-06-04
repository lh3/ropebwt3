#include "rb3priv.h"

rld_t *rb3_enc_plain2rld(int64_t len, const uint8_t *bwt)
{
	int64_t i, i0;
	rld_t *e;
	rlditr_t ei;
	e = rld_init(6, 3);
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
