#include <string.h>
#include "bwtl.h"
#include "kalloc.h"
#include "libsais16.h" // for libsais16()
#include "io.h" // for rb3_nt6_table[]

static uint32_t bwtl_cnt_table[256];

void rb3_bwtl_init(void)
{
	int i, c;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (c = 0; c != 4; ++c)
			x |= (((i&3) == c) + ((i>>2&3) == c) + ((i>>4&3) == c) + (i>>6 == c)) << (c<<3);
		bwtl_cnt_table[i] = x;
	}
}

rb3_bwtl_t *rb3_bwtl_gen(void *km, int len, const uint8_t *seq)
{
	rb3_bwtl_t *b;
	int32_t i;
	b = Kcalloc(km, rb3_bwtl_t, 1); // allocate long-term memory first to reduce memory fragmentation
	b->km = km;
	b->seq_len = len;
	b->bwt_size = (len + 15) / 16;
	b->bwt = Kcalloc(km, uint32_t, b->bwt_size);
	b->n_occ = (len + 16) / 16 * 4;
	b->occ = Kcalloc(km, int32_t, b->n_occ);

	{ // calculate b->bwt
		uint8_t *s;
		uint16_t *s16;
		s16 = Kcalloc(km, uint16_t, len);
		for (i = 0; i < len; ++i) s16[i] = rb3_nt6_table[seq[i]];
		b->sa = Kcalloc(km, int32_t, len + 1);
		libsais16(s16, &b->sa[1], len, 0, 0);
		b->sa[0] = len;
		s = Kcalloc(km, uint8_t, len + 1);
		for (i = 0; i <= len; ++i) {
			if (b->sa[i] == 0) b->primary = i;
			else s[i] = s16[b->sa[i] - 1] - 1;
		}
		kfree(km, s16);
		for (i = b->primary; i < len; ++i) s[i] = s[i + 1];
		for (i = 0; i < len; ++i)
			b->bwt[i>>4] |= s[i] << ((15 - (i&15)) << 1);
		kfree(km, s);
	}
	{ // calculate b->occ
		int32_t c[4];
		memset(c, 0, 16);
		for (i = 0; i < len; ++i) {
			if (i % 16 == 0)
				memcpy(b->occ + (i/16) * 4, c, 16);
			++c[rb3_bwtl_B0(b, i)];
		}
		if (i % 16 == 0)
			memcpy(b->occ + (i/16) * 4, c, 16);
		memcpy(&b->acc[1], c, 16);
		b->acc[0] = 1;
		for (i = 1; i < 5; ++i) b->acc[i] += b->acc[i-1];
	}
	return b;
}

void rb3_bwtl_rank1a(const rb3_bwtl_t *bwt, int32_t k, int32_t cnt[4])
{
	uint32_t x, b;
	if (k > bwt->primary) --k; // because $ is not in bwt
	memcpy(cnt, bwt->occ + (k>>4<<2), 16);
	if (k % 16 == 0) return;
	--k;
	b = bwt->bwt[k>>4] & ~((1U<<((~k&15)<<1)) - 1);
	x = bwtl_cnt_table[b&0xff] + bwtl_cnt_table[b>>8&0xff] + bwtl_cnt_table[b>>16&0xff] + bwtl_cnt_table[b>>24];
	x -= 15 - (k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

void rb3_bwtl_rank2a(const rb3_bwtl_t *bwt, int32_t k, int32_t l, int32_t cntk[4], int32_t cntl[4])
{
	rb3_bwtl_rank1a(bwt, k, cntk);
	rb3_bwtl_rank1a(bwt, l, cntl);
}

void rb3_bwtl_destroy(rb3_bwtl_t *bwt)
{
	kfree(bwt->km, bwt->occ);
	kfree(bwt->km, bwt->bwt);
	kfree(bwt->km, bwt->sa);
	kfree(bwt->km, bwt);
}
