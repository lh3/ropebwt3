#include <stdio.h>
#include "rld0.h"
int main(int argc, char *argv[]) {
	if (argc < 2) return 1;
	rld_t *e = rld_restore(argv[1]);
	rlditr_t ei; // iterator
	rld_itr_init(e, &ei, 0);
	uint64_t ok[6];
	rld_rank1a(e, 100000, ok);
	rld_destroy(e);
	return 0;
}
