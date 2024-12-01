#ifndef BRE_H
#define BRE_H

#include <stdint.h>

#define BRE_ERR_TRUNCATE  (-1)
#define BRE_ERR_INCONSIS  (-2)

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint8_t b_per_sym; // bytes per symbol
	uint8_t b_per_run; // bytes per run
	uint8_t atype;     // alphabet type
	uint8_t mtype;     // type of multi-string BWT
	int64_t asize;     // alphabet size
	int64_t l_aux;     // length of auxiliary data
	uint8_t *aux;      // auxiliary data
} bre_hdr_t;

struct bre_file_s;
typedef struct bre_file_s bre_file_t;

typedef enum { BRE_AT_UNKNOWN=0, BRE_AT_ASCII, BRE_AT_DNA6, BRE_AT_DNA16 } bre_atype_t;

extern int bre_verbose;

/**
 * Open a BRE file for reading
 *
 * @param fn   file name. NULL or "-" for stdin
 *
 * @return file handler. NULL for error
 */
bre_file_t *bre_open_read(const char *fn);

// the first four magic bytes have been read
bre_file_t *bre_open_no_magic(FILE *fp);

/**
 * Read a run
 *
 * In a BRE file, a run longer than (1<<b_per_run)-1 may be split into separate
 * records. This function combines them. As a result, _c_ returned by the
 * function is different upon each call.
 *
 * @param f    file handler
 * @param c    symbol (out). -1 if nothing is read
 *
 * @return run length if positive. 0 if nothing is read. -1 for error
 */
int64_t bre_read(bre_file_t *f, int64_t *c);

/**
 * Fill a BRE header structure
 *
 * @param h          header to fill (out). bre_hdr_t::aux is discarded (be careful of memory leak!)
 * @param at         alphabet type (atype)
 * @param b_per_run  bytes per run
 *
 * @return always 0 at the moment
 */
int bre_hdr_init(bre_hdr_t *h, bre_atype_t at, int32_t b_per_run);

/**
 * Create a BRE file
 *
 * @param fn   file name
 * @param hdr  BRE header
 *
 * @return file handker, NULL for error
 */
bre_file_t *bre_open_write(const char *fn, const bre_hdr_t *hdr);

/**
 * Write a run
 *
 * This function accumulates run length if _c_ is the same in consecutive
 * function calls. It also splits a run longer than (1<<b_per_run)-1 into
 * multiple records in the output file. It is ok to write a long run through
 * multiple function calls. It is also ok to write a long run. This function
 * automatically handles both cases.
 *
 * Remember to call bre_close() at the end; otherwise some runs will not be
 * written to the output file.
 *
 * @param f    file handler
 * @param c    symbol
 * @param l    run length
 *
 * @return 0 for success and -1 for error
 */
int bre_write(bre_file_t *f, int64_t c, int64_t l);

/**
 * Close a file handler
 *
 * @param f    file handler
 */
void bre_close(bre_file_t *f);

int bre_error(const bre_file_t *f);
int64_t bre_n_rec(const bre_file_t *f);
int64_t bre_n_sym(const bre_file_t *f);
int64_t bre_n_run(const bre_file_t *f);
const bre_hdr_t *bre_get_hdr(const bre_file_t *f);

#ifdef __cplusplus
}
#endif

#endif
