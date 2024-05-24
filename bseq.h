#ifndef MM_BSEQ_H
#define MM_BSEQ_H

#include <stdint.h>

struct bseq_file_s;
typedef struct bseq_file_s bseq_file_t;

typedef struct {
	int l_seq, rid;
	char *name, *seq;
} bseq1_t;

typedef struct {
	int l_seq, rid, l_name;
	uint32_t seq_s, seq_e;
	uint32_t name_s, name_e;
	//char *name, *seq;
} bseq_info_t;

typedef struct {
	char *name, *seq;
	bseq_info_t *bseq_info;
	uint32_t tot_bp, tot_name;
} bseqs_t;

bseq_file_t *bseq_open(const char *fn);
void bseq_close(bseq_file_t *fp);
bseq1_t *bseq_read(bseq_file_t *fp, uint64_t chunk_size, int *n_);
bseqs_t bseqs_read(bseq_file_t *fp, uint64_t chunk_size, int *n_);
int bseq_eof(bseq_file_t *fp);

#endif
