#ifndef RSEQ_H
#define RSEQ_H
#include "bseq.cuh"
#include "minimap.cuh"

typedef struct {
	int tbatch_size, n_processed, keep_name;
	bseq_file_t *fp;
	uint64_t ibatch_size, n_read;
	mm_idx_t *mi;
} share_t;

typedef struct {
    int n_seq;
	bseq1_t *seq;
	mm128_v a;
} step_t;

__host__ void* read_seq(share_t *shared);
#endif
