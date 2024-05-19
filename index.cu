#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "minimap.cuh"
#include "khash.cuh"
#include "read_seq.cuh"

__host__ mm_idx_t* mm_idx_init(int w, int k, int b)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b;
	mi->max_occ = UINT32_MAX;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	return mi;
}

__host__ void *mm_idx_gen(bseq_file_t *fp, int w, int k, int b, int tbatch_size, uint64_t ibatch_size, int keep_name)
{
	share_t pl;
	memset(&pl, 0, sizeof(share_t));
	pl.tbatch_size = tbatch_size;
	pl.keep_name = keep_name;
	pl.ibatch_size = ibatch_size;
	pl.fp = fp;
	if (pl.fp == 0) return 0;
	pl.mi = mm_idx_init(w, k, b);

    read_seq(&pl);

	
}