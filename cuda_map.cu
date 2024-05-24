#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bseq.h"
#include "kvec.h"
#include "minimap.h"
#include "sdust.h"
void mm_mapopt_init(mm_mapopt_t *opt)
{
	opt->radius = 500;
	opt->max_gap = 10000;
	opt->min_cnt = 4;
	opt->min_match = 40;
	opt->sdust_thres = 0;
	opt->flag = MM_F_WITH_REP;
	opt->merge_frac = .5;
}

/****************************
 * Find approxiate mappings *
 ****************************/

struct mm_tbuf_s { // per-thread buffer
	mm128_v mini; // query minimizers
	mm128_v coef; // Hough transform coefficient
	mm128_v intv; // intervals on sorted coef
	uint32_v reg2mini;
	uint32_v rep_aux;
	sdust_buf_t *sdb;
	// the following are for computing LIS
	uint32_t n, m;
	uint64_t *a;
	size_t *b, *p;
	// final output
	kvec_t(mm_reg1_t) reg;
};

typedef struct {
	int batch_size, n_processed;
	const mm_mapopt_t *opt;
	bseq_file_t *fp;
	const mm_idx_t *mi;
} streams_t;

typedef struct {
	const streams_t *p;
    int n_seq;
	bseq1_t *seq;
	int *n_reg;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
} step_t;


__host__ step_t* stream_read_seq(streams_t *st) {
	int i, j;
    streams_t *p = st;
        step_t *s;
    s = (step_t*)calloc(1, sizeof(step_t));
    s->seq = bseq_read(p->fp, p->batch_size, &s->n_seq);
    if (s->seq) {
        s->p = p;
        for (i = 0; i < s->n_seq; ++i)
            s->seq[i].rid = p->n_processed++;
        //s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
        //for (i = 0; i < p->n_threads; ++i)
        //    s->buf[i] = mm_tbuf_init();
        //s->n_reg = (int*)calloc(s->n_seq, sizeof(int));
        //->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t*));
        return s;
    } else free(s);
}

int cuda_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int ibatch_size)
{
	streams_t pl;
	memset(&pl, 0, sizeof(streams_t));
	pl.fp = bseq_open(fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.mi = idx;
	pl.batch_size = ibatch_size;
	
    step_t *s = stream_read_seq(&pl);

	bseq_close(pl.fp);
	return 0;
}