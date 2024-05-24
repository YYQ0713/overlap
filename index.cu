#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "minimap.h"
#include "kvec.h"
#include "khash.h"

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

//sketch
__device__ unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

__device__ __forceinline__ void hash64(uint64_t key, uint64_t mask, uint64_t *hash_value)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	*hash_value = key;
}

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param p      minimizers; p->a[i].x is the 2k-bit hash value;
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
__device__ void mm_sketch(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos;
	mm128_t *buf, min = { UINT64_MAX, UINT64_MAX };

	assert(len > 0 && w > 0 && k > 0);
	//buf = (mm128_t*)alloca(w * 16);
	cudaMalloc(&buf, w * 16);
	memset(buf, 0xff, w * 16);
	//cudaMemset(buf, 0xff, w * 16);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			if (++l >= k)
				hash64(kmer[z], mask, &info.x), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		} else l = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		/*
		if (l == w + k - 1) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k) kv_push(mm128_t, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1) kv_push(mm128_t, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, *p, buf[j]);
			}
		}
		*/
		if (++buf_pos == w) buf_pos = 0;
	}
	//if (min.x != UINT64_MAX)
		//kv_push(mm128_t, *p, min);
}



mm_idx_t *mm_idx_init(int w, int k, int b)
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

void mm_idx_destroy(mm_idx_t *mi)
{
	int i;
	if (mi == 0) return;
	for (i = 0; i < 1<<mi->b; ++i) {
		free(mi->B[i].p);
		free(mi->B[i].a.a);
		kh_destroy(idx, (idxhash_t*)mi->B[i].h);
	}
	free(mi->B);
	if (mi->name)
		for (i = 0; i < mi->n; ++i) free(mi->name[i]);
	free(mi->len); free(mi->name);
	free(mi);
}

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
	int mask = (1<<mi->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) {
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

uint32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return UINT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

void mm_idx_set_max_occ(mm_idx_t *mi, float f)
{
	mi->freq_thres = f;
	mi->max_occ = mm_idx_cal_max_occ(mi, f);
}

/*********************************
 * Sort and generate hash tables *
 *********************************/

static void worker_post(void *g, long i, int tid)
{
	int j, start_a, start_p, n, n_keys;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>mi->b<<1, &absent);
			assert(absent && j - start_a == n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == start_p);

	// deallocate and clear b->a
	free(b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}
 
static void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}

/******************
 * Generate index *
 ******************/

#include <string.h>
#include <zlib.h>
#include "bseq.h"

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int tbatch_size, n_processed, keep_name;
	bseq_file_t *fp;
	uint64_t ibatch_size, n_read;
	mm_idx_t *mi;
} pipeline_t;

typedef struct {
    int n_seq;
	bseq1_t *seq;
	mm128_v a;
} step_t;

typedef struct {
	int n_processed, keep_name;
	bseq_file_t *fp;
	uint64_t ibatch_size, n_read;
	mm_idx_t *mi;
} cu_shared_t;

static void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &mi->B[a[i].x&mask].a;
		kv_push(mm128_t, *p, a[i]);
	}
}

/*
static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->n_read > p->ibatch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = bseq_read(p->fp, p->tbatch_size, &s->n_seq);
		if (s->seq) {
			uint32_t old_m = p->mi->n, m, n;
			assert((uint64_t)p->n_processed + s->n_seq <= INT32_MAX);
			m = n = p->mi->n + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m) {
				if (p->keep_name)
					p->mi->name = (char**)realloc(p->mi->name, m * sizeof(char*));
				p->mi->len = (int*)realloc(p->mi->len, m * sizeof(int));
			}
			for (i = 0; i < s->n_seq; ++i) {
				if (p->keep_name) {
					assert(strlen(s->seq[i].name) <= 254);
					p->mi->name[p->mi->n] = strdup(s->seq[i].name);
				}
				p->mi->len[p->mi->n++] = s->seq[i].l_seq;
				s->seq[i].rid = p->n_processed++;
				p->n_read += s->seq[i].l_seq;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			bseq1_t *t = &s->seq[i];
			mm_sketch(t->seq, t->l_seq, p->mi->w, p->mi->k, t->rid, &s->a);
			free(t->seq); free(t->name);
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
        step_t *s = (step_t*)in;
		mm_idx_add(p->mi, s->a.n, s->a.a);
		free(s->a.a); free(s);
	}
    return 0;
}
*/

/*
mm_idx_t *mm_idx_gen(bseq_file_t *fp, int w, int k, int b, int tbatch_size, int n_threads, uint64_t ibatch_size, int keep_name)
{
	pipeline_t pl;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.tbatch_size = tbatch_size;
	pl.keep_name = keep_name;
	pl.ibatch_size = ibatch_size;
	pl.fp = fp;
	if (pl.fp == 0) return 0;
	pl.mi = mm_idx_init(w, k, b);

	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_post(pl.mi, n_threads);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	return pl.mi;
}

mm_idx_t *mm_idx_build(const char *fn, int w, int k, int n_threads) // a simpler interface
{
	bseq_file_t *fp;
	mm_idx_t *mi;
	fp = bseq_open(fn);
	if (fp == 0) return 0;
	mi = mm_idx_gen(fp, w, k, MM_IDX_DEF_B, 1<<18, n_threads, UINT64_MAX, 1);
	mm_idx_set_max_occ(mi, 0.001);
	bseq_close(fp);
	return mi;
}
*/

/*************
 * index I/O *
 *************/

#define MM_IDX_MAGIC "MMI\1"

void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint32_t x[6];
	int i;
	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n, x[4] = mi->name? 1 : 0, x[5] = mi->max_occ;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 6, fp);
	fwrite(&mi->freq_thres, sizeof(float), 1, fp);
	fwrite(mi->len, 4, mi->n, fp);
	if (mi->name) {
		for (i = 0; i < mi->n; ++i) {
			uint8_t l;
			l = strlen(mi->name[i]);
			fwrite(&l, 1, 1, fp);
			fwrite(mi->name[i], 1, l, fp);
		}
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
}

mm_idx_t *mm_idx_load(FILE *fp)
{
	int i;
	char magic[4];
	uint32_t x[6];
	mm_idx_t *mi;
	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 6, fp) != 6) return 0;
	mi = mm_idx_init(x[0], x[1], x[2]);
	mi->n = x[3], mi->max_occ = x[5];
	fread(&mi->freq_thres, sizeof(float), 1, fp);
	mi->len = (int32_t*)malloc(mi->n * 4);
	fread(mi->len, 4, mi->n, fp);
	if (x[4]) { // has names
		mi->name = (char**)calloc(mi->n, sizeof(char*));
		for (i = 0; i < mi->n; ++i) {
			uint8_t l;
			fread(&l, 1, 1, fp);
			mi->name[i] = (char*)malloc(l + 1);
			fread(mi->name[i], 1, l, fp);
			mi->name[i][l] = 0;
		}
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		fread(&b->n, 4, 1, fp);
		b->p = (uint64_t*)malloc(b->n * 8);
		fread(b->p, 8, b->n, fp);
		fread(&size, 4, 1, fp);
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
			fread(x, 8, 2, fp);
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}
	return mi;
}

/**************
 * CUDA INDEX *
 **************/
__device__ size_t gpu_strlen(const char* str) {
    size_t len = 0;
    while (str[len] != '\0') {
        len++;
    }
    return len;
}

__device__ char* gpu_strdup(const char* src) {
    
    size_t len = 0;
    while (src[len] != '\0') {
        len++;
    }
    len++;

    char* dst;
    cudaMalloc(&dst, len * sizeof(char));

    for (size_t i = 0; i < len; i++) {
        dst[i] = src[i];
    }

    return dst;
}

__global__ void seq_init(stepcu_t *s, mm_idx_t *mi) {
	uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
	uint32_t n_threads = blockDim.x * gridDim.x;
	int count;
	for(uint32_t i = idx; i < s->n_seq; i += n_threads) {
		//count = gpu_strlen(s->seq[i].name);
		mi->name[i] = &s->seq.name[s->seq.bseq_info[i].name_s];
		mi->len[i] = s->seq.bseq_info[i].l_seq;
		s->seq.bseq_info[i].rid = i;
	}
}

__global__ void compute_sketch(stepcu_t *d_s, mm_idx_t *d_m, mm128_v *vec) {
	uint32_t idx = threadIdx.x + blockIdx.x * blockDim.x;
	uint32_t n_threads = blockDim.x * gridDim.x;
	stepcu_t *s = d_s;
	
	for (int i = idx; i < s->n_seq; i += n_threads) {
		bseq_info_t *t = &d_s->seq.bseq_info[i];
		mm_sketch(&d_s->seq.seq[t->name_s], t->l_seq, d_m->w, d_m->k, t->rid, &vec[i]);
	}
}

void cu_index_gen(bseq_file_t *fp, int w, int k, int b, uint64_t ibatch_size, int keep_name) {
	cu_shared_t h_cush;
	cu_shared_t *d_cush;
	mm_idx_t *d_mi;
	memset(&h_cush, 0, sizeof(cu_shared_t));
	h_cush.keep_name = keep_name;
	h_cush.ibatch_size = ibatch_size;
	h_cush.fp = fp;
	if (h_cush.fp == 0) exit(0);
	mm_idx_t *h_mi = mm_idx_init(w, k, b);

	int i;
	stepcu_t *h_s;
	if (h_cush.n_read > h_cush.ibatch_size) exit(0);
	h_s = (stepcu_t*)calloc(1, sizeof(stepcu_t));
	h_s->seq = bseqs_read(h_cush.fp, h_cush.ibatch_size, &h_s->n_seq);
	//printf("test\n");
	if (h_s->seq.seq) {
		uint32_t m = h_s->n_seq;
		kroundup32(m);
		stepcu_t *d_s = copy_step_t_to_gpu(h_s);
		print_stepcu<<<1, 1>>>(d_s);
		CUDA_CHECK(cudaMalloc((void**)&d_cush, sizeof(cu_shared_t)));
		CUDA_CHECK(cudaMemcpy(d_cush, &h_cush, sizeof(cu_shared_t), cudaMemcpyHostToDevice));
		if(h_cush.keep_name) {
			CUDA_CHECK(cudaMalloc((void**)&d_mi, sizeof(mm_idx_t)));
			CUDA_CHECK(cudaMemcpy(d_mi, h_mi, sizeof(mm_idx_t), cudaMemcpyHostToDevice));
			char **d_name;
    		CUDA_CHECK(cudaMalloc((void**)&d_name, m * sizeof(char*)));
			CUDA_CHECK(cudaMemcpy(&(d_mi->name), &d_name, sizeof(char**), cudaMemcpyHostToDevice));
		}
		int32_t *d_len;
    	CUDA_CHECK(cudaMalloc((void**)&d_len, m * sizeof(int32_t)));
		CUDA_CHECK(cudaMemcpy(&(d_mi->len), &d_len, sizeof(int32_t*), cudaMemcpyHostToDevice));

		mm_idx_bucket_t *d_B;
    	CUDA_CHECK(cudaMalloc((void**)&d_B, (1 << b) * sizeof(mm_idx_bucket_t)));
		CUDA_CHECK(cudaMemcpy(&(d_mi->B), &d_B, sizeof(mm_idx_bucket_t*), cudaMemcpyHostToDevice));

		mm128_v *gpu_vecinfo;
		CUDA_CHECK(cudaMalloc(&gpu_vecinfo, h_s->n_seq * sizeof(mm128_v)));

		dim3 block(32);
		seq_init<<<(h_s->n_seq + block.x - 1) / block.x, block>>>(d_s, d_mi);
		cudaDeviceSynchronize();
		compute_sketch<<<(h_s->n_seq + block.x - 1) / block.x, block>>>(d_s, d_mi, gpu_vecinfo);
	}






	//dont forget
	//CUDA_CHECK(cudaFree(d_cush));
	//CUDA_CHECK(cudaFree(d_len));
	//CUDA_CHECK(cudaFree(d_mi));
	//CUDA_CHECK(cudaFree(d_s));
	//CUDA_CHECK(cudaFree(d_name));
}
