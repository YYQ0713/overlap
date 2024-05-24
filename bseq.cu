#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "bseq.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

extern unsigned char seq_nt4_table[256];

struct bseq_file_s {
	int is_eof;
	gzFile fp;
	kseq_t *ks;
};

bseq_file_t *bseq_open(const char *fn)
{
	bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (f == 0) return 0;
	fp = (bseq_file_t*)calloc(1, sizeof(bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void bseq_close(bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

/*
bseq1_t *bseq_read(bseq_file_t *fp, uint64_t chunk_size, int *n_)
{
	int size = 0, m, n;
	bseq1_t *seqs;
	kseq_t *ks = fp->ks;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		bseq1_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = (bseq1_t*)realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->name = strdup(ks->name.s);
		s->seq = strdup(ks->seq.s);
		s->l_seq = ks->seq.l;
		size += seqs[n++].l_seq;
		if (size >= chunk_size) break;
	}
	if (n == 0) fp->is_eof = 1;
	*n_ = n;
	return seqs;
}
*/


bseqs_t bseqs_read(bseq_file_t *fp, uint64_t chunk_size, int *n_)
{
	uint32_t total_bp = 0, total_name = 0, m, n;
	uint32_t total_m_bp, total_m_name;
	uint32_t name_len, bp_len; 
	bseqs_t seqs;
	memset(&seqs, 0, sizeof(bseqs_t));
	kseq_t *ks = fp->ks;
	m = n = 0;
	total_m_bp = total_m_name = 0;

	while (kseq_read(ks) >= 0) {
		bseq_info_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs.bseq_info = (bseq_info_t*)realloc(seqs.bseq_info, m * sizeof(bseq_info_t));
            if (!seqs.bseq_info) {
                // Handle memory allocation failure
                fprintf(stderr, "Info Memory allocation failed\n");
				free(seqs.bseq_info); free(seqs.name); free(seqs.seq);
                exit(1);
            }
		}
		s = &seqs.bseq_info[n];

		name_len = strlen(ks->name.s) + 1;
		while (total_name + name_len >= total_m_name) {
			total_m_name = total_m_name ? total_m_name << 1 : 1024;
		}
		seqs.name = (char*)realloc(seqs.name, total_m_name);
		if (!seqs.name) {
			// Handle memory allocation failure
			fprintf(stderr, "Name Memory allocation failed\n");
			free(seqs.bseq_info); free(seqs.name); free(seqs.seq);
			exit(1);
		}
		memcpy(seqs.name + total_name, ks->name.s, name_len);
		s->l_name = name_len - 1;
		s->name_s = total_name;
		s->name_e = total_name + name_len;
		total_name += name_len;

		bp_len = ks->seq.l + 1;
		while (total_bp + bp_len >= total_m_bp) {
			total_m_bp = total_m_bp ? total_m_bp << 1 : 1048576;
		}
		seqs.seq = (char*)realloc(seqs.seq, total_m_bp);
		if (!seqs.seq) {
			// Handle memory allocation failure
			fprintf(stderr, "Seqs Memory allocation failed\n");
			free(seqs.bseq_info); free(seqs.name); free(seqs.seq);
			exit(1);
		}
		memcpy(seqs.seq + total_bp, ks->seq.s, bp_len);
		s->l_seq = ks->seq.l;
		s->seq_s = total_bp;
		s->seq_e = total_bp + bp_len;
		total_bp += bp_len;
		n++;
		
		
		if (total_bp >= chunk_size) break;
	}
	if (n == 0) fp->is_eof = 1;
	*n_ = n;
	seqs.tot_name = total_name;
	seqs.tot_bp = total_bp;
	return seqs;
}

int bseq_eof(bseq_file_t *fp)
{
	return fp->is_eof;
}
