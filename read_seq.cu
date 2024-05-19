#include "read_seq.cuh"
#include "bseq.cuh"
#include "minimap.cuh"
#include "khash.cuh"
#include <assert.h>



__host__ void* read_seq(share_t *shared) {
    int i;
    share_t *p = shared;
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
}