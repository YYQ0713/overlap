#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "minimap.cuh"


#define MM_VERSION "20240519"

void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

int main(int argc, char *argv[]) {
    mm_mapopt_t opt;
	int i, c, k = 15, w = -1, b = MM_IDX_DEF_B, n_threads = 3, keep_name = 1, is_idx = 0;
	int tbatch_size = 100000000;
	uint64_t ibatch_size = 4000000000ULL;
	float f = 0.001;
	bseq_file_t *fp = 0;
	char *fnw = 0;
	FILE *fpr = 0, *fpw = 0;




    fp = bseq_open(argv[optind]);
    for (;;) {
		mm_idx_t *mi = 0;
        if (!bseq_eof(fp))
			mm_idx_gen(fp, w, k, b, tbatch_size, ibatch_size, keep_name);

        
    }

}
