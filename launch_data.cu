#include <cuda_runtime.h>
#include <iostream>
#include <cstring>
#include <stdio.h>
#include "minimap.h"

// 分配并拷贝step_t结构体到GPU
__host__ stepcu_t* copy_step_t_to_gpu(const stepcu_t *cpu_step) {
    stepcu_t *gpu_step;
    CUDA_CHECK(cudaMalloc(&gpu_step, sizeof(stepcu_t)));

    CUDA_CHECK(cudaMemcpy(gpu_step, cpu_step, sizeof(stepcu_t), cudaMemcpyHostToDevice));

    if (cpu_step->seq.bseq_info != nullptr) {
        bseq_info_t *gpu_info;
        CUDA_CHECK(cudaMalloc((void**)&gpu_info, cpu_step->n_seq * sizeof(bseq_info_t)));
        CUDA_CHECK(cudaMemcpy(gpu_info, cpu_step->seq.bseq_info, cpu_step->n_seq * sizeof(bseq_info_t), cudaMemcpyHostToDevice));

        bseq_info_t **gpu_info_ptr;
        CUDA_CHECK(cudaMalloc(&gpu_info_ptr, sizeof(bseq_info_t*)));
        CUDA_CHECK(cudaMemcpy(gpu_info_ptr, &gpu_info, sizeof(bseq_info_t*), cudaMemcpyHostToDevice));

        CUDA_CHECK(cudaMemcpy(&(gpu_step->seq.bseq_info), &gpu_info, sizeof(bseq_info_t*), cudaMemcpyHostToDevice));
    }
    
    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
    if (cpu_step->seq.name != nullptr) {
        char *gpu_name;
        //size_t name_len = strlen(cpu_step->seq.name) + 1;
        size_t name_len = cpu_step->seq.tot_name;
        CUDA_CHECK(cudaMalloc(&gpu_name, name_len));
        CUDA_CHECK(cudaMemcpy(gpu_name, cpu_step->seq.name, name_len, cudaMemcpyHostToDevice));

        char **gpu_name_ptr;
        CUDA_CHECK(cudaMalloc(&gpu_name_ptr, sizeof(char*)));
        CUDA_CHECK(cudaMemcpy(gpu_name_ptr, &gpu_name, sizeof(char*), cudaMemcpyHostToDevice));

        CUDA_CHECK(cudaMemcpy(&gpu_step->seq.name, &gpu_name, sizeof(char*), cudaMemcpyHostToDevice));
    }

    if (cpu_step->seq.seq != nullptr) {
        char *gpu_seq_data;
        size_t seq_len = cpu_step->seq.tot_bp;
        CUDA_CHECK(cudaMalloc(&gpu_seq_data, seq_len));
        CUDA_CHECK(cudaMemcpy(gpu_seq_data, cpu_step->seq.seq, seq_len, cudaMemcpyHostToDevice));

        char **gpu_seq_ptr;
        CUDA_CHECK(cudaMalloc(&gpu_seq_ptr, sizeof(char*)));
        CUDA_CHECK(cudaMemcpy(gpu_seq_ptr, &gpu_seq_data, sizeof(char*), cudaMemcpyHostToDevice));

        CUDA_CHECK(cudaMemcpy(&gpu_step->seq.seq, &gpu_seq_data, sizeof(char*), cudaMemcpyHostToDevice));
    }
    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
    return gpu_step;
}

__global__ void print_stepcu(stepcu_t *s) {
    //printf("%s\n", s->seq.name + s->seq.bseq_info[1].name_s);
    printf("%s\n", &(s->seq.seq[s->seq.bseq_info[s->n_seq - 2].seq_s]));
}

__global__ void print_idxcu(mm_idx_t *mi) {
    printf("%s\n", mi->name[0]);
}
