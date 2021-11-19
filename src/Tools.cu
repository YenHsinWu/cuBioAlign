#include "Tools.h"

namespace BioAlign{
    ////////////////////
    // Tool functions //
    ////////////////////

    void LowerNode(Node *nd_ptr){
        char *seq_device, *seq_host;
        int sz = nd_ptr->Len() * sizeof(char);

        int thread_num = 512, block_num = (sz / thread_num) + 1;

        seq_host = (char*)malloc(sz + 1);

        cudaMalloc(&seq_device, sz);
        cudaMemcpy(seq_device, nd_ptr->Sequence(), sz, cudaMemcpyHostToDevice);

        Lowercase<<<block_num, thread_num>>>(seq_device, sz);

        cudaMemcpy(seq_host, seq_device, sz, cudaMemcpyDeviceToHost);

        nd_ptr->Update(seq_host);
    }

    void UpperNode(Node *nd_ptr){
        char *seq_device, *seq_host;
        int sz = nd_ptr->Len() * sizeof(char);

        int thread_num = 512, block_num = (sz / thread_num) + 1;

        seq_host = (char*)malloc(sz + 1);

        cudaMalloc(&seq_device, sz);
        cudaMemcpy(seq_device, nd_ptr->Sequence(), sz, cudaMemcpyHostToDevice);

        Uppercase<<<block_num, thread_num>>>(seq_device, sz);

        cudaMemcpy(seq_host, seq_device, sz, cudaMemcpyDeviceToHost);

        nd_ptr->Update(seq_host);
    }
};