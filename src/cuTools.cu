#include "cuTools.cuh"

namespace BioAlign{
    /////////////////////////
    // CUDA Tool functions //
    /////////////////////////

    __global__ void Lowercase(char *str, int n){
        int idx = blockDim.x * blockIdx.x + threadIdx.x;

        if(idx < n){
            if(str[idx] >= 65 || str[idx] <= 90)
                str[idx] += 32;
        }
    }

    __global__ void Uppercase(char *str, int n){
        int idx = blockDim.x * blockIdx.x + threadIdx.x;

        if(idx < n){
            if(str[idx] >= 97 || str[idx] <= 122)
                str[idx] -= 32;
        }
    }
};