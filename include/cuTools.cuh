#ifndef __CUTOOLS_H__
#define __CUTOOLS_H__

#include <cuda_runtime.h>

namespace BioAlign{
    /////////////////////////
    // CUDA Tool functions //
    /////////////////////////

    __global__ void Lowercase(char*, int);
    __global__ void Uppercase(char*, int);

    __global__ void LongestSubstringLen(int*, char*, char*, int, int);
};

#endif