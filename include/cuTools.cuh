#ifndef __CUTOOLS_H__
#define __CUTOOLS_H__

#include <stdio.h>
#include <string.h>
#include <cuda_runtime.h>

namespace BioAlign{
    /////////////////////////
    // CUDA Tool functions //
    /////////////////////////
    __global__ void VectorHalfAdd(float*, float*, float*, int);

    __global__ void Lowercase(char*, int);
    __global__ void Uppercase(char*, int);

    __global__ void LongestSubstringLen(int*, char*, char*, int, int);

    __device__ void SliceCString(char*, char*, int, int);
    __device__ bool CompareCStrings(char*, char*, int);
    __global__ void KmersFrequency(float*, char*, char*, char*, char*, int, int, int);

    __global__ void MinPos(double*, int, int, int*, double*);
};

#endif