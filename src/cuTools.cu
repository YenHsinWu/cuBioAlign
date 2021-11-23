#include "cuTools.cuh"

namespace BioAlign{
    /////////////////////////
    // CUDA Tool functions //
    /////////////////////////

    __global__ void Lowercase(char *str, int n){
        int idx = blockDim.x * blockIdx.x + threadIdx.x;

        if(idx < n){
            if(str[idx] >= 65 && str[idx] <= 90)
                str[idx] += 32;
        }
    }

    __global__ void Uppercase(char *str, int n){
        int idx = blockDim.x * blockIdx.x + threadIdx.x;

        if(idx < n){
            if(str[idx] >= 97 && str[idx] <= 122)
                str[idx] -= 32;
        }
    }

    __global__ void LongestSubstringLen(int *substr_lens, char *str_a, char *str_b, int len_a, int len_b){
        int idx = blockDim.x * blockIdx.x + threadIdx.x;
        int a, b;
        int current_len;

        if(idx < len_a){
            a = idx;
            current_len = 0;
            substr_lens[idx] = 0;
            b = 0;
            while(b < len_b){
                if(str_a[a] == str_b[b]){
                    current_len ++;
                    a ++;
                    b ++;

                    if(a >= len_a){
                        if(current_len > substr_lens[idx])
                            substr_lens[idx] = current_len;
                        break;
                    }
                }
                else{
                    if(current_len == 0)
                        b ++;
                    else{
                        if(current_len > substr_lens[idx]){
                            substr_lens[idx] = current_len;
                            current_len = 0;
                            a = idx;
                            b ++;
                        }
                        else{
                            current_len = 0;
                            a = idx;
                            b ++;
                        }
                    }
                }
            }
        }
    }
};