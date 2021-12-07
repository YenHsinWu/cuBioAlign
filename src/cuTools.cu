#include "cuTools.cuh"

namespace BioAlign{
    /////////////////////////
    // CUDA Tool functions //
    /////////////////////////

    __global__ void VectorHalfAdd(float *vec_a, float *vec_b, float *vec_c, int n){
        int idx = blockDim.x * blockIdx.x + threadIdx.x;

        if(idx < n)
            vec_c[idx] = (vec_a[idx] + vec_b[idx]) * 0.5;
    }

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

    __device__ char* SliceCString(char *source, int start, int end){
        char *target = (char*)malloc((end - start + 1) * sizeof(char));
        for(int i = start; i < end; i ++)
            target[i - start] = source[i];
        target[end] = '\0';

        return target;
    }

    __device__ bool CompareCStrings(char *str1, char *str2, int len){
        for(int i = 0; i < len; i ++){
            if(str1[i] != str2[i])
                return false;
        }

        return true;
    }

    __global__ void KmersFrequency(float* freq, char* seq, char* kmers, char *subseq, char *subkmers, int k, int kmers_size, int seq_len){
        int idx = blockDim.x * blockIdx.x + threadIdx.x;
        int cnt_idx = threadIdx.x;
        int tmp = 0;
        __shared__ int cnt[512];

        if(idx < kmers_size){
            subkmers = SliceCString(kmers, k * idx, k * (idx + 1));
            for(int i = 0; i < (seq_len - k + 1); i ++){
                subseq = SliceCString(seq, i, i + k);
                if(CompareCStrings(subseq, subkmers, k))
                    freq[idx] += 1;
            }
        }

        while(idx < kmers_size){
            tmp += freq[idx];
            idx += blockDim.x;
        }
        cnt[cnt_idx] = tmp;
        __syncthreads();

        int b_idx = blockDim.x / 2;
        while(b_idx != 0){
            if(cnt_idx < b_idx)
                cnt[cnt_idx] += cnt[cnt_idx + b_idx];
            __syncthreads();

            b_idx /= 2;
        }

        if(cnt_idx == 0){
            for(int i = 0; i < kmers_size; i ++)
                freq[i] /= cnt[0]; 
        }
    }

    __global__ void MinPos(double *elems, int row_num, int col_num, int *min_idxs, double* mins){
        int r_idx = blockDim.x * blockIdx.x + threadIdx.x;
        double minimum = 10e8;
        int min_idx;

        if(r_idx < row_num){
            for(int c_idx = r_idx + 1; c_idx < col_num; c_idx ++){
                if(elems[r_idx * row_num + c_idx] < minimum){
                    minimum = elems[r_idx * row_num + c_idx];
                    min_idx = c_idx;
                }
            }
        }

        min_idxs[r_idx] = min_idx;
        mins[r_idx] = minimum;
    }
};