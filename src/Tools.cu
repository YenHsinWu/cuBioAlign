#include "Tools.h"

namespace BioAlign{
    ////////////////////
    // Tool functions //
    ////////////////////

    void LowerNode(Node *nd_ptr){
        char *seq_device, *seq_host;
        int sz = nd_ptr->Len() * sizeof(char);

        int thread_num = 512, block_num = (sz / thread_num) + 1;

        seq_host = (char*)malloc(sz);

        cudaMalloc(&seq_device, sz);
        cudaMemcpy(seq_device, nd_ptr->Sequence(), sz, cudaMemcpyHostToDevice);

        Lowercase<<<block_num, thread_num>>>(seq_device, sz);

        cudaMemcpy(seq_host, seq_device, sz, cudaMemcpyDeviceToHost);
        
        nd_ptr->Update(seq_host);

        free(seq_host);
        cudaFree(seq_device);
    }

    void UpperNode(Node *nd_ptr){
        char *seq_device, *seq_host;
        int sz = nd_ptr->Len() * sizeof(char);

        int thread_num = 512, block_num = (sz / thread_num) + 1;

        seq_host = (char*)malloc(sz);

        cudaMalloc(&seq_device, sz);
        cudaMemcpy(seq_device, nd_ptr->Sequence(), sz, cudaMemcpyHostToDevice);

        Uppercase<<<block_num, thread_num>>>(seq_device, sz);

        cudaMemcpy(seq_host, seq_device, sz, cudaMemcpyDeviceToHost);
        
        nd_ptr->Update(seq_host);

        free(seq_host);
        cudaFree(seq_device);
    }

    void ReadFasta(std::string fname, std::vector<Node*>& nodes_ptr){
        try{
            std::size_t pos = fname.find_last_of('.');
            std::string file_extension = fname.substr(pos + 1);

            if((file_extension.compare("fa") != 0) && (file_extension.compare("fasta") != 0))
                throw "In ReadFasta(std::string, std::vector<Node*>&) : File type not supported.";

            std::ifstream input_file;
            std::string line;

            input_file.open(fname);
            if(!input_file)
                throw "In ReadFasta(std::string, std::vector<Node*>&) : Cannot open file. Maybe wrong file name is given.";

            while(!input_file.eof()){
                std::getline(input_file, line);
                if(line[0] == '>')
                    std::getline(input_file, line);

                nodes_ptr.push_back(new Node(line));
            }

            input_file.close();
        }
        catch(const char* error_message){
            std::cout << error_message << std::endl;
        }
    }

    void GenerateLeaves(int n, std::vector<Vertex*>& vtxs_ptr){
        for(int i = 0; i < n; i ++)
            vtxs_ptr.push_back(new Vertex(std::to_string(i), 0.0));
    }

    double ACSDistance(Node *ndptr_a, Node *ndptr_b){
        int n = ndptr_a->Len(), m = ndptr_b->Len(), sum;
        int *lens_ab, *lens_aa, *lens_ba, *lens_bb, *lens_device;
        char *seq_a_device, *seq_b_device;
        double result, l_ab, l_aa, l_ba, l_bb, d_ab, d_ba;

        int thread_num = 512, block_num;

        lens_ab = (int*)malloc(n * sizeof(int));
        lens_aa = (int*)malloc(n * sizeof(int));
        lens_ba = (int*)malloc(m * sizeof(int));
        lens_bb = (int*)malloc(m * sizeof(int));

        cudaMalloc(&lens_device, n * sizeof(int));
        cudaMalloc(&seq_a_device, n * sizeof(char));
        cudaMalloc(&seq_b_device, m * sizeof(char));

        cudaMemcpy(seq_a_device, ndptr_a->Sequence(), n * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(seq_b_device, ndptr_b->Sequence(), m * sizeof(char), cudaMemcpyHostToDevice);

        block_num = (n / thread_num) + 1;
        LongestSubstringLen<<<block_num, thread_num>>>(lens_device, seq_a_device, seq_b_device, n, m);
        cudaMemcpy(lens_ab, lens_device, n * sizeof(int), cudaMemcpyDeviceToHost);

        block_num = (m / thread_num) + 1;
        LongestSubstringLen<<<block_num, thread_num>>>(lens_device, seq_b_device, seq_a_device, m, n);
        cudaMemcpy(lens_ba, lens_device, m * sizeof(int), cudaMemcpyDeviceToHost);

        for(int i = 0; i < n; i ++)
            lens_aa[n - i - 1] = i + 1;

        for(int i = 0; i < m; i ++)
            lens_bb[m - i - 1] = i + 1;

        sum = 0;
        for(int i = 0; i < n; i ++)
            sum += lens_ab[i];
        l_ab = (double)sum / n;

        sum = 0;
        for(int i = 0; i < n; i ++)
            sum += lens_aa[i];
        l_aa = (double)sum / n;

        sum = 0;
        for(int i = 0; i < m; i ++)
            sum += lens_ba[i];
        l_ba = (double)sum / m;

        sum = 0;
        for(int i = 0; i < m; i ++)
            sum += lens_bb[i];
        l_bb = (double)sum / m;

        d_ab = (log(m) / l_ab) - (log(n) / l_aa);
        d_ba = (log(n) / l_ba) - (log(m) / l_bb);

        result = (d_ab + d_ba) * 0.5;

        free(lens_ab);
        free(lens_aa);
        free(lens_ba);
        free(lens_bb);
        cudaFree(lens_device);
        cudaFree(seq_a_device);
        cudaFree(seq_b_device);

        return result;
    }

    void InitKmers(char *kmers, int k){
        std::vector<std::string> init_kmers = {"A", "C", "G", "T"};
        std::vector<std::string> next_gene = {"A", "C", "G", "T"};
        std::vector<std::string> kmers_tmp;
        std::string string_kmers;
        int cnt = 0;

        while(cnt < (k - 1)){
            for(int i = 0; i < init_kmers.size(); i ++){
                for(int j = 0; j < next_gene.size(); j ++){
                    kmers_tmp.push_back(init_kmers[i] + next_gene[j]);
                }
            }
            init_kmers = kmers_tmp;
            kmers_tmp.clear();
            cnt ++;
        }

        for(int i = 0; i < init_kmers.size(); i ++){
            string_kmers += init_kmers[i];
        }

        strncpy(kmers, string_kmers.c_str(), string_kmers.size());
        kmers[string_kmers.size()] = '\0';
    }

    double FFPDistance(Node *ndptr_a, Node *ndptr_b, int k){
        int kmers_size = 1 << (2 * k);
        float result, entropy_a = 0, entropy_b = 0;
        float *freq_a, *freq_b, *freq_a_device, *freq_b_device;
        float *freq_m, *freq_m_device;
        char *kmers_host, *kmers_device;
        char *seq_a_device, *seq_b_device, *subseq_device, *subkmers_device;

        int thread_num = 512, block_num = (kmers_size / thread_num) + 1;

        //UpperNode(ndptr_a);
        //UpperNode(ndptr_b);

        freq_a = (float*)malloc(kmers_size * sizeof(float));
        freq_b = (float*)malloc(kmers_size * sizeof(float));
        freq_m = (float*)malloc(kmers_size * sizeof(float));

        kmers_host = (char*)malloc(k * kmers_size * sizeof(char) + 1);

        cudaMalloc(&seq_a_device, ndptr_a->Len() * sizeof(char));
        cudaMalloc(&seq_b_device, ndptr_b->Len() * sizeof(char));
        cudaMalloc(&subseq_device, k * sizeof(char) + 1);

        cudaMalloc(&kmers_device, k * kmers_size * sizeof(char) + 1);
        cudaMalloc(&subkmers_device, k * sizeof(char) + 1);
        cudaMalloc(&freq_a_device, kmers_size * sizeof(float));
        cudaMalloc(&freq_b_device, kmers_size * sizeof(float));
        cudaMalloc(&freq_m_device, kmers_size * sizeof(float));

        for(int i = 0; i < kmers_size; i ++){
            freq_a[i] = 0;
            freq_b[i] = 0;
        }

        cudaMemcpy(seq_a_device, ndptr_a->Sequence(), ndptr_a->Len() * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(seq_b_device, ndptr_b->Sequence(), ndptr_b->Len() * sizeof(char), cudaMemcpyHostToDevice);

        InitKmers(kmers_host, k);
        cudaMemcpy(kmers_device, kmers_host, k * kmers_size * sizeof(char) + 1, cudaMemcpyHostToDevice);

        cudaMemcpy(freq_a_device, freq_a, kmers_size * sizeof(float), cudaMemcpyHostToDevice);
        KmersFrequency<<<block_num, thread_num>>>(freq_a_device, seq_a_device, kmers_device, subseq_device, subkmers_device, k, kmers_size, ndptr_a->Len());
        cudaMemcpy(freq_a, freq_a_device, kmers_size * sizeof(float), cudaMemcpyDeviceToHost);

        cudaMemcpy(freq_b_device, freq_b, kmers_size * sizeof(float), cudaMemcpyHostToDevice);
        KmersFrequency<<<block_num, thread_num>>>(freq_b_device, seq_b_device, kmers_device, subseq_device, subkmers_device, k, kmers_size, ndptr_b->Len());
        cudaMemcpy(freq_b, freq_b_device, kmers_size * sizeof(float), cudaMemcpyDeviceToHost);

        VectorHalfAdd<<<block_num, thread_num>>>(freq_a_device, freq_b_device, freq_m_device, kmers_size);
        cudaMemcpy(freq_m, freq_m_device, kmers_size * sizeof(float), cudaMemcpyDeviceToHost);

        for(int i = 0; i < kmers_size; i ++){
            if(freq_a[i] > 10e-8 && freq_m[i] > 10e-8)
                entropy_a += freq_a[i] * log(freq_m[i] / freq_a[i]);
            if(freq_b[i] > 10e-8 && freq_m[i] > 10e-8)
                entropy_b += freq_b[i] * log(freq_m[i] / freq_b[i]);
        }

        result = -0.5 * (entropy_a + entropy_b);

        free(freq_m); free(freq_a); free(freq_b);
        cudaFree(seq_a_device);
        cudaFree(seq_b_device);
        cudaFree(subseq_device);
        cudaFree(kmers_device);
        cudaFree(subkmers_device);
        cudaFree(freq_a_device);
        cudaFree(freq_b_device);
        cudaFree(freq_m_device);

        return result;
    }

    void FindMinPosition(double* mat_elems, int row_num, int col_num, int pos[]){
        int *min_idxs_host, *min_idxs_device;
        double *mins_host, *mins_device;
        double *elems_device;
        double min = 10e8;
        int min_ridx;

        int thread_num = 512, block_num = (row_num / thread_num) + 1;

        min_idxs_host = (int*)malloc(row_num * sizeof(int));
        mins_host = (double*)malloc(row_num * sizeof(double));
        
        cudaMalloc(&min_idxs_device, row_num * sizeof(int));
        cudaMalloc(&mins_device, row_num * sizeof(double));
        cudaMalloc(&elems_device, row_num * col_num * sizeof(double));

        cudaMemcpy(elems_device, mat_elems, row_num * col_num * sizeof(double), cudaMemcpyHostToDevice);

        MinPos<<<block_num, thread_num>>>(elems_device, row_num, col_num, min_idxs_device, mins_device);

        cudaMemcpy(min_idxs_host, min_idxs_device, row_num * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(mins_host, mins_device, row_num * sizeof(double), cudaMemcpyDeviceToHost);

        for(int i = 0; i < row_num; i ++){
            if(mins_host[i] < min){
                min = mins_host[i];
                min_ridx = i;
            }
        }

        pos[0] = min_ridx;
        pos[1] = min_idxs_host[min_ridx];
    }


};