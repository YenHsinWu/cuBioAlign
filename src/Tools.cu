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

        free(seq_host);
        cudaFree(seq_device);
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

    double ACSDistance(Node *ndptr_a, Node *ndptr_b){
        int n = ndptr_a->Len(), m = ndptr_b->Len(), sum;
        int *lens_ab, *lens_aa, *lens_ba, *lens_bb, *lens_device;
        char *str_a_device, *str_b_device;
        double result, l_ab, l_aa, l_ba, l_bb, d_ab, d_ba;

        int thread_num = 512, block_num;

        lens_ab = (int*)malloc(n * sizeof(int));
        lens_aa = (int*)malloc(n * sizeof(int));
        lens_ba = (int*)malloc(m * sizeof(int));
        lens_bb = (int*)malloc(m * sizeof(int));

        cudaMalloc(&lens_device, n * sizeof(int));
        cudaMalloc(&str_a_device, n * sizeof(char));
        cudaMalloc(&str_b_device, m * sizeof(char));

        cudaMemcpy(str_a_device, ndptr_a->Sequence(), n * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(str_b_device, ndptr_b->Sequence(), m * sizeof(char), cudaMemcpyHostToDevice);

        block_num = (n / thread_num) + 1;
        LongestSubstringLen<<<block_num, thread_num>>>(lens_device, str_a_device, str_b_device, n, m);
        cudaMemcpy(lens_ab, lens_device, n * sizeof(int), cudaMemcpyDeviceToHost);

        block_num = (m / thread_num) + 1;
        LongestSubstringLen<<<block_num, thread_num>>>(lens_device, str_b_device, str_a_device, m, n);
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
        cudaFree(str_a_device);
        cudaFree(str_b_device);

        return result;
    }
};