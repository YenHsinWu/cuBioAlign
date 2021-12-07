# cuBioAlign
A bioinformatics library based on C++ and CUDA

## Compile
Please compile the codes with the following command :   
`nvcc -o main ./src/Node.cpp ./src/DistMatrix.cpp ./src/Vertex.cpp ./src/Tools.cu ./src/cuTools.cu ./main.cpp -I./include`  
Make sure you have nvidia cuda compiler in your machine.