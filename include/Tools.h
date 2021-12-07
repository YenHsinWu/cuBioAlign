#ifndef __TOOLS_H__
#define __TOOLS_H__

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <malloc.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#include "Node.h"
#include "Vertex.h"
#include "cuTools.cuh"

namespace BioAlign{
    ////////////////////
    // Tool functions //
    ////////////////////

    void LowerNode(Node*);
    void UpperNode(Node*);

    void ReadFasta(std::string, std::vector<Node*>&);

    double ACSDistance(Node*, Node*);

    void InitKmers(char*, int);
    double FFPDistance(Node*, Node*, int);

    double FindMin(double*, int, int);
    void CreateWeightTree(std::vector<Vertex*>&);
};

#endif