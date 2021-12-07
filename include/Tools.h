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
    void GenerateLeaves(int, std::vector<Vertex*>&);

    void FreeNodes(std::vector<Node*>&);
    void FreeVertices(std::vector<Vertex*>&);

    double ACSDistance(Node*, Node*);

    void InitKmers(char*, int);
    double FFPDistance(Node*, Node*, int);

    void FindMinPosition(double*, int, int, int[]);
    void CreateWeightTree(std::vector<Vertex*>&);
};

#endif