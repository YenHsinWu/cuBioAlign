#ifndef __DISTMATRIX_H__
#define __DISTMATRIX_H__

#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include "Node.h"
#include "Tools.h"

namespace BioAlign{
    //////////////////////
    // DistMatrix class //
    //////////////////////

    class DistMatrix{
        private:
            int m_row, m_col;
            double *m_dists;

        public:
            DistMatrix();
            DistMatrix(int, int);
            DistMatrix(const DistMatrix&);
            ~DistMatrix();

            int GetRowNum() const;
            int GetColNum() const;

            double* GetElements() const;
            double GetElement(int, int) const;
            double* GetRow(int) const;
            double* GetCol(int) const;
            void SetElement(int, int, double);

            void ACS(const std::vector<Node*>&);

            void WriteToFile(std::string);

            friend std::ostream& operator<< (std::ostream&, const DistMatrix&);
    };
};

#endif