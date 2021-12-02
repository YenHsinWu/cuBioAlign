#ifndef __DISTMATRIX_H__
#define __DISTMATRIX_H__

#include <stdlib.h>
#include <string.h>

#include <iostream>

namespace BioAlign{
    ///////////////////
    // DistMatrix class //
    ///////////////////

    class DistMatrix{
        private:
            int m_row, m_col;
            double *m_dists;

        public:
            DistMatrix();
            DistMatrix(int, int);
            DistMatrix(const DistMatrix&);
            ~DistMatrix();

            friend std::ostream& operator<< (std::ostream&, const DistMatrix&);
    };
};

#endif