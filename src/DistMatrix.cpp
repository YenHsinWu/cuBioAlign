#include "DistMatrix.h"

namespace BioAlign{
    //////////////////////
    // DistMatrix class //
    //////////////////////

    DistMatrix::DistMatrix():m_row(0), m_col(0){
        m_dists = nullptr;
    }

    DistMatrix::DistMatrix(int row, int col):m_row(row), m_col(col){
        m_dists = (double*)calloc(m_row * m_col, sizeof(double));
    }

    DistMatrix::DistMatrix(const DistMatrix& mat){
        m_row = mat.m_row;
        m_col = mat.m_col;

        free(m_dists);
        m_dists = (double*)malloc(m_row * m_col * sizeof(double));
        memcpy(m_dists, mat.m_dists, m_row * m_col);
    }

    DistMatrix::~DistMatrix(){
        free(m_dists);

        m_row = 0;
        m_col = 0;
        m_dists = nullptr;
    }

    std::ostream& operator<< (std::ostream& os, const DistMatrix& mat){
        for(int i = 0; i < mat.m_row; i ++){
            for(int j = 0; j < mat.m_col; j ++)
                os << mat.m_dists[i * mat.m_row + j] << ' ';
            os << std::endl;
        }
        os << "Rows : " << mat.m_row << ", Cols : " << mat.m_col << std::endl;

        return os;
    }
};