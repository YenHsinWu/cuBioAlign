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
};