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
        m_dists = (double*)malloc(mat.m_row * mat.m_col * sizeof(double));
        memcpy(m_dists, mat.m_dists, mat.m_row * mat.m_col);
    }

    DistMatrix::~DistMatrix(){
        free(m_dists);

        m_row = 0;
        m_col = 0;
        m_dists = nullptr;
    }

    double DistMatrix::GetElement(int r, int c) const{
        return m_dists[r * m_row + c];
    }

    void DistMatrix::SetElement(int r, int c, double val){
        m_dists[r * m_row + c] = val;
    }

    std::ostream& operator<< (std::ostream& os, const DistMatrix& mat){
        for(int i = -1; i < mat.m_row; i ++){
            for(int j = -1; j < mat.m_col; j ++){
                if(i == -1 && j == -1)
                    os << "      ";
                else if(i == -1 && j != -1)
                    os << std::setw(6) << j << ' ';
                else if(i != -1 && j == -1)
                    os << std::setw(4) << i << "   ";
                else
                    os << std::fixed << std::setw(6) << std::setprecision(2) << mat.m_dists[i * mat.m_row + j] << ' ';
            }
            os << std::endl;
        }
        os << "Rows : " << mat.m_row << ", Cols : " << mat.m_col << std::endl;

        return os;
    }
};