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

    double DistMatrix::GetElement(int r, int c) const{
        return m_dists[r * m_row + c];
    }

    void DistMatrix::SetElement(int r, int c, double val){
        m_dists[r * m_row + c] = val;
    }

    double* DistMatrix::GetRow(int r) const{
        double *row;
        row = (double*)alloca(m_row * sizeof(double));

        for(int c = 0; c < m_col; c ++)
            row[c] = m_dists[r * m_row + c];

        return row;
    }

    double* DistMatrix::GetCol(int c) const{
        double *col;
        col = (double*)alloca(m_col * sizeof(double));

        for(int r = 0; r < m_row; r ++)
            col[r] = m_dists[r * m_row + c];

        return col;
    }

    void DistMatrix::ACS(const std::vector<Node*>& nodes_ptr){
        m_row = nodes_ptr.size();
        m_col = nodes_ptr.size();

        free(m_dists);
        m_dists = (double*)calloc(m_row * m_col,  sizeof(double));

        double dist;

        for(int i = 0; i < m_row; i ++){
            for(int j = i + 1; j < m_col; j ++){
                dist = ACSDistance(nodes_ptr[i], nodes_ptr[j]);
                this->SetElement(i, j, dist);
                this->SetElement(j, i, dist);
            }
        }
    }

    void DistMatrix::WriteToFile(std::string fname){
        std::ofstream output_file;
        output_file.open(fname);

        for(int i = -1; i < m_row; i ++){
            for(int j = -1; j < m_col; j ++){
                if(i == -1 && j == -1)
                    output_file << "      ";
                else if(i == -1 && j != -1)
                    output_file << std::setw(6) << j << ' ';
                else if(i != -1 && j == -1)
                    output_file << std::setw(4) << i << "   ";
                else
                    output_file << std::fixed << std::setw(6) << std::setprecision(2) << m_dists[i * m_row + j] << ' ';
            }
            output_file << std::endl;
        }
        output_file << "Rows : " << m_row << ", Cols : " << m_col << std::endl;

        output_file.close();
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