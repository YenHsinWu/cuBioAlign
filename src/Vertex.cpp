#include "Vertex.h"

namespace BioAlign{
    //////////////////
    // Vertex class //
    //////////////////

    Vertex::Vertex():m_weight(0){
        m_name.clear();
        m_left = nullptr;
        m_right = nullptr;
    }

    Vertex::Vertex(std::string name, double w):m_weight(w){
        m_name = name;
        m_left = nullptr;
        m_right = nullptr;
    }

    Vertex::Vertex(Vertex *left, Vertex *right, double w):m_weight(w){
        m_name = '(' + left->m_name + ',' + right->m_name + ')';
        m_left = left;
        m_right = right;
    }

    Vertex::Vertex(const Vertex& vtx){
        m_name = vtx.m_name;
        m_weight = vtx.m_weight;
        m_left = vtx.m_left;
        m_right = vtx.m_right;
    }

    Vertex::~Vertex(){
        m_name.clear();
        m_weight = 0;
        m_left = nullptr;
        m_right = nullptr;
    }

    std::ostream& operator<< (std::ostream& os, const Vertex& vtx){
        os << "Vertex name : " << vtx.m_name << ", Weight : " << vtx.m_weight << std::endl;
        return os;
    }
};