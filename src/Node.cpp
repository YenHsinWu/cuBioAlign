#include "Node.h"

namespace BioAlign{
    ////////////////
    // Node class //
    ////////////////

    Node::Node():m_size(0){
        m_sequence = nullptr;
    }

    Node::Node(std::string seq){
        m_size = seq.size();

        m_sequence = (char*)malloc(m_size * sizeof(char) + 1);
        std::strncpy(m_sequence, seq.c_str(), m_size);
        m_sequence[m_size] = '\0';
    }

    Node::Node(const char* seq, int sz){
        m_size = sz;

        m_sequence = (char*)malloc(m_size * sizeof(char) + 1);
        std::strncpy(m_sequence, seq, m_size);
        m_sequence[m_size] = '\0';
    }

    Node::Node(const Node& nd){
        m_size = nd.m_size;

        m_sequence = (char*)malloc(m_size * sizeof(char) + 1);
        std::strncpy(m_sequence, nd.m_sequence, m_size);
        m_sequence[m_size] = '\0';
    }

    Node::~Node(){
        free(m_sequence);

        m_size = 0;
        m_sequence = nullptr;
    }

    const char* Node::Sequence() const{
        return m_sequence;
    }

    int Node::Len() const{
        return m_size;
    }

    char Node::At(int idx) const{
        char result = '\0';

        try{
            if(idx < 0 || idx >= m_size)
                throw "In Node::At(int) : Index out of range.";
            result = m_sequence[idx];
        }
        catch(const char* error_message){
            std::cout << error_message << std::endl;
        }

        return result;
    }

    std::ostream& operator<< (std::ostream& os, const Node nd){
        os << "Sequence : " << nd.m_sequence << ", Length = " << nd.m_size;
        return os
    }

    bool Node::operator== (const Node& other){
        if(m_size != other.m_size)
            return false;

        if(std::strncmp(m_sequence, other.m_sequence, m_size) == 0)
            return true;
        return false;
    }

    bool Node::operator!= (const Node& other){
        return !(*this == other);
    }

    void Node::operator= (const Node& other){
        m_size = other.m_size;

        free(m_sequence);
        m_sequence = (char*)malloc(m_size * sizeof(char) + 1);
        std::strncpy(m_sequence, other.m_sequence, m_size);
        m_sequence[m_size] = '\0';
    }

    char& Node::operator[] (int idx){
        try{
            if(idx < 0 || idx > m_size)
                throw "In Node::operator[] (int) : Index out of range.";
        }
        catch(const char* error_message){
            std::cout << error_message << std::endl;
        }

        return m_sequence[idx];
    }

    const char& Node::operator[] (int idx) const{
        try{
            if(idx < 0 || idx > m_size)
                throw "In Node::operator[] (int) : Index out of range.";
        }
        catch(const char* error_message){
            std::cout << error_message << std::endl;
        }

        return m_sequence[idx];
    }
};