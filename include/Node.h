#ifndef __CUBIOALIGN_H__
#define __CUBIOALIGN_H__

#include <stdlib.h>
#include <cstring>

#include <iostream>
#include <fstream>
#include <string>

namespace BioAlign{
    ////////////////
    // Node class //
    ////////////////

    class Node{
        private:
            int m_size;
            char *m_sequence;

        public:
            Node();
            Node(std::string);
            Node(const char*, int);
            Node(const Node&);
            ~Node();

            const char* Sequence() const;
            int Len() const;
            char At(int) const;

            friend std::ostream& operator<< (std::ostream&, const Node&);
            bool operator== (const Node&);
            bool operator!= (const Node&);
            void operator= (const Node&);
            char& operator[] (int);
            const char& operator[] (int) const;
    };
};

#endif