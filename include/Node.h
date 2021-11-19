#ifndef __NODE_H__
#define __NODE_H__

#include <stdlib.h>
#include <string.h>

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

            void Update(char*);
            void Update(std::string);

            friend std::ostream& operator<< (std::ostream&, const Node&);
            bool operator== (const Node&);
            bool operator!= (const Node&);
            void operator= (const Node&);
            char& operator[] (int);
            const char& operator[] (int) const;
    };
};

#endif