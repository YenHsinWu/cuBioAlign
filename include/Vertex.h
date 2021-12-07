#ifndef __VERTEX_H__
#define __VERTEX_H__

#include <iostream>
#include <string>

namespace BioAlign{
    //////////////////
    // Vertex class //
    //////////////////

    class Vertex{
        private:
            std::string m_name;
            double m_weight;
            Vertex *m_left;
            Vertex *m_right;

        public:
            Vertex();
            Vertex(std::string, double);
            Vertex(Vertex*, Vertex*, double);
            Vertex(const Vertex&);
            ~Vertex();

            friend std::ostream& operator<< (std::ostream&, const Vertex&);
    };
};

#endif