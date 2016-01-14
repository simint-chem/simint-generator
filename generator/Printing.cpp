#include <iostream>
#include "generator/Printing.hpp"


void PrintQuartetSet(const QuartetSet & q, const std::string & title)
{
    std::cout << title << ": " << q.size() << "\n";
    for(const auto & it : q)
        std::cout << "    " << it << "\n";
    std::cout << "\n";
}


void PrintDoubletSet(const DoubletSet & d, const std::string & title)
{
    std::cout << title << ": " << d.size() << "\n";
    for(const auto & it : d)
        std::cout << "    " << it << "\n";
    std::cout << "\n";
}

void PrintGaussianSet(const GaussianSet & g, const std::string & title)
{
    std::cout << title << ": " << g.size() << "\n";
    for(const auto & it : g)
        std::cout << "    " << it << "\n";
    std::cout << "\n";
}

