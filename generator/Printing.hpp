/*! \file
 *
 * \brief Miscellaneous functions for printing information
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef SIMINT_GUARD_GENERATOR__PRINTING_HPP_
#define SIMINT_GUARD_GENERATOR__PRINTING_HPP_

#include <iostream>
#include "generator/Types.hpp"


static const std::string indent1(4, ' ');
static const std::string indent2(8, ' ');
static const std::string indent3(12, ' ');
static const std::string indent4(16, ' ');
static const std::string indent5(20, ' ');
static const std::string indent6(24, ' ');
static const std::string indent7(28, ' ');
static const std::string indent8(32, ' ');
static const std::string indent9(36, ' ');
static const std::string indent10(40, ' ');


inline void PrintQuartetSet(const QuartetSet & q, const std::string & title)
{
    std::cout << title << ": " << q.size() << "\n";
    for(const auto & it : q)
        std::cout << "    " << it << "\n";
    std::cout << "\n";
}


inline void PrintDoubletSet(const DoubletSet & d, const std::string & title)
{
    std::cout << title << ": " << d.size() << "\n";
    for(const auto & it : d)
        std::cout << "    " << it << "\n";
    std::cout << "\n";
}

inline void PrintGaussianSet(const GaussianSet & g, const std::string & title)
{
    std::cout << title << ": " << g.size() << "\n";
    for(const auto & it : g)
        std::cout << "    " << it << "\n";
    std::cout << "\n";
}

#endif
