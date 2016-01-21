#ifndef SIMINT_GUARD_GENERATOR__PRINTING_HPP_
#define SIMINT_GUARD_GENERATOR__PRINTING_HPP_

#include "generator/Types.hpp"

static const std::string indent1(4, ' ');
static const std::string indent2(8, ' ');
static const std::string indent3(12, ' ');
static const std::string indent4(16, ' ');
static const std::string indent5(20, ' ');
static const std::string indent6(24, ' ');
static const std::string indent7(28, ' ');
static const std::string indent8(32, ' ');

void PrintDoubletSet(const DoubletSet & d, const std::string & title);
void PrintQuartetSet(const QuartetSet & q, const std::string & title);
void PrintGaussianSet(const GaussianSet & g, const std::string & title);

#endif
