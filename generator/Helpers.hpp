#ifndef HELPERS_HPP
#define HELPERS_HPP

#include "generator/Classes.hpp"

QuartetSet GenerateInitialQuartetTargets(std::array<int, 4> amlst, bool initial);
DoubletSet GenerateInitialDoubletTargets(std::array<int, 2> amlst, DoubletType type, bool initial);

void PruneRight(QuartetSet & qs, DoubletType type);

void PrintQuartetSet(const QuartetSet & q, const std::string & title);

void PrintQuartetSet_Arr(const QuartetSet & q, const std::string & title);

int GaussianOrder(const std::array<int, 4> & ijk);

#endif
