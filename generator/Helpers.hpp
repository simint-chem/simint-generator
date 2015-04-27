#ifndef HELPERS_HPP
#define HELPERS_HPP

#include "generator/Classes.hpp"

QuartetSet GenerateInitialTargets(std::array<int, 4> amlst, bool initial);

void PruneRight(QuartetSet & qs, DoubletType type);

void PrintQuartetSet(const QuartetSet & q, const std::string & title);

void PrintQuartetSet_Arr(const QuartetSet & q, const std::string & title);

int GaussianOrder(const std::array<int, 4> & ijk);

#endif
