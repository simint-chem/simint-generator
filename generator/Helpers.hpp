#ifndef HELPERS_HPP
#define HELPERS_HPP

#include "generator/Classes.hpp"
#include "generator/Options.hpp"



QuartetSet GenerateInitialQuartetTargets(QAM amlst);
DoubletSet GenerateInitialDoubletTargets(DAM amlst, DoubletType type);
int GaussianOrder(const QAM & ijk);
GaussianSet AllGaussiansForAM(int am);


#endif
