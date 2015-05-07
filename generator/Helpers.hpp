#ifndef HELPERS_HPP
#define HELPERS_HPP

#include "generator/Classes.hpp"

QuartetSet GenerateInitialQuartetTargets(QAMList amlst, bool initial);
DoubletSet GenerateInitialDoubletTargets(DAMList amlst, DoubletType type, bool initial);

void PruneRight(DoubletSet & ds);
void PruneRight(QuartetSet & qs, DoubletType type);
void PruneET(QuartetSet & qs);

void PrintDoubletSet(const DoubletSet & d, const std::string & title);
void PrintQuartetSet(const QuartetSet & q, const std::string & title);

int GaussianOrder(const QAMList & ijk);


GaussianSet AllGaussiansForAM(int am);


#endif
