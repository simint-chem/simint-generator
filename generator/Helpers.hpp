#ifndef HELPERS_HPP
#define HELPERS_HPP

#include "generator/Classes.hpp"


static const std::string indent1(4, ' ');
static const std::string indent2(8, ' ');
static const std::string indent3(12, ' ');
static const std::string indent4(16, ' ');
static const std::string indent5(20, ' ');
static const std::string indent6(24, ' ');
static const std::string indent7(28, ' ');
static const std::string indent8(32, ' ');


QuartetSet GenerateInitialQuartetTargets(QAM amlst, bool initial);
DoubletSet GenerateInitialDoubletTargets(DAM amlst, DoubletType type, bool initial);

void PrintDoubletSet(const DoubletSet & d, const std::string & title);
void PrintQuartetSet(const QuartetSet & q, const std::string & title);
void PrintGaussianSet(const GaussianSet & g, const std::string & title);
void PrintGaussianMap(const GaussianMap & g, const std::string & title);

int GaussianOrder(const QAM & ijk);


GaussianSet AllGaussiansForAM(int am);


#endif
