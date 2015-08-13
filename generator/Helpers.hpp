#ifndef HELPERS_HPP
#define HELPERS_HPP

#include "generator/Classes.hpp"
#include "generator/Options.hpp"


static const std::string indent1(4, ' ');
static const std::string indent2(8, ' ');
static const std::string indent3(12, ' ');
static const std::string indent4(16, ' ');
static const std::string indent5(20, ' ');
static const std::string indent6(24, ' ');
static const std::string indent7(28, ' ');
static const std::string indent8(32, ' ');

std::string GetNextArg(int & i, int argc, char ** argv);
int GetIArg(int & i, int argc, char ** argv);

// versions std::vector<std::string>
std::string GetNextArg(size_t & i, const std::vector<std::string> & opt);
int GetIArg(size_t & i, const std::vector<std::string> & opt); 

QuartetSet GenerateInitialQuartetTargets(QAM amlst);
DoubletSet GenerateInitialDoubletTargets(DAM amlst, DoubletType type);

void PrintDoubletSet(const DoubletSet & d, const std::string & title);
void PrintQuartetSet(const QuartetSet & q, const std::string & title);
void PrintGaussianSet(const GaussianSet & g, const std::string & title);

int GaussianOrder(const QAM & ijk);


GaussianSet AllGaussiansForAM(int am);

OptionsMap DefaultOptions(void);
std::vector<std::string> ParseCommonOptions(OptionsMap & options, int argc, char ** argv);


#endif
