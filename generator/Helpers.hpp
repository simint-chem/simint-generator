#ifndef HELPERS_HPP
#define HELPERS_HPP

#include "generator/Classes.hpp"
#include "generator/Options.hpp"



std::string GetNextArg(int & i, int argc, char ** argv);
int GetIArg(int & i, int argc, char ** argv);

// versions std::vector<std::string>
std::string GetNextArg(size_t & i, const std::vector<std::string> & opt);
int GetIArg(size_t & i, const std::vector<std::string> & opt); 

QuartetSet GenerateInitialQuartetTargets(QAM amlst);
DoubletSet GenerateInitialDoubletTargets(DAM amlst, DoubletType type);


int GaussianOrder(const QAM & ijk);


GaussianSet AllGaussiansForAM(int am);

OptionMap DefaultOptions(void);
std::vector<std::string> ParseCommonOptions(OptionMap & options, int argc, char ** argv);


#endif
