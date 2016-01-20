#ifndef SIMINT_GUARD_GENERATOR__COMMANDLINE_HPP_
#define SIMINT_GUARD_GENERATOR__COMMANDLINE_HPP_

#include <vector>
#include "generator/Options.hpp"

std::string GetNextArg(int & i, int argc, char ** argv);

int GetIArg(int & i, int argc, char ** argv);

std::string GetNextArg(size_t & i, const std::vector<std::string> & opt);

int GetIArg(size_t & i, const std::vector<std::string> & opt);

std::vector<std::string> ParseCommonOptions(OptionMap & options, int argc, char ** argv);


// A helper macro
#define CMDLINE_ASSERT(val, desc) if( !(val) ) { std::cout << "\n" << (desc) << "\n\n"; return 1; }


#endif
