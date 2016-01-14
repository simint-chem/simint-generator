#include <vector>
#include "generator/Options.hpp"

std::string GetNextArg(int & i, int argc, char ** argv);

int GetIArg(int & i, int argc, char ** argv);

std::string GetNextArg(size_t & i, const std::vector<std::string> & opt);

int GetIArg(size_t & i, const std::vector<std::string> & opt);

std::vector<std::string> ParseCommonOptions(OptionMap & options, int argc, char ** argv);

