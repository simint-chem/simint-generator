#ifndef SIMINT_GUARD_GENERATOR__NAMING_HPP_
#define SIMINT_GUARD_GENERATOR__NAMING_HPP_

#include <string>

#include "generator/Classes.hpp"

std::string ArrVarName(const QAM & am, const std::string & prefix = "");
std::string ArrVarName(int am1, int am2, const std::string & ketstr, const std::string & prefix = "");
std::string ArrVarName(const std::string & brastr, int am3, int am4, const std::string & prefix = "");

std::string HRRVarName(const QAM & am);
std::string HRRVarName(int am1, int am2, const std::string & ketstr);
std::string HRRVarName(const std::string & brastr, int am3, int am4);

std::string PrimVarName(const QAM & am);
std::string PrimPtrName(const QAM & am);

#endif
