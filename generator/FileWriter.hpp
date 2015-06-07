#ifndef WRITER_HPP
#define WRITER_HPP

#include <map>
#include <string>
#include <ostream>

#include "generator/Options.hpp"

class BoysGen;
class VRR_Algorithm_Base;
class ET_Algorithm_Base;
class HRR_Algorithm_Base;

typedef std::array<int, 4> QAM;

void WriteFile(std::ostream & os,
               const QAM & am,
               const std::string & prefix,
               const OptionsMap & options,
               const BoysGen & bg,
               VRR_Algorithm_Base & vrralgo,
               ET_Algorithm_Base & etalgo,
               HRR_Algorithm_Base & hrralgo);

#endif
