#ifndef WRITER_HPP
#define WRITER_HPP

#include <map>
#include <string>
#include <ostream>

class BoysGen;
class VRR_Algorithm_Base;
class ET_Algorithm_Base;
class HRR_Algorithm_Base;

typedef std::array<int, 4> QAMList;


// options
typedef std::map<int, int> OptionsMap;
#define OPTION_FLATPRIM 1


void Writer_Looped(std::ostream & os,
                   const QAMList & am,
                   const std::string & nameappend,
                   const OptionsMap & options,
                   const BoysGen & bg,
                   VRR_Algorithm_Base & vrralgo,
                   ET_Algorithm_Base & etalgo,
                   HRR_Algorithm_Base & hrralgo);

#endif
