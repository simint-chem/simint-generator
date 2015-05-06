#ifndef WRITER_HPP
#define WRITER_HPP

#include "generator/Classes.hpp"

void Writer_Unrolled(std::ostream & os,
                     const QAMList & am,
                     const std::string & nameappend,
                     const BoysMap & bm,
                     const HRRQuartetStepList & hrrsteps);

void Writer_Looped(std::ostream & os,
                   const QAMList & am,
                   const std::string & nameappend,
                   const BoysMap & bm,
                   VRR_Algorithm_Base & vrralgo,
                   const HRRBraKetStepList & hrrsteps);

#endif
