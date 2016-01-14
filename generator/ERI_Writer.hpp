#ifndef ERIWRITER_HPP
#define ERIWRITER_HPP

#include <ostream>
#include <array>

class BoysGen;
class VRR_Writer;
class ET_Writer;
class HRR_Writer;
class ERIGeneratorInfo;

typedef std::array<int, 4> QAM;

void WriteFile(std::ostream & os,
               ERIGeneratorInfo & info,
               const BoysGen & bg,
               const VRR_Writer & vrr_writer,
               const ET_Writer & et_writer,
               const HRR_Writer & hrr_writer);

void WriteFile_Permute(std::ostream & os, ERIGeneratorInfo & info);

#endif
