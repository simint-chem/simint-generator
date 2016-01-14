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

void WriteERIFile(std::ostream & os,
               const QAM & am,
               const BoysGen & bg,
               const VRR_Writer & vrr_writer,
               const ET_Writer & et_writer,
               const HRR_Writer & hrr_writer,
               const ERIGeneratorInfo & info);

void WriteERIFile_Permute(std::ostream & os);

#endif
