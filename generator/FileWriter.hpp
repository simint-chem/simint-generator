#ifndef WRITER_HPP
#define WRITER_HPP

#include <map>
#include <string>
#include <ostream>

#include "generator/Options.hpp"

class BoysGen;
class VRR_Writer;
class ET_Writer;
class HRR_Writer;

typedef std::array<int, 4> QAM;

void WriteFile(std::ostream & os,
               const QAM & am,
               const BoysGen & bg,
               const VRR_Writer & vrr_writer,
               const ET_Writer & et_writer,
               const HRR_Writer & hrr_writer);

void WriteFile_Permute(std::ostream & os);

#endif
