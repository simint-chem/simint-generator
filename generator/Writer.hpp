#ifndef WRITER_HPP
#define WRITER_HPP

#include <ostream>
#include <array>
#include <string>

#include "generator/Classes.hpp"

void Write_Generic(std::ostream & os,
                   const std::array<int, 4> & am,
                   const std::string & nameappend,
                   const BoysMap & bm,
                   const VRRInfo & vrinfo,
                   const ETInfo & etinfo,
                   const HRRInfo & hrrinfo);

#endif
