#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_VRR_Writer.hpp"
#include "generator/ostei/OSTEI_HRR_Writer.hpp"
#include "generator/ostei/OSTEI_Writer.hpp"


///////////////////////////
// OSTEI_Writer Base Class //
///////////////////////////
OSTEI_Writer::OSTEI_Writer(std::ostream & os,
                       std::ostream & osh,
                       const OSTEI_GeneratorInfo & info,
                       const OSTEI_VRR_Writer & vrr_writer,
                       const OSTEI_HRR_Writer & hrr_writer)
   : os_(os), osh_(osh), info_(info), vinfo_(info_.GetVectorInfo()),
     vrr_writer_(vrr_writer), hrr_writer_(hrr_writer)
{ }

