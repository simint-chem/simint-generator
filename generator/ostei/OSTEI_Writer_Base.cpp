#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_VRR_Writer.hpp"
#include "generator/ostei/OSTEI_HRR_Writer.hpp"
#include "generator/ostei/OSTEI_Writer.hpp"
#include "generator/Printing.hpp"
#include "generator/Naming.hpp"


///////////////////////////
// OSTEI_Writer Base Class //
///////////////////////////
OSTEI_Writer_Base::OSTEI_Writer_Base(std::ostream & os,
                                     std::ostream & osh,
                                     const OSTEI_GeneratorInfo & info,
                                     const OSTEI_VRR_Writer & vrr_writer,
                                     const OSTEI_HRR_Writer & hrr_writer)
       : os_(os), osh_(osh), info_(info),
         vrr_writer_(vrr_writer), hrr_writer_(hrr_writer)
{ }


void OSTEI_Writer_Base::DeclarePrimPointers(void) const
{
    for(const auto & it : info_.GetContQ())
    {
        os_ << indent4 << "double * restrict " << PrimPtrName(it)
            << " = " << ArrVarName(it) << " + abcd * " << NCART(it) << ";\n";
    }

    os_ << "\n\n";
}


