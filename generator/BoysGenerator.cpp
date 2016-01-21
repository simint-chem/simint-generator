/*! \file
 *
 * \brief Base class for generating the Boys function
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#include "generator/BoysGenerator.hpp"
#include "generator/GeneratorInfoBase.hpp"


////////////////
// Base class
////////////////
BoysGenerator::BoysGenerator(const GeneratorInfoBase & info)
    : info_(info), vinfo_(info.GetVectorInfo())
{ }
        

ConstantMap BoysGenerator::GetConstants(void) const
{
    return ConstantMap();
}
