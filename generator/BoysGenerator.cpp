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
