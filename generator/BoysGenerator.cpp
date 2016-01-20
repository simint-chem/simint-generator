#include "generator/BoysGenerator.hpp"
#include "generator/ERIGeneratorInfo.hpp"


////////////////
// Base class
////////////////
BoysGenerator::BoysGenerator(const ERIGeneratorInfo & info)
    : info_(info), vinfo_(info.GetVectorInfo())
{ }
        

ConstantMap BoysGenerator::GetConstants(void) const
{
    return ConstantMap();
}
