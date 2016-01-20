#include "generator/Boys.hpp"
#include "generator/ERIGeneratorInfo.hpp"


////////////////
// Base class
////////////////
BoysGen::BoysGen(const ERIGeneratorInfo & info)
    : info_(info), vinfo_(info.GetVectorInfo())
{ }
        

ConstantMap BoysGen::GetConstants(void) const
{
    return ConstantMap();
}
