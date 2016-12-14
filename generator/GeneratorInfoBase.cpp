/*! \file
 *
 * \brief Holds information about the requested generation (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#include <algorithm>

#include "generator/GeneratorInfoBase.hpp"


/*! \brief Returns the lowercase verion of c and replaces commas with spaces
 */
static char TransformPred(char c)
{
    if(c == ',')
        return ' ';
    else
        return (char)(tolower(c));
}


IncludeSet GeneratorInfoBase::ConvertCPUFlags(std::string cpuflagsstr)
{
    StringSet cpuflags;

    // separate by spaces rather than commas
    // and convert to lower case
    std::transform(cpuflagsstr.begin(), cpuflagsstr.end(), cpuflagsstr.begin(), TransformPred);

    std::stringstream ss(cpuflagsstr);

    while(ss)
    {
        std::string s;
        ss >> s;

        if(s.size())
            cpuflags.insert(s);
    }

    return cpuflags;
}


GeneratorInfoBase::GeneratorInfoBase(QAM finalam,
                                     int deriv,
                                     Compiler compiler,
                                     const std::string & cpuflagsstr,
                                     const OptionMap & options)
    : finalam_(finalam),
      deriv_(deriv),
      compiler_(compiler),
      cpuflags_(ConvertCPUFlags(cpuflagsstr)),
      options_(options),
      scalar_(false)
{
    if(HasCPUFlag("avx512"))
    {
        bool has_fma = HasCPUFlag("fma");
        vector_ = std::unique_ptr<VectorInfo>(new BasicIntelSIMDVector(512, has_fma, false));
    }
    else if(HasCPUFlag("avx"))
    {
        bool has_fma = HasCPUFlag("fma");
        vector_ = std::unique_ptr<VectorInfo>(new BasicIntelSIMDVector(256, has_fma, true));
    }
    else if(HasCPUFlag("sse2"))
        vector_ = std::unique_ptr<VectorInfo>(new BasicIntelSIMDVector(128, false, true));
    else
    {
        scalar_ = true;
        vector_ = std::unique_ptr<VectorInfo>(new ScalarVector);
    }
}

