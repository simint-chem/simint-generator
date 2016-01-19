#include <algorithm>

#include "generator/ERIGeneratorInfo.hpp"

static char TransformPred(char c)
{
    if(c == ',')
        return ' ';
    else
        return (char)(tolower(c));
}


std::set<std::string> ERIGeneratorInfo::ConvertCPUFlags(std::string cpuflagsstr)
{
    std::set<std::string> cpuflags;

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


ERIGeneratorInfo::ERIGeneratorInfo(QAM finalam,
                                   Compiler compiler,
                                   const std::string & cpuflagsstr,
                                   const OptionMap & options)
    : finalam_(finalam),
      compiler_(compiler),
      cpuflags_(ConvertCPUFlags(cpuflagsstr)),
      options_(options),
      scalar_(false)
{
    AddInclude("\"vectorization/vectorization.h\"");

    if(HasCPUFlag("avx"))
    {
        vector_ = std::unique_ptr<VectorInfo>(new BasicIntelSIMDVector(256));
        AddInclude("\"vectorization/intrinsics_avx.h\"");
    }
    else if(HasCPUFlag("sse2"))
    {
        vector_ = std::unique_ptr<VectorInfo>(new BasicIntelSIMDVector(128));
        AddInclude("\"vectorization/intrinsics_sse.h\"");
    }
    else
    {
        scalar_ = true;
        vector_ = std::unique_ptr<VectorInfo>(new ScalarVector);
    }
}

