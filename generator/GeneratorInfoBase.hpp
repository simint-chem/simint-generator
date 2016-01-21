#ifndef SIMINT_GUARD_GENERATOR__GENERATORINFOBASE_HPP_
#define SIMINT_GUARD_GENERATOR__GENERATORINFOBASE_HPP_

#include <memory>
#include "generator/Types.hpp"
#include "generator/VectorInfo.hpp"
#include "generator/Ncart.hpp"
#include "generator/Options.hpp"




class GeneratorInfoBase
{
public:    
    GeneratorInfoBase(QAM finalam,
                      Compiler compiler,
                      const std::string & cpuflagsstr,
                      const OptionMap & options);

    QAM FinalAM(void) const
    {
        return finalam_;
    }

    int GetOption(Option opt) const
    {
        return options_.at(opt);
    }

    int L(void) const
    {
        return finalam_[0] + finalam_[1] + finalam_[2] + finalam_[3];
    }

    bool IsFinalAM(const QAM & am) const
    {
        return am == finalam_;
    }

    bool HasCPUFlag(const std::string & flag) const
    {
        return cpuflags_.count(flag);
    }

    bool HasFMA(void) const
    {
        return HasCPUFlag("fma");
    }

    const VectorInfo & GetVectorInfo(void) const
    {
        return *(vector_);
    }

    bool Scalar(void) const
    {
        return scalar_;
    }

    bool Vectorized(void) const
    {
        return !Scalar();
    }

    IncludeSet GetIncludes(void) const
    {
        return includes_;
    }


private:
    QAM finalam_;

    Compiler compiler_;
    std::set<std::string> cpuflags_;
    OptionMap options_;
    bool scalar_;

    std::unique_ptr<VectorInfo> vector_;

    std::set<std::string> includes_;

    static std::set<std::string> ConvertCPUFlags(std::string cpuflagsstr);
};



#endif
