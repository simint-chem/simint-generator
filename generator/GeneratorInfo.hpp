#ifndef GENERATORINFO_HPP
#define GENERATORINFO_HPP

#include "generator/Classes.hpp"

enum class Options
{
    StackMem,
    InlineVRR,
    InlineET,
    InlineHRR,
    Scalar,
    NoSingleET,
    NoET
};


enum class Compiler
{
    Intel,
    GCC;
};

typedef std::map<Options, int> OptionsMap;




class GeneratorInfo
{
public:    

    QAM FinalAM(void) const
    {
        return finalam_;
    }

    int L(void) const
    {
        return finalam_[0] + finalam_[1] + finalam_[2] + finalam_[3];
    }

    bool IsFinalAM(const QAM & am) const
    {
        return am == finalam_;
    }

    
    bool HasBraVRR(void) const
    {
        return ((finalam_[0] + finalam_[1] + finalam_[2] + finalam_[3]) > 0);
    }

    bool HasKetVRR(void) const
    {
        return false; // for now
    }

    bool HasVRR(void) const
    {
        return ((finalam_[0] + finalam_[1] + finalam_[2] + finalam_[3]) > 0);
    }

    bool HasVRR_I(void) const
    {
        return HasVRR() && finalam_[0] >= finalam_[1];
    }

    bool HasVRR_J(void) const
    {
        return HasVRR() && finalam_[1] > finalam_[0];
    }

    bool HasVRR_K(void) const
    {
        return HasVRR() && finalam_[2] >= finalam_[3];
    }

    bool HasVRR_L(void) const
    {
        return HasVRR() && finalam_[3] > finalam_[2];
    }

    bool HasET(void) const
    {
        return (!GetOption(Options::NoET)) && (finalam_[2]+finalam_[3] > 0);
    }

    bool HasHRR(void) const
    {
        return (HasBraHRR() || HasKetHRR());
    }

    bool HasBraHRR(void) const
    {
        return ( (finalam_[0] > 0) && (finalam_[1] > 0) );
    }

    bool HasBraHRR_I(void) const // going from J -> I
    {
        return HasBraHRR() && (finalam_[1] > finalam_[0]);
    }

    bool HasBraHRR_J(void) const // going from I -> J
    {
        return HasBraHRR() && (finalam_[0] >= finalam_[1]);
    }

    bool HasKetHRR(void) const
    {
        return ( (finalam_[2] > 0) && (finalam_[3] > 0) );
    }

    bool HasKetHRR_K(void) const // going from L -> K
    {
        return HasKetHRR() && (finalam_[3] > finalam_[2]);
    }

    bool HasKetHRR_L(void) const // going from K -> L 
    {
        return HasKetHRR() && (finalam_[2] >= finalam_[3]);
    }

    bool HasCPUFlag(const std::string & flag) const
    {
        return cpuflags_.count(flag);
    } 

    bool Scalar(void) const
    {
        return scalar_;
    }

    bool Vectorized(void) const
    {
        return !scalar_;
    }

    bool HasFMA(void) const
    {
        return HasCPUFlag("fma");
    }

    const VectorInfo & GetVectorInfo(void) const
    {
        return *(vector_);
    }



private:
    QAM finalam_;

    Compiler comp_;
    std::set<std::string> cpuflags_;
    OptionsMap options_;
    bool scalar_;

    std::unique_ptr<VectorInfo> vector_;

    static std::set<std::string> ConvertCPUFlags(std::string cpuflagsstr);

};



#endif
