#ifndef GENERATORINFO_HPP
#define GENERATORINFO_HPP

#include <memory>
#include "generator/Classes.hpp"
#include "generator/VectorInfo.hpp"
#include "generator/Ncart.hpp"
#include "generator/Options.hpp"




class ERIGeneratorInfo
{
public:    
    ERIGeneratorInfo(QAM finalam,
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

/*    
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
        return (!GetOption(Option::NoET)) && (finalam_[2]+finalam_[3] > 0);
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
*/
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


    //////////////////////
    // Constants
    //////////////////////
    void AddNamedConstant(const std::string & name, const std::string & val)
    {
        constants_.emplace(name, val);
    }

    void AddIntegerConstant(int i)
    {
        AddNamedConstant(StringBuilder("const_", i), std::to_string(i));
    }

    std::map<std::string, std::string> GetConstants(void) const
    {
        return constants_;
    }


    
    //////////////////////////////////
    // Memory for contracted quartets
    //////////////////////////////////
    void SetContQ(const QAMSet & q)
    {
        contq_ = q;
        nelements_ = 0;
        for(const auto & it : contq_)
        {
            if(it != finalam_)
                nelements_ += NCART(it);
        }

        memory_ = nelements_ * sizeof(double);
    }

    QAMSet GetContQ(void) const
    {
        return contq_;
    }

    size_t ContMemoryReq(void) const
    {
        return memory_;
    }
    
    bool IsContQ(const QAM & am) const
    {
        return contq_.count(am);
    }

    bool UseStack(void) const
    {
        return memory_ <= static_cast<size_t>(GetOption(Option::StackMem));
    }

    bool UseHeap(void) const
    {
        return !UseStack();
    }


    ////////////////////
    // Include files
    ////////////////////
    void AddInclude(const std::string & inc)
    {
        includes_.insert(inc);
    }

    void AddIncludes(const std::set<std::string> & incs)
    {
        for(const auto & it : incs)
            includes_.insert(it);
    }

    std::set<std::string> GetIncludes(void) const
    {
        return includes_;
    }


private:
    QAM finalam_;

    Compiler comp_;
    std::set<std::string> cpuflags_;
    OptionMap options_;
    bool scalar_;

    std::unique_ptr<VectorInfo> vector_;

    static std::set<std::string> ConvertCPUFlags(std::string cpuflagsstr);

    ////////////////////
    // Include files
    ////////////////////
    std::set<std::string> includes_;

    /////////////
    // Constants
    /////////////
    std::map<std::string, std::string> constants_;

    //////////////////////////////////
    // Memory for contracted quartets
    //////////////////////////////////
    QAMSet contq_;
    size_t memory_;
    size_t nelements_;

};



#endif
