#ifndef SIMINT_GUARD_ERI__ERIGENERATORINFO_HPP_
#define SIMINT_GUARD_ERI__ERIGENERATORINFO_HPP_

#include <memory>
#include "generator/Types.hpp"
#include "generator/GeneratorInfoBase.hpp"
#include "generator/VectorInfo.hpp"
#include "generator/Ncart.hpp"
#include "generator/Options.hpp"




class ERIGeneratorInfo : public GeneratorInfoBase
{
public:    
    using GeneratorInfoBase::GeneratorInfoBase;



    //////////////////////////////////
    // Memory for contracted quartets
    //////////////////////////////////
    void SetContQ(const QAMSet & q)
    {
        contq_ = q;
        nelements_ = 0;
        for(const auto & it : contq_)
        {
            if(it != FinalAM())
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


private:
    //////////////////////////////////
    // Memory for contracted quartets
    //////////////////////////////////
    QAMSet contq_;
    size_t memory_;
    size_t nelements_;

};



#endif
