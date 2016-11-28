#pragma once


#include "generator/Types.hpp"
#include "generator/GeneratorInfoBase.hpp"
#include "generator/Options.hpp"




class OSTEI_GeneratorInfo : public GeneratorInfoBase
{
public:    

    using GeneratorInfoBase::GeneratorInfoBase;

    //////////////////////////////////
    // Memory for contracted quartets
    //////////////////////////////////
    void SetContQ(const TaggedQAMSet & q)
    {
        contq_ = q;
        nelements_ = 0;
        for(const auto & it : contq_)
        {
            if(it.qam != FinalAM())
                nelements_ += NCART(it.qam);
        }

        memory_ = nelements_ * sizeof(double);
    }

    void SetContQ(const QAMSet & q)
    {
        TaggedQAMSet q2;
        for(const auto & it : q)
            q2.insert({it, ""});

        SetContQ(q2);
    }

    TaggedQAMSet GetContQ(void) const
    {
        return contq_;
    }

    size_t ContMemoryReq(void) const
    {
        return memory_;
    }

    size_t ContNElements(void) const
    {
        return nelements_;
    }
    
    bool IsContQ(const QAM & am) const
    {
        for(const auto & it : contq_)
            if(it.qam == am)
                return true;
        return false;
    }

    bool IsContQ(const TaggedQAM & am) const
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
    TaggedQAMSet contq_;
    size_t memory_;
    size_t nelements_;

};




