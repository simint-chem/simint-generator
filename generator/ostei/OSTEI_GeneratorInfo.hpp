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
    void SetContQ(const QAMSet & q)
    {
        contq_ = q;
        cont_nelements_ = 0;
        for(const auto & it : contq_)
        {
            if(it != FinalAM())
                cont_nelements_ += NCART(it.qam);
        }

        cont_memory_ = cont_nelements_ * sizeof(double);
    }

    QAMSet GetContQ(void) const
    {
        return contq_;
    }

    size_t ContMemoryReq(void) const
    {
        return cont_memory_;
    }

    size_t ContNElements(void) const
    {
        return cont_nelements_;
    }

    size_t PrimNElements(void) const
    {
        if(PrimUseHeap())
            return prim_nelements_;
        else
            return 0;
    }
    
    bool IsContQ(const QAM & am) const
    {
        return contq_.count(am);
    }

    void SetPrimNElements(size_t nelements)
    {
        prim_nelements_ = nelements;
    }

    bool PrimUseStack(void) const
    {
        // I have found heap to be slightly faster, possibly
        // due to better cache behavior. Feel free to test
        // in the future
        return true;
    }

    bool PrimUseHeap(void) const
    {
        return !(PrimUseStack());
    }

    void SetDeriv1_MissingCenter(int center)
    {
        deriv1_missing_center_ = center;
    }

    int Deriv1_MissingCenter(void) const
    {
        return deriv1_missing_center_;
    }

    bool IsUnique(void) const
    {
        QAM am = FinalAM();

        if(am[0] < am[1])
            return false;
        if(am[2] < am[3])
            return false;
        if( (am[0] + am[1]) < (am[2] + am[3]) )
            return false;
        if( (am[0] + am[1]) == (am[2] + am[3]) && (am[0] < am[2]) ) 
            return false;
        return true;
    }

private:
    //////////////////////////////////
    // Memory for contracted quartets
    //////////////////////////////////
    QAMSet contq_;
    size_t cont_memory_;
    size_t cont_nelements_;
    size_t prim_nelements_;

    int deriv1_missing_center_;

};




