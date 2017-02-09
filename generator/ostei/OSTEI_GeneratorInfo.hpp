#pragma once


#include "generator/Types.hpp"
#include "generator/GeneratorInfoBase.hpp"

class OSTEI_GeneratorInfo : public GeneratorInfoBase
{
public:    

    using GeneratorInfoBase::GeneratorInfoBase;

    bool UseStack(void) const
    {
        // I have found heap to be slightly faster, possibly
        // due to better cache behavior. Feel free to test
        // in the future
        return false;
    }

    bool UseHeap(void) const
    {
        return !(UseStack());
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
    QAMSet bcontq_;
    QAMSet primq_;
    QAMSet contq_;

    size_t bcont_nelements_;
    size_t cont_nelements_;
    size_t prim_nelements_;

    int deriv1_missing_center_;

};




