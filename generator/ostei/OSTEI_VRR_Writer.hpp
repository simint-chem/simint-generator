#pragma once


#include <iostream>

#include "generator/Types.hpp"
#include "generator/ostei/OSTEI_VRR_Algorithm_Base.hpp"

// foward declare
class OSTEI_GeneratorInfo;

class OSTEI_VRR_Writer
{   
    public:
        OSTEI_VRR_Writer(const OSTEI_VRR_Algorithm_Base & vrr_algo, const OSTEI_GeneratorInfo & info,
                         int start_external = 0, int start_general = 0);

        const OSTEI_VRR_Algorithm_Base & Algo(void) const;

        virtual ConstantMap GetConstants(void) const;
        virtual void WriteVRR(std::ostream & os) const;
        virtual void WriteVRRFile(std::ostream & os, std::ostream & osh) const;

    private:
        const OSTEI_VRR_Algorithm_Base & vrr_algo_;
        const OSTEI_GeneratorInfo & info_;

        int start_external_;
        int start_general_;

        virtual void WriteVRR_Inline_(std::ostream & os, QAM am) const;
        virtual void WriteVRR_External_(std::ostream & os, QAM am) const;
        virtual void WriteVRR_General_(std::ostream & os, QAM am) const;
        void WriteVRRSteps_(std::ostream & os, QAM qam, const VRR_StepSet & vs, const std::string & num_n) const;

};

