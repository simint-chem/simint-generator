#pragma once


#include <iostream>

#include "generator/Types.hpp"
#include "generator/ostei/OSTEI_HRR_Algorithm_Base.hpp"

// foward declare
class OSTEI_GeneratorInfo;



class OSTEI_HRR_Writer
{   
    public:
        OSTEI_HRR_Writer(const OSTEI_HRR_Algorithm_Base & hrr_algo, const OSTEI_GeneratorInfo & info,
                         int start_external = 0, int start_general = 0);

        bool HasHRR(void) const;
        bool HasBraHRR(void) const;
        bool HasKetHRR(void) const;


        virtual ConstantMap GetConstants(void) const;
        virtual void WriteHRR(std::ostream & os) const;
        virtual void WriteHRRFile(std::ostream & of, std::ostream & ofh) const;


    protected:
        const OSTEI_HRR_Algorithm_Base & hrr_algo_; 
        const OSTEI_GeneratorInfo & info_;

        int start_external_;
        int start_general_;

        void WriteBraSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                            const std::string & ncart_ket, const std::string & ketstr) const;
        void WriteKetSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                            const std::string & ncart_ket, const std::string & brastr) const;

        void WriteHRR_Bra_Inline_(std::ostream & os, QAM am) const;
        void WriteHRR_Ket_Inline_(std::ostream & os, QAM am) const;
        void WriteHRR_Bra_External_(std::ostream & os, QAM am) const;
        void WriteHRR_Ket_External_(std::ostream & os, QAM am) const;
        void WriteHRR_Bra_General_(std::ostream & os, QAM am) const;
        void WriteHRR_Ket_General_(std::ostream & os, QAM am) const;

        std::string HRRBraStepVar_(const Doublet & d, const std::string & ncart_ket, const std::string & ketstr) const;
        std::string HRRKetStepVar_(const Doublet & d, const std::string & brastr) const;
};

