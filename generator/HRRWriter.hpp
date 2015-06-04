#ifndef HRRWRITER_HPP
#define HRRWRITER_HPP

#include <iostream>

#include "generator/Classes.hpp"

class WriterBase;


class HRRWriter
{   
    public:
        HRRWriter(const HRRBraKetStepList & hrrsteps, const QAMList & finalam);

        void WriteHRR(std::ostream & os, const WriterBase & base) const;

        void WriteHRRFile(std::ostream & os, const WriterBase & base) const;
        void WriteHRRHeaderFile(std::ostream & os, const WriterBase & base) const;

        bool HasHRR(void) const;
        bool HasBraHRR(void) const;
        bool HasKetHRR(void) const;

        const DoubletSetMap & TopBras(void) const;
        const DoubletSetMap & TopKets(void) const;
        const QuartetSet & TopQuartets(void) const;

    private:
        HRRBraKetStepList hrrsteps_;
        DoubletSetMap hrrtopbras_, hrrtopkets_;
        QuartetSet hrrtopquartets_;

        DAMListSet brahrr_ptrs_;
        
        void WriteHRRInline_(std::ostream & os, const WriterBase & base) const;
        void WriteHRRExternal_(std::ostream & os, const WriterBase & base) const;

        void WriteHRRBraSteps_(std::ostream & os, const WriterBase & base, const std::string & ncart_ket) const;

        std::string HRRBraStepArrVar_(const Doublet & d, const std::string & ncart_ket, bool istarget, const WriterBase & base) const;
        std::string HRRKetStepArrVar_(const Doublet & d, const DAMList & braam, bool istarget, const WriterBase & base) const;
};

#endif
