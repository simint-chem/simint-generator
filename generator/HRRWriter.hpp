#ifndef HRRWRITER_HPP
#define HRRWRITER_HPP

#include <iostream>

#include "generator/Classes.hpp"

class WriterBase;


class HRRWriter
{   
    public:
        HRRWriter(const HRRBraKetStepList & hrrsteps, const QAMList & finalam);

        void WriteHRRInline(std::ostream & os, const WriterBase & base) const;

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

        std::string HRRBraStepArrVar(const Doublet & d, int ketam, bool istarget, const WriterBase & base) const;
        std::string HRRKetStepArrVar(const Doublet & d, const DAMList & braam, bool istarget, const WriterBase & base) const;
        std::string HRRBraStepString(const HRRDoubletStep & hrr, int ketam, const WriterBase & base) const;
        std::string HRRKetStepString(const HRRDoubletStep & hrr, const DAMList & braam, const WriterBase & base) const;
};

#endif
