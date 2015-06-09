#ifndef HRRWRITER_HPP
#define HRRWRITER_HPP

#include <iostream>

#include "generator/Classes.hpp"

class WriterBase;
class HRR_Algorithm_Base;


class HRR_Writer
{   
    public:
        HRR_Writer(const HRR_Algorithm_Base & hrr_algo);

        void WriteHRR(std::ostream & os, const WriterBase & base) const;

        void WriteIncludes(std::ostream & os, const WriterBase & base) const;
        void WriteHRRFile(std::ostream & ofb, std::ostream & ofk, const WriterBase & base) const;
        void WriteHRRHeaderFile(std::ostream & os, const WriterBase & base) const;

    private:
        HRRBraKetStepList hrrsteps_;
        std::pair<DoubletSet, DoubletSet>  brakettop_;
        std::pair<DAMSet, DAMSet> brakettopam_;
        QAMSet topquartetam_;
        

        void WriteHRRInline_(std::ostream & os, const WriterBase & base) const;
        void WriteHRRExternal_(std::ostream & os, const WriterBase & base) const;

        void WriteBraSteps_(std::ostream & os, const WriterBase & base, const std::string & ncart_ket, const std::string & ketstr) const;
        void WriteKetSteps_(std::ostream & os, const WriterBase & base, const std::string & ncart_ket, const std::string & brastr) const;

        std::string HRRBraStepVar_(const Doublet & d, const std::string & ncart_ket, const std::string & ketstr, bool istarget, const WriterBase & base) const;
        std::string HRRKetStepVar_(const Doublet & d, const std::string & ncart_bra, const std::string & brastr, bool istarget, const WriterBase & base) const;
};

#endif
