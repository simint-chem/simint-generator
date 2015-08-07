#ifndef HRRWRITER_HPP
#define HRRWRITER_HPP

#include <iostream>

#include "generator/Classes.hpp"


class HRR_Algorithm_Base;


class HRR_Writer
{   
    public:
        HRR_Writer(const HRR_Algorithm_Base & hrr_algo);

        void WriteHRR(std::ostream & os) const;

        void AddConstants(void) const;
        void WriteHRRFile(std::ostream & ofb, std::ostream & ofk) const;
        void WriteHRRHeaderFile(std::ostream & os) const;

    private:
        const HRR_Algorithm_Base & hrr_algo_; 

        void WriteHRRInline_(std::ostream & os) const;
        void WriteHRRExternal_(std::ostream & os) const;

        void WriteBraSteps_(std::ostream & os, const std::string & ncart_ket, const std::string & ketstr) const;
        void WriteKetSteps_(std::ostream & os, const std::string & ncart_ket, const std::string & brastr) const;

        std::string HRRBraStepVar_(const Doublet & d, const std::string & ncart_ket, const std::string & ketstr, bool istarget) const;
        std::string HRRKetStepVar_(const Doublet & d, const std::string & ncart_bra, const std::string & brastr, bool istarget) const;
};

#endif
