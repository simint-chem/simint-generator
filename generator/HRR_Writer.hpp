#ifndef HRRWRITER_HPP
#define HRRWRITER_HPP

#include <iostream>

#include "generator/Classes.hpp"
#include "generator/WriterBase.hpp"
#include "generator/HRR_Algorithm_Base.hpp"

// foward declare
class ERIGeneratorInfo;



class HRR_Writer : public WriterBase
{   
    public:
        HRR_Writer(const HRR_Algorithm_Base & hrr_algo);

        bool HasHRR(void) const;
        bool HasBraHRR(void) const;
        bool HasKetHRR(void) const;


        virtual void AddConstants(ERIGeneratorInfo & info) const = 0;
        virtual void WriteHRR(std::ostream & os, const ERIGeneratorInfo & info) const = 0;
        virtual void WriteHRRFile(std::ostream & ofb, std::ostream & ofk, std::ostream & ofh, const ERIGeneratorInfo & info) const = 0;


    protected:
        const HRR_Algorithm_Base & hrr_algo_; 

        void WriteBraSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                            const std::string & ncart_ket, const std::string & ketstr) const;
        void WriteKetSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                            const std::string & ncart_ket, const std::string & brastr) const;

        std::string HRRBraStepVar_(const Doublet & d, const std::string & ncart_ket, const std::string & ketstr) const;
        std::string HRRKetStepVar_(const Doublet & d, const std::string & brastr) const;
};



class HRR_Writer_Inline : public HRR_Writer, public IsInlineRR
{
    public:
        HRR_Writer_Inline(const HRR_Algorithm_Base & hrr_algo);

        virtual void AddConstants(ERIGeneratorInfo & info) const;
        virtual void WriteHRR(std::ostream & os, const ERIGeneratorInfo & info) const;
        virtual void WriteHRRFile(std::ostream & ofb, std::ostream & ofk, std::ostream & ofh, const ERIGeneratorInfo & info) const;

};


class HRR_Writer_External : public HRR_Writer, public IsExternalRR
{
    public:
        HRR_Writer_External(const HRR_Algorithm_Base & hrr_algo);

        virtual void AddConstants(ERIGeneratorInfo & info) const;
        virtual void WriteHRR(std::ostream & os, const ERIGeneratorInfo & info) const;
        virtual void WriteHRRFile(std::ostream & ofb, std::ostream & ofk, std::ostream & ofh, const ERIGeneratorInfo & info) const;

};

#endif
