#ifndef SIMINT_GUARD_GENERATOR__HRR_WRITER_HPP_
#define SIMINT_GUARD_GENERATOR__HRR_WRITER_HPP_

#include <iostream>

#include "generator/Classes.hpp"
#include "generator/WriterBase.hpp"
#include "generator/ERI_HRR_Algorithm_Base.hpp"

// foward declare
class ERIGeneratorInfo;



class HRR_Writer : public WriterBase
{   
    public:
        HRR_Writer(const HRR_Algorithm_Base & hrr_algo, const ERIGeneratorInfo & info);

        bool HasHRR(void) const;
        bool HasBraHRR(void) const;
        bool HasKetHRR(void) const;


        virtual ConstantMap GetConstants(void) const;
        virtual void WriteHRR(std::ostream & os) const = 0;
        virtual void WriteHRRFile(std::ostream & ofb, std::ostream & ofk, std::ostream & ofh) const = 0;


    protected:
        const HRR_Algorithm_Base & hrr_algo_; 
        const ERIGeneratorInfo & info_;
        const VectorInfo & vinfo_;

        void WriteBraSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                            const std::string & ncart_ket, const std::string & ketstr) const;
        void WriteKetSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                            const std::string & ncart_ket, const std::string & brastr) const;

        std::string HRRBraStepVar_(const Doublet & d, const std::string & ncart_ket, const std::string & ketstr) const;
        std::string HRRKetStepVar_(const Doublet & d, const std::string & brastr) const;
};



class HRR_Writer_Inline : public HRR_Writer
{
    public:
        using HRR_Writer::HRR_Writer;

        bool IsInline(void) const { return true; }
        bool IsExternal(void) const { return false; }

        virtual void WriteHRR(std::ostream & os) const;
        virtual void WriteHRRFile(std::ostream & ofb, std::ostream & ofk, std::ostream & ofh) const;

};


class HRR_Writer_External : public HRR_Writer
{
    public:
        using HRR_Writer::HRR_Writer;

        bool IsInline(void) const { return false; }
        bool IsExternal(void) const { return true; }

        virtual void WriteHRR(std::ostream & os) const;
        virtual void WriteHRRFile(std::ostream & ofb, std::ostream & ofk, std::ostream & ofh) const;

};

#endif
