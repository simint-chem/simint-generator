#ifndef SIMINT_GUARD_ERI__ERI_HRR_WRITER_HPP_
#define SIMINT_GUARD_ERI__ERI_HRR_WRITER_HPP_

#include <iostream>

#include "generator/Types.hpp"
#include "generator/WriterBase.hpp"
#include "generator/eri/ERI_HRR_Algorithm_Base.hpp"

// foward declare
class ERIGeneratorInfo;



class ERI_HRR_Writer : public WriterBase
{   
    public:
        ERI_HRR_Writer(const ERI_HRR_Algorithm_Base & hrr_algo, const ERIGeneratorInfo & info);

        bool HasHRR(void) const;
        bool HasBraHRR(void) const;
        bool HasKetHRR(void) const;


        virtual ConstantMap GetConstants(void) const;
        virtual void WriteHRR(std::ostream & os) const = 0;
        virtual void WriteHRRFile(std::ostream & of, std::ostream & ofh) const = 0;


    protected:
        const ERI_HRR_Algorithm_Base & hrr_algo_; 
        const ERIGeneratorInfo & info_;
        const VectorInfo & vinfo_;

        void WriteBraSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                            const std::string & ncart_ket, const std::string & ketstr) const;
        void WriteKetSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                            const std::string & ncart_ket, const std::string & brastr) const;

        std::string HRRBraStepVar_(const Doublet & d, const std::string & ncart_ket, const std::string & ketstr) const;
        std::string HRRKetStepVar_(const Doublet & d, const std::string & brastr) const;
};



class ERI_HRR_Writer_Inline : public ERI_HRR_Writer
{
    public:
        using ERI_HRR_Writer::ERI_HRR_Writer;

        bool IsInline(void) const { return true; }
        bool IsExternal(void) const { return false; }

        virtual void WriteHRR(std::ostream & os) const;
        virtual void WriteHRRFile(std::ostream & of, std::ostream & ofh) const;

};


class ERI_HRR_Writer_External : public ERI_HRR_Writer
{
    public:
        using ERI_HRR_Writer::ERI_HRR_Writer;

        bool IsInline(void) const { return false; }
        bool IsExternal(void) const { return true; }

        virtual void WriteHRR(std::ostream & os) const;
        virtual void WriteHRRFile(std::ostream & of, std::ostream & ofh) const;

};

#endif
