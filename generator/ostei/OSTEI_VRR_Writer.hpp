#pragma once


#include <iostream>

#include "generator/Types.hpp"
#include "generator/WriterBase.hpp"
#include "generator/ostei/OSTEI_VRR_Algorithm_Base.hpp"

// foward declare
class OSTEI_GeneratorInfo;



class OSTEI_VRR_Writer : public WriterBase
{   
    public:
        OSTEI_VRR_Writer(const OSTEI_VRR_Algorithm_Base & vrr_algo, const OSTEI_GeneratorInfo & info);

        bool HasVRR(void) const;
        bool HasBraVRR(void) const;
        bool HasKetVRR(void) const;
        bool HasVRR_I(void) const;
        bool HasVRR_J(void) const;
        bool HasVRR_K(void) const;
        bool HasVRR_L(void) const;

        void DeclarePrimArrays(std::ostream & os) const;
        void DeclarePrimPointers(std::ostream & os) const;


        virtual ConstantMap GetConstants(void) const;
        virtual void WriteVRR(std::ostream & os) const = 0;
        virtual void WriteVRRFile(std::ostream & os, std::ostream & osh) const = 0;

    protected:
        const OSTEI_VRR_Algorithm_Base & vrr_algo_;
        const OSTEI_GeneratorInfo & info_;
        const VectorInfo & vinfo_;

        void WriteVRRSteps_(std::ostream & os, QAM qam, const VRR_StepSet & vs, const std::string & num_n) const;
};



class OSTEI_VRR_Writer_Inline : public OSTEI_VRR_Writer
{
    public:
        using OSTEI_VRR_Writer::OSTEI_VRR_Writer;

        bool IsInline(void) const { return true; }
        bool IsExternal(void) const { return false; }

        virtual ConstantMap GetConstants(void) const;
        virtual void WriteVRR(std::ostream & os) const;
        virtual void WriteVRRFile(std::ostream & os, std::ostream & osh) const;
};


class OSTEI_VRR_Writer_External : public OSTEI_VRR_Writer
{
    public:
        using OSTEI_VRR_Writer::OSTEI_VRR_Writer;

        bool IsInline(void) const { return false; }
        bool IsExternal(void) const { return true; }

        virtual void WriteVRR(std::ostream & os) const;
        virtual void WriteVRRFile(std::ostream & os, std::ostream & osh) const;
};




