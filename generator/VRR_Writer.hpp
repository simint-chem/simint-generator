#ifndef VRRWRITER_HPP
#define VRRWRITER_HPP

#include <iostream>

#include "generator/Classes.hpp"
#include "generator/WriterBase.hpp"
#include "generator/VRR_Algorithm_Base.hpp"

// foward declare
class ERIGeneratorInfo;



class VRR_Writer : public WriterBase
{   
    public:
        VRR_Writer(const VRR_Algorithm_Base & vrr_algo, const ERIGeneratorInfo & info);

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
        const VRR_Algorithm_Base & vrr_algo_;
        const ERIGeneratorInfo & info_;
        const VectorInfo & vinfo_;

        void WriteVRRSteps_(std::ostream & os, QAM qam, const VRR_StepSet & vs, const std::string & num_n) const;
};



class VRR_Writer_Inline : public VRR_Writer
{
    public:
        using VRR_Writer::VRR_Writer;

        bool IsInline(void) const { return true; }
        bool IsExternal(void) const { return false; }

        virtual ConstantMap GetConstants(void) const;
        virtual void WriteVRR(std::ostream & os) const;
        virtual void WriteVRRFile(std::ostream & os, std::ostream & osh) const;
};


class VRR_Writer_External : public VRR_Writer
{
    public:
        using VRR_Writer::VRR_Writer;

        bool IsInline(void) const { return false; }
        bool IsExternal(void) const { return true; }

        virtual void WriteVRR(std::ostream & os) const;
        virtual void WriteVRRFile(std::ostream & os, std::ostream & osh) const;
};



#endif
