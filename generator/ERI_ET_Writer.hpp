#ifndef SIMINT_GUARD_GENERATOR__ERI_ET_WRITER_HPP_
#define SIMINT_GUARD_GENERATOR__ERI_ET_WRITER_HPP_

#include <iostream>

#include "generator/Classes.hpp"
#include "generator/WriterBase.hpp"
#include "generator/ERI_ET_Algorithm_Base.hpp"

// foward declare
class ERIGeneratorInfo;



class ERI_ET_Writer : public WriterBase
{   
    public:
        ERI_ET_Writer(const ERI_ET_Algorithm_Base & et_algo, const ERIGeneratorInfo & info); 

        bool HasET(void) const;
        bool HasBraET(void) const;
        bool HasKetET(void) const;

        void DeclarePrimArrays(std::ostream & os) const;
        void DeclarePrimPointers(std::ostream & os) const;

        virtual ConstantMap GetConstants(void) const;
        virtual void WriteET(std::ostream & os) const = 0;
        virtual void WriteETFile(std::ostream & os, std::ostream & osh) const = 0;


    protected:
        const ERI_ET_Algorithm_Base & et_algo_;
        const ERIGeneratorInfo & info_;
        const VectorInfo & vinfo_;

        std::string ETStepString_(const ETStep & et) const;
        std::string ETStepVar_(const Quartet & q) const;

};



class ERI_ET_Writer_Inline : public ERI_ET_Writer
{
    public:
        using ERI_ET_Writer::ERI_ET_Writer;

        bool IsInline(void) const { return true; }
        bool IsExternal(void) const { return false; }

        virtual ConstantMap GetConstants(void) const;
        virtual void WriteET(std::ostream & os) const;
        virtual void WriteETFile(std::ostream & os, std::ostream & osh) const;
};


class ERI_ET_Writer_External : public ERI_ET_Writer
{
    public:
        using ERI_ET_Writer::ERI_ET_Writer;

        bool IsInline(void) const { return false; }
        bool IsExternal(void) const { return true; }

        virtual void WriteET(std::ostream & os) const;
        virtual void WriteETFile(std::ostream & os, std::ostream & osh) const;
};




#endif
