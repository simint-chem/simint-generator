#pragma once


#include <iostream>

#include "generator/Types.hpp"
#include "generator/ostei/OSTEI_ET_Algorithm_Base.hpp"

// foward declare
class OSTEI_GeneratorInfo;



class OSTEI_ET_Writer
{   
    public:
        OSTEI_ET_Writer(const OSTEI_ET_Algorithm_Base & et_algo, const OSTEI_GeneratorInfo & info); 

        bool HasET(void) const;
        bool HasBraET(void) const;
        bool HasKetET(void) const;

        void DeclarePrimArrays(std::ostream & os) const;
        void DeclarePrimPointers(std::ostream & os) const;

        virtual ConstantMap GetConstants(void) const;
        virtual void WriteET(std::ostream & os) const = 0;
        virtual void WriteETFile(std::ostream & os, std::ostream & osh) const = 0;


    protected:
        const OSTEI_ET_Algorithm_Base & et_algo_;
        const OSTEI_GeneratorInfo & info_;
        const VectorInfo & vinfo_;

        std::string ETStepString_(const ETStep & et) const;
        std::string ETStepVar_(const Quartet & q) const;

};



class OSTEI_ET_Writer_Inline : public OSTEI_ET_Writer
{
    public:
        using OSTEI_ET_Writer::OSTEI_ET_Writer;

        virtual ConstantMap GetConstants(void) const;
        virtual void WriteET(std::ostream & os) const;
        virtual void WriteETFile(std::ostream & os, std::ostream & osh) const;
};


class OSTEI_ET_Writer_External : public OSTEI_ET_Writer
{
    public:
        using OSTEI_ET_Writer::OSTEI_ET_Writer;

        virtual void WriteET(std::ostream & os) const;
        virtual void WriteETFile(std::ostream & os, std::ostream & osh) const;
};





