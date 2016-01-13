#ifndef ETWRITER_HPP
#define ETWRITER_HPP

#include <iostream>

#include "generator/Classes.hpp"
#include "generator/ET_Algorithm_Base.hpp"

// foward declare
class ERIGeneratorInfo;



class ET_Writer
{   
    public:
        ET_Writer(const ET_Algorithm_Base & et_algo); 

        bool HasET(void) const;
        bool HasBraET(void) const;
        bool HasKetET(void) const;

        void DeclarePrimArrays(std::ostream & os, const ERIGeneratorInfo & info) const;
        void DeclarePrimPointers(std::ostream & os, const ERIGeneratorInfo & info) const;

        virtual void AddConstants(ERIGeneratorInfo & info) const = 0;
        virtual void WriteET(std::ostream & os, const ERIGeneratorInfo & info) const = 0;
        virtual void WriteETFile(std::ostream & os, std::ostream & osh, const ERIGeneratorInfo & info) const = 0;


    protected:
        const ET_Algorithm_Base & et_algo_;

        std::string ETStepString_(const ETStep & et, const ERIGeneratorInfo & info) const;
        std::string ETStepVar_(const Quartet & q) const;

};



class ET_Writer_Inline : public ET_Writer
{
    public:
        virtual void AddConstants(ERIGeneratorInfo & info) const;
        virtual void WriteET(std::ostream & os, const ERIGeneratorInfo & info) const;
        virtual void WriteETFile(std::ostream & os, std::ostream & osh, const ERIGeneratorInfo & info) const;
};


class ET_Writer_External : public ET_Writer
{
    public:
        virtual void AddConstants(ERIGeneratorInfo & info) const;
        virtual void WriteET(std::ostream & os, const ERIGeneratorInfo & info) const;
        virtual void WriteETFile(std::ostream & os, std::ostream & osh, const ERIGeneratorInfo & info) const;
};




#endif
