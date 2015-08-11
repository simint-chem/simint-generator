#ifndef ETWRITER_HPP
#define ETWRITER_HPP

#include <iostream>

#include "generator/Classes.hpp"


class ET_Algorithm_Base;

class ET_Writer
{   
    public:
        ET_Writer(const ET_Algorithm_Base & et_algo); 

        void WriteET(std::ostream & os) const;

        void AddConstants(void) const;
        void DeclarePrimArrays(std::ostream & os) const;
        void DeclarePrimPointers(std::ostream & os) const;

        void WriteETFile(std::ostream & os, std::ostream & osh) const;

        bool HasBraET(void) const;
        bool HasKetET(void) const;

    private:
        const ET_Algorithm_Base & et_algo_;

        std::string ETStepString_(const ETStep & et) const;
        std::string ETStepVar_(const Quartet & q) const;

        void WriteETInline_(std::ostream & os) const;
        void WriteETExternal_(std::ostream & os) const;
};

#endif
