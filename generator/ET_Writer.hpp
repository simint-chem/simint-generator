#ifndef ETWRITER_HPP
#define ETWRITER_HPP

#include <iostream>

#include "generator/Classes.hpp"

class WriterBase;


class ET_Writer
{   
    public:
        ET_Writer(const ETStepList & etsl); 

        bool HasET(void) const;

        void DeclarePointers(std::ostream & os, const WriterBase & base) const;
        void DeclarePrimArrays(std::ostream & os, const WriterBase & base) const;

        void WriteETInline(std::ostream & os, const WriterBase & base) const;

    private:
        ETStepList etsl_;

        static std::string ETStepString(const ETStep & et, const WriterBase & base);
        static std::string ETStepVar(const Quartet & q, const WriterBase & base);
};

#endif
