#ifndef ETWRITER_HPP
#define ETWRITER_HPP

#include <iostream>

#include "generator/Classes.hpp"


class ET_Algorithm_Base;

class ET_Writer
{   
    public:
        ET_Writer(const ET_Algorithm_Base & et_algo); 

        void WriteIncludes(std::ostream & os) const;
        void AddConstants(void) const;
        void DeclarePrimArrays(std::ostream & os) const;
        void DeclarePrimPointers(std::ostream & os) const;

        void WriteETInline(std::ostream & os) const;

    private:
        ETStepList etsl_;
        QAMSet etint_;

        std::set<std::string> et_i_; // gets multiplied by one_over_2q

        static std::string ETStepString(const ETStep & et);
        static std::string ETStepVar(const Quartet & q);
};

#endif
