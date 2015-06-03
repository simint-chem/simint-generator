#ifndef VRRWRITER_HPP
#define VRRWRITER_HPP

#include <iostream>
#include <utility>

#include "generator/Classes.hpp"

class WriterBase;


class VRRWriter
{   
    public:
        VRRWriter(const VRRMap & vrrmap, const VRRReqMap & vrrreqmap);

        void WriteVRR(std::ostream & os, const WriterBase & base) const;

        void WriteIncludes(std::ostream & os, const WriterBase & base) const;
        void DeclarePointers(std::ostream & os, const WriterBase & base) const;
        void DeclareAuxArrays(std::ostream & os, const WriterBase & base) const;

        bool HasVRR(void) const;

    private:
        VRRMap vrrmap_;
        VRRReqMap vrrreqmap_;

        void WriteVRRInline_(std::ostream & os, const WriterBase & base) const;
        void WriteVRRExternal_(std::ostream & os, const WriterBase & base) const;

        void WriteVRRFile_(std::ostream & os, const WriterBase & base) const;

        void WriteVRRSteps_(std::ostream & os, const WriterBase & base) const;
        void WriteVRRSteps_(std::ostream & os, const WriterBase & base, const GaussianSet & greq, const std::string & num_m) const;
        void WriteVRRSteps_(std::ostream & os, const WriterBase & base, const GaussianSet & greq, int num_m) const;

};

#endif
