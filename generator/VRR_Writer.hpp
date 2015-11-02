#ifndef VRRWRITER_HPP
#define VRRWRITER_HPP

#include <iostream>
#include <utility>

#include "generator/Classes.hpp"
#include "generator/VRR_Algorithm_Base.hpp"

class VRR_Writer
{   
    public:
        VRR_Writer(const VRR_Algorithm_Base & vrr_algo);

        void WriteVRR(std::ostream & os) const;

        void AddConstants(void) const;
        void DeclarePrimArrays(std::ostream & os) const;
        void DeclarePrimPointers(std::ostream & os) const;

        void WriteVRRFile(std::ostream & os, std::ostream & osh) const;

        bool HasBraVRR(void) const;
        bool HasKetVRR(void) const;

    private:
        const VRR_Algorithm_Base & vrr_algo_;

        void WriteVRRInline_(std::ostream & os) const;
        void WriteVRRExternal_(std::ostream & os) const;

        void WriteVRRSteps_(std::ostream & os) const;
        void WriteVRRSteps_(std::ostream & os, QAM qam, const VRR_StepSet & vs, const std::string & num_n) const;
};

#endif
