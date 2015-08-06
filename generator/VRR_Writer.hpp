#ifndef VRRWRITER_HPP
#define VRRWRITER_HPP

#include <iostream>
#include <utility>

#include "generator/Classes.hpp"


class VRR_Algorithm_Base;

class VRR_Writer
{   
    public:
        VRR_Writer(const VRR_Algorithm_Base & vrr_algo_);

        void WriteVRR(std::ostream & os) const;

        void AddConstants(void) const;
        void DeclarePrimArrays(std::ostream & os) const;
        void DeclarePrimPointers(std::ostream & os) const;

        void WriteVRRFile(std::ostream & os) const;
        void WriteVRRHeaderFile(std::ostream & os) const;

    private:
        VRRMap vrrmap_;
        int maxFm_;
        VRRMReq vrrmreq_;

        std::map<QAM, QAMSet> qamreq_; // quartets required for a particular QAM
        std::map<QAM, std::set<std::string>> varreq_; // other variables required for a particular QAM
        std::set<std::string> allvarreq_;

        int maxint_;

        void WriteVRRInline_(std::ostream & os) const;
        void WriteVRRExternal_(std::ostream & os) const;

        void WriteVRRSteps_(std::ostream & os) const;
        void WriteVRRSteps_(std::ostream & os, QAM qam, const VRRStepList & vs, const std::string & num_n) const;
};

#endif
