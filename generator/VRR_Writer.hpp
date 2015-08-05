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

        void WriteIncludes(std::ostream & os) const;
        void AddConstants(void) const;
        void DeclarePrimArrays(std::ostream & os) const;
        void DeclarePrimPointers(std::ostream & os) const;

        void WriteVRRFile(std::ostream & os) const;
        void WriteVRRHeaderFile(std::ostream & os) const;

    private:
        VRRMap vrrmap_;
        int maxFm_;
        VRRMReq vrrmreq_;

        // all the vrr_i parameters needed
        // (multiplied by 1/2p or other similar factors)
        std::set<int> vrr_bra_i_, vrr_bra_j_, vrr_bra_k_, vrr_bra_l_;
        std::set<int> vrr_ket_i_, vrr_ket_j_, vrr_ket_k_, vrr_ket_l_;
        std::set<int> vrr_2pq_;

        void WriteVRRInline_(std::ostream & os) const;
        void WriteVRRExternal_(std::ostream & os) const;

        void WriteVRRSteps_(std::ostream & os) const;
        void WriteVRRSteps_(std::ostream & os, QAM qam, const VRRStepList & vs, const std::string & num_n) const;
};

#endif
