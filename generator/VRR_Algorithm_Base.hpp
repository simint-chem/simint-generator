#ifndef VRR_ALGORITHM_BASE_HPP
#define VRR_ALGORITHM_BASE_HPP

#include "generator/Classes.hpp"

class VRR_Algorithm_Base
{
    public:
        void Create(const QuartetSet & q);

        VRRMap GetVRRMap(void) const;
        int GetMaxFm(void) const;
        VRRMReq GetVRRMReq(void) const;

        virtual ~VRR_Algorithm_Base() = default; 

    private:
        // VRRMap maps a AM quartet to its steps
        VRRMap vrrmap_;

        // Maximum m value needed for a quartet
        VRRMReq vrrmreq_;

        QuartetSet top_;

        void PruneQuartets_(QuartetSet & q) const;
        virtual VRRStep VRRStep_(const Quartet & q) = 0;
};


#endif
