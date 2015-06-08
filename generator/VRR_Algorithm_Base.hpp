#ifndef VRR_ALGORITHM_BASE_HPP
#define VRR_ALGORITHM_BASE_HPP

#include "generator/Classes.hpp"

class VRR_Algorithm_Base
{
    public:
        // this will create a map for all possible
        // components, but only some may
        // be needed (although I haven't seen such a case...)
        std::pair<VRRMap, GaussianMap> CreateAllMaps(const GaussianMap & greq);

        virtual ~VRR_Algorithm_Base() = default; 

    private:
        virtual VRRMap CreateVRRMap_(int am) = 0;
};


#endif
