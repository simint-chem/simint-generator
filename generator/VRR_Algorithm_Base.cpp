#include <iostream>

#include "generator/VRR_Algorithm_Base.hpp"

using namespace std;

std::pair<VRRMap, GaussianMap> VRR_Algorithm_Base::CreateAllMaps(const GaussianMap & greq)
{
    // holds the requirements for each am
    GaussianMap vrm = greq;

    // max am
    int maxam = vrm.rbegin()->first;

    // holds the VRR steps
    VRRMap vm;

    for(int i = 0; i <= maxam; i++)
    {
        VRRMap vm2 = CreateVRRMap(i);
        vm.insert(vm2.begin(), vm2.end());
    }

    // recurse down from the end
    for(int i = maxam-1; i >= 0; i--)
    {
        // get the set from the previous am
        GaussianSet prev = vrm[i+1];
        
        // for each one, look up its requirements
        for(const auto & it : prev)
        {
            XYZStep s = vm[it];

            // add what it needs to this step
            Gaussian g1 = it.StepDown(s, 1);
            Gaussian g2 = it.StepDown(s, 2);

            // if these are valid
            if(g1)
                vrm[g1.am()].insert(g1);
            if(g2)
                vrm[g2.am()].insert(g2);
        }
    }

    return {vm, vrm}; 
}

