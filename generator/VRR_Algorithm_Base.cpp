#include <iostream>

#include "generator/Helpers.hpp"
#include "generator/VRR_Algorithm_Base.hpp"

using namespace std;

VRRMap VRR_Algorithm_Base::GetVRRMap(void) const
{
    return vrrmap_;
}


GaussianMap VRR_Algorithm_Base::GetAMReq(void) const
{
    return amreq_;
}

void VRR_Algorithm_Base::CreateAllMaps(const GaussianMap & greq)
{
    // holds the requirements for each am
    // so store what we initially want
    amreq_ = greq;
    PrintGaussianMap(greq, "Initial VRR");

    // max am is the last entry of the map
    int maxam = amreq_.rbegin()->first;

    // generate the steps for all am
    for(int i = 0; i <= maxam; i++)
    {
        VRRMap vm2 = CreateVRRMap_(i);
        vrrmap_.insert(vm2.begin(), vm2.end());
    }

    // recurse down from the end
    for(int i = maxam-1; i >= 0; i--)
    {
        // get the set from the previous am
        GaussianSet prev = amreq_[i+1];
        
        // for each one, look up its requirements
        for(const auto & it : prev)
        {
            XYZStep s = vrrmap_[it];

            // add what it needs to this step
            Gaussian g1 = it.StepDown(s, 1);
            Gaussian g2 = it.StepDown(s, 2);

            // if these are valid
            if(g1)
                amreq_[g1.am()].insert(g1);
            if(g2)
                amreq_[g2.am()].insert(g2);
        }
    }

    cout << "\n\n";
}

