#include <iostream>

#include "generator/Helpers.hpp"
#include "generator/VRR_Algorithm_Base.hpp"

using namespace std;

VRRMap VRR_Algorithm_Base::GetVRRMap(void) const
{
    return vrrmap_;
}

int VRR_Algorithm_Base::GetMaxFm(void) const
{
    return vrrmreq_.at({0, 0, 0, 0});
}

VRRMReq VRR_Algorithm_Base::GetVRRMReq(void) const
{
    return vrrmreq_;
}


void VRR_Algorithm_Base::PruneQuartets_(QuartetSet & q) const
{
    QuartetSet qnew;
 
    for(const auto & it : q)
    {
        if(it && it.am() > 0)
            qnew.insert(it);
    }

    q = qnew;
}


void VRR_Algorithm_Base::Create(const QuartetSet & q)
{
    // holds the requirements for each am
    // so store what we initially want
    PrintQuartetSet(q, "Initial VRR");

    VRRMap newmap;

    // solved quartets
    QuartetSet solvedquartets;

    // unsolved quartets
    QuartetSet targets = q;
    PruneQuartets_(targets);

    // add max m for initial targets
    for(const auto & it : q)
    {
        QAM qam = it.amlist();
        if( (vrrmreq_.count(qam) == 0) || (it.m > vrrmreq_.at(qam)))
            vrrmreq_[qam] = it.m;
    }


    while(targets.size())
    {
        QuartetSet newtargets;

        for(auto it = targets.rbegin(); it != targets.rend(); ++it)
        {
            VRRStep vs = VRRStep_(*it);
            QAM qam = it->amlist();

            newmap[qam].push_back(vs);

            if((vrrmreq_.count(qam) == 0) || (it->m > vrrmreq_.at(qam)))
                vrrmreq_[qam] = it->m;
            
            // add new targets and m values 
            for(const auto & it2 : vs.src)
            {
                if(it2)
                {
                    if(solvedquartets.count(it2) == 0)
                        newtargets.insert(it2);

                    QAM qam2 = it2.amlist();
                    if( (vrrmreq_.count(qam2) == 0) || (it2.m > vrrmreq_.at(qam2)))
                        vrrmreq_[qam2] = it2.m;
                }
            }

            solvedquartets.insert(*it);
        }

        PruneQuartets_(newtargets);

        targets = newtargets;
    } 


    // remove all steps where the target m value > 0
    // This is called from within a loop over m
    // so it is handled by index math
    for(auto & it : newmap)
    {
        VRRStepList vs = it.second;
        for(auto & it2 : vs)
        {
            if(it2.target.m == 0)
                vrrmap_[it.first].push_back(it2);
        } 
    }

    // and empty step list for s s s s
    vrrmap_[{0,0,0,0}] = VRRStepList();

}

