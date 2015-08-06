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

std::map<QAM, std::set<std::string>> VRR_Algorithm_Base::GetVarReq(void) const
{
    return varreq_;
}

std::map<QAM, QAMSet> VRR_Algorithm_Base::GetQAMReq(void) const
{
    return qamreq_;
}


std::set<std::string> VRR_Algorithm_Base::GetAllVarReq(void) const
{
    std::set<std::string> str;
    for(const auto & it : varreq_)
        str.insert(it.second.begin(), it.second.end());
    return str;
}


int VRR_Algorithm_Base::GetMaxInt(void) const
{
    return maxint_;
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
    // (should be 0)
    for(const auto & it : q)
    {
        QAM qam = it.amlist();
        if( (vrrmreq_.count(qam) == 0) || (it.m > vrrmreq_.at(qam)))
            vrrmreq_[qam] = it.m;
    }


    while(targets.size())
    {
        QuartetSet newtargets;

        for(const auto & it : targets)
        {
            VRRStep vs = VRRStep_(it);
            QAM qam = it.amlist();

            // set the max m value
            if((vrrmreq_.count(qam) == 0) || (it.m > vrrmreq_.at(qam)))
                vrrmreq_[qam] = it.m;

            // fill in the ijkl members
            int istep = XYZStepToIdx(vs.xyz);
            vs.ijkl = {0, 0, 0, 0};
            if(vs.type == DoubletType::BRA)
            {
                if(vs.src[2] || vs.src[3])
                    vs.ijkl[0] = it.bra.left.ijk[istep]-1;
                if(vs.src[4] || vs.src[5])
                    vs.ijkl[1] = it.bra.right.ijk[istep];
                if(vs.src[6])
                    vs.ijkl[2] = it.ket.left.ijk[istep];
                if(vs.src[7])
                    vs.ijkl[3] = it.ket.right.ijk[istep];
            }
            else
            {
                if(vs.src[2] || vs.src[3])
                    vs.ijkl[0] = it.bra.left.ijk[istep];
                if(vs.src[4] || vs.src[5])
                    vs.ijkl[1] = it.bra.right.ijk[istep];
                if(vs.src[6])
                    vs.ijkl[2] = it.ket.left.ijk[istep]-1;
                if(vs.src[7])
                    vs.ijkl[3] = it.ket.right.ijk[istep];
            }
             
            
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

            newmap[qam].push_back(vs);
            solvedquartets.insert(it);
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



    // determine requirements for all QAM
    // and the maximum integer constant
    maxint_ = 0;
    for(const auto & it : vrrmap_) // for all QAM
    {
        QAM qam = it.first;

        for(const auto & its : it.second)  // for all steps for this QAM
        {
            for(const auto & it3 : its.src)
            {
                if(it3)
                    qamreq_[qam].insert(it3.amlist());
            }

            std::string sstep = XYZStepToStr(its.xyz);

            maxint_ = std::max(its.ijkl[0], maxint_);
            maxint_ = std::max(its.ijkl[1], maxint_);
            maxint_ = std::max(its.ijkl[2], maxint_);
            maxint_ = std::max(its.ijkl[3], maxint_);

            if(its.type == DoubletType::BRA)
            {
                if(its.src[0])
                    varreq_[qam].insert(std::string("P_PA_") + sstep);

                if(its.src[1])
                    varreq_[qam].insert(std::string("aop_PQ_") + sstep);

                if(its.src[2] || its.src[3])
                {
                    varreq_[qam].insert(std::string("vrr_") + std::to_string(its.ijkl[0]) + "_over_2p");
                    varreq_[qam].insert(std::string("a_over_p"));
                }

                if(its.src[4] || its.src[5])
                {
                    varreq_[qam].insert(std::string("vrr_") + std::to_string(its.ijkl[1]) + "_over_2p");
                    varreq_[qam].insert(std::string("a_over_p"));
                }

                if(its.src[6])
                    varreq_[qam].insert(std::string("vrr_") + std::to_string(its.ijkl[2]) + "_over_2pq");

                if(its.src[7])
                    varreq_[qam].insert(std::string("vrr_") + std::to_string(its.ijkl[3]) + "_over_2pq");
            }
            else
            {
                if(its.src[0])
                    varreq_[qam].insert(std::string("Q_PA_") + sstep);

                if(its.src[1])
                    varreq_[qam].insert(std::string("aoq_PQ_") + sstep);

                if(its.src[2] || its.src[3])
                {
                    varreq_[qam].insert(std::string("vrr_") + std::to_string(its.ijkl[2]) + "_over_2q");
                    varreq_[qam].insert(std::string("a_over_q"));
                }

                if(its.src[4] || its.src[5])
                {
                    varreq_[qam].insert(std::string("vrr_") + std::to_string(its.ijkl[3]) + "_over_2q");
                    varreq_[qam].insert(std::string("a_over_p"));
                }

                if(its.src[6])
                    varreq_[qam].insert(std::string("vrr_") + std::to_string(its.ijkl[0]) + "_over_2pq");

                if(its.src[7])
                    varreq_[qam].insert(std::string("vrr_") + std::to_string(its.ijkl[1]) + "_over_2pq");
            }
        }
    }



}

