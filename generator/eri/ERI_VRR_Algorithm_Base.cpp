/*! \file
 * 
 * \brief Base class for ERI Vertical Recurrence Relation steps (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <iostream>
#include <algorithm>
#include "generator/Printing.hpp"
#include "generator/eri/ERI_VRR_Algorithm_Base.hpp"


ERI_VRR_Algorithm_Base::ERI_VRR_Algorithm_Base(const OptionMap & options)
    : options_(options)
{
}

int ERI_VRR_Algorithm_Base::GetOption(Option opt) const
{
    return options_.at(opt);
}

QAMList ERI_VRR_Algorithm_Base::GetAMOrder(void) const
{
    return amorder_; 
}

RRStepType ERI_VRR_Algorithm_Base::GetRRStep(QAM am) const
{
    // should all be the same?
    return vrrmap_.at(am).begin()->type;
}

QAMSet ERI_VRR_Algorithm_Base::GetAllAM(void) const
{
    return allam_;
}

int ERI_VRR_Algorithm_Base::GetMReq(QAM am) const
{
    return vrrmreq_max_.at(am);
}

VRR_StepSet ERI_VRR_Algorithm_Base::GetSteps(QAM am) const
{
    return vrrmap_.at(am);
}


QAMSet ERI_VRR_Algorithm_Base::GetAMReq(QAM am) const
{
    return qamreq_.at(am);
}


IntSet ERI_VRR_Algorithm_Base::GetIntReq_2p(QAM am) const
{
    if(qamint_2p_.count(am) == 0)
        return IntSet();
    else
        return qamint_2p_.at(am);
}


IntSet ERI_VRR_Algorithm_Base::GetIntReq_2q(QAM am) const
{
    if(qamint_2q_.count(am) == 0)
        return IntSet();
    else
        return qamint_2q_.at(am);
}


IntSet ERI_VRR_Algorithm_Base::GetIntReq_2pq(QAM am) const
{
    if(qamint_2pq_.count(am) == 0)
        return IntSet();
    else
        return qamint_2pq_.at(am);
}

IntSet ERI_VRR_Algorithm_Base::GetAllInt_2p(void) const
{
    IntSet iset;
    for(const auto & it : qamint_2p_)
        iset.insert(it.second.begin(), it.second.end());
    return iset;
}

IntSet ERI_VRR_Algorithm_Base::GetAllInt_2q(void) const
{
    IntSet iset;
    for(const auto & it : qamint_2q_)
        iset.insert(it.second.begin(), it.second.end());
    return iset;
}

IntSet ERI_VRR_Algorithm_Base::GetAllInt_2pq(void) const
{
    IntSet iset;
    for(const auto & it : qamint_2pq_)
        iset.insert(it.second.begin(), it.second.end());
    return iset;
}


StringSet ERI_VRR_Algorithm_Base::GetVarReq(QAM am) const
{
    return varreq_.at(am);
}


StringSet ERI_VRR_Algorithm_Base::GetAllVarReq(void) const
{
    StringSet sset;
    for(const auto & it : varreq_)
        sset.insert(it.second.begin(), it.second.end());
    return sset;
}



int ERI_VRR_Algorithm_Base::GetMaxInt(void) const
{
    return maxint_;
}



void ERI_VRR_Algorithm_Base::PruneQuartets_(QuartetSet & q) const
{
    QuartetSet qnew;
 
    for(const auto & it : q)
    {
        if(it && it.am() > 0)
            qnew.insert(it);
    }

    q = qnew;
}

QuartetSet ERI_VRR_Algorithm_Base::PruneQuartets_(QuartetSet & q, RRStepType rrstep) const
{
    // this function is a mess

    // remove ssss
    PruneQuartets_(q);
 
    // now remove any where the rrstep index is 0
    QuartetSet qnew;
    QuartetSet qret;

    int stepidx = static_cast<int>(rrstep);

    for(const auto & it : q)
    {
        if(it && it.amlist()[stepidx] > 0)
            qnew.insert(it);
        else
            qret.insert(it);
    }

    q = qnew;
    return qret;
}


bool ERI_VRR_Algorithm_Base::HasVRR(void) const
{
    return HasBraVRR() || HasKetVRR();
}

bool ERI_VRR_Algorithm_Base::HasVRROfType(RRStepType steptype) const
{
    for(const auto & it : vrrmap_)
    {
        for(const auto & it2 : it.second)
        {
            if(it2.type == steptype)
                return true;
        }
    }

    return false;
}

bool ERI_VRR_Algorithm_Base::HasBraVRR(void) const
{
    return HasVRROfType(RRStepType::I) || HasVRROfType(RRStepType::J);
}

bool ERI_VRR_Algorithm_Base::HasKetVRR(void) const
{
    return HasVRROfType(RRStepType::K) || HasVRROfType(RRStepType::L);
}

bool ERI_VRR_Algorithm_Base::HasVRR_I(void) const
{
    return HasVRROfType(RRStepType::I);
}

bool ERI_VRR_Algorithm_Base::HasVRR_J(void) const
{
    return HasVRROfType(RRStepType::J);
}

bool ERI_VRR_Algorithm_Base::HasVRR_K(void) const
{
    return HasVRROfType(RRStepType::K);
}

bool ERI_VRR_Algorithm_Base::HasVRR_L(void) const
{
    return HasVRROfType(RRStepType::L);
}


void ERI_VRR_Algorithm_Base::AMOrder_AddWithDependencies_(QAMList & order, QAM am) const
{
    // skip if it was already done somewhere
    if(std::find(order.begin(), order.end(), am) != order.end())
        return;

    // skip if it's not done by VRR
    if(allam_.count(am) == 0 || am == QAM{0,0,0,0})
        return;

    // get requirements
    QAMSet req = GetAMReq(am);

    for(const auto & it : req)
        AMOrder_AddWithDependencies_(order, it);

    order.push_back(am);
}

VRRStep ERI_VRR_Algorithm_Base::Create_Single_(Quartet q, RRStepType rrstep)
{
    VRRStep vs = VRRStep_(q, rrstep);

    // fill in the ijkl members
    int istep = XYZStepToIdx(vs.xyz);
    vs.ijkl = {0, 0, 0, 0};
    if(rrstep == RRStepType::I)
    {
        if(vs.src[2] || vs.src[3])
            vs.ijkl[0] = q.bra.left.ijk[istep]-1;
        if(vs.src[4] || vs.src[5])
            vs.ijkl[1] = q.bra.right.ijk[istep];
        if(vs.src[6])
            vs.ijkl[2] = q.ket.left.ijk[istep];
        if(vs.src[7])
            vs.ijkl[3] = q.ket.right.ijk[istep];
    }
    else if(rrstep == RRStepType::J)
    {
        if(vs.src[2] || vs.src[3])
            vs.ijkl[0] = q.bra.left.ijk[istep];
        if(vs.src[4] || vs.src[5])
            vs.ijkl[1] = q.bra.right.ijk[istep]-1;
        if(vs.src[6])
            vs.ijkl[2] = q.ket.left.ijk[istep];
        if(vs.src[7])
            vs.ijkl[3] = q.ket.right.ijk[istep];
    }
    else if(rrstep == RRStepType::K)
    {
        if(vs.src[2] || vs.src[3])
            vs.ijkl[2] = q.ket.left.ijk[istep]-1;
        if(vs.src[4] || vs.src[5])
            vs.ijkl[3] = q.ket.right.ijk[istep];
        if(vs.src[6])
            vs.ijkl[0] = q.bra.left.ijk[istep];
        if(vs.src[7])
            vs.ijkl[1] = q.bra.right.ijk[istep];
    }
    else
    {
        if(vs.src[2] || vs.src[3])
            vs.ijkl[2] = q.ket.left.ijk[istep];
        if(vs.src[4] || vs.src[5])
            vs.ijkl[3] = q.ket.right.ijk[istep]-1;
        if(vs.src[6])
            vs.ijkl[0] = q.bra.left.ijk[istep];
        if(vs.src[7])
            vs.ijkl[1] = q.bra.right.ijk[istep];
    }

    return vs;
}



void ERI_VRR_Algorithm_Base::Create_WithOrder(const QuartetSet & q, IdxOrder idx_order)
{
    // holds the requirements for each am
    // so store what we initially want
    PrintQuartetSet(q, "Initial VRR");

    VRR_StepMap newmap;

    // solved quartets
    QuartetSet solvedquartets;

    // unsolved quartets
    QuartetSet targets = q;
    PruneQuartets_(targets);

    for(int i = 0; i < 4; i++)
    {
        // check if we should stop
        if(idx_order[i] < 0)
            break;

        RRStepType rrstep = static_cast<RRStepType>(idx_order[i]);

        QuartetSet my_targets = targets;

        // set targets to be everything we aren't doing 
        targets = PruneQuartets_(my_targets, rrstep);

        while(my_targets.size())
        {
            QuartetSet newtargets;

            for(const auto & it : my_targets)
            {
                QAM qam = it.amlist();

                VRRStep vs = Create_Single_(it, rrstep);

                // add new targets
                for(const auto & it2 : vs.src)
                {
                    if(it2)
                    {
                        if(solvedquartets.count(it2) == 0)
                            newtargets.insert(it2);
                    }
                }

                newmap[qam].insert(vs);
                solvedquartets.insert(it);
            }

            // now get what I need to do
            my_targets = newtargets;
            QuartetSet not_for_me = PruneQuartets_(my_targets, rrstep);

            // add anything that is not my responsibility to the
            // main targets list. These might not be a part of
            // that list right now. 
            targets.insert(not_for_me.begin(),
                           not_for_me.end());
        }
    } 

    // save the solved quartets
    for(const auto & it : solvedquartets)
        allam_.insert(it.amlist());
    allam_.insert({0,0,0,0});

    // and we must do at least (s s | s s)^0
    vrrmreq_max_[{0,0,0,0}] = 0;

    // Handle the needed m values
    // For now, add all to the zero m value
    // but store the max m values needed
    for(auto & it : newmap)
    {
        VRR_StepSet vs = it.second;

        // notice it2 is a copy
        for(auto it2 : vs)
        {
            int m = it2.target.m;

            vrrmreq_max_[it2.target.amlist()] = std::max(m, vrrmreq_max_[it2.target.amlist()]);

            it2.target.m = 0;
            for(auto & mit : it2.src)
            {
                vrrmreq_max_[mit.amlist()] = std::max(mit.m, vrrmreq_max_[mit.amlist()]);
                mit.m -= m;
            }
   
            vrrmap_[it.first].insert(it2);
        } 
    }


    // and empty step list for s s s s
    vrrmap_[{0,0,0,0}] = VRR_StepSet();

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

            if(its.type == RRStepType::I || its.type == RRStepType::J)
            {
                if(its.src[0])
                {
                    if(its.type == RRStepType::I)
                        varreq_[qam].insert(std::string("P_PA_") + sstep);
                    else
                        varreq_[qam].insert(std::string("P_PB_") + sstep);
                }

                if(its.src[1])
                    varreq_[qam].insert(std::string("aop_PQ_") + sstep);

                if(its.src[2] || its.src[3])
                {
                    varreq_[qam].insert("one_over_2p");
                    varreq_[qam].insert("a_over_p");
                    qamint_2p_[qam].insert(its.ijkl[0]);
                }

                if(its.src[4] || its.src[5])
                {
                    varreq_[qam].insert("one_over_2p");
                    varreq_[qam].insert("a_over_p");
                    qamint_2p_[qam].insert(its.ijkl[1]);
                }

                if(its.src[6])
                {
                    varreq_[qam].insert("one_over_2pq");
                    qamint_2pq_[qam].insert(its.ijkl[2]);
                }

                if(its.src[7])
                {
                    varreq_[qam].insert("one_over_2pq");
                    qamint_2pq_[qam].insert(its.ijkl[3]);
                }
            }
            else
            {
                if(its.src[0])
                {
                    if(its.type == RRStepType::K)
                        varreq_[qam].insert(std::string("Q_PA_") + sstep);
                    else
                        varreq_[qam].insert(std::string("Q_PB_") + sstep);

                }

                if(its.src[1])
                    varreq_[qam].insert(std::string("aoq_PQ_") + sstep);

                if(its.src[2] || its.src[3])
                {
                    varreq_[qam].insert("one_over_2q");
                    varreq_[qam].insert("a_over_q");
                    qamint_2q_[qam].insert(its.ijkl[2]);
                }

                if(its.src[4] || its.src[5])
                {
                    varreq_[qam].insert("one_over_2q");
                    varreq_[qam].insert("a_over_q");
                    qamint_2q_[qam].insert(its.ijkl[3]);
                }

                if(its.src[6])
                {
                    varreq_[qam].insert("one_over_2pq");
                    qamint_2pq_[qam].insert(its.ijkl[0]);
                }

                if(its.src[7])
                {
                    varreq_[qam].insert("one_over_2pq");
                    qamint_2pq_[qam].insert(its.ijkl[1]);
                }
            }
        }
    }


    // order the AM
    for(const auto & it : q) 
        AMOrder_AddWithDependencies_(amorder_, it.amlist());

} 


void ERI_VRR_Algorithm_Base::Create(const QuartetSet & q)
{
    // Find where the lowest AM is first. This is our initial RR
    QAM qam_max{0,0,0,0};
    for(const Quartet & iq : q)
    {
        QAM aml = iq.amlist();
        for(int i = 0; i < 4; i++)
            qam_max[i] = std::max(qam_max[i], aml[i]);
    }


    // reversed, so we end up with a more "traditional" order
    // (ie, form dsss, then dsps).
    IdxOrder idx_order{3, 2, 1, 0}; 

    std::stable_sort(idx_order.begin(), idx_order.end(),
                     [ & qam_max ] ( int i, int j ) { return qam_max[i] < qam_max[j]; });

    Create_WithOrder(q, idx_order);
}


void ERI_VRR_Algorithm_Base::Create_WithOrder(QAM q, IdxOrder idx_order)
{
    Create_WithOrder(GenerateInitialQuartetTargets(q), idx_order);
}

void ERI_VRR_Algorithm_Base::Create(QAM q)
{
    Create(GenerateInitialQuartetTargets(q));
}


