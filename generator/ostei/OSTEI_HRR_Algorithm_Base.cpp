/*! \file
 * 
 * \brief Base class for OSTEI Horizontal Recurrence Relation steps (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#include <iostream>
#include <algorithm>
#include "generator/Printing.hpp"
#include "generator/ostei/OSTEI_HRR_Algorithm_Base.hpp"


OSTEI_HRR_Algorithm_Base::OSTEI_HRR_Algorithm_Base(const OptionMap & options)
    : options_(options)
{
}

int OSTEI_HRR_Algorithm_Base::GetOption(Option opt) const
{
    return options_.at(opt);
}

void OSTEI_HRR_Algorithm_Base::PruneDoublets_(DoubletSet & d, DoubletSet & pruned, RRStepType steptype)
{
    DoubletSet dnew;

    for(auto & it : d)
    {
        if(steptype == RRStepType::J || steptype == RRStepType::L)
        {
            if(it.right && it.right.am() != 0)
                dnew.insert(it);
            else if(it)
                pruned.insert(it);
        }
        else
        {
            if(it.left && it.left.am() != 0)
                dnew.insert(it);
            else if(it)
                pruned.insert(it);
        }
    }

    d = dnew; 
}

void OSTEI_HRR_Algorithm_Base::AMOrder_AddWithDependencies_(DAMList & order, DAM am, DoubletType type) const
{
    // skip if already done
    if(std::find(order.begin(), order.end(), am) != order.end())
        return;

    // skip if it's not done by HRR
    if(type == DoubletType::BRA)
    {
        if(allbraam_.count(am) == 0)
            return;
    }
    else
    {
        if(allketam_.count(am) == 0)
            return;
    }

    // get requirements
    DAMSet dam;
    if(type == DoubletType::BRA)
        dam = GetBraAMReq(am);
    else
        dam = GetKetAMReq(am);

    for(const auto & it : dam)
        AMOrder_AddWithDependencies_(order, it, type);

    order.push_back(am);
}



void OSTEI_HRR_Algorithm_Base::HRRDoubletLoop_(HRRDoubletStepList & hrrlist,
                                         const DoubletSet & inittargets,
                                         DoubletSet & solveddoublets,
                                         DoubletSet & pruned,
                                         RRStepType steptype)
{
    DoubletSet targets = inittargets;

    while(targets.size())
    {
        DoubletSet newtargets;

        for(auto it = targets.rbegin(); it != targets.rend(); ++it)
        {
            // skip if done alread
            // this can happen with some of the more complex trees
            if(solveddoublets.count(*it) > 0)
                continue;

            HRRDoubletStep hrr = this->DoubletStep_(*it, steptype);
            hrrlist.push_back(hrr);

            if(solveddoublets.count(hrr.src[0]) == 0)
                newtargets.insert(hrr.src[0]);
            if(solveddoublets.count(hrr.src[1]) == 0)
                newtargets.insert(hrr.src[1]);
            
            solveddoublets.insert(*it);
        }

        //cout << "Generated " << newtargets.size() << " new targets\n";

        PruneDoublets_(newtargets, pruned, steptype);

        //cout << "After pruning: " << newtargets.size() << " new targets\n";
        //for(const auto & it : newtargets)
        //    std::cout << "    " << it << "\n";

        targets = newtargets;
    } 

}



void OSTEI_HRR_Algorithm_Base::Create(QAM am)
{
    finalam_ = am;

    // First, we need a list of doublet steps for the bra

    // holds all the 'solved' doublets
    DoubletSet solvedbras;

    // generate initial targets
    DoubletSet initbras = GenerateInitialDoubletTargets({am[0],am[1]}, DoubletType::BRA);
    //PrintDoubletSet(initbras, "Initial Targets");

    // What direction are we going?
    RRStepType steptype = RRStepType::J; // Moving I -> J
    if(am[0] < am[1])
        steptype = RRStepType::I;  // Moving J -> I

    // Initial bra targets
    DoubletSet targets = initbras;
    PruneDoublets_(targets, bratop_, steptype);
    //PrintDoubletSet(targets, "Initial bra targets");

    // Solve the bra part
    HRRDoubletStepList brasteps;
    HRRDoubletLoop_(brasteps, targets, solvedbras, bratop_, steptype);
    std::reverse(brasteps.begin(), brasteps.end());

    // store in map by AM
    // and store requirements
    for(const auto & it : brasteps)
    {
        brasteps_[it.target.amlist()].push_back(it);
        for(const auto & it2 : it.src)
        {
            if(it2)
                brareq_[it.target.amlist()].insert(it2.amlist());
        }
    }


    //std::cout << "\n\n";
    //std::cout << "--------------------------------------------------------------------------------\n";
    //std::cout << "BRA HRR step done.\n";
    //std::cout << "--------------------------------------------------------------------------------\n";
    //for(const auto & it : brasteps_)
    //for(const auto & it2 : it.second)
    //    std::cout << it2 << "\n";
    

    // now do kets
    // What direction are we going?
    steptype = RRStepType::L;  // Moving K->L
    if(am[2] < am[3])
        steptype = RRStepType::K; // Moving L->K

    // we only need the bras from the original targets
    DoubletSet solvedkets;
    DoubletSet initkets = GenerateInitialDoubletTargets({am[2],am[3]}, DoubletType::KET);
    targets = initkets;
    PruneDoublets_(targets, kettop_, steptype);

    std::cout << "\n\n";
    //PrintDoubletSet(targets, "Initial ket targets");

    // Solve the ket part
    HRRDoubletStepList ketsteps;
    HRRDoubletLoop_(ketsteps, targets, solvedkets, kettop_, steptype);
    std::reverse(ketsteps.begin(), ketsteps.end());

    // store in map by AM
    // and store requirements
    for(const auto & it : ketsteps)
    {
        ketsteps_[it.target.amlist()].push_back(it);
        for(const auto & it2 : it.src)
            ketreq_[it.target.amlist()].insert(it2.amlist());
    }
    

    //std::cout << "\n\n";
    //std::cout << "--------------------------------------------------------------------------------\n";
    //std::cout << "KET HRR step done.\n";
    //std::cout << "--------------------------------------------------------------------------------\n";
    //for(const auto & it : ketsteps_)
    //for(const auto & it2 : it.second)
    //    std::cout << it2 << "\n";

    std::cout << "\n\n";

    // fill in top info
    for(const auto & it : bratop_)
        bratopam_.insert(it.amlist());
    for(const auto & it : kettop_)
        kettopam_.insert(it.amlist());

    // form top quartets
    // all top bra/ket combo
    for(const auto & it : bratop_)
    for(const auto & it2 : kettop_)
        topquartets_.insert({it, it2, 0});

    for(const auto & it : topquartets_)
        topqam_.insert(it.amlist());


    // solved info
    for(const auto & it : solvedbras)
        allbraam_.insert(it.amlist());
    for(const auto & it : solvedkets)
        allketam_.insert(it.amlist());


    // Put am in the order of calculation
    AMOrder_AddWithDependencies_(braamorder_, {am[0], am[1]}, DoubletType::BRA);
    AMOrder_AddWithDependencies_(ketamorder_, {am[2], am[3]}, DoubletType::KET);

    //PrintDoubletSet(bratop_, "HRR top level bras");
    //PrintDoubletSet(kettop_, "HRR top level kets");
    //PrintQuartetSet(topquartets_, "HRR top quartets");
}


DAMSet OSTEI_HRR_Algorithm_Base::TopBraAM(void) const
{
    return bratopam_;
}

DAMSet OSTEI_HRR_Algorithm_Base::TopKetAM(void) const
{
    return kettopam_;
}

DoubletSet OSTEI_HRR_Algorithm_Base::TopBraDoublets(void) const
{
    return bratop_;
}

DoubletSet OSTEI_HRR_Algorithm_Base::TopKetDoublets(void) const
{
    return kettop_;
}

QAMSet OSTEI_HRR_Algorithm_Base::TopAM(void) const
{
    QAMSet qs;
    for(const auto & it1 : TopBraAM())
    for(const auto & it2 : TopKetAM())
        qs.insert({it1[0], it1[1], it2[0], it2[1]});
    return qs;
}

QuartetSet OSTEI_HRR_Algorithm_Base::TopQuartets(void) const
{
    return topquartets_;
}

RRStepType OSTEI_HRR_Algorithm_Base::GetBraRRStep(DAM dam) const
{
    return brasteps_.at(dam).begin()->type;
}

RRStepType OSTEI_HRR_Algorithm_Base::GetKetRRStep(DAM dam) const
{
    return ketsteps_.at(dam).begin()->type;
}

DAMSet OSTEI_HRR_Algorithm_Base::GetBraAMReq(DAM am) const
{
    return brareq_.at(am);
}

DAMSet OSTEI_HRR_Algorithm_Base::GetKetAMReq(DAM am) const
{
    return ketreq_.at(am);
}


HRRDoubletStepList OSTEI_HRR_Algorithm_Base::GetBraSteps(DAM am) const
{
    return brasteps_.at(am);
}

HRRDoubletStepList OSTEI_HRR_Algorithm_Base::GetKetSteps(DAM am) const
{
    return ketsteps_.at(am);
}


DAMList OSTEI_HRR_Algorithm_Base::GetBraAMOrder(void) const
{
    return braamorder_;
}

DAMList OSTEI_HRR_Algorithm_Base::GetKetAMOrder(void) const
{
    return ketamorder_;
}

bool OSTEI_HRR_Algorithm_Base::HasHRR(void) const
{
    return HasBraHRR() || HasKetHRR();
}

bool OSTEI_HRR_Algorithm_Base::HasBraHRR(void) const
{
    return brasteps_.size();
}

bool OSTEI_HRR_Algorithm_Base::HasKetHRR(void) const
{
    return ketsteps_.size();
}

QAMList OSTEI_HRR_Algorithm_Base::GenerateAMReq(QAM am, RRStepType rrstep) const
{
    QAMList req;
    switch(rrstep)
    {
        case RRStepType::I:
            req.push_back({am[0]-1, am[1]+1, am[2]  , am[3]  });
            req.push_back({am[0]-1, am[1]  , am[2]  , am[3]  });
            break;
        case RRStepType::J:
            req.push_back({am[0]+1, am[1]-1, am[2]  , am[3]  });
            req.push_back({am[0]  , am[1]-1, am[2]  , am[3]  });
            break;
        case RRStepType::K:
            req.push_back({am[0]  , am[1]  , am[2]-1, am[3]+1});
            req.push_back({am[0]  , am[1]  , am[2]-1, am[3]  });
            break;
        case RRStepType::L:
            req.push_back({am[0]  , am[1]  , am[2]+1, am[3]-1});
            req.push_back({am[0]  , am[1]  , am[2]  , am[3]-1});
            break;
    }

    return req;
}

