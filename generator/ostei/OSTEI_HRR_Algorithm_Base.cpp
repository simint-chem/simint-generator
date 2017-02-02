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

void OSTEI_HRR_Algorithm_Base::AMOrder_AddWithDependencies_(QAMList & order,
                                                            QAMSet & topqam,
                                                            std::map<QAM, DoubletType> & smap,
                                                            QAMSet am) const
{
    QAMSet qamnew;
    QAMSet amdone;

    for(auto & iam : am)
    {
        DAM braam = {iam[0], iam[1]};
        DAM ketam = {iam[2], iam[3]};

        // skip if already done
        if(std::find(order.begin(), order.end(), iam) != order.end())
            continue;

        // skip if it's not done by HRR
        if(allbraam_.count(braam) == 0 && allketam_.count(ketam) == 0)
        {
            topqam.insert(iam);
            continue;
        }


        // prefer ket first
        if(allketam_.count(ketam))
        {
            // tag gets re-added here
            smap[iam] = DoubletType::KET;
            DAMSet dam = GetKetAMReq(ketam);
            for(const auto & it : dam)
                qamnew.insert({iam[0], iam[1], it[0], it[1], iam.tag});
        }
        else
        {
            // do bra
            // tag gets re-added here
            smap[iam] = DoubletType::BRA;
            DAMSet dam = GetBraAMReq(braam);
            for(const auto & it : dam)
                qamnew.insert({it[0], it[1], iam[2], iam[3], iam.tag});
        }

        amdone.insert(iam);
    }

    if(qamnew.size())
        AMOrder_AddWithDependencies_(order, topqam, smap, qamnew);

    // check for duplicates since they may have been added by a dependency
    for(auto & iam : amdone)
    {
        if(std::find(order.begin(), order.end(), iam) == order.end())
            order.push_back(iam);
    }
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
    Create(std::set<QAM>{am});
}

void OSTEI_HRR_Algorithm_Base::Create(std::set<QAM> am)
{
    // What is the am we base the direction on?
    QAM highest_am = *(am.rbegin());

    // First, we need a list of doublet steps for the ket
    // What direction are we going?
    RRStepType ketsteptype = RRStepType::L;  // Moving K->L
    if(highest_am[2] < highest_am[3])
        ketsteptype = RRStepType::K; // Moving L->K

    // we need the kets from the original targets
    DoubletSet initkets, solvedkets, prunedket;

    for(auto & iam : am)
    {
        // strip of the tag here
        auto add_targets = GenerateDoubletTargets({iam[2],iam[3]}, DoubletType::KET);
        initkets.insert(add_targets.begin(), add_targets.end());
    }

    DoubletSet targets = initkets;
    PruneDoublets_(targets, prunedket, ketsteptype);

    //std::cout << "\n\n";
    //PrintDoubletSet(targets, "Initial ket targets");

    // Solve the ket part
    HRRDoubletStepList ketsteps;
    HRRDoubletLoop_(ketsteps, targets, solvedkets, prunedket, ketsteptype);
    std::reverse(ketsteps.begin(), ketsteps.end());

    // store in map by AM
    // and store requirements
    for(const auto & it : ketsteps)
    {
        ketsteps_[it.target.amlist()].push_back(it);
        for(const auto & it2 : it.src)
            ketreq_[it.target.amlist()].insert(it2.amlist());
    }
    

    std::cout << "\n\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    std::cout << "KET HRR step done.\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    for(const auto & it : ketsteps_)
    for(const auto & it2 : it.second)
        std::cout << it2 << "\n";


    // Now for the bras
    // What direction are we going?
    RRStepType brasteptype = RRStepType::J; // Moving I -> J
    if(highest_am[0] < highest_am[1])
        brasteptype = RRStepType::I;  // Moving J -> I

    // we need the bras from the original targets
    DoubletSet initbras, solvedbras, prunedbra;
    for(auto & iam : am)
    {
        // this strips off the tag
        auto add_targets = GenerateDoubletTargets({iam[0],iam[1]}, DoubletType::BRA);
        initbras.insert(add_targets.begin(), add_targets.end());
    }


    // Initial bra targets
    targets = initbras;
    PruneDoublets_(targets, prunedbra, brasteptype);

    //std::cout << "\n\n";
    //PrintDoubletSet(targets, "Initial bra targets");

    // Solve the bra part
    HRRDoubletStepList brasteps;
    HRRDoubletLoop_(brasteps, targets, solvedbras, prunedbra, brasteptype);
    std::reverse(brasteps.begin(), brasteps.end());

    // store in map by AM
    // and store requirements
    for(const auto & it : brasteps)
    {
        brasteps_[it.target.amlist()].push_back(it);
        for(const auto & it2 : it.src)
            brareq_[it.target.amlist()].insert(it2.amlist());
    }


    std::cout << "\n\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    std::cout << "BRA HRR step done.\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    for(const auto & it : brasteps_)
    for(const auto & it2 : it.second)
        std::cout << it2 << "\n";
    

    // solved info
    for(const auto & it : solvedbras)
        allbraam_.insert(it.amlist());
    for(const auto & it : solvedkets)
        allketam_.insert(it.amlist());

    // Put am in the order of calculation
    AMOrder_AddWithDependencies_(amorder_, topqam_, stepmap_, am);

    // form the top quartets
    for(auto & it : topqam_)
    {
        auto amquartets = GenerateQuartetTargets(it);
        topquartets_.insert(amquartets.begin(), amquartets.end());
    }

    PrintQuartetSet(topquartets_, "HRR top quartets");
}

QAMSet OSTEI_HRR_Algorithm_Base::TopAM(void) const
{
    return topqam_;
}

QuartetSet OSTEI_HRR_Algorithm_Base::TopQuartets(void) const
{
    return topquartets_;
}

RRStepType OSTEI_HRR_Algorithm_Base::GetBraRRStep(DAM dam) const
{
    return brasteps_.at(dam.notag()).begin()->type;
}

RRStepType OSTEI_HRR_Algorithm_Base::GetKetRRStep(DAM dam) const
{
    return ketsteps_.at(dam.notag()).begin()->type;
}

DAMSet OSTEI_HRR_Algorithm_Base::GetBraAMReq(DAM am) const
{
    return brareq_.at(am.notag());
}

DAMSet OSTEI_HRR_Algorithm_Base::GetKetAMReq(DAM am) const
{
    return ketreq_.at(am.notag());
}

HRRDoubletStepList OSTEI_HRR_Algorithm_Base::GetBraSteps(DAM am) const
{
    return brasteps_.at(am.notag());
}

HRRDoubletStepList OSTEI_HRR_Algorithm_Base::GetKetSteps(DAM am) const
{
    return ketsteps_.at(am.notag());
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

DoubletType OSTEI_HRR_Algorithm_Base::GetDoubletStep(QAM am) const
{
    return stepmap_.at(am); 
}

QAMList OSTEI_HRR_Algorithm_Base::GetAMOrder(void) const
{
    return amorder_;
}

QAMList OSTEI_HRR_Algorithm_Base::GenerateAMReq(QAM am, RRStepType rrstep) const
{
    QAMList req;
    switch(rrstep)
    {
        case RRStepType::I:
            req.push_back({am[0]-1, am[1]+1, am[2]  , am[3]  , am.tag });
            req.push_back({am[0]-1, am[1]  , am[2]  , am[3]  , am.tag });
            break;
        case RRStepType::J:
            req.push_back({am[0]+1, am[1]-1, am[2]  , am[3]  , am.tag });
            req.push_back({am[0]  , am[1]-1, am[2]  , am[3]  , am.tag });
            break;
        case RRStepType::K:
            req.push_back({am[0]  , am[1]  , am[2]-1, am[3]+1, am.tag});
            req.push_back({am[0]  , am[1]  , am[2]-1, am[3]  , am.tag});
            break;
        case RRStepType::L:
            req.push_back({am[0]  , am[1]  , am[2]+1, am[3]-1, am.tag});
            req.push_back({am[0]  , am[1]  , am[2]  , am[3]-1, am.tag});
            break;
    }

    return req;
}

