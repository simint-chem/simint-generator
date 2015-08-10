#include <iostream>
#include <algorithm>

#include "generator/Helpers.hpp"
#include "generator/ET_Algorithm_Base.hpp"

using namespace std;

void ET_Algorithm_Base::PruneQuartets_(QuartetSet & qs, QuartetSet & pruned)
{
    QuartetSet qsnew;

    // by this point, only kets of the form | X 0 ) should be here
    // so prune | s s )
    for(auto & it : qs)
    {
        if(it && it.ket.left.am() != 0)
            qsnew.insert(it);
        else if(it)
            pruned.insert(it);
    }

    qs = qsnew; 
}


void ET_Algorithm_Base::AMOrder_AddWithDependencies_(QAMList & order, QAM am)
{
    // skip if it was already done somewhere
    if(std::find(order.begin(), order.end(), am) != order.end()) 
        return;

    // skip if it's not done by ET
    if(allqam_.count(am) == 0)
        return;

    // get requirements
    QAMSet req = GetAMReq(am);

    for(const auto & it : req)
        AMOrder_AddWithDependencies_(order, it);

    order.push_back(am);
}


void ET_Algorithm_Base::ETStepLoop_(ETStepList & etsl,
                                   const QuartetSet & inittargets,
                                   QuartetSet & solvedquartets,
                                   QuartetSet & pruned)
{
    QuartetSet targets = inittargets;

    // steps, sorted by their am
    std::map<QAM, ETStepList> etslmap;

    while(targets.size())
    {
        QuartetSet newtargets;

        for(auto it = targets.begin(); it != targets.end(); ++it)
        {
            // skip if done already
            // this can happen with some of the more complex trees
            if(solvedquartets.count(*it) > 0)
                continue;

            ETStep ets = this->ETStep_(*it);

            // add to the list
            etsl.push_back(ets);

            for(const auto & it2 : ets.src)
            {
                if(it2 && solvedquartets.count(it2) == 0)
                    newtargets.insert(it2);
            }
       
            // insert copy with no flags 
            solvedquartets.insert(*it);
        }

        cout << "Generated " << newtargets.size() << " new targets\n";

        PruneQuartets_(newtargets, pruned);

        //cout << "After pruning: " << newtargets.size() << " new targets\n";
        //for(const auto & it : newtargets)
        //    cout << "    " << it << "\n";

        targets = newtargets;
    } 
}


void ET_Algorithm_Base::Create(QAM am)
{
    Create(GenerateInitialQuartetTargets(am));
}


void ET_Algorithm_Base::Create(const QuartetSet & inittargets)
{
    // holds all the 'solved' quartets
    QuartetSet solvedquartets;

    // generate initial targets
    QuartetSet targets = inittargets;

    PruneQuartets_(targets, ettop_);
    PrintQuartetSet(targets, "Initial ET Targets");

    // Solve
    ETStepList etstep;
    ETStepLoop_(etstep, targets, solvedquartets, ettop_);

    // reverse the steps
    // (now done in ETStepLoop_)
    //std::reverse(etsteps.begin(), etsteps.end());

    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << "ET step done. Solution is " << etsteps_.size() << " steps\n";
    cout << "--------------------------------------------------------------------------------\n";
    //for(auto & it : etsteps_)
    //    cout << it << "\n";

    cout << "\n\n";

    // store top level stuff
    for(const auto & it : ettop_)
        ettopam_.insert(it.amlist());

    // sort the step by AM and store the requirements
    for(const auto & it : etstep)
    {
        QAM qam = it.target.amlist();
        etsteps_[qam].push_back(it);

        for(const auto & it2 : it.src)
        {
            if(it2)
                etreq_[qam].insert(it2.amlist());
        }

        int stepidx = XYZStepToIdx(it.xyz);
        int ival = it.target.bra.left.ijk[stepidx];
        int kval = (it.target.ket.left.ijk[stepidx]-1);

        // ok to add duplicates. They are stored as a map
        if(ival > 0)
            intreq_[qam].insert(ival);
        if(kval > 0)
            intreq_[qam].insert(kval);
    }

    // store all intreq
    for(const auto & it : intreq_)
        allintreq_.insert(it.second.begin(), it.second.end());

    // store all solved qam
    for(const auto & it : solvedquartets)
        allqam_.insert(it.amlist());


    // determine the proper order to do these in
    for(const auto & q : inittargets)
        AMOrder_AddWithDependencies_(amorder_, q.amlist()); 

    PrintQuartetSet(ettop_, "Top level ET");
}


ETStepList ET_Algorithm_Base::GetSteps(QAM am) const
{
    return etsteps_.at(am);
}

QAMList ET_Algorithm_Base::GetAMOrder(void) const
{
    return amorder_;
}

QAMSet ET_Algorithm_Base::GetAMReq(QAM am) const
{
    return etreq_.at(am);
}

QuartetSet ET_Algorithm_Base::TopQuartets(void) const
{
    return ettop_;
}

QAMSet ET_Algorithm_Base::TopAM(void) const
{
    return ettopam_;
}

IntSet ET_Algorithm_Base::GetIntReq(QAM am) const
{
    return intreq_.at(am);
}

IntSet ET_Algorithm_Base::GetAllInt(void) const
{
    return allintreq_;
}


