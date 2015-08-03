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

void ET_Algorithm_Base::ETAddWithDependencies_(std::vector<QAM> & amorder, QAM am)
{
    // skip if it was already done somewhere
    if(std::find(amorder.begin(), amorder.end(), am) != amorder.end()) 
        return;

    // also skip if this is a VRR result
    if(am[2] == 0)
        return;
 
    // 4 possible dependencies
    std::vector<QAM> deps;  

    if(am[0] >= 0 && am[2] >= 1)
    {
        deps.push_back({am[0],   0, am[2]-1, 0});
        deps.push_back({am[0]+1, 0, am[2]-1, 0});
    }
    if(am[0] >= 0 && am[2] >= 2)
    {
        deps.push_back({am[0]  , 0, am[2]-2, 0});
    }
    if(am[0] >= 1 && am[2] >= 1)
    {
        deps.push_back({am[0]-1, 0, am[2]-1, 0});
    }


    // Yo dawg, I heard your dependencies might have dependencies
    for(const auto & it : deps)
        ETAddWithDependencies_(amorder, it);

    // add am to amorder
    // If it exists, then it was a dependency of itself. Uh oh
    if(std::find(amorder.begin(), amorder.end(), am) != amorder.end()) 
        throw std::runtime_error("Circular dependency in electron transfer tree");

    amorder.push_back(am);
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

            // add to the map
            etslmap[ets.target.amlist()].push_back(ets);

            if(solvedquartets.count(ets.src[0]) == 0)
                newtargets.insert(ets.src[0]);

            if(solvedquartets.count(ets.src[1]) == 0)
                newtargets.insert(ets.src[1]);

            if(solvedquartets.count(ets.src[2]) == 0)
                newtargets.insert(ets.src[2]);

            if(solvedquartets.count(ets.src[3]) == 0)
                newtargets.insert(ets.src[3]);
       
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

    ///////////////////////////////////////////////////
    // we need to do a dependency graph type of thing
    // This only really comes about when an initial
    // quartet depends on another initial quartet
    // But I might as well do the whole thing, just in case
    std::vector<QAM> amorder;

    // start with the initial targets
    // (well, their am)
    std::set<QAM> initamlist;
    for(const auto & it : inittargets)
        initamlist.insert(it.amlist());

    // add the dependencies for this am
    for(const auto & it : initamlist)
        ETAddWithDependencies_(amorder, it);

    // now create the final step list
    //cout << "AMORDER: " << amorder.size() << "\n";
    //for(const auto & ait : amorder)
    //    cout << "   " << ait[0] << " , " << ait[1] << " , " << ait[2] << " , " << ait[3] << "\n";

    //cout << "\nMAP: " << etslmap.size() << "\n";
    //for(const auto & mit : etslmap)
    //   cout << "   " << mit.first[0] << " , " << mit.first[1] << " , " << mit.first[2] << " , " << mit.first[3] << " Steps: " << mit.second.size() << "\n";

    cout << "\n";

    if(amorder.size() != etslmap.size())
        throw runtime_error("Error - inconsistent map and ordering sizes");

    for(const auto & ait : amorder)
        etsl.insert(etsl.end(), etslmap.at(ait).begin(), etslmap.at(ait).end());
}



void ET_Algorithm_Base::Create_ETStepList(const QuartetSet & inittargets)
{
    // holds all the 'solved' quartets
    QuartetSet solvedquartets;

    // generate initial targets
    QuartetSet targets = inittargets;
    PruneQuartets_(targets, ettop_);
    PrintQuartetSet(targets, "Initial ET Targets");

    // Solve
    ETStepLoop_(etsteps_, targets, solvedquartets, ettop_);

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
    {
        ettopam_.insert(it.amlist());
        ettopgauss_[it.bra.left.am()].insert(it.bra.left);
    }

    PrintQuartetSet(ettop_, "Top level ET");
    PrintGaussianMap(ettopgauss_, "ET Top level gaussian map");
}

QAMSet ET_Algorithm_Base::TopQAM(void) const
{
    return ettopam_;
}

QuartetSet ET_Algorithm_Base::TopQuartets(void) const
{
    return ettop_;
}


GaussianMap ET_Algorithm_Base::TopGaussians(void) const
{
    return ettopgauss_;
}

ETStepList ET_Algorithm_Base::ETSteps(void) const
{
    return etsteps_;
}
