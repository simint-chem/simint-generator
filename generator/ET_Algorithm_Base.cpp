#include <iostream>

#include "generator/ET_Algorithm_Base.hpp"

using namespace std;

void ET_Algorithm_Base::ETAddWithDependencies(std::vector<QAM> & amorder, QAM am)
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
        ETAddWithDependencies(amorder, it);

    // add am to amorder
    // If it exists, then it was a dependency of itself. Uh oh
    if(std::find(amorder.begin(), amorder.end(), am) != amorder.end()) 
        throw std::runtime_error("Circular dependency in electron transfer tree");

    amorder.push_back(am);
}


void ET_Algorithm_Base::ETStepLoop(ETStepList & etsl,
                                   const QuartetSet & inittargets,
                                   QuartetSet & solvedquartets)
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
            if(solvedquartets.count(it) > 0)
                continue;

            ETStep ets = this->etstep(*it);

            // add to the map
            etslmap[ets.target.amlist()].push_back(ets);

            if(solvedquartets.count(ets.src1) == 0)
                newtargets.insert(ets.src1);

            if(solvedquartets.count(ets.src2) == 0)
                newtargets.insert(ets.src2);

            if(solvedquartets.count(ets.src3) == 0)
                newtargets.insert(ets.src3);

            if(solvedquartets.count(ets.src4) == 0)
                newtargets.insert(ets.src4);
       
            // insert copy with no flags 
            solvedquartets.insert(it);
        }

        cout << "Generated " << newtargets.size() << " new targets\n";

        PruneET(newtargets);

        cout << "After pruning: " << newtargets.size() << " new targets\n";
        for(const auto & it : newtargets)
            cout << "    " << it << "\n";

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
        ETAddWithDependencies(amorder, it);

    // now create the final step list
    cout << "AMORDER: " << amorder.size() << "\n";
    for(const auto & ait : amorder)
        cout << "   " << ait[0] << " , " << ait[1] << " , " << ait[2] << " , " << ait[3] << "\n";

    cout << "\nMAP: " << etslmap.size() << "\n";
    for(const auto & mit : etslmap)
        cout << "   " << mit.first[0] << " , " << mit.first[1] << " , " << mit.first[2] << " , " << mit.first[3] << " Steps: " << mit.second.size() << "\n";
    cout << "\n";

    if(amorder.size() != etslmap.size())
        throw runtime_error("Error - inconsistent map and ordering sizes");

    for(const auto & ait : amorder)
        etsl.insert(etsl.end(), etslmap.at(ait).begin(), etslmap.at(ait).end());
}



ETStepList ET_Algorithm_Base::Create_ETStepList(const QuartetSet & inittargets)
{
    ETStepList etsl;

    // holds all the 'solved' quartets
    QuartetSet solvedquartets;

    // generate initial targets
    QuartetSet targets = inittargets;
    PruneET(targets);
    PrintQuartetSet(targets, "Initial ET Targets");

    // Solve
    ETStepLoop(etsl, targets, solvedquartets);

    // reverse the steps
    // (now done in ETStepLoop)
    //std::reverse(etsl.begin(), etsl.end());

    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << "ET step done. Solution is " << etsl.size() << " steps\n";
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : etsl)
        cout << it << "\n";

    cout << "\n\n";

    // mark top level reqs
    for(const auto & it : etsl)
    {
        if(it.src1 && solvedquartets.count(it.src1) == 0)
            ettop_.insert(it.src1);
        if(it.src2 && solvedquartets.count(it.src2) == 0)
            ettop_.insert(it.src2);
        if(it.src3 && solvedquartets.count(it.src3) == 0)
            ettop_.insert(it.src3);
        if(it.src4 && solvedquartets.count(it.src4) == 0)
            ettop_.insert(it.src4);
    } 

    for(const auto & it : ettop_)
        ettopam_.insert(it.amlist());

    return etsl;
}

QAMSet ET_Algorithm_Base::TopAM(Void) const
{
    return ettopam_;
}
