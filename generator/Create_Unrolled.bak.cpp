#include <iostream>
#include <memory>

#include "generator/Classes.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Boys.hpp"
#include "generator/Writer.hpp"
#include "generator/Helpers.hpp"

using namespace std;



HRRQuartetStepList Create_QuartetStepListsUnrolled(QAMList amlist,
                                                   std::unique_ptr<HRR_Algorithm_Base> & hrralgo,
                                                   std::ostream & out)
{
    // all HRR steps
    HRRQuartetStepList hrrlist;

    // holds all the 'solved' quartets
    QuartetSet solvedquartets;

    // generate initial targets
    QuartetSet inittargets = GenerateInitialQuartetTargets(amlist, true);
    PrintQuartetSet(inittargets, "Initial Targets");

    // Inital bra targets
    QuartetSet targets = inittargets;
    PruneRight(targets, DoubletType::BRA);
    PrintQuartetSet(targets, "Inital bra targets");

    // Solve the bra part
    HRRLoop(hrrlist, targets, solvedquartets, DoubletType::BRA, hrralgo);

    // now do kets
    // targets are src1 and src2 of hrrlist, pruned on the ket side
    // and include any init targets that weren't solved
    targets.clear();
    for(const auto & it : hrrlist)
    {
        if(solvedquartets.count(it.src1) == 0)
            targets.insert(it.src1);
        if(solvedquartets.count(it.src2) == 0)
            targets.insert(it.src2);
    }
    
    if(targets.size() == 0)
        targets = inittargets; // might be some there? ie  ( ss | ps )

    PruneRight(targets, DoubletType::KET);

    cout << "\n\n";
    PrintQuartetSet(targets, "Initial ket targets");

    // Solve the ket part
    HRRLoop(hrrlist, targets, solvedquartets, DoubletType::KET, hrralgo);

    // reverse the steps
    std::reverse(hrrlist.begin(), hrrlist.end());

    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << "HRR step done. Solution is " << hrrlist.size() << " steps\n";
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : hrrlist)
        cout << it << "\n";

    cout << "\n\n";

    return hrrlist;
}
