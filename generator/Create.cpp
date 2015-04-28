#include <iostream>
#include <memory>

#include "generator/Classes.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Boys.hpp"
#include "generator/Writer.hpp"
#include "generator/Helpers.hpp"

using namespace std;


static void HRRDoubletLoop(HRRDoubletStepList & hrrlist,
                           const DoubletSet & inittargets,
                           DoubletSet & solveddoublets,
                           std::unique_ptr<HRR_Algorithm_Base> & hrralgo)
{
    DoubletSet targets = inittargets;

    while(targets.size())
    {
        DoubletSet newtargets;

        for(auto it = targets.rbegin(); it != targets.rend(); ++it)
        {
            HRRDoubletStep hrr = hrralgo->doubletstep(*it);
            hrrlist.push_back(hrr);

            if(solveddoublets.count(hrr.src1) == 0)
                newtargets.insert(hrr.src1);
            if(solveddoublets.count(hrr.src2) == 0)
                newtargets.insert(hrr.src2);
            
            solveddoublets.insert(*it);
        }

        cout << "Generated " << newtargets.size() << " new targets\n";

        PruneRight(newtargets);
        cout << "After pruning: " << newtargets.size() << " new targets\n";
        for(const auto & it : newtargets)
            cout << "    " << it << "\n";

        targets = newtargets;
    } 
}


static void HRRQuartetLoop(HRRQuartetStepList & hrrlist,
                           const QuartetSet & inittargets,
                           QuartetSet & solvedquartets,
                           DoubletType type,
                           std::unique_ptr<HRR_Algorithm_Base> & hrralgo)
{
    QuartetSet targets = inittargets;

    while(targets.size())
    {
        QuartetSet newtargets;

        for(auto it = targets.rbegin(); it != targets.rend(); ++it)
        {
            HRRQuartetStep hrr = hrralgo->quartetstep(*it, type);
            hrrlist.push_back(hrr);

            if(solvedquartets.count(hrr.src1) == 0)
                newtargets.insert(hrr.src1);
            if(solvedquartets.count(hrr.src2) == 0)
                newtargets.insert(hrr.src2);
            
            solvedquartets.insert(*it);
        }

        cout << "Generated " << newtargets.size() << " new targets\n";

        PruneRight(newtargets, type);
        cout << "After pruning: " << newtargets.size() << " new targets\n";
        for(const auto & it : newtargets)
            cout << "    " << it << "\n";

        targets = newtargets;
    } 
}
               


HRRBraKetStepList Create_DoubletStepLists(QAMList amlist, std::unique_ptr<HRR_Algorithm_Base> & hrralgo)
{
    // First, we need a list of doublet steps for the bra
    HRRDoubletStepList bralist;

    // holds all the 'solved' doublets
    DoubletSet solvedbras;

    // generate initial targets
    DoubletSet initbras = GenerateInitialDoubletTargets({amlist[0],amlist[1]}, DoubletType::BRA, true);
    PrintDoubletSet(initbras, "Initial Targets");

    // Inital bra targets
    DoubletSet targets = initbras;
    PruneRight(targets);
    PrintDoubletSet(targets, "Inital bra targets");

    // Solve the bra part
    HRRDoubletLoop(bralist, targets, solvedbras, hrralgo);
    std::reverse(bralist.begin(), bralist.end());

    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << "BRA HRR step done. Solution is " << bralist.size() << " steps\n";
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : bralist)
        cout << it << "\n";
    

    // now do kets
    // we only need the bras from the original targets
    HRRDoubletStepList ketlist;
    DoubletSet solvedkets;
    DoubletSet initkets = GenerateInitialDoubletTargets({amlist[2],amlist[3]}, DoubletType::KET, true);
    targets = initkets;
    PruneRight(targets);

    cout << "\n\n";
    PrintDoubletSet(targets, "Initial ket targets");

    // Solve the ket part
    HRRDoubletLoop(ketlist, targets, solvedkets, hrralgo);
    std::reverse(ketlist.begin(), ketlist.end());

    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << "KET HRR step done. Solution is " << ketlist.size() << " steps\n";
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : ketlist)
        cout << it << "\n";

    cout << "\n\n";

    // mark top level reqs
    for(auto & it : bralist)
    {
        if(solvedbras.count(it.src1) == 0)
            it.src1.flags |= DOUBLET_HRRTOPLEVEL;
        if(solvedbras.count(it.src2) == 0)
            it.src2.flags |= DOUBLET_HRRTOPLEVEL;
    } 
    for(auto & it : ketlist)
    {
        if(solvedkets.count(it.src1) == 0)
            it.src1.flags |= DOUBLET_HRRTOPLEVEL;
        if(solvedkets.count(it.src2) == 0)
            it.src2.flags |= DOUBLET_HRRTOPLEVEL;
    } 

    return HRRBraKetStepList(bralist, ketlist);
}


HRRQuartetStepList Create_QuartetStepList(QAMList amlist,
                                          std::unique_ptr<HRR_Algorithm_Base> & hrralgo)
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
    HRRQuartetLoop(hrrlist, targets, solvedquartets, DoubletType::BRA, hrralgo);

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
    HRRQuartetLoop(hrrlist, targets, solvedquartets, DoubletType::KET, hrralgo);

    // reverse the steps
    std::reverse(hrrlist.begin(), hrrlist.end());

    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << "HRR step done. Solution is " << hrrlist.size() << " steps\n";
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : hrrlist)
        cout << it << "\n";

    // mark top level reqs
    for(auto & it : hrrlist)
    {
        if(solvedquartets.count(it.src1) == 0)
            it.src1.flags |= QUARTET_HRRTOPLEVEL;
        if(solvedquartets.count(it.src2) == 0)
            it.src2.flags |= QUARTET_HRRTOPLEVEL;
    } 

    cout << "\n\n";

    return hrrlist;
}
