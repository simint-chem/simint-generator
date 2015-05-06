#include <iostream>
#include <algorithm>

#include "generator/AlgorithmBase.hpp"
#include "generator/Helpers.hpp"

using namespace std;

void HRR_Algorithm_Base::HRRDoubletLoop(HRRDoubletStepList & hrrlist,
                                        const DoubletSet & inittargets,
                                        DoubletSet & solveddoublets)
{
    DoubletSet targets = inittargets;

    while(targets.size())
    {
        DoubletSet newtargets;

        for(auto it = targets.rbegin(); it != targets.rend(); ++it)
        {
            HRRDoubletStep hrr = this->doubletstep(*it);
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


void HRR_Algorithm_Base::HRRQuartetLoop(HRRQuartetStepList & hrrlist,
                                        const QuartetSet & inittargets,
                                        QuartetSet & solvedquartets,
                                        DoubletType type)
{
    QuartetSet targets = inittargets;

    while(targets.size())
    {
        QuartetSet newtargets;

        for(auto it = targets.rbegin(); it != targets.rend(); ++it)
        {
            HRRQuartetStep hrr = this->quartetstep(*it, type);
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
               

HRRQuartetStep HRR_Algorithm_Base::quartetstep(const Quartet & target, DoubletType steptype)
{
    // Call doublet step, then create the Quartet 

    Doublet d = target.get(steptype);

    // call doubletstep
    HRRDoubletStep hds = this->doubletstep(d);

    // Create the HRR step
    // (this preserves the flags of target
    if(steptype == DoubletType::BRA)
    {
        HRRQuartetStep hrr{target, 
                          {hds.src1, target.ket, target.m, 0},   // src1 quartet
                          {hds.src2, target.ket, target.m, 0},   // src2 quartet
                          steptype,                                 // bra or ket being stepped
                          hds.xyz};                                 // cartesian direction being stepped
        return hrr;
    }
    else
    {
        HRRQuartetStep hrr{target, 
                          {target.bra, hds.src1, target.m, 0},   // src1 quartet
                          {target.bra, hds.src2, target.m, 0},   // src2 quartet
                          steptype,                                 // bra or ket being stepped
                          hds.xyz};                                 // cartesian direction being stepped
        return hrr;
    }
}


std::pair<VRRMap, VRRReqMap> VRR_Algorithm_Base::CreateAllMaps(const GaussianSet & greq)
{
    // holds the requirements for each am
    VRRReqMap vrm;

    // find the requested gaussians for each am
    for(const auto & it : greq)
        vrm[it.am()].insert(it); 

    // max am
    int maxam = vrm.rbegin()->first;

    // holds the VRR steps
    VRRMap vm;

    for(int i = 0; i <= maxam; i++)
    {
        VRRMap vm2 = CreateVRRMap(i);
        vm.insert(vm2.begin(), vm2.end());
    }

    // recurse down from the end
    for(int i = maxam-1; i > 0; i--)
    {
        // get the set from the previous am
        GaussianSet prev = vrm[i+1];
        
        // for each one, look up its requirements
        for(const auto & it : prev)
        {
            XYZStep s = vm[it];

            // add what it needs to this step
            Gaussian g1 = it.StepDown(s, 1);
            Gaussian g2 = it.StepDown(s, 2);

            // if these are valid
            if(g1)
                vrm[g1.am()].insert(g1);
            if(g2)
                vrm[g2.am()].insert(g2);
        }
    }

    return {vm, vrm}; 
}



HRRBraKetStepList HRR_Algorithm_Base::Create_DoubletStepLists(QAMList amlist)
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
    HRRDoubletLoop(bralist, targets, solvedbras);
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
    HRRDoubletLoop(ketlist, targets, solvedkets);
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


HRRQuartetStepList HRR_Algorithm_Base::Create_QuartetStepList(QAMList amlist)
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
    HRRQuartetLoop(hrrlist, targets, solvedquartets, DoubletType::BRA);

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
    HRRQuartetLoop(hrrlist, targets, solvedquartets, DoubletType::KET);

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





