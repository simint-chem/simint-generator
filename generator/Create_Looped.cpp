#include <iostream>
#include <memory>

#include "generator/Classes.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Boys.hpp"
#include "generator/Writer.hpp"
#include "generator/Helpers.hpp"

using namespace std;



static void HRRLoop(HRRDoubletStepList & hrrlist,
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
               

void Create_Looped(std::array<int, 4> amlist,
                   std::unique_ptr<HRR_Algorithm_Base> & hrralgo,
                   std::ostream & out)
{
    // read information about the boys function
    BoysMap bm = ReadBoysFitInfo("/home/ben/programming/simint/generator/dat");

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
    HRRLoop(bralist, targets, solvedbras, hrralgo);

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
    HRRLoop(ketlist, targets, solvedkets, hrralgo);

    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << "KET HRR step done. Solution is " << ketlist.size() << " steps\n";
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : ketlist)
        cout << it << "\n";


    // find top-level requirements at the end of HRR
    // these are src1 and src2 that are not solved
    // also set this in the hrrstep members
    DoubletSet topbras, topkets;

    for(auto & it : bralist)
    {
        if(solvedbras.count(it.src1) == 0)
        {
            it.src1.flags |= DOUBLET_HRRTOPLEVEL;
            topbras.insert(it.src1);
        }
        else
            it.src1.flags &= ~(DOUBLET_HRRTOPLEVEL);
            
        if(solvedbras.count(it.src2) == 0)
        {
            it.src2.flags |= DOUBLET_HRRTOPLEVEL;
            topbras.insert(it.src2);
        }
        else
            it.src2.flags &= ~(DOUBLET_HRRTOPLEVEL);
    }

    if(topbras.size() == 0)
        topbras.insert({DoubletType::BRA, {0,0,0}, {0,0,0}, DOUBLET_HRRTOPLEVEL});

    for(auto & it : ketlist)
    {
        if(solvedkets.count(it.src1) == 0)
        {
            it.src1.flags |= DOUBLET_HRRTOPLEVEL;
            topkets.insert(it.src1);
        }
        else
            it.src1.flags &= ~(DOUBLET_HRRTOPLEVEL);
            
        if(solvedkets.count(it.src2) == 0)
        {
            it.src2.flags |= DOUBLET_HRRTOPLEVEL;
            topkets.insert(it.src2);
        }
        else
            it.src2.flags &= ~(DOUBLET_HRRTOPLEVEL);
    }

    if(topkets.size() == 0)
        topkets.insert({DoubletType::KET, {0,0,0}, {0,0,0}, DOUBLET_HRRTOPLEVEL});

    ShellQuartetSet topreq;
    for(const auto & it : topbras)
    {
        for(const auto & it2 : topkets)
        {
            topreq.insert({it, it2});
        }
    }

    // reverse the steps
    std::reverse(bralist.begin(), bralist.end());
    std::reverse(ketlist.begin(), ketlist.end());

    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << " HRR TOP LEVEL REQ\n";
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : topreq)
        cout << "    " << it << "\n"; 

    
    cout << "\n\n";

    // write out

    HRRDoubletStepInfo hrrinfo{bralist, ketlist, topreq};

    Writer_Looped(out,
                  amlist,
                  "FOcombined",
                  bm,
                  VRRInfo{2, {}},
                  ETInfo(),
                  hrrinfo);

}
