#include <iostream>
#include <memory>

#include "generator/Classes.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Boys.hpp"
#include "generator/Writer.hpp"
#include "generator/Helpers.hpp"

using namespace std;


static void HRRLoop(HRRQuartetStepList & hrrlist,
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
               

void Create_Unrolled(std::array<int, 4> amlist,
                     std::unique_ptr<HRR_Algorithm_Base> & hrralgo,
                     std::ostream & out)
{
    // read information about the boys function
    BoysMap bm = ReadBoysFitInfo("/home/ben/programming/simint/generator/dat");

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
    PrintQuartetSet(targets, "Initial ket targets");

    // Solve the ket part
    HRRLoop(hrrlist, targets, solvedquartets, DoubletType::KET, hrralgo);

    // find top-level requirements at the end of HRR
    // these are src1 and src2 that are not solved
    // also set this in the hrrstep members
    QuartetSet topreq;

    for(auto & it : hrrlist)
    {
        if(solvedquartets.count(it.src1) == 0)
        {
            topreq.insert(it.src1);
            it.src1.flags |= QUARTET_HRRTOPLEVEL;
        }
        else
            it.src1.flags &= ~(QUARTET_HRRTOPLEVEL);
            
        if(solvedquartets.count(it.src2) == 0)
        {
            topreq.insert(it.src2);
            it.src2.flags |= QUARTET_HRRTOPLEVEL;
        }
        else
            it.src2.flags &= ~(QUARTET_HRRTOPLEVEL);
    }


    // reverse the steps
    std::reverse(hrrlist.begin(), hrrlist.end());

    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << "HRR step done. Solution is " << hrrlist.size() << " steps\n";
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : hrrlist)
        cout << it << "\n";

    /*
    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << " CODE\n";
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : hrrlist)
        cout << it.code_line() << "\n";
    */

    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << " HRR TOP LEVEL REQ\n";
    cout << "--------------------------------------------------------------------------------\n";
    PrintQuartetSet(topreq, "Top level req");

    
    cout << "\n\n";

    // write out
    HRRQuartetStepInfo hrrinfo{hrrlist, topreq};

    Writer_Unrolled(out,
                    amlist,
                    "FOcombined",
                    bm,
                    VRRInfo{2, {}},
                    ETInfo(),
                    hrrinfo);
}
