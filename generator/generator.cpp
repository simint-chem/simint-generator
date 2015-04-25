#include <iostream>
#include <memory>
#include <stdexcept>

#include "generator/Classes.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Helpers.hpp"
#include "generator/Boys.hpp"

using namespace std;

void PruneRight(QuartetSet & qs, DoubletType type)
{
    QuartetSet qsnew;

    for(auto & it : qs)
    {
        if(it.get(type).right.am() != 0)
            qsnew.insert(it);
    }

    qs = qsnew; 
}


void HRRLoop(HRRList & hrrlist,
             const QuartetSet & inittargets,
             QuartetSet & solvedquartets,
             DoubletType type,
             std::unique_ptr<HRR_Algorithm_Base> & hrralgo)
{
    QuartetSet targets = inittargets;

    while(targets.size())
    {
        QuartetSet newtargets;

        for(const auto & it : targets)
        {
            HRRStep hrr = hrralgo->step(it, type);
            hrrlist.push_back(hrr);

            if(solvedquartets.count(hrr.src1) == 0)
                newtargets.insert(hrr.src1);
            if(solvedquartets.count(hrr.src2) == 0)
                newtargets.insert(hrr.src2);
            
            solvedquartets.insert(it);
        }

        cout << "Generated " << newtargets.size() << " new targets\n";

        PruneRight(newtargets, type);
        cout << "After pruning: " << newtargets.size() << " new targets\n";
        for(const auto & it : newtargets)
            cout << "    " << it << "\n";

        targets = newtargets;
    } 
}
               


int main(void)
{
    try {

    // read information about the boys function
    BoysMap bm = ReadBoysInfo("/home/ben/programming/simint/generator/dat");



    // generate initial targets
    QuartetSet inittargets = GenerateInitialTargets({0,0,0,1});

    cout << "Initial targets: " << inittargets.size() << "\n";
    for(auto & it : inittargets)
        cout << "    " << it << "\n";
    cout << "\n";


    QuartetSet solvedquartets;

    std::unique_ptr<HRR_Algorithm_Base> hrralgo(new Makowski_HRR);
    HRRList hrrlist;

    // bras
    QuartetSet targets = inittargets;
    PruneRight(targets, DoubletType::BRA);

    cout << "BRA Initial targets: " << targets.size() << "\n";
    for(auto & it : targets)
        cout << "    " << it << "\n";
    cout << "\n";

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


    cout << "KET Initial targets: " << targets.size() << "\n";
    for(auto & it : targets)
        cout << "    " << it << "\n";
    cout << "\n";

    HRRLoop(hrrlist, targets, solvedquartets, DoubletType::KET, hrralgo);

    // find top-level requirements
    // these are src1 and src2 that are not solved
    QuartetSet topreq;

    for(const auto & it : hrrlist)
    {
        if(solvedquartets.count(it.src1) == 0)
            topreq.insert(it.src1);
        if(solvedquartets.count(it.src2) == 0)
            topreq.insert(it.src2);
    }



    // reverse
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
    cout << " TOP LEVEL REQ\n";
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : topreq)
        cout << "    " << it << "\n"; 

    
    cout << "\n\n";


    }
    catch(std::exception & ex)
    {
        cout << "\n\n";
        cout << "Caught exception\n";
        cout << "What = " << ex.what() << "\n\n";
        return 100;
    }
    return 0;
}
