#include <iostream>
#include <memory>

#include "Classes.h"
#include "Algorithms.h"
#include "Helpers.h"

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
    // generate initial targets
    QuartetSet targets = GenerateInitialTargets({0,2,0,0});

    cout << "Initial targets: " << targets.size() << "\n";
    for(auto & it : targets)
        cout << "    " << it << "\n";
    cout << "\n";


    QuartetSet solvedquartets;

    std::unique_ptr<HRR_Algorithm_Base> hrralgo(new Makowski_HRR);
    HRRList hrrlist;

    // bras
    HRRLoop(hrrlist, targets, solvedquartets, DoubletType::BRA, hrralgo);

    // now do kets
    // targets are what has been "solved" by the bra step, pruned on the ket side
    targets = solvedquartets;
    PruneRight(targets, DoubletType::KET);
    HRRLoop(hrrlist, targets, solvedquartets, DoubletType::KET, hrralgo);

    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << "HRR step done. Solution is " << hrrlist.size() << " steps\n";

    // reverse
    std::reverse(hrrlist.begin(), hrrlist.end());
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : hrrlist)
        cout << it << "\n";
    cout << "\n\n";
    cout << "--------------------------------------------------------------------------------\n";
    cout << " CODE\n";
    cout << "--------------------------------------------------------------------------------\n";
    for(auto & it : hrrlist)
        cout << it.code_line() << "\n";
    cout << "--------------------------------------------------------------------------------\n";

    

    return 0;
}
