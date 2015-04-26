#include <iostream>
#include <memory>
#include <stdexcept>

#include "generator/Classes.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Boys.hpp"
#include "generator/Writer.hpp"

using namespace std;


QuartetSet GenerateInitialTargets(std::array<int, 4> amlst)
{
    QuartetSet qs;
    int nam1 = ((amlst[0] + 1) * (amlst[0] + 2)) / 2;
    int nam2 = ((amlst[1] + 1) * (amlst[1] + 2)) / 2;
    int nam3 = ((amlst[2] + 1) * (amlst[2] + 2)) / 2;
    int nam4 = ((amlst[3] + 1) * (amlst[3] + 2)) / 2;

    Gaussian cur1 = Gaussian{amlst[0], 0, 0};
    int ijkl = 0;
    for(int i = 0; i < nam1; i++)
    {
        Gaussian cur2 = Gaussian{amlst[1], 0, 0};
        for(int j = 0; j < nam2; j++)
        {
            Doublet bra{DoubletType::BRA, cur1, cur2};
            Gaussian cur3 = Gaussian{amlst[2], 0, 0};
            for(int k = 0; k < nam3; k++)
            {
                Gaussian cur4 = Gaussian{amlst[3], 0, 0};
                for(int l = 0; l < nam4; l++)
                {
                    Doublet ket{DoubletType::KET, cur3, cur4};
                    qs.insert(Quartet{bra, ket, 0, QUARTET_INITIAL, ijkl});
                    cur4.Iterate();
                    ijkl++;
                }

                cur3.Iterate();
            }

            cur2.Iterate(); 
        }

        cur1.Iterate();
    }

    return qs;
}




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


void HRRLoop(HRRStepList & hrrlist,
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
            HRRStep hrr = hrralgo->step(*it, type);
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
               

void PrintQuartetSet(const QuartetSet & q,
                     const std::string & title)
{
    cout << title << ": " << q.size() << "\n";
    for(auto & it : q)
        cout << "    " << it << "\n";
    cout << "\n";
}


int main(void)
{
    try {

    // read information about the boys function
    BoysMap bm = ReadBoysFitInfo("/home/ben/programming/simint/generator/dat");

    // all HRR steps
    HRRStepList hrrlist;

    // holds all the 'solved' quartets
    QuartetSet solvedquartets;

    // The algorithm we are using
    std::unique_ptr<HRR_Algorithm_Base> hrralgo(new Makowski_HRR);


    // generate initial targets
    QuartetSet inittargets = GenerateInitialTargets({1,1,1,1});
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
    for(auto & it : topreq)
        cout << "    " << it << "\n"; 

    
    cout << "\n\n";

    // write out
    HRRInfo hrrinfo{hrrlist, topreq};

    Write_Generic(cout,
                  {1,1,1,1},
                  "FOcombined",
                  bm,
                  VRRInfo{2, {}},
                  ETInfo(),
                  hrrinfo);


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
