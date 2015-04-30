#include <iostream>
#include <memory>
#include <stdexcept>

#include "generator/Classes.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Boys.hpp"
#include "generator/Create.hpp"
#include "generator/Writer.hpp"

using namespace std;


int main(int argc, char ** argv)
{
    try {

    if(argc != 5)
    {
        printf("Need 4 arguments\n");
        return 1;
    }

    QAMList amlist{
                    atoi(argv[1]),
                    atoi(argv[2]),
                    atoi(argv[3]),
                    atoi(argv[4])
                  };

    const int L = amlist[0] + amlist[1] + amlist[2] + amlist[3] + 1;

    // Read in the boys map
    BoysMap bm = ReadBoysFitInfo("/home/ben/programming/simint/generator/dat");

    // algorithms used
    //std::unique_ptr<HRR_Algorithm_Base> hrralgo(new Makowski_HRR);
    std::unique_ptr<VRR_Algorithm_Base> vrralgo(new Makowski_VRR);
    VRRMap vm = vrralgo->CreateAllMaps(L);

    for(const auto & it : vm)
        cout << "    " << it.first << "  <-  " << it.second << "\n";

    /*
    // Create the quartet list
    HRRQuartetStepList hqsl = Create_QuartetStepList(amlist, hrralgo);

    // Write it out
    Writer_Unrolled(cout, amlist, "FOCombined", bm, hqsl);
    */

    // Create the step lists
    //HRRBraKetStepList bksl = Create_DoubletStepLists(amlist, hrralgo);

    // Write it out
    //Writer_Looped(cout, amlist, "FOCombined", bm, bksl);


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
