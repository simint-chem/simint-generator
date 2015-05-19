#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>

#include "generator/Classes.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Boys.hpp"
#include "generator/Writer.hpp"

using namespace std;


int main(int argc, char ** argv)
{
    try {

    if(argc != 6)
    {
        std::cout << "usage: generator am1 am2 am3 am4 outfile\n";
        return 1;
    }

    QAMList amlist{
                    atoi(argv[1]),
                    atoi(argv[2]),
                    atoi(argv[3]),
                    atoi(argv[4])
                  };

    // Read in the boys map
    BoysMap bm = ReadBoysFitInfo("/home/ben/programming/simint/generator/dat");

    // algorithms used
    std::unique_ptr<HRR_Algorithm_Base> hrralgo(new Makowski_HRR);
    std::unique_ptr<VRR_Algorithm_Base> vrralgo(new Makowski_VRR);
    std::unique_ptr<ET_Algorithm_Base> etalgo(new Makowski_ET);

    /*
    // Create the quartet list
    HRRQuartetStepList hqsl = hrralgo->Create_QuartetStepList(amlist);

    // Write it out
    Writer_Unrolled(cout, amlist, "FOcombined", bm, hqsl);
    */

    // Create/Write

    std::ofstream of(argv[5]);
    if(!of.is_open())
    {
        std::cout << "Cannot open file: " << argv[5] << "\n";
        return 2; 
    }

    Writer_Looped(of, amlist, "FO", bm, *vrralgo, *etalgo, *hrralgo);


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
