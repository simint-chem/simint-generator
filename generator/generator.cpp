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

    if(argc != 8)
    {
        std::cout << "usage: generator boystype prefix am1 am2 am3 am4 outfile\n";
        return 1;
    }

    std::string boystype(argv[1]);
    std::string prefix(argv[2]);

    QAMList amlist{
                    atoi(argv[3]),
                    atoi(argv[4]),
                    atoi(argv[5]),
                    atoi(argv[6])
                  };

    std::string fpath(argv[7]);

    // Read in the boys map
    std::unique_ptr<BoysGen> bg;

    if(boystype == "FO")
        bg = std::unique_ptr<BoysGen>(new BoysFO("/home/ben/programming/simint/generator/dat"));
    else if(boystype == "split")
        bg = std::unique_ptr<BoysGen>(new BoysSplit());
    else
    {
        std::cout << "Unknown boys type \"" << boystype << "\"\n";
        return 3;
    }

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

    std::ofstream of(fpath);
    if(!of.is_open())
    {
        std::cout << "Cannot open file: " << fpath << "\n";
        return 2; 
    }

    Writer_Looped(of, amlist, prefix, *bg, *vrralgo, *etalgo, *hrralgo);


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
