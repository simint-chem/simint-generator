#include <stdexcept>
#include <fstream>

#include "generator/CommandLine.hpp"
#include "generator/StringBuilder.hpp"

#include "generator/ostei/Algorithms.hpp"
#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_VRR_Writer.hpp"
#include "generator/ostei/OSTEI_HRR_Writer.hpp"
#include "generator/ostei/OSTEI_Writer.hpp"
#include "generator/Printing.hpp"

int main(int argc, char ** argv)
{
    try {


    // other stuff
    std::string fpath;
    std::string hpath;
    QAM finalam{0,0,0,0};

    bool finalamset = false;

    // parse command line
    OptionMap options = DefaultOptions();
    std::vector<std::string> otheropt = ParseCommonOptions(options, argc, argv);

    // parse specific options
    size_t iarg = 0;
    while(iarg < otheropt.size())
    {
        std::string argstr(GetNextArg(iarg, otheropt));
        if(argstr == "-o")
            fpath = GetNextArg(iarg, otheropt);
        else if(argstr == "-oh")
            hpath = GetNextArg(iarg, otheropt);
        else if(argstr == "-q")
        {
            finalam[0] = GetIArg(iarg, otheropt);   
            finalam[1] = GetIArg(iarg, otheropt);   
            finalam[2] = GetIArg(iarg, otheropt);   
            finalam[3] = GetIArg(iarg, otheropt);   
            finalamset = true;
        }
        else
        {
            std::cout << "\n\n";
            std::cout << "--------------------------------\n";
            std::cout << "Unknown argument: " << argstr << "\n";
            std::cout << "--------------------------------\n";
            return 1; 
        } 
    }


    // check for required options
    CMDLINE_ASSERT( fpath != "", "output source file path (-o) required" )
    CMDLINE_ASSERT( hpath != "", "output header file path (-oh) required" )
    CMDLINE_ASSERT( finalamset == true, "AM quartet (-q) required" )


    // open the output file
    // actually create
    std::cout << "Generating source file " << fpath << "\n";
    std::cout << "Appending to header file " << hpath << "\n";

    std::ofstream of(fpath);
    if(!of.is_open())
        throw std::runtime_error(StringBuilder("Cannot open file: ", fpath, "\n"));

    std::ofstream ofh(hpath, std::ofstream::app);
    if(!ofh.is_open())
        throw std::runtime_error(StringBuilder("Cannot open file: ", hpath, "\n"));
    

    // Information for this OSTEI
    OSTEI_GeneratorInfo info(finalam, 1, options);


    //////////////////////////////////////////////////////////////
    //! \todo Are there special permutations for derivatives?
    //////////////////////////////////////////////////////////////

    // algorithms used
    Makowski_HRR hrralgo(info);
    Makowski_VRR vrralgo(info);

    // Working backwards, I need:
    // 1.) HRR Steps
    //     We need the different increments, decrements for derivatives
    std::set<QAM> needed_am;

    // We only need to do three centers. The fourth is free
    // Find which one to skip
    int max = *(std::max_element(finalam.begin(), finalam.end()));
    
    int missing_center;
    if(finalam[0] == max)
        missing_center = 0;
    else if(finalam[1] == max)
        missing_center = 1;
    else if(finalam[2] == max)
        missing_center = 2;
    else
        missing_center = 3;
    
    const char * dir[4] = {"2a", "2b", "2c", "2d"};
    for(int i = 0; i < 4; i++)
    {
        if(i == missing_center)
            continue;

        QAM amtmp_p(finalam.qam, dir[i]);
        QAM amtmp_m(finalam);
        amtmp_p.qam[i]++;
        amtmp_m.qam[i]--;

        needed_am.insert(amtmp_p);
        if(ValidQAM(amtmp_m))
            needed_am.insert(amtmp_m);
    }

    info.SetDeriv1_MissingCenter(missing_center);

    hrralgo.Create(needed_am);
    OSTEI_HRR_Writer hrr_writer(hrralgo, info,
                                options[Option::ExternalHRR],
                                options[Option::GeneralHRR]);

    QuartetSet topquartets = hrralgo.TopQuartets();

    // 2.) VRR Steps
    vrralgo.Create(topquartets);
    OSTEI_VRR_Writer vrr_writer(vrralgo, info,
                                options[Option::ExternalVRR],
                                options[Option::GeneralVRR]);


    // Create the OSTEI_Writer and write the file
    OSTEIDeriv1_Writer ostei_deriv1_writer(of, ofh, info, vrr_writer, hrr_writer);
    ostei_deriv1_writer.WriteFile();



    }
    catch(std::exception & ex)
    {
        std::cout << "\n\n";
        std::cout << "Caught exception\n";
        std::cout << "What = " << ex.what() << "\n\n";
        return 1;
    }
    return 0;
}
