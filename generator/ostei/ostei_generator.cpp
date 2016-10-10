#include <stdexcept>
#include <fstream>

#include "generator/CommandLine.hpp"

#include "generator/ostei/Algorithms.hpp"
#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_VRR_Writer.hpp"
#include "generator/ostei/OSTEI_ET_Writer.hpp"
#include "generator/ostei/OSTEI_HRR_Writer.hpp"
#include "generator/ostei/OSTEI_Writer.hpp"


int main(int argc, char ** argv)
{
    try {


    // other stuff
    std::string fpath;
    std::string hpath;
    std::string cpuflags;
    QAM finalam;

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
        else if(argstr == "-c")
            cpuflags = GetNextArg(iarg, otheropt);
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
    OSTEI_GeneratorInfo info(finalam,
                          Compiler::Intel,
                          cpuflags, options);


    //////////////////////////////////////////////////////////////
    //! \todo We are doing all this work even if it is a special
    //        permutation
    //////////////////////////////////////////////////////////////

    // algorithms used
    std::unique_ptr<OSTEI_HRR_Algorithm_Base> hrralgo(new Makowski_HRR(options));
    std::unique_ptr<OSTEI_VRR_Algorithm_Base> vrralgo(new Makowski_VRR(options));
    std::unique_ptr<OSTEI_ET_Algorithm_Base> etalgo(new Makowski_ET(options));

    // Writers
    std::unique_ptr<OSTEI_HRR_Writer> hrr_writer;
    std::unique_ptr<OSTEI_VRR_Writer> vrr_writer;
    std::unique_ptr<OSTEI_ET_Writer> et_writer;

    // Working backwards, I need:
    // 1.) HRR Steps
    hrralgo->Create(finalam);
    if(options[Option::InlineHRR])
        hrr_writer = std::unique_ptr<OSTEI_HRR_Writer>(new OSTEI_HRR_Writer_Inline(*hrralgo, info));
    else
        hrr_writer = std::unique_ptr<OSTEI_HRR_Writer>(new OSTEI_HRR_Writer_External(*hrralgo, info));

    QuartetSet topquartets = hrralgo->TopQuartets();

    // 2.) ET steps
    //     with the HRR top level stuff as the initial targets
    if(!options[Option::NoET])
    {
        DoubletType etdirection = DoubletType::KET;
        if((finalam[2]+finalam[3]) > (finalam[0]+finalam[1]))
            etdirection = DoubletType::BRA;

        etalgo->Create(topquartets, etdirection);
        topquartets = etalgo->TopQuartets();
    }


    if(options[Option::InlineET])
        et_writer = std::unique_ptr<OSTEI_ET_Writer>(new OSTEI_ET_Writer_Inline(*etalgo, info));
    else
        et_writer = std::unique_ptr<OSTEI_ET_Writer>(new OSTEI_ET_Writer_External(*etalgo, info));


    // 3.) VRR Steps
    // requirements for vrr are the top level stuff from ET
    vrralgo->Create(topquartets);

    if(options[Option::InlineVRR])
        vrr_writer = std::unique_ptr<OSTEI_VRR_Writer>(new OSTEI_VRR_Writer_Inline(*vrralgo, info));
    else
        vrr_writer = std::unique_ptr<OSTEI_VRR_Writer>(new OSTEI_VRR_Writer_External(*vrralgo, info));


    // set the contracted quartets
    info.SetContQ(hrralgo->TopAM());

    // Create the OSTEI_Writer and write the file
    OSTEI_Writer_Basic ostei_writer(of, ofh, info, *vrr_writer, *et_writer, *hrr_writer);
    ostei_writer.WriteFile();


    // For information
    std::cout << "\nCONTWORK SIZE: " << info.ContNElements() << "  " << info.ContMemoryReq() << "\n";

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
