#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>

#include "generator/CommandLine.hpp"
#include "generator/Algorithms.hpp"

#include "generator/ERIGeneratorInfo.hpp"
#include "generator/Boys.hpp"
#include "generator/VRR_Writer.hpp"
#include "generator/ET_Writer.hpp"
#include "generator/HRR_Writer.hpp"
#include "generator/ERI_Writer.hpp"


int main(int argc, char ** argv)
{
    try {

    OptionMap options = DefaultOptions();

    // other stuff
    std::string boystype;
    std::string fpath;
    std::string cpuflags;
    std::string datdir;
    QAM finalam;

    bool finalamset = false;

    // parse command line
    std::vector<std::string> otheropt = ParseCommonOptions(options, argc, argv);

    // parse specific options
    size_t iarg = 0;
    while(iarg < otheropt.size())
    {
        std::string argstr(GetNextArg(iarg, otheropt));
        if(argstr == "-o")
            fpath = GetNextArg(iarg, otheropt);
        else if(argstr == "-c")
            cpuflags = GetNextArg(iarg, otheropt);
        else if(argstr == "-d")
            datdir = GetNextArg(iarg, otheropt);
        else if(argstr == "-b")
            boystype = GetNextArg(iarg, otheropt);
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
    CMDLINE_ASSERT( boystype != "", "Boys type (-b) required" )
    CMDLINE_ASSERT( fpath != "", "output path (-o) required" )
    CMDLINE_ASSERT( datdir != "", "dat directory (-d) required" )
    CMDLINE_ASSERT( finalamset == true, "AM quartet (-q) required" )
    CMDLINE_ASSERT( cpuflags != "", "CPU flags (-c) required" )


    // open the output file
    std::ofstream of(fpath);
    if(!of.is_open())
        throw std::runtime_error(StringBuilder("Cannot open file: ", fpath, "\n"));
    

    // Information for this ERI
    ERIGeneratorInfo info(finalam,
                          Compiler::Intel,
                          cpuflags, options);


    // is this one a special permutation?
    if( ( (finalam[0] == 0 && finalam[1] > 0)  && ( finalam[2] == 0 || finalam[3] == 0 ) ) ||
        ( (finalam[2] == 0 && finalam[3] > 0)  && ( finalam[0] == 0 || finalam[1] == 0 ) ) )
        WriteFile_Permute(of, info);
    else
    {
        // Read in the boys map
        std::unique_ptr<BoysGen> bg;

        if(boystype == "FO")
            bg = std::unique_ptr<BoysGen>(new BoysFO(datdir));
        else if(boystype == "split")
            bg = std::unique_ptr<BoysGen>(new BoysSplit());
        else
        {
            std::cout << "Unknown boys type \"" << boystype << "\"\n";
            return 3;
        }

        // algorithms used
        std::unique_ptr<HRR_Algorithm_Base> hrralgo(new Makowski_HRR(options));
        std::unique_ptr<VRR_Algorithm_Base> vrralgo(new Makowski_VRR(options));
        std::unique_ptr<ET_Algorithm_Base> etalgo(new Makowski_ET(options));

        // Writers
        std::unique_ptr<HRR_Writer> hrr_writer;
        std::unique_ptr<VRR_Writer> vrr_writer;
        std::unique_ptr<ET_Writer> et_writer;

        // Working backwards, I need:
        // 1.) HRR Steps
        hrralgo->Create(finalam);
        if(options[Option::InlineHRR])
            hrr_writer = std::unique_ptr<HRR_Writer>(new HRR_Writer_Inline(*hrralgo));
        else
            hrr_writer = std::unique_ptr<HRR_Writer>(new HRR_Writer_External(*hrralgo));


        // 2.) ET steps
        //     with the HRR top level stuff as the initial targets
        DoubletType etdirection = DoubletType::KET;
        if((finalam[2]+finalam[3]) > (finalam[0]+finalam[1]))
            etdirection = DoubletType::BRA;


        etalgo->Create(hrralgo->TopQuartets(), etdirection);

        if(options[Option::InlineET])
            et_writer = std::unique_ptr<ET_Writer>(new ET_Writer_Inline(*etalgo));
        else
            et_writer = std::unique_ptr<ET_Writer>(new ET_Writer_External(*etalgo));

        // 3.) VRR Steps
        // requirements for vrr are the top level stuff from ET
        vrralgo->Create(etalgo->TopQuartets());

        if(options[Option::InlineVRR])
            vrr_writer = std::unique_ptr<VRR_Writer>(new VRR_Writer_Inline(*vrralgo));
        else
            vrr_writer = std::unique_ptr<VRR_Writer>(new VRR_Writer_External(*vrralgo));

        // set the contracted quartets
        info.SetContQ(hrralgo->TopAM());

        // print out some info
        std::cout << "MEMORY (per shell quartet): " << info.ContMemoryReq() << "\n";

        WriteFile(of, info, *bg, *vrr_writer, *et_writer, *hrr_writer);
    }
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
