#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>

#include "generator/Helpers.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Boys.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/VRR_Writer.hpp"
#include "generator/ET_Writer.hpp"
#include "generator/HRR_Writer.hpp"
#include "generator/FileWriter.hpp"

using namespace std;


int main(int argc, char ** argv)
{
    try {

    OptionsMap options = DefaultOptions();

    // other stuff
    std::string boystype;
    std::string fpath;
    std::string cpuinfofile;
    std::string datdir;
    QAM amlist;

    bool amlistset = false;

    // parse command line
    int i = 1;
    while(i < argc)
    {
        std::string argstr(GetNextArg(i, argc, argv));
        if(argstr == "-ve")
            options[OPTION_INLINEVRR] = 0;
        else if(argstr == "-he")
            options[OPTION_INLINEHRR] = 0;
        else if(argstr == "-s")
            options[OPTION_STACKMEM] = GetIArg(i, argc, argv);
        else if(argstr == "-d")
            datdir = GetNextArg(i, argc, argv);

        else if(argstr == "-c")
            cpuinfofile = GetNextArg(i, argc, argv);

        else if(argstr == "-i")
            options[OPTION_INTRINSICS] = 1;
        else if(argstr == "-S")
            options[OPTION_SCALAR] = 1;

        else if(argstr == "-q")
        {
            amlist[0] = GetIArg(i, argc, argv);   
            amlist[1] = GetIArg(i, argc, argv);   
            amlist[2] = GetIArg(i, argc, argv);   
            amlist[3] = GetIArg(i, argc, argv);   
            amlistset = true;
        }
        else if(argstr == "-b")
            boystype = GetNextArg(i, argc, argv);
        else if(argstr == "-o")
            fpath = GetNextArg(i, argc, argv);
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
    if(boystype == "")
    {
        std::cout << "\nBoys type (-b) required\n\n";
        return 2;
    }

    if(fpath == "")
    {
        std::cout << "\noutput path (-o) required\n\n";
        return 2;
    }

    if(datdir == "")
    {
        std::cout << "\ndat directory (-d) required\n\n";
        return 2;
    }

    if(amlistset == false)
    {
        std::cout << "\nAM quartet (-q) required\n\n";
        return 2;
    }

    if(cpuinfofile == "")
    {
        std::cout << "\nCPU info file required\n\n";
        return 2;
    }

    if(options[OPTION_INTRINSICS] != 0 && options[OPTION_SCALAR] != 0)
    {
        std::cout << "\nUsing intrinsics with a scalar calculation doesn't make sense...\n\n";
        return 2;
    }


    // open the output file
    std::ofstream of(fpath);
    if(!of.is_open())
    {
        std::cout << "Cannot open file: " << fpath << "\n";
        return 2; 
    }
    

    // Read in the boys map
    std::unique_ptr<BoysGen> bg;

    if(boystype == "FO")
        bg = std::unique_ptr<BoysGen>(new BoysFO(datdir));
    else if(boystype == "split")
        bg = std::unique_ptr<BoysGen>(new BoysSplit());
    else if(boystype == "vref")
        bg = std::unique_ptr<BoysGen>(new BoysVRef());
    else
    {
        std::cout << "Unknown boys type \"" << boystype << "\"\n";
        return 3;
    }

    // algorithms used
    std::unique_ptr<HRR_Algorithm_Base> hrralgo(new Makowski_HRR);
    std::unique_ptr<VRR_Algorithm_Base> vrralgo(new Makowski_VRR);
    std::unique_ptr<ET_Algorithm_Base> etalgo(new Makowski_ET);

    // Base writer information
    WriterInfo::Init(options, amlist, cpuinfofile);

    // Working backwards, I need:
    // 1.) HRR Steps
    hrralgo->Create(amlist);
    HRR_Writer hrr_writer(*hrralgo);

    // 2.) ET steps
    //     with the HRR top level stuff as the initial targets
    etalgo->Create(hrralgo->TopAM());
    ET_Writer et_writer(*etalgo);

    // 3.) VRR Steps
    // requirements for vrr are the top level stuff from ET
    vrralgo->Create(etalgo->TopQuartets());
    VRR_Writer vrr_writer(*vrralgo);

    // set the contracted quartets
    WriterInfo::SetContQ(hrralgo->TopAM());

    // print out some info
    std::cout << "MEMORY (per shell quartet): " << WriterInfo::MemoryReq() << "\n";

    WriteFile(of, *bg, vrr_writer, et_writer, hrr_writer);

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
