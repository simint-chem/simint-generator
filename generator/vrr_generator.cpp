#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>

#include "generator/Algorithms.hpp"
#include "generator/Helpers.hpp"
#include "generator/Options.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/VRR_Writer.hpp"

using namespace std;



string GetNextArg(int & i, int argc, char ** argv)
{
    if(i >= argc)
        throw std::runtime_error("Error - no more arguments!");

    return argv[i++];
}

int GetIArg(int & i, int argc, char ** argv)
{   
    std::string str = GetNextArg(i, argc, argv);
    try {

        return stoi(str);
    }
    catch(...)
    {
        std::stringstream ss;
        ss << "Cannot convert to int: " << str;
        throw std::runtime_error(ss.str());
    }
}


int main(int argc, char ** argv)
{
    try {

    // default options
    OptionsMap options = DefaultOptions();

    // max L value
    int maxL = 0;

    // other stuff
    std::string fpath;
    std::string cpuinfofile;

    // parse command line
    int i = 1;
    while(i < argc)
    {
        std::string argstr(GetNextArg(i, argc, argv));
        if(argstr == "-L")
            maxL = GetIArg(i, argc, argv);
        else if(argstr == "-o")
            fpath = GetNextArg(i, argc, argv);
        else if(argstr == "-i")
        {
            options[OPTION_INTRINSIC] = 1;
            cpuinfofile = GetNextArg(i, argc, argv);
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

    if(fpath == "")
    {
        std::cout << "\noutput path (-o) required\n\n";
        return 2;
    }

    if(maxL == 0)
    {
        std::cout << "\nMaximum L value (-L) required\n\n";
        return 2;
    }

    if(fpath.back() != '/')
        fpath += '/';


    // The algorithm to use 
    std::unique_ptr<VRR_Algorithm_Base> vrralgo(new Makowski_VRR);

    // different source and header files
    std::string srcpath = fpath + "vrr.c";
    std::string headpath = fpath + "vrr.h";
    
    cout << "Generating source file " << srcpath << "\n";
    cout << "Generating header file " << headpath << "\n";

    std::ofstream of(srcpath);
    if(!of.is_open())
    {
        std::cout << "Cannot open file: " << srcpath << "\n";
        return 2; 
    }

    std::ofstream ofh(headpath);
    if(!ofh.is_open())
    {
        std::cout << "Cannot open file: " << headpath << "\n";
        return 2; 
    }


    // we want all gaussians up to the maximum L value
    GaussianMap vreq;
    for(int i = 0; i <= maxL; i++)
        vreq[i] = AllGaussiansForAM(i);

    // Create the mapping
    vrralgo->CreateAllMaps(vreq);

    // Create the writer and base writer
    WriterInfo::Init(options, "", {0, 0, 0, 0});  // the amlist parameter doesn't matter much here

    // read in cpuflags if needed
    if(options[OPTION_INTRINSIC] != 0)
        WriterInfo::ReadCPUFlags(cpuinfofile); 

    VRR_Writer vrr_writer(*vrralgo);

    // write to the output file
    vrr_writer.WriteVRRFile(of);
    vrr_writer.WriteVRRHeaderFile(ofh);
    cout << "Done!\n";

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
