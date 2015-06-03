#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>

#include "generator/Classes.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Helpers.hpp"
#include "generator/Options.hpp"
#include "generator/WriterBase.hpp"
#include "generator/VRRWriter.hpp"

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
    OptionsMap options;

    // max L value
    int maxL = 0;

    // other stuff
    std::string fpath;

    // parse command line
    int i = 1;
    while(i < argc)
    {
        std::string argstr(GetNextArg(i, argc, argv));
        if(argstr == "-L")
            maxL = GetIArg(i, argc, argv);
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


    // Read in the boys map
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
    ETReqMap vreq;
    for(int i = 0; i <= maxL; i++)
        vreq[i] = AllGaussiansForAM(i);

    // Create the mapping
    std::pair<VRRMap, VRRReqMap> vrrinfo = vrralgo->CreateAllMaps(vreq);

    // Create the writer and base writer
    WriterBase base(options, {0, 0, 0, 0});  // the amlist parameter doesn't matter much here
    VRRWriter vrr_writer(vrrinfo.first, vrrinfo.second);

    // write to the output file
    vrr_writer.WriteVRRFile(of, base);
    vrr_writer.WriteVRRHeaderFile(ofh, base);
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
