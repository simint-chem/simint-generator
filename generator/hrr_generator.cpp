#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>

#include "generator/Helpers.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Options.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/HRR_Writer.hpp"

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
            options[OPTION_INTRINSICS] = 1;
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




    // different source and header files
    std::string srcpath_bra = fpath + "hrr_bra.c";
    std::string srcpath_ket = fpath + "hrr_ket.c";
    std::string headpath = fpath + "hrr.h";
    
    cout << "Generating bra source file " << srcpath_bra << "\n";
    cout << "Generating ket source file " << srcpath_ket << "\n";
    cout << "Generating header file " << headpath << "\n";

    std::ofstream ofb(srcpath_bra);
    if(!ofb.is_open())
    {
        std::cout << "Cannot open file: " << srcpath_bra << "\n";
        return 2; 
    }

    std::ofstream ofk(srcpath_ket);
    if(!ofk.is_open())
    {
        std::cout << "Cannot open file: " << srcpath_ket << "\n";
        return 2; 
    }

    std::ofstream ofh(headpath);
    if(!ofh.is_open())
    {
        std::cout << "Cannot open file: " << headpath << "\n";
        return 2; 
    }

    // header guard and include file
    ofh << "#ifndef HRR_H\n";
    ofh << "#define HRR_H\n";

    ofh << "\n";
    ofh << "#include \"vectorization/vectorization.h\"\n";
    ofh << "\n\n";

    // include files for sources
    ofb << "\n";
    ofb << "#include \"vectorization/vectorization.h\"\n";
    ofb << "\n\n";

    ofk << "\n";
    ofk << "#include \"vectorization/vectorization.h\"\n";
    ofk << "\n\n";

    // we want all doublets up to L
    for(int i = 1; i <= maxL; i++)
    for(int j = 1; j <= i; j++)
    {
        // The algorithm to use
        std::unique_ptr<HRR_Algorithm_Base> hrralgo(new Makowski_HRR);

        // we can do both bras/kets in the same loop iteration
        QAM am{i, j, i, j};
        WriterInfo::Init(options, "", am);

        // read in cpuflags if needed
        if(options[OPTION_INTRINSICS] != 0)
            WriterInfo::ReadCPUFlags(cpuinfofile); 

        hrralgo->Create_DoubletStepLists(am);
        HRR_Writer hrr_writer(*hrralgo);

        // write to the output file (appending)
        hrr_writer.WriteHRRFile(ofb, ofk);
        hrr_writer.WriteHRRHeaderFile(ofh);
    }

    ofh << "#endif\n";

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
