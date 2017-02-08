#include <stdexcept>
#include <fstream>

#include "generator/CommandLine.hpp"
#include "generator/StringBuilder.hpp"

#include "generator/ostei/Algorithms.hpp"
#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_VRR_Writer.hpp"


int main(int argc, char ** argv)
{
    try {

    // quartet we are looking for
    QAM finalam{0,0,0,0};

    // other stuff
    std::string fpath;
    std::string hpath;

    bool finalamset = false;

    // parse command line
    OptionMap options = DefaultOptions();
    std::vector<std::string> otheropt = ParseCommonOptions(options, argc, argv);

    IdxOrder order{0, -1, -1, -1};

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
        else if(argstr == "-center_j")
            order = IdxOrder{1, -1, -1, -1};
        else if(argstr == "-center_k")
            order = IdxOrder{2, -1, -1, -1};
        else if(argstr == "-center_l")
            order = IdxOrder{3, -1, -1, -1};
        else
        {
            std::cout << "\n\n";
            std::cout << "--------------------------------\n";
            std::cout << "Unknown argument: " << argstr << "\n";
            std::cout << "--------------------------------\n";
            return 1; 
        } 
    }

    CMDLINE_ASSERT( fpath != "", "output source file path (-o) required" )
    CMDLINE_ASSERT( hpath != "", "output header file path (-oh) required" )
    CMDLINE_ASSERT( finalamset == true, "AM quartet (-q) required" )


    // actually create
    std::cout << "Generating source file " << fpath << "\n";
    std::cout << "Appending to header file " << hpath << "\n";

    std::ofstream of(fpath);
    if(!of.is_open())
        throw std::runtime_error(StringBuilder("Cannot open file: ", fpath, "\n"));

    std::ofstream ofh(hpath, std::ofstream::app);
    if(!ofh.is_open())
        throw std::runtime_error(StringBuilder("Cannot open file: ", hpath, "\n"));

    OSTEI_GeneratorInfo info(finalam, 0, options);

    // output to source file
    // Include the main header file
    of << "\n#include \"simint/ostei/gen/ostei_generated.h\"\n";

    // create the algorithm and write the file
    Makowski_VRR vrralgo(info);
    vrralgo.Create_WithOrder(finalam, order);
    OSTEI_VRR_Writer vrr_writer(vrralgo, info);
    vrr_writer.WriteVRRFile(of, ofh);

    } // close try block
    catch(std::exception & ex)
    {
        std::cout << "\n\n";
        std::cout << "Caught exception\n";
        std::cout << "What = " << ex.what() << "\n\n";
        return 1;
    }
    return 0;
}
