#include <stdexcept>
#include <fstream>

#include "generator/CommandLine.hpp"
#include "generator/StringBuilder.hpp"

#include "generator/ostei/Algorithms.hpp"
#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_HRR_Writer.hpp"


int main(int argc, char ** argv)
{
    try {

    // Doublet we are calculating
    DAM finalam{0,0};  // initialize to prevent compiler warnings

    // other stuff
    std::string fpath;
    std::string hpath;

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


    CMDLINE_ASSERT( fpath != "", "output path (-o) required" )
    CMDLINE_ASSERT( hpath != "", "output header file path (-oh) required" )
    CMDLINE_ASSERT( finalamset == true, "AM doublet (-q) required" )

    // We can do bras and kets at the same time
    QAM am{finalam[0], finalam[1], finalam[0], finalam[1]};

    std::ofstream of(fpath);
    if(!of.is_open())
        throw std::runtime_error(StringBuilder("Cannot open file: ", fpath, "\n"));

    std::ofstream ofh(hpath, std::ofstream::app);
    if(!ofh.is_open())
        throw std::runtime_error(StringBuilder("Cannot open file: ", hpath, "\n"));

    // Write out
    of << "\n#include \"simint/ostei/gen/ostei_generated.h\"\n";

    OSTEI_GeneratorInfo info(am, 0, options);

    // The algorithm to use
    Makowski_HRR hrralgo(options);
    hrralgo.Create(am);

    OSTEI_HRR_Writer hrr_writer(hrralgo, info);
    hrr_writer.WriteHRRFile(of, ofh);

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
