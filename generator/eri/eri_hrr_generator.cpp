#include <stdexcept>
#include <fstream>

#include "generator/CommandLine.hpp"

#include "generator/eri/Algorithms.hpp"
#include "generator/eri/ERIGeneratorInfo.hpp"
#include "generator/eri/ERI_HRR_Writer.hpp"


int main(int argc, char ** argv)
{
    try {

    // Doublet we are calculating
    DAM finalam;

    // other stuff
    std::string fpath;
    std::string hpath;
    std::string cpuflags;

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
    of << "\n#include \"simint/eri/eri.h\"\n";

    ERIGeneratorInfo info(am, Compiler::Intel, cpuflags, options);

    // The algorithm to use
    std::unique_ptr<ERI_HRR_Algorithm_Base> hrralgo(new Makowski_HRR(options));

    hrralgo->Create(am);

    std::unique_ptr<ERI_HRR_Writer> hrr_writer(new ERI_HRR_Writer_External(*hrralgo, info));
    hrr_writer->WriteHRRFile(of, ofh);

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
