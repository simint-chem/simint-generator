#include <stdexcept>
#include <fstream>

#include "generator/CommandLine.hpp"

#include "generator/eri/Algorithms.hpp"
#include "generator/eri/ERIGeneratorInfo.hpp"
#include "generator/eri/ERI_VRR_Writer.hpp"


int main(int argc, char ** argv)
{
    try {

    // quartet we are looking for
    QAM finalam;

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

    ERIGeneratorInfo info(finalam, Compiler::Intel, cpuflags, options);

    // output to source file
    // Include the main eri file
    of << "#include \"simint/eri/eri.h\"\n";

    // create the algorithm
    std::unique_ptr<ERI_VRR_Algorithm_Base> vrralgo(new Makowski_VRR(options));
    vrralgo->Create(finalam);

    // create the writer and write it
    std::unique_ptr<ERI_VRR_Writer> vrr_writer(new ERI_VRR_Writer_External(*vrralgo, info));
    vrr_writer->WriteVRRFile(of, ofh);
    

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
