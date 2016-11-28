#include <stdexcept>
#include <fstream>

#include "generator/CommandLine.hpp"

#include "generator/ostei/Algorithms.hpp"
#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_ET_Writer.hpp"


int main(int argc, char ** argv)
{
    try {

    // quartet we are looking for
    QAM finalam;

    // other stuff
    std::string fpath;
    std::string hpath;
    std::string cpuflags;

    bool finalamset = 0;

    bool direction_bra = false;

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
        else if(argstr == "-b")
            direction_bra = true;
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

    OSTEI_GeneratorInfo info(finalam, 0, Compiler::Intel, cpuflags, options);

    of << "#include \"simint/ostei/ostei.h\"\n";

    std::unique_ptr<OSTEI_ET_Algorithm_Base> etalgo(new Makowski_ET(options));

    if(direction_bra)
        etalgo->Create(finalam, DoubletType::BRA);
    else
        etalgo->Create(finalam, DoubletType::KET);

    std::unique_ptr<OSTEI_ET_Writer> et_writer(new OSTEI_ET_Writer_External(*etalgo, info));
    et_writer->WriteETFile(of, ofh);
    

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
