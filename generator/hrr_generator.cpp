#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>

#include "generator/CommandLine.hpp"
#include "generator/Algorithms.hpp"

#include "generator/ERIGeneratorInfo.hpp"
#include "generator/HRR_Writer.hpp"


static
void CreateHRR(QAM am,
               const std::string & brapath, const std::string & ketpath,
               const std::string & cpuflags,
               const OptionMap & options, std::ostream & ofh)
{
        std::ofstream ofb(brapath);
        std::ofstream ofk(ketpath);

        if(!ofb.is_open())
            throw std::runtime_error(StringBuilder("Cannot open file: ", brapath, "\n"));
        if(!ofk.is_open())
            throw std::runtime_error(StringBuilder("Cannot open file: ", ketpath, "\n"));

        ERIGeneratorInfo info(am, Compiler::Intel, cpuflags, options);

        ofb << "\n#include \"eri/eri.h\"\n\n\n";
        ofk << "\n#include \"eri/eri.h\"\n\n\n";

        // disable this diagnostic (for now)
        ofb << "#ifdef __INTEL_COMPILER\n";
        ofb << "    #pragma warning(disable:2620)\n";
        ofb << "#endif\n\n\n";

        ofk << "#ifdef __INTEL_COMPILER\n";
        ofk << "    #pragma warning(disable:2620)\n";
        ofk << "#endif\n\n\n";

        // The algorithm to use
        std::unique_ptr<HRR_Algorithm_Base> hrralgo(new Makowski_HRR(options));

        hrralgo->Create(am);

        std::unique_ptr<HRR_Writer> hrr_writer(new HRR_Writer_External(*hrralgo));
        hrr_writer->WriteHRRFile(ofb, ofk, ofh, info);
}


int main(int argc, char ** argv)
{
    try {

    // default options
    OptionMap options = DefaultOptions();

    // max L value
    int maxL = 0;

    // other stuff
    std::string fpath;
    std::string cpuflags;

    // parse command line
    std::vector<std::string> otheropt = ParseCommonOptions(options, argc, argv);

    // parse specific options
    size_t iarg = 0;
    while(iarg < otheropt.size())
    {
        std::string argstr(GetNextArg(iarg, otheropt));
        if(argstr == "-L")
            maxL = GetIArg(iarg, otheropt);
        else if(argstr == "-o")
            fpath = GetNextArg(iarg, otheropt);
        else if(argstr == "-c")
            cpuflags = GetNextArg(iarg, otheropt);
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
    CMDLINE_ASSERT( maxL > 0, "Maximum L value (-L) greater than 0 required")

    if(fpath.back() != '/')
        fpath += '/';



    // different source and header files
    std::string headpath = fpath + "hrr.h";
    
    std::cout << "Generating header file " << headpath << "\n";


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

    std::vector<std::pair<int, int>> togen;
    // we want all doublets up to L
    for(int i = 1; i <= maxL; i++)
    for(int j = 1; j <= i; j++)
        togen.push_back({i,j});


    // we also need some more
    for(int i = maxL+1; i < 2*maxL; i++)
    for(int j = 1; j <= (2*maxL-i); j++)
        togen.push_back({i,j});


 

    // create left->right
    for(const auto & it : togen)
    {
        int i = it.first;
        int j = it.second;

        // create left->right
        std::string brapath1 = StringBuilder(fpath, "hrr_bra_J_", amchar[i], "_", amchar[j], ".c");
        std::string ketpath1 = StringBuilder(fpath, "hrr_ket_L_", amchar[i], "_", amchar[j], ".c");

        // create right->left
        std::string brapath2 = StringBuilder(fpath, "hrr_bra_I_", amchar[j], "_", amchar[i], ".c");
        std::string ketpath2 = StringBuilder(fpath, "hrr_ket_K_", amchar[j], "_", amchar[i], ".c");

        // Write out
        CreateHRR({i, j, i, j}, brapath1, ketpath1, cpuflags, options, ofh);
        CreateHRR({j, i, j, i}, brapath2, ketpath2, cpuflags, options, ofh);
    }


    ofh << "#endif\n";

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
