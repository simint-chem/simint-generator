#include <stdexcept>
#include <fstream>

#include "generator/CommandLine.hpp"

#include "generator/eri/Algorithms.hpp"
#include "generator/eri/ERIGeneratorInfo.hpp"
#include "generator/eri/ERI_ET_Writer.hpp"


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
    std::string headpath = fpath + "et.h";
    
    std::cout << "Generating header file " << headpath << "\n";


    std::ofstream ofh(headpath);
    if(!ofh.is_open())
        throw std::runtime_error(StringBuilder("Cannot open file: ", headpath, "\n"));

    // start the header file
    ofh << "#ifndef ET__H\n";
    ofh << "#define ET__H\n";
    ofh << "\n";
    ofh << "#include \"eri/eri.h\"\n\n\n";


    // only do this if we are actually doing ET
    // (we always need to create the header file though)
    if(!options[Option::NoET])
    {
        // we want all gaussians up to the maximum L value
        // First, bra -> ket
        for(int i = 0; i <= maxL; i++)
        for(int j = 1; j <= maxL; j++)
        {
            std::string srcpath = StringBuilder(fpath, "et_ket_", amchar[i], "_s_", amchar[j], "_s.c");
            std::cout << "Generating source file " << srcpath << "\n";
            std::ofstream of(srcpath);

            if(!of.is_open())
                throw std::runtime_error(StringBuilder("Cannot open file: ", srcpath, "\n"));

            QAM am{i, 0, j, 0};
            ERIGeneratorInfo info(am, Compiler::Intel, cpuflags, options);

            of << "#include \"eri/eri.h\"\n";

            std::unique_ptr<ERI_ET_Algorithm_Base> etalgo(new Makowski_ET(options));
            etalgo->Create(am, DoubletType::KET);

            std::unique_ptr<ERI_ET_Writer> et_writer(new ERI_ET_Writer_External(*etalgo, info));
            et_writer->WriteETFile(of, ofh);
        }


        // Now, ket->bra
        for(int i = 1; i <= maxL; i++)
        for(int j = 0; j <= maxL; j++)
        {
            std::string srcpath = StringBuilder(fpath, "et_bra_", amchar[i], "_s_", amchar[j], "_s.c");
            std::cout << "Generating source file " << srcpath << "\n";
            std::ofstream of(srcpath);

            if(!of.is_open())
                throw std::runtime_error(StringBuilder("Cannot open file: ", srcpath, "\n"));

            QAM am{i, 0, j, 0};
            ERIGeneratorInfo info(am, Compiler::Intel, cpuflags, options);

            of << "#include \"eri/eri.h\"\n";

            std::unique_ptr<ERI_ET_Algorithm_Base> etalgo(new Makowski_ET(options));
            etalgo->Create(am, DoubletType::BRA);

            std::unique_ptr<ERI_ET_Writer> et_writer(new ERI_ET_Writer_External(*etalgo, info));
            et_writer->WriteETFile(of, ofh);
        }
    }

    ofh << "\n";
    ofh << "#endif\n\n";
    

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
