#include <stdexcept>
#include <fstream>

#include "generator/CommandLine.hpp"

#include "generator/eri/Algorithms.hpp"
#include "generator/eri/ERIGeneratorInfo.hpp"
#include "generator/eri/ERI_VRR_Writer.hpp"


static
void CreateVRR(QAM am, const std::string & path, const std::string & cpuflags,
               const OptionMap & options, std::ostream & ofh)
{
    std::cout << "Generating source file " << path << "\n";
    std::ofstream of(path);
    if(!of.is_open())
        throw std::runtime_error(StringBuilder("Cannot open file: ", path, "\n"));

    ERIGeneratorInfo info(am, Compiler::Intel, cpuflags, options);

    // output to source file
    // Include the main eri file
    of << "#include \"eri/eri.h\"\n";

    // create the algorithm
    std::unique_ptr<ERI_VRR_Algorithm_Base> vrralgo(new Makowski_VRR(options));
    vrralgo->Create(am);

    // create the writer and write it
    std::unique_ptr<ERI_VRR_Writer> vrr_writer(new ERI_VRR_Writer_External(*vrralgo, info));
    vrr_writer->WriteVRRFile(of, ofh);
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
    std::string headpath = fpath + "vrr.h";
    
    std::cout << "Generating header file " << headpath << "\n";

    std::ofstream ofh(headpath);
    if(!ofh.is_open())
        throw std::runtime_error(StringBuilder("Cannot open file: ", headpath, "\n"));

    // start the header file
    ofh << "#ifndef VRR__H\n";
    ofh << "#define VRR__H\n";
    ofh << "\n";
    ofh << "#include \"eri/eri.h\"\n";
    ofh << "\n";

    // init once here to get the includes
    ERIGeneratorInfo info({maxL, 0, 0, 0},
                          Compiler::Intel,
                          cpuflags, options);


    // write out the includes
    for(const auto & it : info.GetIncludes())
        ofh << "#include " << it << "\n";
    ofh << "\n";


    // we want all VRR up to the maximum L value
    for(int i = 1; i <= maxL; i++)
    {
        std::string srcpath1 = StringBuilder(fpath, "vrr_", amchar[i], "_s_s_s.c");
        std::string srcpath2 = StringBuilder(fpath, "vrr_s_", amchar[i], "_s_s.c");
        std::string srcpath3 = StringBuilder(fpath, "vrr_s_s_", amchar[i], "_s.c");
        std::string srcpath4 = StringBuilder(fpath, "vrr_s_s_s_", amchar[i], ".c");

        CreateVRR({i, 0, 0, 0}, srcpath1, cpuflags, options, ofh);
        CreateVRR({0, i, 0, 0}, srcpath2, cpuflags, options, ofh);
        CreateVRR({0, 0, i, 0}, srcpath3, cpuflags, options, ofh);
        CreateVRR({0, 0, 0, i}, srcpath4, cpuflags, options, ofh);
    }


    if(options[Option::NoET] > 0)
    {
        // we have to generate ( X s | X s ) VRR relations
        // (and sXsX, and XssX, and sXXs)
        for(int i = 1; i <= maxL; i++)
        for(int j = 1; j <= maxL; j++)
        {
            std::string srcpath1 = StringBuilder(fpath, "vrr_", amchar[i], "_s_", amchar[j], "_s.c");
            std::string srcpath2 = StringBuilder(fpath, "vrr_s_", amchar[i], "_s_", amchar[j], ".c");
            std::string srcpath3 = StringBuilder(fpath, "vrr_", amchar[i], "_s_s_", amchar[j], ".c");
            std::string srcpath4 = StringBuilder(fpath, "vrr_s_", amchar[i], "_", amchar[j], "_s.c");

            CreateVRR({i, 0, j, 0}, srcpath1, cpuflags, options, ofh);
            CreateVRR({0, i, 0, j}, srcpath2, cpuflags, options, ofh);
            CreateVRR({i, 0, 0, j}, srcpath3, cpuflags, options, ofh);
            CreateVRR({0, i, j, 0}, srcpath4, cpuflags, options, ofh);
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
