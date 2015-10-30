#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>

#include "generator/Algorithms.hpp"
#include "generator/Helpers.hpp"
#include "generator/Options.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/VRR_Writer.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    try {

    // default options
    OptionsMap options = DefaultOptions();

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

    if(cpuflags == "")
    {
        std::cout << "\nCPU flags required\n\n";
        return 2;
    }

    if(fpath.back() != '/')
        fpath += '/';



    // different source and header files
    std::string headpath = fpath + "vrr.h";
    
    cout << "Generating header file " << headpath << "\n";

    std::ofstream ofh(headpath);
    if(!ofh.is_open())
    {
        std::cout << "Cannot open file: " << headpath << "\n";
        return 2; 
    }

    // start the header file
    ofh << "#ifndef VRR__H\n";
    ofh << "#define VRR__H\n";
    ofh << "\n";
    ofh << "#include \"eri/eri.h\"\n";
    ofh << "\n";

    // init once here to get the includes
    WriterInfo::Init(options, {maxL, 0, 0, 0}, cpuflags);
    WriterInfo::WriteIncludes(ofh);


    // we want all gaussians up to the maximum L value
    for(int i = 1; i <= maxL; i++)
    {
        std::stringstream ss1, ss2, ss3, ss4;
        ss1 << fpath << "vrr_" << amchar[i] << "_s_s_s.c";
        ss2 << fpath << "vrr_s_" << amchar[i] << "_s_s.c";
        ss3 << fpath << "vrr_s_s_" << amchar[i] << "_s.c";
        ss4 << fpath << "vrr_s_s_s_" << amchar[i] << "_.c";

        std::string srcpath1 = ss1.str();
        std::string srcpath2 = ss2.str();
        std::string srcpath3 = ss3.str();
        std::string srcpath4 = ss4.str();
        cout << "Generating source file " << srcpath1 << "\n";
        cout << "Generating source file " << srcpath2 << "\n";
        cout << "Generating source file " << srcpath3 << "\n";
        cout << "Generating source file " << srcpath4 << "\n";

        std::ofstream of1(srcpath1);
        std::ofstream of2(srcpath2);
        std::ofstream of3(srcpath3);
        std::ofstream of4(srcpath4);

        if(!of1.is_open())
        {
            std::cout << "Cannot open file: " << srcpath1 << "\n";
            return 2; 
        }

        if(!of2.is_open())
        {
            std::cout << "Cannot open file: " << srcpath2 << "\n";
            return 2; 
        }

        if(!of3.is_open())
        {
            std::cout << "Cannot open file: " << srcpath3 << "\n";
            return 2; 
        }

        if(!of4.is_open())
        {
            std::cout << "Cannot open file: " << srcpath4 << "\n";
            return 2; 
        }


        // output to source file
        of1 << "#include \"eri/eri.h\"\n";
        of2 << "#include \"eri/eri.h\"\n";
        of3 << "#include \"eri/eri.h\"\n";
        of4 << "#include \"eri/eri.h\"\n";

        // The algorithm to use 
        std::unique_ptr<VRR_Algorithm_Base> vrralgo1(new Makowski_VRR(options));
        std::unique_ptr<VRR_Algorithm_Base> vrralgo2(new Makowski_VRR(options));
        std::unique_ptr<VRR_Algorithm_Base> vrralgo3(new Makowski_VRR(options));
        std::unique_ptr<VRR_Algorithm_Base> vrralgo4(new Makowski_VRR(options));

        QAM am1{i, 0, 0, 0};
        QAM am2{0, i, 0, 0};
        QAM am3{0, 0, i, 0};
        QAM am4{0, 0, 0, i};


        // Create the mapping
        WriterInfo::Init(options, am1, cpuflags);
        vrralgo1->Create(am1);
        VRR_Writer vrr_writer1(*vrralgo1);
        vrr_writer1.WriteVRRFile(of1, ofh);

        WriterInfo::Init(options, am2, cpuflags);
        vrralgo2->Create(am2);
        VRR_Writer vrr_writer2(*vrralgo2);
        vrr_writer2.WriteVRRFile(of2, ofh);

        WriterInfo::Init(options, am3, cpuflags);
        vrralgo3->Create(am3);
        VRR_Writer vrr_writer3(*vrralgo3);
        vrr_writer3.WriteVRRFile(of3, ofh);

        WriterInfo::Init(options, am4, cpuflags);
        vrralgo4->Create(am4);
        VRR_Writer vrr_writer4(*vrralgo4);
        vrr_writer4.WriteVRRFile(of4, ofh);

        cout << "Done!\n";

    }

    if(options[OPTION_NOET] > 0)
    {
        // we have to generate ( X s | X s ) VRR relations
        for(int i = 1; i <= maxL; i++)
        for(int j = 1; j <= maxL; j++)
        {
            std::stringstream ss1, ss2;
            ss1 << fpath << "vrr_" << amchar[i] << "_s_" << amchar[j] << "_s.c";
            ss2 << fpath << "vrr_s_" << amchar[i] << "_s_" << amchar[j] << ".c";

            std::string srcpath1 = ss1.str();
            std::string srcpath2 = ss2.str();
            cout << "Generating source file " << srcpath1 << "\n";
            cout << "Generating source file " << srcpath2 << "\n";

            std::ofstream of1(srcpath1);
            std::ofstream of2(srcpath2);

            if(!of1.is_open())
            {
                std::cout << "Cannot open file: " << srcpath1 << "\n";
                return 2; 
            }

            if(!of2.is_open())
            {
                std::cout << "Cannot open file: " << srcpath2 << "\n";
                return 2; 
            }

            // output to source file
            of1 << "#include \"eri/eri.h\"\n";
            of2 << "#include \"eri/eri.h\"\n";

            // The algorithm to use 
            std::unique_ptr<VRR_Algorithm_Base> vrralgo1(new Makowski_VRR(options));
            std::unique_ptr<VRR_Algorithm_Base> vrralgo2(new Makowski_VRR(options));

            QAM am1{i, 0, j, 0};
            QAM am2{0, i, 0, j};

            WriterInfo::Init(options, am1, cpuflags);
            vrralgo1->Create(am1);
            VRR_Writer vrr_writer1(*vrralgo1);
            vrr_writer1.WriteVRRFile(of1, ofh);

            WriterInfo::Init(options, am2, cpuflags);
            vrralgo2->Create(am2);
            VRR_Writer vrr_writer2(*vrralgo2);
            vrr_writer2.WriteVRRFile(of2, ofh);

            cout << "Done!\n";

        }
    }

    ofh << "\n";
    ofh << "#endif\n\n";
    

    } // close try block
    catch(std::exception & ex)
    {
        cout << "\n\n";
        cout << "Caught exception\n";
        cout << "What = " << ex.what() << "\n\n";
        return 100;
    }
    return 0;
}
