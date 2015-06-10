#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>

#include "generator/Classes.hpp"
#include "generator/Algorithms.hpp"
#include "generator/Boys.hpp"
#include "generator/FileWriter.hpp"

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
    OptionsMap options;
    options[OPTION_STACKMEM] = 0;
    options[OPTION_INLINEVRR] = 1;
    options[OPTION_INLINEHRR] = 1;
    options[OPTION_PERMUTE] = 0;

    // other stuff
    std::string prefix;
    std::string boystype;
    std::string fpath;
    QAM amlist;
    bool amlistset = false;

    // parse command line
    int i = 1;
    while(i < argc)
    {
        std::string argstr(GetNextArg(i, argc, argv));
        if(argstr == "-ve")
            options[OPTION_INLINEVRR] = 0;
        else if(argstr == "-he")
            options[OPTION_INLINEHRR] = 0;
        else if(argstr == "-s")
            options[OPTION_STACKMEM] = GetIArg(i, argc, argv);
        else if(argstr == "-P")
            options[OPTION_PERMUTE] = 1;
        else if(argstr == "-p")
            prefix = GetNextArg(i, argc, argv);

        else if(argstr == "-q")
        {
            amlist[0] = GetIArg(i, argc, argv);   
            amlist[1] = GetIArg(i, argc, argv);   
            amlist[2] = GetIArg(i, argc, argv);   
            amlist[3] = GetIArg(i, argc, argv);   
            amlistset = true;
        }
        else if(argstr == "-b")
            boystype = GetNextArg(i, argc, argv);
        else if(argstr == "-o")
            fpath = GetNextArg(i, argc, argv);
        else
        {
            std::cout << "\n\n";
            std::cout << "--------------------------------\n";
            std::cout << "Unknown argument: " << argstr << "\n";
            std::cout << "--------------------------------\n";
            return 1; 
        } 
    }

    // check for required options
    if(prefix == "")
    {
        std::cout << "\nprefix (-p) required\n\n";
        return 2;
    }

    if(boystype == "")
    {
        std::cout << "\nBoys type (-b) required\n\n";
        return 2;
    }

    if(fpath == "")
    {
        std::cout << "\noutput path (-o) required\n\n";
        return 2;
    }

    if(amlistset == false)
    {
        std::cout << "\nAM quartet (-q) required\n\n";
        return 2;
    }

    if(options[OPTION_PERMUTE] != 0)
    {
        std::cout << "\nPermutation of integrals not implemented\n\n";
        return 10;
    }

    // Read in the boys map
    std::unique_ptr<BoysGen> bg;

    if(boystype == "FO")
        bg = std::unique_ptr<BoysGen>(new BoysFO("/home/ben/programming/simint/generator/dat"));
    else if(boystype == "split")
        bg = std::unique_ptr<BoysGen>(new BoysSplit());
    else if(boystype == "vref")
        bg = std::unique_ptr<BoysGen>(new BoysVRef());
    else
    {
        std::cout << "Unknown boys type \"" << boystype << "\"\n";
        return 3;
    }


    // algorithms used
    std::unique_ptr<HRR_Algorithm_Base> hrralgo(new Makowski_HRR);
    std::unique_ptr<VRR_Algorithm_Base> vrralgo(new Makowski_VRR);
    std::unique_ptr<ET_Algorithm_Base> etalgo(new Makowski_ET);

    std::ofstream of(fpath);
    if(!of.is_open())
    {
        std::cout << "Cannot open file: " << fpath << "\n";
        return 2; 
    }

    WriteFile(of, amlist, prefix, options, *bg, *vrralgo, *etalgo, *hrralgo);

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
