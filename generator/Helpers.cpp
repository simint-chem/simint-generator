#include <iostream>
#include <map>

#include "generator/Classes.hpp"
#include "generator/Options.hpp"
#include "generator/Helpers.hpp"

using std::cout;

static std::map<ExpList, int> ordermap_;


// For options parsing
std::string GetNextArg(int & i, int argc, char ** argv)
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


std::string GetNextArg(size_t & i, const std::vector<std::string> & opt)
{
    if(i >= opt.size())
        throw std::runtime_error("Error - no more arguments!");

    return opt[i++];
}


int GetIArg(size_t & i, const std::vector<std::string> & opt)
{
    std::string str = GetNextArg(i, opt);
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


QuartetSet GenerateInitialQuartetTargets(QAM amlst)
{
    QuartetSet qs;
    int nam1 = ((amlst[0] + 1) * (amlst[0] + 2)) / 2;
    int nam2 = ((amlst[1] + 1) * (amlst[1] + 2)) / 2;
    int nam3 = ((amlst[2] + 1) * (amlst[2] + 2)) / 2;
    int nam4 = ((amlst[3] + 1) * (amlst[3] + 2)) / 2;

    Gaussian cur1 = Gaussian{amlst[0], 0, 0};
    for(int i = 0; i < nam1; i++)
    {
        Gaussian cur2 = Gaussian{amlst[1], 0, 0};
        for(int j = 0; j < nam2; j++)
        {
            Doublet bra{DoubletType::BRA, cur1, cur2};
            Gaussian cur3 = Gaussian{amlst[2], 0, 0};
            for(int k = 0; k < nam3; k++)
            {
                Gaussian cur4 = Gaussian{amlst[3], 0, 0};
                for(int l = 0; l < nam4; l++)
                {
                    Doublet ket{DoubletType::KET, cur3, cur4};
                    qs.insert(Quartet{bra, ket, 0});
                    cur4.Iterate();
                } 

                cur3.Iterate();
            }

            cur2.Iterate(); 
        }       
                
        cur1.Iterate();
    }
                
    return qs;
}

DoubletSet GenerateInitialDoubletTargets(DAM amlst, DoubletType type)
{
    DoubletSet ds;
    int nam1 = ((amlst[0] + 1) * (amlst[0] + 2)) / 2;
    int nam2 = ((amlst[1] + 1) * (amlst[1] + 2)) / 2;

    Gaussian cur1 = Gaussian{amlst[0], 0, 0};
    for(int i = 0; i < nam1; i++)
    {
        Gaussian cur2 = Gaussian{amlst[1], 0, 0};
        for(int j = 0; j < nam2; j++)
        {
            ds.insert(Doublet{type, cur1, cur2});
            cur2.Iterate();
        }
        cur1.Iterate();
    }
                
    return ds;
}



void PrintQuartetSet(const QuartetSet & q, const std::string & title)
{
    cout << title << ": " << q.size() << "\n";
    for(const auto & it : q)
        cout << "    " << it << "\n";
    cout << "\n";
}


void PrintDoubletSet(const DoubletSet & d, const std::string & title)
{
    cout << title << ": " << d.size() << "\n";
    for(const auto & it : d)
        cout << "    " << it << "\n";
    cout << "\n";
}

void PrintGaussianSet(const GaussianSet & g, const std::string & title)
{
    cout << title << ": " << g.size() << "\n";
    for(const auto & it : g)
        cout << "    " << it << "\n";
    cout << "\n";
}

int GaussianOrder(const ExpList & ijk)
{
    if(!ordermap_.size())
    {
        // generate the order map
        for(int i = 0; i <= 32; i++)
        {
            int n = 0;
            Gaussian g{i, 0, 0};
            ordermap_[g.ijk] = n++;
            
            while(g.Iterate())
                ordermap_[g.ijk] = n++;
        }
    }

    return ordermap_.at(ijk);
}


GaussianSet AllGaussiansForAM(int am)
{
    GaussianSet gs;
    Gaussian g{am, 0, 0};
    do {
        gs.insert(g);
    } while(g.Iterate());

    return gs;
}


OptionsMap DefaultOptions(void)
{
    // default options
    OptionsMap options;
    options[OPTION_STACKMEM] = 0;
    options[OPTION_INLINEVRR] = 1;
    options[OPTION_INLINEET] = 1;
    options[OPTION_INLINEHRR] = 1;
    options[OPTION_INTRINSICS] = 0;
    options[OPTION_SCALAR] = 0;

    options[OPTION_NOSINGLEET] = 0;
    options[OPTION_NOET] = 1;
    return options;
}


std::vector<std::string> ParseCommonOptions(OptionsMap & options, int argc, char ** argv)
{
    std::vector<std::string> ret;

    int i = 1;
    while(i < argc)
    {
        std::string argstr(GetNextArg(i, argc, argv));
        if(argstr == "-et")
            options[OPTION_NOET] = 0;
        if(argstr == "-etvrr")
            options[OPTION_NOSINGLEET] = 1;
        if(argstr == "-ve")
            options[OPTION_INLINEVRR] = 0;
        else if(argstr == "-ee")
            options[OPTION_INLINEET] = 0;
        else if(argstr == "-he")
            options[OPTION_INLINEHRR] = 0;
        else if(argstr == "-s")
            options[OPTION_STACKMEM] = GetIArg(i, argc, argv);
        else if(argstr == "-i")
            options[OPTION_INTRINSICS] = 1;
        else if(argstr == "-S")
            options[OPTION_SCALAR] = 1;
        else
            ret.push_back(argstr);
    }

    return ret;
}


