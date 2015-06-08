#include <iostream>
#include <map>

#include "generator/Classes.hpp"

using std::cout;

static std::map<ExpList, int> ordermap_;

QuartetSet GenerateInitialQuartetTargets(QAM amlst, bool initial)
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

DoubletSet GenerateInitialDoubletTargets(DAM amlst, DoubletType type, bool initial)
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

void PrintGaussianMap(const GaussianMap & g, const std::string & title)
{
    cout << title << ": " << g.size() << "\n";
    for(const auto & it : g)
    {
        cout << "  AM: " << it.first << "\n";
        for(const auto & it2 : it.second)
            cout << "      " << it2 << "\n"; 
        cout << "\n";
    }
    
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
