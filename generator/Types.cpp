/*! \file
 *
 * \brief Common types used in the generator, plus some helper functions (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#include <iostream>
#include <map>
#include <algorithm>

#include "generator/Types.hpp"
#include "generator/Options.hpp"


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



int GaussianOrder(const ExpList & ijk)
{
    int am = ijk[0] + ijk[1] + ijk[2];

    const std::vector<ExpList> & v = gorder_map.at(am);

    auto it = std::find(v.begin(), v.end(), ijk);
    if(it == v.end())
        throw std::runtime_error("Gaussian not found in gorder_map");

    return std::distance(v.begin(), it);
}


GaussianSet AllGaussiansForAM(int am)
{
    GaussianSet gs;

    const std::vector<ExpList> & v = gorder_map.at(am);

    for(const auto & it : v)
        gs.insert({it[0], it[1], it[2]});

    return gs;
}

bool ValidQAM(QAM am)
{
    return (am[0] >= 0 &&
            am[1] >= 0 &&
            am[2] >= 0 &&
            am[3] >= 0);
}

