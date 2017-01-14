#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <fstream>
#include <cmath>
#include <array>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "test/Common.hpp"

bool ValidGaussian(const std::array<int, 3> & g)
{
    return (g[0] >= 0 && g[1] >= 0 && g[2] >= 0); 
}


bool IterateGaussian(std::array<int, 3> & g)
{
    int am = g[0] + g[1] + g[2];

    if(g[2] == am)  // at the end
        return false;

    if(g[2] < (am - g[0]))
        g = {{ g[0],   g[1]-1,      g[2]+1 }};
    else
        g = {{ g[0]-1, am-g[0]+1,   0      }};

    return ValidGaussian(g);
}



GaussianVec CopyGaussianVec(const GaussianVec & v)
{
    GaussianVec copy;
    copy.reserve(v.size());
    for(const auto & it : v)
    {
        simint_shell sh;
        simint_initialize_shell(&sh);
        simint_copy_shell(&it, &sh);
        copy.push_back(sh);
    }
    return copy;
}



void FreeGaussianVec(GaussianVec & agv)
{
    for(auto & it : agv)
        simint_free_shell(&it);
    agv.clear();
}



ShellMap CopyShellMap(const ShellMap & m)
{
    ShellMap copy;
    for(auto & it : m)
        copy.insert({it.first, CopyGaussianVec(it.second)});
    return copy;
}



void FreeShellMap(ShellMap & m)
{
    for(auto & it : m)
        FreeGaussianVec(it.second);
    m.clear();
}



ShellMap ReadBasis(const std::string & file)
{
    std::map<char, int> ammap = {
                                         { 'S',     0 },
                                         { 'P',     1 },
                                         { 'D',     2 },
                                         { 'F',     3 },
                                         { 'G',     4 },
                                         { 'H',     5 },
                                         { 'I',     6 },
                                         { 'J',     7 },
                                         { 'K',     8 },
                                         { 'L',     9 },
                                         { 'M',    10 },
                                         { 'N',    11 },
                                         { 'O',    12 },
                                         { 'Q',    13 },
                                         { 'R',    14 },
                                         { 'T',    15 },
                                         { 'U',    16 },
                                         { 'V',    17 },
                                         { 'W',    18 },
                                         { 'X',    19 },
                                         { 'Y',    20 },
                                         { 'Z',    21 },
                                         { 'A',    22 },
                                         { 'B',    23 },
                                         { 'C',    24 },
                                         { 'E',    25 }
                                 };

    std::ifstream f(file.c_str());
    if(!f.is_open())
        throw std::runtime_error(std::string("Error opening file ") + file);

    f.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);

    ShellMap shellmap;

    int natom;
    f >> natom;


    for(int i = 0; i < natom; i++)
    {
        std::string sym; // not really used
        int nshell, nallprim, nallprimg;  // nprim, nprimg not really used
        double x, y, z;

        f >> sym >> nshell >> nallprim >> nallprimg >> x >> y >> z;

        for(int j = 0; j < nshell; j++)
        {
            std::string type;
            int nprim, ngen;
            f >> type >> nprim >> ngen;


            // allocate the shell
            // set the angular momentum
            // (circularly loop throught the type to handle sp, spd, etc)
            auto itg = type.begin();
            GaussianVec gvec(ngen);
            for(auto & it : gvec)
            {
                simint_initialize_shell(&it);
                simint_allocate_shell(nprim, &it);
                it.am = ammap[*itg];
                it.nprim = nprim;
                it.x = x;
                it.y = y;
                it.z = z;

                itg++;
                if(itg == type.end())
                    itg = type.begin();
            }

            // loop through general contractions
            for(int k = 0; k < nprim; k++)
            {
                double alpha;
                f >> alpha;

                for(int l = 0; l < ngen; l++)
                {
                    gvec[l].alpha[k] = alpha;
                    f >> gvec[l].coef[k];
                }

            }

            // add shell to the vector for its am
            for(auto & it : gvec)
                shellmap[it.am].push_back(it);
        }
    }

    return shellmap;
}



std::pair<int, int> FindMaxParams(const ShellMap & m)
{
    int maxam = 0;
    int maxel = 0;

    for(auto & it : m)
    {
        // it.first = am
        // it.second = shells with am of it.first

        if(it.first > maxam)
            maxam = it.first;

        const int ncart = NCART(it.first);
        const int nshell = it.second.size(); // nshell with this AM
        const int n = ncart * nshell;
        if(n > maxel)
            maxel = n;
    }

    return {maxam, maxel};
}







void Chop(double * const restrict calc, int ncalc)
{
    for(int i = 0; i < ncalc; i++)
    {
        if(fabs(calc[i]) < RND_ZERO)
            calc[i] = 0.0;
    }
}



std::pair<double, double> CalcError(double const * const restrict calc, double const * const restrict ref, int ncalc)
{
    double maxerr = 0;
    double maxrelerr = 0;

    for(int i = 0; i < ncalc; i++)
    {
        const double r = ref[i];
        const double diff = fabs(calc[i] - r);
        const double rel = (fabs(ref[i]) > RND_ZERO ? fabs(diff / r) : 0.0);
        if(diff > maxerr)
            maxerr = diff;
        if(rel > maxrelerr)
            maxrelerr = rel;
    }

    return std::pair<double, double>(maxerr, maxrelerr);
}

    


void PrintTimingHeader(void)
{
    // Timing header
    printf("%13s %12s %12s %16s %16s %16s %16s %16s %16s\n",
                           "Quartet", "NCont", "NPrim",
                           "Ticks(Fill)", "Ticks(Ints)", "FullTicks(Ints)",
                           "Time(Fill)", "Time(Ints)", "FullTime(Ints)");
}


void PrintAMTimingInfo(int i, int j, int k, int l, size_t nshell1234, size_t nprim1234, const TimeContrib & info)
{
        printf("( %d %d | %d %d ) %12lu %12lu %16llu %16llu %16llu %16llu %16llu %16llu\n",
                                                                      i, j, k, l,
                                                                      nshell1234, nprim1234,
                                                                      info.ticks_shell_pair.load(),
                                                                      info.ticks_integrals.load(),
                                                                      info.fullticks_integrals.load(),
                                                                      info.time_shell_pair.load(),
                                                                      info.time_integrals.load(),
                                                                      info.fulltime_integrals.load());
}

bool UniqueQuartet(int i, int j, int k, int l)
{
    if(i < j)
        return false;
    if(k < l)
        return false;
    if( (i + j) < (k + l) )
        return false;
    if( (i + j) == (k + l) && (i < k) ) 
        return false;
    return true;
}



