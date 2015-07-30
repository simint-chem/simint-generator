#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <fstream>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "eri/eri.h"

#include "test/common.hpp"
#include "test/valeev.hpp"

static erifunc funcs[MAXAM+1][MAXAM+1][MAXAM+1][MAXAM+1];


/////////////////////////
// For RefIntegralReader
/////////////////////////
RefIntegralReader::RefIntegralReader(const std::string & basfile)
         : file_( (basfile + ".ref").c_str(), std::ifstream::in)
{
    if(!file_.is_open())
        throw std::runtime_error(std::string("Error opening file ") + basfile + ".ref");
    file_.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
}

void RefIntegralReader::ReadNext(double * out, int nsize)
{
    file_.read(reinterpret_cast<char *>(out), nsize * sizeof(double));
}



void RefIntegralReader::Reset(void)
{
    file_.seekg(0, std::ifstream::beg);
}



bool IsValidGaussian(const std::array<int, 3> & g)
{
    return (g[0] >= 0 && g[1] >= 0 && g[2] >= 0);
}



bool IterateGaussian(std::array<int, 3> & g)
{
    int am = g[0] + g[1] + g[2];

    if(g[2] == am)  // at the end
        return false;

    if(g[2] < (am - g[0]))
        g = {g[0],   g[1]-1,      g[2]+1 };
    else
        g = {g[0]-1, am-g[0]+1, 0        };

    return IsValidGaussian(g);
}



bool ValidQuartet(std::array<int, 4> am)
{
    if(am[0] < am[1])
        return false;
    if(am[2] < am[3])
        return false;
    if( (am[0] + am[1]) < (am[2] + am[3]) )
        return false;
    if(am[0] < am[2])
        return false;
    return true;
}



bool ValidQuartet(int i, int j, int k, int l)
{
    return ValidQuartet({i, j, k, l});
}


GaussianVec CopyGaussianVec(const GaussianVec & v)
{
    GaussianVec copy;
    copy.reserve(v.size());
    for(const auto & it : v)
        copy.push_back(copy_gaussian_shell(it));
    return copy;
}



void FreeGaussianVec(GaussianVec & agv)
{
    for(auto & it : agv)
        free_gaussian_shell(it);
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
                allocate_gaussian_shell(nprim, &it);
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



std::array<int, 3> FindMapMaxParams(const ShellMap & m)
{
    int maxnprim = 0;
    int maxam = 0;
    int maxel = 0;
    for(auto & it : m)
    {
        const int nca = NCART(it.first);
        const int nsh = it.second.size();
        const int n = nca * nsh;
        if(n > maxel)
            maxel = n;

        if(it.first > maxam)
            maxam = it.first;

        for(auto & it2 : it.second)
        {
            if(it2.nprim > maxnprim)
                maxnprim = it2.nprim;
        }
    }

    return {maxam, maxnprim, maxel*maxel*maxel*maxel};
}



void ValeevIntegrals(const GaussianVec & gv1, const GaussianVec & gv2,
                     const GaussianVec & gv3, const GaussianVec & gv4,
                     double * const integrals, bool normalize)
{
    int inorm = (normalize ? 1 : 0);
    const gaussian_shell * A = gv1.data();
    const gaussian_shell * B = gv2.data();
    const gaussian_shell * C = gv3.data();
    const gaussian_shell * D = gv4.data();

    const int nshell1 = gv1.size();
    const int nshell2 = gv2.size();
    const int nshell3 = gv3.size();
    const int nshell4 = gv4.size();

    const int am1 = A[0].am;
    const int am2 = B[0].am;
    const int am3 = C[0].am;
    const int am4 = D[0].am;

    const int n234 = nshell2 * nshell3 * nshell4 * NCART(am1) * NCART(am2) * NCART(am3) * NCART(am4);

    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for(int i = 0; i < nshell1; i++)
    {
        const int idxstart = i * n234;
        int idx = 0;

        for(int j = 0; j < nshell2; j++)
        for(int k = 0; k < nshell3; k++)
        for(int l = 0; l < nshell4; l++)
        {
            long double vA[3] = { A[i].x, A[i].y, A[i].z };
            long double vB[3] = { B[j].x, B[j].y, B[j].z };
            long double vC[3] = { C[k].x, C[k].y, C[k].z };
            long double vD[3] = { D[l].x, D[l].y, D[l].z };

            std::array<int, 3> g1{am1, 0, 0};
            do
            {
                std::array<int, 3> g2{am2, 0, 0};
                do
                {
                    std::array<int, 3> g3{am3, 0, 0};
                    do
                    {
                        std::array<int, 3> g4{am4, 0, 0};
                        do
                        {
                            double myint = 0.0;

                            for(int m = 0; m < A[i].nprim; m++)
                            for(int n = 0; n < B[j].nprim; n++)
                            for(int o = 0; o < C[k].nprim; o++)
                            for(int p = 0; p < D[l].nprim; p++)
                            {

                                double val = (double)Valeev_eri(g1[0], g1[1], g1[2], A[i].alpha[m], vA,
                                                                g2[0], g2[1], g2[2], B[j].alpha[n], vB,
                                                                g3[0], g3[1], g3[2], C[k].alpha[o], vC,
                                                                g4[0], g4[1], g4[2], D[l].alpha[p], vD, inorm);
                                myint += val * A[i].coef[m] * B[j].coef[n] * C[k].coef[o] * D[l].coef[p];

/*
                                printf("IDX: %d\n", idx);
                                printf("VAL: %8.3e\n", val);
                                printf("%8.3e %8.3e %8.3e %8.3e\n", vA[0], vA[1], vA[2], A[i].alpha[m]);
                                printf("%8.3e %8.3e %8.3e %8.3e\n", vB[0], vB[1], vB[2], B[j].alpha[n]);
                                printf("%8.3e %8.3e %8.3e %8.3e\n", vC[0], vC[1], vC[2], C[k].alpha[o]);
                                printf("%8.3e %8.3e %8.3e %8.3e\n", vD[0], vD[1], vD[2], D[l].alpha[p]);
                                printf("%8.3e %8.3e %8.3e %8.3e\n", A[i].coef[m], B[j].coef[n], C[k].coef[o], D[l].coef[p]);
                                std::cout << "i,j,k,l,idx = ";
                                std::cout << g1[0] << "," << g1[1] << "," << g1[2] << "  ";
                                std::cout << g2[0] << "," << g2[1] << "," << g2[2] << "  ";
                                std::cout << g3[0] << "," << g3[1] << "," << g3[2] << "  ";
                                std::cout << g4[0] << "," << g4[1] << "," << g4[2] << "  ";
                                std::cout << idx << "  int = " << integrals[idx] << "\n";
*/


                            }

                            integrals[idxstart + idx] = myint;
                            idx++;

                        } while(IterateGaussian(g4));
                    } while(IterateGaussian(g3));
                } while(IterateGaussian(g2));
            } while(IterateGaussian(g1));
        }


    }
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



void Init_Test(void)
{
    for(int i = 0; i <= MAXAM; i++)
    for(int j = 0; j <= MAXAM; j++)
    for(int k = 0; k <= MAXAM; k++)
    for(int l = 0; l <= MAXAM; l++)
    {
        funcs[i][j][k][l] = eri_notyetimplemented;
    }


    funcs[0][0][0][0] = eri_s_s_s_s;
    #if MAXAM > 0
    funcs[1][0][0][0] = eri_p_s_s_s;
    funcs[1][0][1][0] = eri_p_s_p_s;
    funcs[1][1][0][0] = eri_p_p_s_s;
    funcs[1][1][1][0] = eri_p_p_p_s;
    funcs[1][1][1][1] = eri_p_p_p_p;
    #endif
    #if MAXAM > 1
    funcs[2][0][0][0] = eri_d_s_s_s;
    funcs[2][0][1][0] = eri_d_s_p_s;
    funcs[2][0][1][1] = eri_d_s_p_p;
    funcs[2][0][2][0] = eri_d_s_d_s;
    funcs[2][1][0][0] = eri_d_p_s_s;
    funcs[2][1][1][0] = eri_d_p_p_s;
    funcs[2][1][1][1] = eri_d_p_p_p;
    funcs[2][1][2][0] = eri_d_p_d_s;
    funcs[2][1][2][1] = eri_d_p_d_p;
    funcs[2][2][0][0] = eri_d_d_s_s;
    funcs[2][2][1][0] = eri_d_d_p_s;
    funcs[2][2][1][1] = eri_d_d_p_p;
    funcs[2][2][2][0] = eri_d_d_d_s;
    funcs[2][2][2][1] = eri_d_d_d_p;
    funcs[2][2][2][2] = eri_d_d_d_d;
    #endif
    #if MAXAM > 2
    funcs[3][0][0][0] = eri_f_s_s_s;
    funcs[3][0][1][0] = eri_f_s_p_s;
    funcs[3][0][1][1] = eri_f_s_p_p;
    funcs[3][0][2][0] = eri_f_s_d_s;
    funcs[3][0][2][1] = eri_f_s_d_p;
    funcs[3][0][3][0] = eri_f_s_f_s;
    funcs[3][1][0][0] = eri_f_p_s_s;
    funcs[3][1][1][0] = eri_f_p_p_s;
    funcs[3][1][1][1] = eri_f_p_p_p;
    funcs[3][1][2][0] = eri_f_p_d_s;
    funcs[3][1][2][1] = eri_f_p_d_p;
    funcs[3][1][2][2] = eri_f_p_d_d;
    funcs[3][1][3][0] = eri_f_p_f_s;
    funcs[3][1][3][1] = eri_f_p_f_p;
    funcs[3][2][0][0] = eri_f_d_s_s;
    funcs[3][2][1][0] = eri_f_d_p_s;
    funcs[3][2][1][1] = eri_f_d_p_p;
    funcs[3][2][2][0] = eri_f_d_d_s;
    funcs[3][2][2][1] = eri_f_d_d_p;
    funcs[3][2][2][2] = eri_f_d_d_d;
    funcs[3][2][3][0] = eri_f_d_f_s;
    funcs[3][2][3][1] = eri_f_d_f_p;
    funcs[3][2][3][2] = eri_f_d_f_d;
    funcs[3][3][0][0] = eri_f_f_s_s;
    funcs[3][3][1][0] = eri_f_f_p_s;
    funcs[3][3][1][1] = eri_f_f_p_p;
    funcs[3][3][2][0] = eri_f_f_d_s;
    funcs[3][3][2][1] = eri_f_f_d_p;
    funcs[3][3][2][2] = eri_f_f_d_d;
    funcs[3][3][3][0] = eri_f_f_f_s;
    funcs[3][3][3][1] = eri_f_f_f_p;
    funcs[3][3][3][2] = eri_f_f_f_d;
    funcs[3][3][3][3] = eri_f_f_f_f;
    #endif

}


// Disable diagnostic about unused parameter "dummy"
#ifdef __INTEL_COMPILER
  #pragma warning(push)
  #pragma warning(disable:869)
#endif
int eri_notyetimplemented(struct multishell_pair const P,
                          struct multishell_pair const Q,
                          double * const restrict dummy)
{
    printf("****************************\n");
    printf("*** NOT YET IMPLEMENTED! ***\n");
    printf("***  ( %2d %2d | %2d %2d )   ***\n", P.am1, P.am2, Q.am1, Q.am2);
    printf("****************************\n");
    exit(1);
    return 0;
}
#ifdef __INTEL_COMPILER
  #pragma warning(pop)
#endif



TimerInfo Integral(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict integrals)
{
    unsigned long long ticks0, ticks1;

    CLOCK(ticks0);
    funcs[P.am1][P.am2][Q.am1][Q.am2](P, Q, integrals);
    CLOCK(ticks1);
    return {ticks1 - ticks0, (ticks1 - ticks0)/(1.0e9*PROC_CYCLES_PER_SECOND)};
}

