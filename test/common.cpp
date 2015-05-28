#include <iostream>
#include <cstdlib>
#include <ctime>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <cmath>
    
#ifdef _OPENMP
#include <omp.h>
#endif

#include "eri/shell.h"
#include "eri/eri.h"

#include "test/common.hpp"
#include "test/valeev.hpp"
#include "test/erd_interface.hpp"
#include "test/cppvectorization.hpp"

static erifunc funcs_FO[MAXAM+1][MAXAM+1][MAXAM+1][MAXAM+1];
static erifunc funcs_vref[MAXAM+1][MAXAM+1][MAXAM+1][MAXAM+1];
static eriflatfunc funcs_vref_flat[MAXAM+1][MAXAM+1][MAXAM+1][MAXAM+1];

static void AppendSlash(std::string & str)
{
    if(str[str.size()-1] != '/')
        str = str + '/';
}

/*
struct gaussian_shell
random_shell(int nprim, int am)
{
    struct gaussian_shell G;
    allocate_gaussian_shell(nprim, &G);

    G.am = am;
    G.x = 2.0 * MAX_COORD * rand() / ((double)RAND_MAX) - MAX_COORD;
    G.y = 2.0 * MAX_COORD * rand() / ((double)RAND_MAX) - MAX_COORD;
    G.z = 2.0 * MAX_COORD * rand() / ((double)RAND_MAX) - MAX_COORD;

    G.nprim = nprim;

    for(int i = 0; i < nprim; i++)
    {
        G.alpha[i] = MAX_EXP * rand() / ((double)RAND_MAX);
        G.coef[i] = MAX_COEF * rand() / ((double)RAND_MAX);
    }

    return G;
}

VecQuartet
CreateRandomQuartets(std::array<int, 4> nshell,
                     std::array<int, 4> nprim,
                     std::array<int, 4> am,
                     bool normalize = true)
{
    VecQuartet arr;

    for(int q = 0; q < 4; q++)
    {
        arr[q].reserve(nshell[q]);

        for(int i = 0; i < nshell[q]; i++)
            arr[q].push_back(random_shell(nprim[q], am[q]));

        if(normalize)
            normalize_gaussian_shells(nshell[q], arr[q].data());
    }

    return arr;
}
*/

VecQuartet ReadQuartets(std::array<int, 4> am,
                        std::string basedir,
                        bool normalize)
{
    AppendSlash(basedir);

    VecQuartet arr;

    std::array<std::string, 4> files{basedir + "1.dat",
                                     basedir + "2.dat",
                                     basedir + "3.dat",
                                     basedir + "4.dat"};

    for(int q = 0; q < 4; q++)
    {
        std::ifstream file(files[q].c_str());
        if(!file.is_open())
            throw std::runtime_error(std::string("Error opening file ") + files[q]);

        file.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);

        int fnshell, fnprim;
        file >> fnshell >> fnprim;

        arr[q].reserve(fnshell);

        for(int i = 0; i < fnshell; i++)
        {
            gaussian_shell G;
            allocate_gaussian_shell(fnprim, &G);
            G.am = am[q];

            G.nprim = fnprim;
            file >> G.x >> G.y >> G.z;

            for(int j = 0; j < fnprim; j++)
                file >> G.alpha[j] >> G.coef[j];

            arr[q].push_back(G);
        }
    }

    if(normalize)
    {
        for(auto & it : arr)
            normalize_gaussian_shells(it.size(), it.data());
    }

    return arr;
}


VecQuartet CopyQuartets(const VecQuartet & orig)
{
    VecQuartet v;
    for(int q = 0; q < 4; q++)
    {
        v[q].reserve(orig[q].size());
        for(const auto & it : orig[q])
            v[q].push_back(copy_gaussian_shell(it));
    }

    return v;
}


void FreeQuartets(VecQuartet & arr)
{
    for(auto & it : arr)
    {
        for(auto & it2 : it)
            free_gaussian_shell(it2);
        it.clear();
    }
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


int ReadValeevIntegrals(std::string basedir,
                        const std::array<int, 4> & am,
                        double * res)
{
    AppendSlash(basedir);

    const char * amchar = "spdfghijklmnoqrtuvwxyzabe";

    uint32_t nsize = 0;

    std::stringstream ss;
    ss << basedir << "ref_" << amchar[am[0]] << "_"
                            << amchar[am[1]] << "_"
                            << amchar[am[2]] << "_"
                            << amchar[am[3]] << ".dat";

    std::ifstream infile(ss.str().c_str(), std::ifstream::binary);
    if(!infile.is_open())
        throw std::runtime_error(std::string("Unable to open ") + ss.str());

    infile.read(reinterpret_cast<char *>(&nsize), sizeof(uint32_t));
    infile.read(reinterpret_cast<char *>(res), nsize * sizeof(double));
    infile.close();

    //std::cout << "Valeev: Read from " << ss.str() << "\n";
    //for(int i = 0; i < nsize; i++)
    //    std::cout << i << " : " << res[i] << "\n";

    return nsize;
}



void ValeevIntegrals(const VecQuartet & gshells,
                     double * integrals, bool normalize)
{
    int inorm = (normalize ? 1 : 0);
    const gaussian_shell * A = gshells[0].data();
    const gaussian_shell * B = gshells[1].data();
    const gaussian_shell * C = gshells[2].data();
    const gaussian_shell * D = gshells[3].data();

    const int nshell1 = gshells[0].size();
    const int nshell2 = gshells[1].size();
    const int nshell3 = gshells[2].size();
    const int nshell4 = gshells[3].size();
    const int n234 = nshell2 * nshell3 * nshell4 * NCART(gshells[0][0].am) * NCART(gshells[1][0].am) * NCART(gshells[2][0].am) * NCART(gshells[3][0].am);

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
            double vA[3] = { A[i].x, A[i].y, A[i].z };
            double vB[3] = { B[j].x, B[j].y, B[j].z };
            double vC[3] = { C[k].x, C[k].y, C[k].z };
            double vD[3] = { D[l].x, D[l].y, D[l].z };

            std::array<int, 3> g1{gshells[0][i].am, 0, 0};
            do
            {
                std::array<int, 3> g2{gshells[1][j].am, 0, 0};
                do
                {
                    std::array<int, 3> g3{gshells[2][k].am, 0, 0};
                    do
                    {
                        std::array<int, 3> g4{gshells[3][l].am, 0, 0};
                        do
                        {
                            double myint = 0.0;

                            for(int m = 0; m < A[i].nprim; m++)
                            for(int n = 0; n < B[j].nprim; n++)
                            for(int o = 0; o < C[k].nprim; o++)
                            for(int p = 0; p < D[l].nprim; p++)
                            {
                                
                                double val = Valeev_eri(g1[0], g1[1], g1[2], A[i].alpha[m], vA,
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


void ERDIntegrals(const VecQuartet & gshells,
                  double * integrals)
{
    const gaussian_shell * A = gshells[0].data();
    const gaussian_shell * B = gshells[1].data();
    const gaussian_shell * C = gshells[2].data();
    const gaussian_shell * D = gshells[3].data();

    const int nshell1 = gshells[0].size();
    const int nshell2 = gshells[1].size();
    const int nshell3 = gshells[2].size();
    const int nshell4 = gshells[3].size();
    const int nshell1234 = nshell1 * nshell2 * nshell3 * nshell4;

    const int am1 = A[0].am;
    const int am2 = B[0].am;
    const int am3 = C[0].am;
    const int am4 = D[0].am;
    const int ncart1234 = NCART(am1) * NCART(am2) * NCART(am3) * NCART(am4);


    const int ncart = nshell1234 * ncart1234; 
    std::fill(integrals, integrals + ncart, 0.0);

    int idx = 0;
    for(int i = 0; i < nshell1; i++)
    for(int j = 0; j < nshell2; j++)
    for(int k = 0; k < nshell3; k++)
    for(int l = 0; l < nshell4; l++)
    {
        ERD_Compute_shell(A[i], B[j], C[k], D[l], integrals + idx);
        idx += ncart1234;
    }
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


std::pair<double, double> CalcError(double const * const calc, double const * const ref, int ncalc)
{
    double maxerr = 0;
    double maxrelerr = 0;

    for(int i = 0; i < ncalc; i++)
    {
        const double r = ref[i];
        const double diff = fabs(calc[i] - r);
        const double rel = fabs(diff / r);
        if(diff > maxerr)
            maxerr = diff;
        if(rel > maxrelerr)
            maxrelerr = rel;
    }
 
    return std::pair<double, double>(maxerr, maxrelerr);   
}


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

int eriflat_notyetimplemented(struct multishell_pair_flat const P,
                              struct multishell_pair_flat const Q,
                              double * const restrict dummy)
{
    printf("****************************\n");
    printf("*** NOT YET IMPLEMENTED! ***\n");
    printf("***  ( %2d %2d | %2d %2d )   ***\n", P.am1, P.am2, Q.am1, Q.am2);
    printf("****************************\n");
    exit(1);
    return 0;
}

void Init_Test(void)
{
    for(int i = 0; i <= MAXAM; i++)
    for(int j = 0; j <= MAXAM; j++)
    for(int k = 0; k <= MAXAM; k++)
    for(int l = 0; l <= MAXAM; l++)
    {
        funcs_FO[i][j][k][l] = eri_notyetimplemented;
        funcs_vref[i][j][k][l] = eri_notyetimplemented;
        funcs_vref_flat[i][j][k][l] = eriflat_notyetimplemented;
    }


    funcs_FO[0][0][0][0] = eri_FO_s_s_s_s;
    #if MAXAM > 0
    funcs_FO[1][0][0][0] = eri_FO_p_s_s_s;
    funcs_FO[1][0][1][0] = eri_FO_p_s_p_s;
    funcs_FO[1][1][0][0] = eri_FO_p_p_s_s;
    funcs_FO[1][1][1][0] = eri_FO_p_p_p_s;
    funcs_FO[1][1][1][1] = eri_FO_p_p_p_p;
    #endif
    #if MAXAM > 1
    funcs_FO[2][0][0][0] = eri_FO_d_s_s_s;
    funcs_FO[2][0][1][0] = eri_FO_d_s_p_s;
    funcs_FO[2][0][1][1] = eri_FO_d_s_p_p;
    funcs_FO[2][0][2][0] = eri_FO_d_s_d_s;
    funcs_FO[2][1][0][0] = eri_FO_d_p_s_s;
    funcs_FO[2][1][1][0] = eri_FO_d_p_p_s;
    funcs_FO[2][1][1][1] = eri_FO_d_p_p_p;
    funcs_FO[2][1][2][0] = eri_FO_d_p_d_s;
    funcs_FO[2][1][2][1] = eri_FO_d_p_d_p;
    funcs_FO[2][2][0][0] = eri_FO_d_d_s_s;
    funcs_FO[2][2][1][0] = eri_FO_d_d_p_s;
    funcs_FO[2][2][1][1] = eri_FO_d_d_p_p;
    funcs_FO[2][2][2][0] = eri_FO_d_d_d_s;
    funcs_FO[2][2][2][1] = eri_FO_d_d_d_p;
    funcs_FO[2][2][2][2] = eri_FO_d_d_d_d;
    #endif

    funcs_vref[0][0][0][0] = eri_vref_s_s_s_s;
    #if MAXAM > 0
    funcs_vref[1][0][0][0] = eri_vref_p_s_s_s;
    funcs_vref[1][0][1][0] = eri_vref_p_s_p_s;
    funcs_vref[1][1][0][0] = eri_vref_p_p_s_s;
    funcs_vref[1][1][1][0] = eri_vref_p_p_p_s;
    funcs_vref[1][1][1][1] = eri_vref_p_p_p_p;
    #endif
    #if MAXAM > 1
    funcs_vref[2][0][0][0] = eri_vref_d_s_s_s;
    funcs_vref[2][0][1][0] = eri_vref_d_s_p_s;
    funcs_vref[2][0][1][1] = eri_vref_d_s_p_p;
    funcs_vref[2][0][2][0] = eri_vref_d_s_d_s;
    funcs_vref[2][1][0][0] = eri_vref_d_p_s_s;
    funcs_vref[2][1][1][0] = eri_vref_d_p_p_s;
    funcs_vref[2][1][1][1] = eri_vref_d_p_p_p;
    funcs_vref[2][1][2][0] = eri_vref_d_p_d_s;
    funcs_vref[2][1][2][1] = eri_vref_d_p_d_p;
    funcs_vref[2][2][0][0] = eri_vref_d_d_s_s;
    funcs_vref[2][2][1][0] = eri_vref_d_d_p_s;
    funcs_vref[2][2][1][1] = eri_vref_d_d_p_p;
    funcs_vref[2][2][2][0] = eri_vref_d_d_d_s;
    funcs_vref[2][2][2][1] = eri_vref_d_d_d_p;
    funcs_vref[2][2][2][2] = eri_vref_d_d_d_d;
    #endif
    #if MAXAM > 2
        funcs_vref[3][0][0][0] = eri_vref_f_s_s_s;
        funcs_vref[3][0][1][0] = eri_vref_f_s_p_s;
        funcs_vref[3][0][1][1] = eri_vref_f_s_p_p;
        funcs_vref[3][0][2][0] = eri_vref_f_s_d_s;
        funcs_vref[3][0][2][1] = eri_vref_f_s_d_p;
        funcs_vref[3][0][3][0] = eri_vref_f_s_f_s;
        funcs_vref[3][1][0][0] = eri_vref_f_p_s_s;
        funcs_vref[3][1][1][0] = eri_vref_f_p_p_s;
        funcs_vref[3][1][1][1] = eri_vref_f_p_p_p;
        funcs_vref[3][1][2][0] = eri_vref_f_p_d_s;
        funcs_vref[3][1][2][1] = eri_vref_f_p_d_p;
        funcs_vref[3][1][2][2] = eri_vref_f_p_d_d;
        funcs_vref[3][1][3][0] = eri_vref_f_p_f_s;
        funcs_vref[3][1][3][1] = eri_vref_f_p_f_p;
        funcs_vref[3][2][0][0] = eri_vref_f_d_s_s;
        funcs_vref[3][2][1][0] = eri_vref_f_d_p_s;
        funcs_vref[3][2][1][1] = eri_vref_f_d_p_p;
        funcs_vref[3][2][2][0] = eri_vref_f_d_d_s;
        funcs_vref[3][2][2][1] = eri_vref_f_d_d_p;
        funcs_vref[3][2][2][2] = eri_vref_f_d_d_d;
        funcs_vref[3][2][3][0] = eri_vref_f_d_f_s;
        funcs_vref[3][2][3][1] = eri_vref_f_d_f_p;
        funcs_vref[3][2][3][2] = eri_vref_f_d_f_d;
        funcs_vref[3][3][0][0] = eri_vref_f_f_s_s;
        funcs_vref[3][3][1][0] = eri_vref_f_f_p_s;
        funcs_vref[3][3][1][1] = eri_vref_f_f_p_p;
        funcs_vref[3][3][2][0] = eri_vref_f_f_d_s;
        funcs_vref[3][3][2][1] = eri_vref_f_f_d_p;
        funcs_vref[3][3][2][2] = eri_vref_f_f_d_d;
        funcs_vref[3][3][3][0] = eri_vref_f_f_f_s;
        funcs_vref[3][3][3][1] = eri_vref_f_f_f_p;
        funcs_vref[3][3][3][2] = eri_vref_f_f_f_d;
        funcs_vref[3][3][3][3] = eri_vref_f_f_f_f;
    #endif

    funcs_vref_flat[0][0][0][0] = eri_vref_flat_s_s_s_s;
    #if MAXAM > 0
    funcs_vref_flat[1][0][0][0] = eri_vref_flat_p_s_s_s;
    funcs_vref_flat[1][0][1][0] = eri_vref_flat_p_s_p_s;
    funcs_vref_flat[1][1][0][0] = eri_vref_flat_p_p_s_s;
    funcs_vref_flat[1][1][1][0] = eri_vref_flat_p_p_p_s;
    funcs_vref_flat[1][1][1][1] = eri_vref_flat_p_p_p_p;
    #endif
    #if MAXAM > 1
    funcs_vref_flat[2][0][0][0] = eri_vref_flat_d_s_s_s;
    funcs_vref_flat[2][0][1][0] = eri_vref_flat_d_s_p_s;
    funcs_vref_flat[2][0][1][1] = eri_vref_flat_d_s_p_p;
    funcs_vref_flat[2][0][2][0] = eri_vref_flat_d_s_d_s;
    funcs_vref_flat[2][1][0][0] = eri_vref_flat_d_p_s_s;
    funcs_vref_flat[2][1][1][0] = eri_vref_flat_d_p_p_s;
    funcs_vref_flat[2][1][1][1] = eri_vref_flat_d_p_p_p;
    funcs_vref_flat[2][1][2][0] = eri_vref_flat_d_p_d_s;
    funcs_vref_flat[2][1][2][1] = eri_vref_flat_d_p_d_p;
    funcs_vref_flat[2][2][0][0] = eri_vref_flat_d_d_s_s;
    funcs_vref_flat[2][2][1][0] = eri_vref_flat_d_d_p_s;
    funcs_vref_flat[2][2][1][1] = eri_vref_flat_d_d_p_p;
    funcs_vref_flat[2][2][2][0] = eri_vref_flat_d_d_d_s;
    funcs_vref_flat[2][2][2][1] = eri_vref_flat_d_d_d_p;
    funcs_vref_flat[2][2][2][2] = eri_vref_flat_d_d_d_d;
    #endif
    #if MAXAM > 2
        funcs_vref_flat[3][0][0][0] = eri_vref_flat_f_s_s_s;
        funcs_vref_flat[3][0][1][0] = eri_vref_flat_f_s_p_s;
        funcs_vref_flat[3][0][1][1] = eri_vref_flat_f_s_p_p;
        funcs_vref_flat[3][0][2][0] = eri_vref_flat_f_s_d_s;
        funcs_vref_flat[3][0][2][1] = eri_vref_flat_f_s_d_p;
        funcs_vref_flat[3][0][3][0] = eri_vref_flat_f_s_f_s;
        funcs_vref_flat[3][1][0][0] = eri_vref_flat_f_p_s_s;
        funcs_vref_flat[3][1][1][0] = eri_vref_flat_f_p_p_s;
        funcs_vref_flat[3][1][1][1] = eri_vref_flat_f_p_p_p;
        funcs_vref_flat[3][1][2][0] = eri_vref_flat_f_p_d_s;
        funcs_vref_flat[3][1][2][1] = eri_vref_flat_f_p_d_p;
        funcs_vref_flat[3][1][2][2] = eri_vref_flat_f_p_d_d;
        funcs_vref_flat[3][1][3][0] = eri_vref_flat_f_p_f_s;
        funcs_vref_flat[3][1][3][1] = eri_vref_flat_f_p_f_p;
        funcs_vref_flat[3][2][0][0] = eri_vref_flat_f_d_s_s;
        funcs_vref_flat[3][2][1][0] = eri_vref_flat_f_d_p_s;
        funcs_vref_flat[3][2][1][1] = eri_vref_flat_f_d_p_p;
        funcs_vref_flat[3][2][2][0] = eri_vref_flat_f_d_d_s;
        funcs_vref_flat[3][2][2][1] = eri_vref_flat_f_d_d_p;
        funcs_vref_flat[3][2][2][2] = eri_vref_flat_f_d_d_d;
        funcs_vref_flat[3][2][3][0] = eri_vref_flat_f_d_f_s;
        funcs_vref_flat[3][2][3][1] = eri_vref_flat_f_d_f_p;
        funcs_vref_flat[3][2][3][2] = eri_vref_flat_f_d_f_d;
        funcs_vref_flat[3][3][0][0] = eri_vref_flat_f_f_s_s;
        funcs_vref_flat[3][3][1][0] = eri_vref_flat_f_f_p_s;
        funcs_vref_flat[3][3][1][1] = eri_vref_flat_f_f_p_p;
        funcs_vref_flat[3][3][2][0] = eri_vref_flat_f_f_d_s;
        funcs_vref_flat[3][3][2][1] = eri_vref_flat_f_f_d_p;
        funcs_vref_flat[3][3][2][2] = eri_vref_flat_f_f_d_d;
        funcs_vref_flat[3][3][3][0] = eri_vref_flat_f_f_f_s;
        funcs_vref_flat[3][3][3][1] = eri_vref_flat_f_f_f_p;
        funcs_vref_flat[3][3][3][2] = eri_vref_flat_f_f_f_d;
        funcs_vref_flat[3][3][3][3] = eri_vref_flat_f_f_f_f;
    #endif
}

int Integral_FO(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict integrals)
{
    return funcs_FO[P.am1][P.am2][Q.am1][Q.am2](P, Q, integrals);
}

int Integral_vref(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict integrals)
{
    return funcs_vref[P.am1][P.am2][Q.am1][Q.am2](P, Q, integrals);
}

int Integral_vref_flat(struct multishell_pair_flat const P, struct multishell_pair_flat const Q, double * const restrict integrals)
{
    return funcs_vref_flat[P.am1][P.am2][Q.am1][Q.am2](P, Q, integrals);
}

