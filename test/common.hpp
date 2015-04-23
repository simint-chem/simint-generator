#include <cstdlib>
#include <ctime>
#include <array>
#include <vector>

#include "eri/shell.h"

#include "valeev.hpp"
#include "erd_interface.hpp"
#include "cppvectorization.hpp"

#define MAX_COORD 0.5
#define MAX_EXP   50.0
#define MAX_COEF  2.0

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

typedef std::vector<gaussian_shell, AlignedAllocator<gaussian_shell>> AlignedGaussianVec;
typedef std::array<AlignedGaussianVec, 4> VecQuartet;


inline
struct gaussian_shell
random_shell(int nprim)
{
    struct gaussian_shell G;
    allocate_gaussian_shell(nprim, &G);

    G.am = 0;
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

inline
void free_random_shell(struct gaussian_shell G)
{
    free_gaussian_shell(G);
}


inline
VecQuartet
CreateRandomQuartets(std::array<int, 4> nshell,
                     std::array<int, 4> nprim,
                     bool normalize = true)
{
    VecQuartet arr;

    for(int q = 0; q < 4; q++)
    {
        arr[q].reserve(nshell[q]);

        for(int i = 0; i < nshell[q]; i++)
            arr[q].push_back(random_shell(nprim[q]));

        if(normalize)
            normalize_gaussian_shells(nshell[q], arr[q].data());
    }

    return arr;
}


inline
void
FreeRandomQuartets(VecQuartet & arr)
{
    for(auto & it : arr)
    {
        for(auto & it2 : it)
            free_random_shell(it2);
        it.clear();
    }
}

inline
void ValeevIntegrals(const VecQuartet & gshells,
                     double * integrals)
{
    const gaussian_shell * A = gshells[0].data();
    const gaussian_shell * B = gshells[1].data();
    const gaussian_shell * C = gshells[2].data();
    const gaussian_shell * D = gshells[3].data();

    int nshell1 = gshells[0].size();
    int nshell2 = gshells[1].size();
    int nshell3 = gshells[2].size();
    int nshell4 = gshells[3].size();

    int idx = 0;
    for(int i = 0; i < nshell1; i++)
    for(int j = 0; j < nshell2; j++)
    for(int k = 0; k < nshell3; k++)
    for(int l = 0; l < nshell4; l++)
    {
        double vA[3] = { A[i].x, A[i].y, A[i].z };
        double vB[3] = { B[j].x, B[j].y, B[j].z };
        double vC[3] = { C[k].x, C[k].y, C[k].z };
        double vD[3] = { D[l].x, D[l].y, D[l].z };

        integrals[idx] = 0.0;

        for(int m = 0; m < A[i].nprim; m++)
        for(int n = 0; n < B[j].nprim; n++)
        for(int o = 0; o < C[k].nprim; o++)
        for(int p = 0; p < D[l].nprim; p++)
        {
            double val = Valeev_eri(0, 0, 0, A[i].alpha[m], vA,
                                    0, 0, 0, B[j].alpha[n], vB,
                                    0, 0, 0, C[k].alpha[o], vC,
                                    0, 0, 0, D[l].alpha[p], vD, 0);
            integrals[idx] += val * A[i].coef[m] * B[j].coef[n] * C[k].coef[o] * D[l].coef[p];
            /*
            printf("IDX: %d\n", idx);
            printf("VAL: %8.3e\n", val);
            printf("%8.3e %8.3e %8.3e %8.3e\n", vA[0], vA[1], vA[2], A[i].alpha[m]);
            printf("%8.3e %8.3e %8.3e %8.3e\n", vB[0], vB[1], vB[2], B[j].alpha[n]);
            printf("%8.3e %8.3e %8.3e %8.3e\n", vC[0], vC[1], vC[2], C[k].alpha[o]);
            printf("%8.3e %8.3e %8.3e %8.3e\n", vD[0], vD[1], vD[2], D[l].alpha[p]);
            printf("%8.3e %8.3e %8.3e %8.3e\n", A[i].coef[m], B[j].coef[n], C[k].coef[o], D[l].coef[p]);
            */
        }

        idx++;
    }
}

void ERDIntegrals(const VecQuartet & gshells,
                  double * integrals)
{
    const gaussian_shell * A = gshells[0].data();
    const gaussian_shell * B = gshells[1].data();
    const gaussian_shell * C = gshells[2].data();
    const gaussian_shell * D = gshells[3].data();

    int nshell1 = gshells[0].size();
    int nshell2 = gshells[1].size();
    int nshell3 = gshells[2].size();
    int nshell4 = gshells[3].size();
   
    ERD_Init(nshell1, A, nshell2, B, nshell3, C, nshell4, D); 

    int ncart = 0;
    for(int i = 0; i < nshell1; i++)
    for(int j = 0; j < nshell2; j++)
    for(int k = 0; k < nshell3; k++)
    for(int l = 0; l < nshell4; l++)
        ncart += NCART(A[i].am)
               * NCART(B[j].am)
               * NCART(C[k].am)
               * NCART(D[l].am);

    std::fill(integrals, integrals + ncart, 0.0);

    int idx = 0;
    for(int i = 0; i < nshell1; i++)
    for(int j = 0; j < nshell2; j++)
    for(int k = 0; k < nshell3; k++)
    for(int l = 0; l < nshell4; l++)
    {
        ERD_Compute_shell(A[i], B[j], C[k], D[l], integrals + idx);
        idx += NCART(A[i].am)
             * NCART(B[j].am)
             * NCART(C[k].am)
             * NCART(D[l].am);
    }
}

