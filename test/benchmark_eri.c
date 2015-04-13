#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "liberd/erd_interface.h"
#include "eri/eri.h"
#include "boys/boys.h"

#define MAX_COORD 0.5
#define MAX_EXP   50.0
#define MAX_COEF  2.0

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

void free_random_shell(struct gaussian_shell G)
{
    free_gaussian_shell(G);
}


int main(int argc, char ** argv)
{
    if(argc != 9)
    {
        printf("Give me 8 arguments! I got %d\n", argc-1);
        return 1;
    }

    const int ntest = 500;

    int nshell1 = atoi(argv[1]);
    int nshell2 = atoi(argv[2]);
    int nshell3 = atoi(argv[3]);
    int nshell4 = atoi(argv[4]);
    int nprim1 = atoi(argv[5]);
    int nprim2 = atoi(argv[6]);
    int nprim3 = atoi(argv[7]);
    int nprim4 = atoi(argv[8]);


    int nshell1234 = nshell1 * nshell2 * nshell3 * nshell4;

    /* Storage of test results */
    double * res_erf         = ALLOC(nshell1234 * sizeof(double));
    double * res_split       = ALLOC(nshell1234 * sizeof(double));
    double * res_taylor      = ALLOC(nshell1234 * sizeof(double));
    double * res_erf_comb    = ALLOC(nshell1234 * sizeof(double));
    double * res_taylor_comb = ALLOC(nshell1234 * sizeof(double));
    double * res_cheby       = ALLOC(nshell1234 * sizeof(double));
    double * res_liberd      = ALLOC(nshell1234 * sizeof(double));

    // no need to round each pair
    int worksize = SIMD_ROUND_DBL(nshell1*nprim1 * nshell2*nprim2 * nshell3*nprim3 * nshell4*nprim4);
    double * intwork1 = ALLOC(3*worksize * sizeof(double));
    double * intwork2 = ALLOC(3*worksize * sizeof(double));

    srand(time(NULL));

    struct gaussian_shell * A = ALLOC(nshell1 * sizeof(struct gaussian_shell));
    struct gaussian_shell * B = ALLOC(nshell2 * sizeof(struct gaussian_shell));
    struct gaussian_shell * C = ALLOC(nshell3 * sizeof(struct gaussian_shell));
    struct gaussian_shell * D = ALLOC(nshell4 * sizeof(struct gaussian_shell));

    // find the maximum possible x value for the boys function
    const double maxR2 = 12.0 * MAX_COORD * MAX_COORD;
    const double max_x = maxR2 * (MAX_EXP*MAX_EXP) / (2.0 * MAX_EXP);
    printf("Maximum parameter to boys: %12.8e\n", max_x);
    //Boys_Init(max_x, 7); // need F0 + 7 for interpolation
    Boys_Init(1e6, BOYS_GRID_MAXN); // need F0 + 7 for interpolation
    printf("Boys grid generated\n");


    printf("%17s %17s %17s %17s %17s %17s %17s\n",
                "res_liberd", "res_erf", "res_split", "res_taylor",
                "res_erf_comb", "res_taylor_comb", "res_cheby");

    for(int n = 0; n < ntest; n++)
    {
        for(int i = 0; i < nshell1; i++)
            A[i] = random_shell(nprim1);

        for(int i = 0; i < nshell2; i++)
            B[i] = random_shell(nprim2);

        for(int i = 0; i < nshell3; i++)
            C[i] = random_shell(nprim3);

        for(int i = 0; i < nshell4; i++)
            D[i] = random_shell(nprim4);

        // normalize the shells
        normalize_gaussian_shells(nshell1, A);
        normalize_gaussian_shells(nshell2, B);
        normalize_gaussian_shells(nshell3, C);

        // Actually calculate
        struct multishell_pair P = create_ss_multishell_pair(nshell1, A, nshell2, B);
        struct multishell_pair Q = create_ss_multishell_pair(nshell3, C, nshell4, D);
        eri_erf_multi_ssss(           P, Q, res_erf,         intwork1, intwork2);
        eri_erf_split_ssss(           P, Q, res_split,       intwork1, intwork2);
        eri_erf_combined_ssss(        P, Q, res_erf_comb,    intwork1, intwork2);
        eri_cheby_ssss(               P, Q, res_cheby,       intwork1, intwork2);
        eri_taylor_ssss(              P, Q, res_taylor,      intwork1, intwork2);
        eri_taylorcombined_ssss(      P, Q, res_taylor_comb, intwork1, intwork2);
        free_multishell_pair(P);
        free_multishell_pair(Q);

        // test with erd
        int idx = 0;
        ERD_Init(nshell1, A, nshell2, B,
                 nshell3, C, nshell4, D);
        for(int i = 0; i < nshell1; i++)
        for(int j = 0; j < nshell2; j++)
        for(int k = 0; k < nshell3; k++)
        for(int l = 0; l < nshell4; l++)
        {
            res_liberd[idx] = 0.0;
            ERD_Compute_shell(A[i], B[j], C[k], D[l], res_liberd + idx);
            idx++;
        }
        ERD_Finalize();


        // print some results
        printf("%17e %17e %17e %17e %17e %17e %17e\n",
                    res_liberd[0], res_erf[0], res_split[0], res_taylor[0],
                    res_erf_comb[0], res_taylor_comb[0], res_cheby[0]);

        // free memory
        for(int i = 0; i < nshell1; i++)
            free_random_shell(A[i]);

        for(int i = 0; i < nshell2; i++)
            free_random_shell(B[i]);

        for(int i = 0; i < nshell3; i++)
            free_random_shell(C[i]);

        for(int i = 0; i < nshell4; i++)
            free_random_shell(D[i]);


    }


    Boys_Finalize();

    FREE(res_liberd);
    FREE(res_erf);
    FREE(res_split);
    FREE(res_taylor);
    FREE(res_erf_comb);
    FREE(res_taylor_comb);
    FREE(res_cheby);
    FREE(A); FREE(B); FREE(C); FREE(D);
    FREE(intwork1); FREE(intwork2);

    return 0;
}
