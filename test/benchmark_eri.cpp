#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "vectorization.h"
#include "eri/eri.h"
#include "boys/boys.h"

#include "erd_interface.hpp"

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
    srand(time(NULL));

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
    double * res_split           = (double *)ALLOC(nshell1234 * sizeof(double));
    double * res_splitcombined   = (double *)ALLOC(nshell1234 * sizeof(double));
    double * res_taylorcombined  = (double *)ALLOC(nshell1234 * sizeof(double));
    double * res_FOcombined      = (double *)ALLOC(nshell1234 * sizeof(double));

    double * res_liberd          = (double *)ALLOC(nshell1234 * sizeof(double));

    // for split
    // no need to round each pair
    int worksize = SIMD_ROUND_DBL(nshell1*nprim1 * nshell2*nprim2 * nshell3*nprim3 * nshell4*nprim4);
    double * intwork1 = (double *)ALLOC(3 * worksize * sizeof(double));
    double * intwork2 = (double *)ALLOC(3 * worksize * sizeof(double));

    // allocate gaussian shell memory
    struct gaussian_shell * A = (struct gaussian_shell *)ALLOC(nshell1 * sizeof(struct gaussian_shell));
    struct gaussian_shell * B = (struct gaussian_shell *)ALLOC(nshell2 * sizeof(struct gaussian_shell));
    struct gaussian_shell * C = (struct gaussian_shell *)ALLOC(nshell3 * sizeof(struct gaussian_shell));
    struct gaussian_shell * D = (struct gaussian_shell *)ALLOC(nshell4 * sizeof(struct gaussian_shell));

    // find the maximum possible x value for the boys function
    const double maxR2 = 12.0 * MAX_COORD * MAX_COORD;
    const double max_x = maxR2 * (MAX_EXP*MAX_EXP) / (2.0 * MAX_EXP);
    printf("Maximum parameter to boys: %12.8e\n", max_x);
    //Boys_Init(max_x, 7); // need F0 + 7 for interpolation
    Boys_Init(1e6, 10); // need F0 + 7 for interpolation
    printf("Boys initialized\n");


    printf("%22s %22s %22s %22s %22s\n",
           "liberd", "split", "splitcombined", "taylorcombined", "FOcombined");

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
        eri_taylorcombined_ssss(  P, Q, res_taylorcombined                      );
        eri_split_ssss(           P, Q, res_split,          intwork1, intwork2  );
        eri_splitcombined_ssss(   P, Q, res_splitcombined                       );
        eri_FOcombined_ssss(      P, Q, res_FOcombined                          );
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
        printf("%22e %22e %22e %22e %22e\n",
                    res_liberd[0], res_split[0], res_splitcombined[0], 
                    res_taylorcombined[0], res_FOcombined[0]);

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
    FREE(res_split);
    FREE(res_splitcombined);
    FREE(res_taylorcombined);
    FREE(res_FOcombined);
    FREE(A); FREE(B); FREE(C); FREE(D);
    FREE(intwork1); FREE(intwork2);

    return 0;
}
