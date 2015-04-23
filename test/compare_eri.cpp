#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "vectorization.h"
#include "liberd/erd_interface.h"
#include "eri/eri.h"
#include "boys/boys.h"

#include "valeev/valeev.h"

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
    double * res_split           = ALLOC(nshell1234 * sizeof(double));
    double * res_splitcombined   = ALLOC(nshell1234 * sizeof(double));
    double * res_taylorcombined  = ALLOC(nshell1234 * sizeof(double));
    double * res_FOcombined      = ALLOC(nshell1234 * sizeof(double));

    double * res_liberd          = ALLOC(nshell1234 * sizeof(double));
    double * res_valeev      = ALLOC(nshell1234 * sizeof(double));

    // for split
    // no need to round each pair
    int worksize = SIMD_ROUND_DBL(nshell1*nprim1 * nshell2*nprim2 * nshell3*nprim3 * nshell4*nprim4);
    double * intwork1 = ALLOC(3 * worksize * sizeof(double));
    double * intwork2 = ALLOC(3 * worksize * sizeof(double));

    // allocate gaussian shell memory
    struct gaussian_shell * A = ALLOC(nshell1 * sizeof(struct gaussian_shell));
    struct gaussian_shell * B = ALLOC(nshell2 * sizeof(struct gaussian_shell));
    struct gaussian_shell * C = ALLOC(nshell3 * sizeof(struct gaussian_shell));
    struct gaussian_shell * D = ALLOC(nshell4 * sizeof(struct gaussian_shell));

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
    normalize_gaussian_shells(nshell4, D);

    // find the maximum possible x value for the boys function
    const double maxR2 = 12.0 * MAX_COORD * MAX_COORD;
    const double max_x = maxR2 * (MAX_EXP*MAX_EXP) / (2.0 * MAX_EXP);
    printf("Maximum parameter to boys: %12.8e\n", max_x);
    Boys_Init(max_x, 7); // need F0 + 7 for interpolation
    printf("Boys grid generated\n");

    Valeev_Init();

    // Actually calculate
    struct multishell_pair P = create_ss_multishell_pair(nshell1, A, nshell2, B);
    struct multishell_pair Q = create_ss_multishell_pair(nshell3, C, nshell4, D);
    eri_taylorcombined_ssss(  P, Q, res_taylorcombined                      );
    eri_split_ssss(           P, Q, res_split,          intwork1, intwork2  );
    eri_splitcombined_ssss(   P, Q, res_splitcombined                       );
    eri_FOcombined_ssss(      P, Q, res_FOcombined                          );
    free_multishell_pair(P);
    free_multishell_pair(Q);


    // test with valeev
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

        res_valeev[idx] = 0.0;

        for(int m = 0; m < A[i].nprim; m++)
        for(int n = 0; n < B[j].nprim; n++)
        for(int o = 0; o < C[k].nprim; o++)
        for(int p = 0; p < D[l].nprim; p++)
        {
            double val = Valeev_eri(0, 0, 0, A[i].alpha[m], vA,
                                    0, 0, 0, B[j].alpha[n], vB,
                                    0, 0, 0, C[k].alpha[o], vC,
                                    0, 0, 0, D[l].alpha[p], vD, 0);
            res_valeev[idx] += val * A[i].coef[m] * B[j].coef[n] * C[k].coef[o] * D[l].coef[p];
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


    // test with erd
    ERD_Init(nshell1, A, nshell2, B,
             nshell3, C, nshell4, D);

    idx = 0;
    for(int i = 0; i < nshell1; i++)
    for(int j = 0; j < nshell2; j++)
    for(int k = 0; k < nshell3; k++)
    for(int l = 0; l < nshell4; l++)
    {
        double erdval[3] = {0,0,0};
        ERD_Compute_shell(A[i], B[j], C[k], D[l], erdval);
        res_liberd[idx] = erdval[0];
        idx++;
    }


    Boys_Finalize();
    Valeev_Finalize();
    ERD_Finalize();


    printf("%22s %22s %22s %22s %22s %22s\n",
           "liberd", "split", "splitcombined", "taylorcombined", "FOcombined", "valeev");

    for(int i = 0; i < nshell1234; i++)
    {
        const double v = res_valeev[i];

        double diff_liberd         = fabs(res_liberd[i]         - v);
        double diff_split          = fabs(res_split[i]          - v);
        double diff_splitcombined  = fabs(res_splitcombined[i]  - v);
        double diff_taylorcombined = fabs(res_taylorcombined[i] - v);
        double diff_FOcombined     = fabs(res_FOcombined[i]     - v);

        printf("%22.4e  %22.4e  %22.4e  %22.4e  %22.4e  %22.4e\n",
                                  res_liberd[i], res_split[i], res_splitcombined[i],
                                  res_taylorcombined[i], res_FOcombined[i],
                                  res_valeev[i]);

        printf("%22.4e  %22.4e  %22.4e  %22.4e  %22.4e\n",
                                  diff_liberd, diff_split, diff_splitcombined,
                                  diff_taylorcombined, diff_FOcombined);

        printf("\n");
    }

    for(int i = 0; i < nshell1; i++)
        free_random_shell(A[i]);

    for(int j = 0; j < nshell2; j++)
        free_random_shell(B[j]);

    for(int k = 0; k < nshell3; k++)
        free_random_shell(C[k]);

    for(int l = 0; l < nshell4; l++)
        free_random_shell(D[l]);

    FREE(res_valeev);
    FREE(res_liberd);
    FREE(res_split);
    FREE(res_splitcombined);
    FREE(res_taylorcombined);
    FREE(res_FOcombined);
    FREE(A); FREE(B); FREE(C); FREE(D);
    FREE(intwork1); FREE(intwork2);

    return 0;
}
