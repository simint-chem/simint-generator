#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "valeev/valeev.h"
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

    int nshell1 = atoi(argv[1]);
    int nshell2 = atoi(argv[2]);
    int nshell3 = atoi(argv[3]);
    int nshell4 = atoi(argv[4]);
    int nprim1 = atoi(argv[5]);
    int nprim2 = atoi(argv[6]);
    int nprim3 = atoi(argv[7]);
    int nprim4 = atoi(argv[8]);


    int nshell1234 = nshell1 * nshell2 * nshell3 * nshell4;

    double * res_f = ALLOC(nshell1234 * sizeof(double));
    double * res_fs = ALLOC(nshell1234 * sizeof(double));
    double * res_ft = ALLOC(nshell1234 * sizeof(double));
    double * res_fc = ALLOC(nshell1234 * sizeof(double));
    double * res_ftc = ALLOC(nshell1234 * sizeof(double));
    double * res_fch = ALLOC(nshell1234 * sizeof(double));
    double * res_e = ALLOC(nshell1234 * sizeof(double));
    double * res_v = ALLOC(nshell1234 * sizeof(double));

    int worksize = SIMD_ROUND_DBL(nshell1*nprim1 * nshell2*nprim2 * nshell3*nprim3 * nshell4*nprim4);
    double * intwork1 = ALLOC(3*worksize * sizeof(double));
    double * intwork2 = ALLOC(3*worksize * sizeof(double));

    srand(time(NULL));

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
    eri_ssss(P, Q, res_f, intwork1, intwork2);
    eri_ssss_cheby(P, Q, res_fch, intwork1, intwork2);
    eri_ssss_split(P, Q, res_fs, intwork1, intwork2);
    eri_ssss_taylor(P, Q, res_ft, intwork1, intwork2);
    eri_ssss_combined(P, Q, res_fc, intwork1, intwork2);
    eri_ssss_taylorcombined(P, Q, res_ftc, intwork1, intwork2);
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

        res_v[idx] = 0.0;

        for(int m = 0; m < A[i].nprim; m++)
        for(int n = 0; n < B[j].nprim; n++)
        for(int o = 0; o < C[k].nprim; o++)
        for(int p = 0; p < D[l].nprim; p++)
        {
            double val = Valeev_eri(0, 0, 0, A[i].alpha[m], vA,
                                    0, 0, 0, B[j].alpha[n], vB,
                                    0, 0, 0, C[k].alpha[o], vC,
                                    0, 0, 0, D[l].alpha[p], vD, 0);
            res_v[idx] += val * A[i].coef[m] * B[j].coef[n] * C[k].coef[o] * D[l].coef[p];
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
        res_e[idx] = 0.0;
        ERD_Compute_shell(A[i], B[j], C[k], D[l], res_e + idx);
        idx++;
    }


    Boys_Finalize();
    Valeev_Finalize();
    ERD_Finalize();


    printf("%11s  %11s  %11s  %11s  %11s  %11s  %11s\n",
                        "diff_e", "diff_f", "diff_fs", "diff_ft", "diff_fc", "diff_ftc", "diff_fch");
    for(int i = 0; i < nshell1234; i++)
    {
        double diff_e = fabs(res_e[i] - res_v[i]);
        double diff_f = fabs(res_f[i] - res_v[i]);
        double diff_fs = fabs(res_fs[i] - res_v[i]);
        double diff_ft = fabs(res_ft[i] - res_v[i]);
        double diff_fc = fabs(res_fc[i] - res_v[i]);
        double diff_fch = fabs(res_fch[i] - res_v[i]);
        double diff_ftc = fabs(res_ftc[i] - res_v[i]);
        //printf("%11.4e  %11.4e  %11.4e  %11.4e  %11.4e  %11.4e  %11.4e  %11.4e\n", res_e[i], res_f[i], res_fs[i], res_ft[i], res_fc[i], res_ftc[i], res_fch[i], res_v[i]);
        printf("%11.4e  %11.4e  %11.4e  %11.4e  %11.4e  %11.4e  %11.4e\n", diff_e, diff_f, diff_fs, diff_ft, diff_fc, diff_ftc, diff_fch);
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

    FREE(res_v); FREE(res_e);
    FREE(res_f); FREE(res_fs);
    FREE(res_ft); FREE(res_fc);
    FREE(res_ftc); FREE(res_fch);
    FREE(A); FREE(B); FREE(C); FREE(D);
    FREE(intwork1); FREE(intwork2);

    return 0;
}
