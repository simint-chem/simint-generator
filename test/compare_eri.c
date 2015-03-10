#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "valeev/valeev.h"
#include "eri/eri.h"
#include "boys/boys.h"

#define MAX_SHELL 10

struct gaussian_shell
random_shell(int nprim)
{
    struct gaussian_shell G;
    G.am = 0;
    G.x = 4.0 * rand() / ((double)RAND_MAX) - 2.0;
    G.y = 4.0 * rand() / ((double)RAND_MAX) - 2.0;
    G.z = 4.0 * rand() / ((double)RAND_MAX) - 2.0;

    G.nprim = nprim;

    G.alpha = ALLOC(nprim * sizeof(double));
    G.coef = ALLOC(nprim * sizeof(double));

    for(int i = 0; i < nprim; i++)
    {
        G.alpha[i] = 10.0 * rand() / ((double)RAND_MAX);
        G.coef[i] = 10.0 * rand() / ((double)RAND_MAX);
    }

    return G;
         
}

void free_random_shell(struct gaussian_shell G)
{
    FREE(G.alpha);
    FREE(G.coef);
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
    int nprim1234 = nprim1 * nprim2 * nprim3 * nprim4;

    double * res1_s = ALLOC(nshell1234 * sizeof(double));
    double * res1_m = ALLOC(nshell1234 * sizeof(double));
    double * res2_s = ALLOC(nshell1234 * sizeof(double));
    double * res2_m = ALLOC(nshell1234 * sizeof(double));
    double * vres = ALLOC(nshell1234 * sizeof(double));
    double * intwork = ALLOC(nprim1234 * nshell1234 * sizeof(double));

    srand(time(NULL));

    struct gaussian_shell A[MAX_SHELL];
    struct gaussian_shell B[MAX_SHELL];
    struct gaussian_shell C[MAX_SHELL];
    struct gaussian_shell D[MAX_SHELL];

    for(int i = 0; i < nshell1; i++)
        A[i] = random_shell(nprim1);

    for(int i = 0; i < nshell2; i++)
        B[i] = random_shell(nprim2);

    for(int i = 0; i < nshell3; i++)
        C[i] = random_shell(nprim3);

    for(int i = 0; i < nshell4; i++)
        D[i] = random_shell(nprim4);


    Boys_Init();
    Valeev_Init();

    struct shell_pair P = create_ss_shell_pair(nshell1, A, nshell2, B);
    struct shell_pair Q = create_ss_shell_pair(nshell3, C, nshell4, D);


    // Actually calculate
    eri_1pair_ssss_multi(nshell1, A, nshell2, B, Q, res1_m, intwork);
    eri_2pair_ssss_multi(P, Q, res2_m, intwork);
    free_shell_pair(P);
    free_shell_pair(Q);


    // do one single run
    P = create_ss_shell_pair(1, A, 1, B);
    Q = create_ss_shell_pair(1, C, 1, D);
    eri_1pair_ssss_single(A[0], B[0], Q, res1_s);
    eri_2pair_ssss_single(P, Q, res2_s);
    free_shell_pair(P);
    free_shell_pair(Q);


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

        vres[idx] = 0.0;

        for(int m = 0; m < A[i].nprim; m++)
        for(int n = 0; n < B[j].nprim; n++)
        for(int o = 0; o < C[k].nprim; o++)
        for(int p = 0; p < D[l].nprim; p++)
        {
            double val = Valeev_eri(0, 0, 0, A[i].alpha[m], vA,
                                    0, 0, 0, B[j].alpha[n], vB,
                                    0, 0, 0, C[k].alpha[o], vC,
                                    0, 0, 0, D[l].alpha[p], vD, 1);
            vres[idx] += val * A[i].coef[m] * B[j].coef[n] * C[k].coef[o] * D[l].coef[p];
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

    Boys_Finalize();
    Valeev_Finalize();


    printf("%11s  %11s  %11s  %11s  %11s  --  %11s  %11s  %11s  %11s\n", "vres", "res1_s", "res1_m", "res2_s", "res2_m", 
                                                                    "diff1_s", "diff1_m", "diff2_s", "diff2_m"); 
    for(int i = 0; i < nshell1234; i++)
    {
        double diff1_s = fabs(res1_s[i] - vres[i]);
        double diff1_m = fabs(res1_m[i] - vres[i]);
        double diff2_s = fabs(res2_s[i] - vres[i]);
        double diff2_m = fabs(res2_m[i] - vres[i]);
        if(diff1_s > 1e-14 || diff1_m > 1e-14 || diff2_s > 1e-14 || diff2_m > 1e-14)
          printf("%11.4e  %11.4e  %11.4e  %11.4e  %11.4e  --  %11.4e  %11.4e  %11.4e  %11.4e\n", 
                                                        vres[i], res1_s[i], res1_m[i], res2_s[i], res2_m[i],
                                                        diff1_s, diff1_m, diff2_s, diff2_m);
    }

    for(int i = 0; i < nshell1; i++)
        free_random_shell(A[i]);

    for(int j = 0; j < nshell2; j++)
        free_random_shell(B[j]);

    for(int k = 0; k < nshell3; k++)
        free_random_shell(C[k]);

    for(int l = 0; l < nshell4; l++)
        free_random_shell(D[l]);

    FREE(res1_s); FREE(res1_m); FREE(res2_s); FREE(res2_m); FREE(vres); FREE(intwork);

    return 0;
}
