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

    double * res1_s = ALLOC(nshell1234 * sizeof(double));
    double * res1_m = ALLOC(nshell1234 * sizeof(double));
    double * res2_s = ALLOC(nshell1234 * sizeof(double));
    double * res2_m = ALLOC(nshell1234 * sizeof(double));

    int wrksize = nprim1*nprim2*SIMD_ROUND_DBL(nprim3*nprim4);
    double * intwork = ALLOC(nshell1234 * wrksize * sizeof(double));

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


    Boys_Finalize();

    printf("%11s  %11s  %11s  %11s\n", "res1_s", "res1_m", "res2_s", "res2_m"); 
    printf("%11.4e  %11.4e  %11.4e  %11.4e\n", res1_s[0], res1_m[0], res2_s[0], res2_m[0]);

    for(int i = 0; i < nshell1; i++)
        free_random_shell(A[i]);

    for(int j = 0; j < nshell2; j++)
        free_random_shell(B[j]);

    for(int k = 0; k < nshell3; k++)
        free_random_shell(C[k]);

    for(int l = 0; l < nshell4; l++)
        free_random_shell(D[l]);

    FREE(res1_s); FREE(res1_m); FREE(res2_s); FREE(res2_m); FREE(intwork);

    return 0;
}
