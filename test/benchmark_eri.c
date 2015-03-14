#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "valeev/valeev.h"
#include "eri/eri.h"
#include "boys/boys.h"

struct gaussian_shell
random_shell(int nprim)
{
    struct gaussian_shell G = allocate_gaussian_shell(nprim);
    G.am = 0;
    //G.x = 4.0 * rand() / ((double)RAND_MAX) - 2.0;
    //G.y = 4.0 * rand() / ((double)RAND_MAX) - 2.0;
    //G.z = 4.0 * rand() / ((double)RAND_MAX) - 2.0;
    G.x = rand() / ((double)RAND_MAX) - 0.5;
    G.y = rand() / ((double)RAND_MAX) - 0.5;
    G.z = rand() / ((double)RAND_MAX) - 0.5;

    G.nprim = nprim;

    for(int i = 0; i < nprim; i++)
    {
        G.alpha[i] = 10.0 * rand() / ((double)RAND_MAX);
        G.coef[i] = 10.0 * rand() / ((double)RAND_MAX);
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

    const int ntest = 50;

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

    int worksize = SIMD_ROUND_DBL(nshell1*nprim1 * nshell2*nprim2 * nshell3*nprim3 * nshell4*nprim4);
    double * intwork1 = ALLOC(3*worksize * sizeof(double));
    double * intwork2 = ALLOC(3*worksize * sizeof(double));

    srand(time(NULL));

    struct gaussian_shell * A = ALLOC(nshell1 * sizeof(struct gaussian_shell));
    struct gaussian_shell * B = ALLOC(nshell2 * sizeof(struct gaussian_shell));
    struct gaussian_shell * C = ALLOC(nshell3 * sizeof(struct gaussian_shell));
    struct gaussian_shell * D = ALLOC(nshell4 * sizeof(struct gaussian_shell));
    
    Boys_Init();
    printf("%11s %11s %11s %11s\n", "res_f", "res_fs", "res_ft", "res_fc"); 

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

        // Actually calculate
        struct shell_pair P = create_ss_shell_pair(nshell1, A, nshell2, B);
        struct shell_pair Q = create_ss_shell_pair(nshell3, C, nshell4, D);
        eri_ssss_flat_combined(P, Q, res_fc, intwork1, intwork2);
        eri_ssss_flat_taylor(P, Q, res_ft, intwork1, intwork2);
        eri_ssss_flat_split(P, Q, res_fs, intwork1, intwork2);
        eri_ssss_flat(P, Q, res_f, intwork1, intwork2);
        free_shell_pair(P);
        free_shell_pair(Q);

        printf("%11e %11e %11e %11e\n", res_f[0], res_fs[0], res_ft[0], res_fc[0]);

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


    FREE(res_f); FREE(res_fs);
    FREE(res_ft); FREE(res_fc);
    FREE(A); FREE(B); FREE(C); FREE(D);
    FREE(intwork1); FREE(intwork2);

    return 0;
}
