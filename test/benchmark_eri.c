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
    G.x = 4.0 * rand() / ((double)RAND_MAX) - 2.0;
    G.y = 4.0 * rand() / ((double)RAND_MAX) - 2.0;
    G.z = 4.0 * rand() / ((double)RAND_MAX) - 2.0;

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


void test_eri_2pair_flat(int na, struct gaussian_shell const * const restrict A,
                         int nb, struct gaussian_shell const * const restrict B,
                         int nc, struct gaussian_shell const * const restrict C,
                         int nd, struct gaussian_shell const * const restrict D,
                         double * const restrict res,
                         double * const restrict work1,
                         double * const restrict work2)
{
    struct shell_pair P = create_ss_shell_pair(na, A, nb, B);
    struct shell_pair Q = create_ss_shell_pair(nc, C, nd, D);
    eri_ssss_flat(P, Q, res, work1, work2);
    free_shell_pair(P);
    free_shell_pair(Q);
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

    double * res2_f = ALLOC(nshell1234 * sizeof(double));

    int worksize = (nshell1*nprim1 * nshell2*nprim2 * nshell3*nprim3 * nshell4*nprim4);
    double * intwork1 = ALLOC(worksize * sizeof(double));
    double * intwork2 = ALLOC(worksize * sizeof(double));

    srand(time(NULL));

    struct gaussian_shell * A = ALLOC(nshell1 * sizeof(struct gaussian_shell));
    struct gaussian_shell * B = ALLOC(nshell2 * sizeof(struct gaussian_shell));
    struct gaussian_shell * C = ALLOC(nshell3 * sizeof(struct gaussian_shell));
    struct gaussian_shell * D = ALLOC(nshell4 * sizeof(struct gaussian_shell));
    
    Boys_Init();
    printf("%11s\n", "res2_f"); 

    for(int i = 0; i < ntest; i++)
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
        test_eri_2pair_flat(nshell1, A, nshell2, B, nshell3, C, nshell4, D, res2_f, intwork1, intwork2);

        printf("%11e\n", res2_f[0]);

        for(int j = 0; j < nshell1; j++)
            free_random_shell(A[j]);

        for(int j = 0; j < nshell2; j++)
            free_random_shell(B[j]);

        for(int j = 0; j < nshell3; j++)
            free_random_shell(C[j]);

        for(int j = 0; j < nshell4; j++)
            free_random_shell(D[j]);

    }


    Boys_Finalize();


    FREE(res2_f);
    FREE(A); FREE(B); FREE(C); FREE(D);
    FREE(intwork1); FREE(intwork2);
    return 0;
}
