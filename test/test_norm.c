#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "eri/shell.h"
#include "vectorization.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_COORD 0.5
#define MAX_EXP   100.0
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
        //G.coef[i] = 1.0;
    }

    return G;
}

void free_random_shell(struct gaussian_shell G)
{
    free_gaussian_shell(G);
}


double compute_overlap(int am,
                       double x1, double y1, double z1, double alpha1,
                       double x2, double y2, double z2, double alpha2)
{
    // only does s for now
    const double X_ab = x1 - x2;    
    const double Y_ab = y1 - y2;    
    const double Z_ab = z1 - z2;    
    const double p = alpha1 + alpha2;
    const double mu = (alpha1 * alpha2) / p;

    const double fac = sqrt(M_PI/p);

    const double Sx = fac * exp(-mu * X_ab);
    const double Sy = fac * exp(-mu * Y_ab);
    const double Sz = fac * exp(-mu * Z_ab);

    return Sx * Sy * Sz;
}
                

void print_overlaps(int nshell, struct gaussian_shell const * const restrict A)
{
    for(int i = 0; i < nshell; i++)
    {
        printf("Shell %d\n", i);
        for(int j = 0; j < A[i].nprim; j++)
        {
            double S = compute_overlap(A[i].am,
                                       A[i].x, A[i].y, A[i].z, A[i].alpha[j],
                                       A[i].x, A[i].y, A[i].z, A[i].alpha[j]);
            S *= A[i].coef[j] * A[i].coef[j];
            printf("Prim %d = %12.8f\n", i, S);
        }

        // do the whole shell
        double S = 0;
        for(int a = 0; a < A[i].nprim; a++)
        for(int b = 0; b < A[i].nprim; b++)
        {
            S += A[i].coef[a] * A[i].coef[b]
               * compute_overlap(A[i].am,
                                 A[i].x, A[i].y, A[i].z, A[i].alpha[a],
                                 A[i].x, A[i].y, A[i].z, A[i].alpha[b]);
        }
        printf("Shell overlap: %12.8f\n", S);
    }

}


int main(int argc, char ** argv)
{
    if(argc != 3)
    {
        printf("Give me 2 argument! I got %d\n", argc-1);
        return 1;
    }

    int nshell = atoi(argv[1]);
    int nprim = atoi(argv[2]);

    //int nshell1234 = nshell1 * nshell2 * nshell3 * nshell4;

    srand(time(NULL));

    struct gaussian_shell * A = ALLOC(nshell * sizeof(struct gaussian_shell));

    for(int i = 0; i < nshell; i++)
        A[i] = random_shell(nprim);


    printf("== Before normalization ==\n");
    print_overlaps(nshell, A);

    printf("\n");
    printf("== After normalization ==\n");
    normalize_gaussian_shells(nshell, A);
    print_overlaps(nshell, A);

    for(int i = 0; i < nshell; i++)
        free_random_shell(A[i]);

    FREE(A);

    return 0;
}
