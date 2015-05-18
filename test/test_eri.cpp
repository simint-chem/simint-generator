#include <cstdio>
#include <cmath>

#include "eri/eri.h"
#include "boys/boys.h"
#include "common.hpp"

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

#define MAXAM 1

typedef int (*erifunc)(struct multishell_pair const, struct multishell_pair const, double * const restrict);


int eri_notyetimplemented(struct multishell_pair const P,
                          struct multishell_pair const Q,
                          double * const restrict dummy)
{
    printf("****************************\n");
    printf("*** NOT YET IMPLEMENTED! ***\n");
    printf("***  ( %2d %2d | %2d %2d )   ***\n", P.am1, P.am2, Q.am1, Q.am2);
    printf("****************************\n");
    exit(1);
    return 0;
}




int main(int argc, char ** argv)
{
    
    erifunc funcs[MAXAM+1][MAXAM+1][MAXAM+1][MAXAM+1];
    for(int i = 0; i <= MAXAM; i++)
    for(int j = 0; j <= MAXAM; j++)
    for(int k = 0; k <= MAXAM; k++)
    for(int l = 0; l <= MAXAM; l++)
        funcs[i][j][k][l] = eri_notyetimplemented; 

    funcs[0][0][0][0] = eri_FOcombined_ssss;
    funcs[1][0][0][0] = eri_FOcombined_psss;
    funcs[1][1][0][0] = eri_FOcombined_ppss;
    funcs[1][0][1][0] = eri_FOcombined_psps;
    funcs[1][1][1][0] = eri_FOcombined_ppps;
    funcs[1][1][1][1] = eri_FOcombined_pppp;


    if(argc != 13)
    {
        printf("Give me 12 arguments! I got %d\n", argc-1);
        return 1;
    }

    srand(time(NULL));

    std::array<int, 4> nshell = {atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4])};
    std::array<int, 4> nprim = {atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8])};
    std::array<int, 4> am = {atoi(argv[9]), atoi(argv[10]), atoi(argv[11]), atoi(argv[12])};


    int nshell1234 = nshell[0] * nshell[1] * nshell[2] * nshell[3];

    int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);   

    /* Storage of test results */
    double * res_FOcombined      = (double *)ALLOC(ncart * nshell1234 * sizeof(double));
    double * res_liberd          = (double *)ALLOC(ncart * nshell1234 * sizeof(double));
    double * res_valeev          = (double *)ALLOC(ncart * nshell1234 * sizeof(double));

    // allocate gaussian shell memory
    VecQuartet gshells(  CreateRandomQuartets(nshell, nprim, am) );

    Boys_Init(0, 7); // need F0 + 7 for interpolation
    Valeev_Init();

    // Actually calculate
    struct multishell_pair P = create_multishell_pair(nshell[0], gshells[0].data(),
                                                      nshell[1], gshells[1].data());
    struct multishell_pair Q = create_multishell_pair(nshell[2], gshells[2].data(),
                                                      nshell[3], gshells[3].data());

    funcs[am[0]][am[1]][am[2]][am[3]](P, Q, res_FOcombined);

    // test with valeev & liberd
    ValeevIntegrals(gshells, res_valeev);
    ERDIntegrals(gshells, res_liberd);


    printf("( %d %d | %d %d )\n", am[0], am[1], am[2], am[3]);
    printf("%22s %22s %22s\n", "liberd", "FOcombined", "valeev");

    int idx = 0;
    for(int i = 0; i < nshell1234; i++)
    {
        for(int j = 0; j < ncart; j++)
        {
            printf("%22.4e  %22.4e  %22.4e\n", res_liberd[idx], res_FOcombined[idx], res_valeev[idx]);

            const double v = res_valeev[idx];
            double diff_liberd         = fabs(res_liberd[idx]         - v);
            double diff_FOcombined     = fabs(res_FOcombined[idx]     - v);
            printf("%22.4e  %22.4e\n", diff_liberd, diff_FOcombined);
            printf("\n");

            idx++;
        }
        printf("\n");
    }

    free_multishell_pair(P);
    free_multishell_pair(Q);

    FreeRandomQuartets(gshells);

    Valeev_Finalize();
    Boys_Finalize();

    FREE(res_valeev);
    FREE(res_liberd);
    FREE(res_FOcombined);

    return 0;
}
