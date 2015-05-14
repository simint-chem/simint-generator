#include <cstdio>
#include <cmath>

#include "eri/eri.h"
#include "boys/boys.h"
#include "common.hpp"



int main(int argc, char ** argv)
{
    if(argc != 9)
    {
        printf("Give me 8 arguments! I got %d\n", argc-1);
        return 1;
    }

    srand(time(NULL));

    std::array<int, 4> nshell = {atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4])};
    std::array<int, 4> nprim = {atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8])};


    int nshell1234 = nshell[0] * nshell[1] * nshell[2] * nshell[3];

    /* Storage of test results */
    double * res_FOcombined      = (double *)ALLOC(81 * nshell1234 * sizeof(double));
    double * res_liberd          = (double *)ALLOC(81 * nshell1234 * sizeof(double));
    double * res_valeev          = (double *)ALLOC(81 * nshell1234 * sizeof(double));

    // allocate gaussian shell memory
    VecQuartet gshells(  CreateRandomQuartets(nshell, nprim, {0,0,0,0}) );

    Boys_Init(0, 7); // need F0 + 7 for interpolation
    Valeev_Init();

    // Actually calculate
    // ( s s | s s )
    struct multishell_pair P = create_multishell_pair(nshell[0], gshells[0].data(),
                                                      nshell[1], gshells[1].data());
    struct multishell_pair Q = create_multishell_pair(nshell[2], gshells[2].data(),
                                                      nshell[3], gshells[3].data());

    eri_FOcombined_ssss(P, Q, res_FOcombined);

    // test with valeev & liberd
    ValeevIntegrals(gshells, res_valeev);
    ERDIntegrals(gshells, res_liberd);


    printf("( s s | s s ) Integrals\n");
    printf("%22s %22s %22s\n",
           "liberd", "FOcombined", "valeev");

    for(int i = 0; i < nshell1234; i++)
    {
        const double v = res_valeev[i];

        double diff_liberd         = fabs(res_liberd[i]         - v);
        double diff_FOcombined     = fabs(res_FOcombined[i]     - v);

        printf("%22.4e  %22.4e  %22.4e\n",
                                  res_liberd[i],
                                  res_FOcombined[i],
                                  res_valeev[i]);

        printf("%22.4e  %22.4e\n",
                                  diff_liberd,
                                  diff_FOcombined);

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
