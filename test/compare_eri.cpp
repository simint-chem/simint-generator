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
    double * res_split           = (double *)ALLOC(nshell1234 * sizeof(double));
    double * res_splitcombined   = (double *)ALLOC(nshell1234 * sizeof(double));
    double * res_taylorcombined  = (double *)ALLOC(nshell1234 * sizeof(double));
    double * res_FOcombined      = (double *)ALLOC(nshell1234 * sizeof(double));

    double * res_liberd          = (double *)ALLOC(nshell1234 * sizeof(double));
    double * res_valeev          = (double *)ALLOC(nshell1234 * sizeof(double));

    // for split
    // no need to round each pair
    int worksize = SIMD_ROUND_DBL(nshell[0]*nprim[0] * nshell[1]*nprim[1] * nshell[2]*nprim[2] * nshell[3]*nprim[3]);

    double * intwork1 = (double *)ALLOC(3 * worksize * sizeof(double));
    double * intwork2 = (double *)ALLOC(3 * worksize * sizeof(double));

    // allocate gaussian shell memory
    VecQuartet gshells(  CreateRandomQuartets(nshell, nprim) );

    // find the maximum possible x value for the boys function
    const double maxR2 = 12.0 * MAX_COORD * MAX_COORD;
    const double max_x = maxR2 * (MAX_EXP*MAX_EXP) / (2.0 * MAX_EXP);
    printf("Maximum parameter to boys: %12.8e\n", max_x);
    Boys_Init(max_x, 7); // need F0 + 7 for interpolation
    Valeev_Init();

    // Actually calculate
    struct multishell_pair P = create_ss_multishell_pair(nshell[0], gshells[0].data(),
                                                         nshell[1], gshells[1].data());
    struct multishell_pair Q = create_ss_multishell_pair(nshell[2], gshells[2].data(),
                                                         nshell[3], gshells[3].data());

    eri_taylorcombined_ssss(  P, Q, res_taylorcombined                      );
    eri_split_ssss(           P, Q, res_split,          intwork1, intwork2  );
    eri_splitcombined_ssss(   P, Q, res_splitcombined                       );
    eri_FOcombined_ssss(      P, Q, res_FOcombined                          );
    free_multishell_pair(P);
    free_multishell_pair(Q);


    // test with valeev & liberd
    ValeevIntegrals(gshells, res_valeev);
    ERDIntegrals(gshells, res_liberd);



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

    FreeRandomQuartets(gshells);

    Valeev_Finalize();
    Boys_Finalize();

    FREE(res_valeev);
    FREE(res_liberd);
    FREE(res_split);
    FREE(res_splitcombined);
    FREE(res_taylorcombined);
    FREE(res_FOcombined);
    FREE(intwork1); FREE(intwork2);

    return 0;
}