#include <cstdio>
#include <cmath>

#include "eri/eri.h"
#include "boys/boys.h"
#include "common.hpp"

#define NTEST 500

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
    double * res_FO              = (double *)ALLOC(nshell1234 * sizeof(double));
    double * res_liberd          = (double *)ALLOC(nshell1234 * sizeof(double));

    printf("%22s %22s\n", "liberd", "FO");

    for(int n = 0; n < NTEST; n++)
    {
        auto gshells(  CreateRandomQuartets(nshell, nprim, {0,0,0,0}) );

        // Actually calculate
        struct multishell_pair P = create_multishell_pair(nshell[0], gshells[0].data(),
                                                             nshell[1], gshells[1].data());
        struct multishell_pair Q = create_multishell_pair(nshell[2], gshells[2].data(),
                                                             nshell[3], gshells[3].data());

        eri_FO_s_s_s_s(           P, Q, res_FO                          );
        free_multishell_pair(P);
        free_multishell_pair(Q);

        // test with liberd
        ERDIntegrals(gshells, res_liberd);

        // print some results
        printf("%22e %22e\n", res_liberd[0], res_FO[0]);

        // free memory
        FreeQuartets(gshells);

    }


    Boys_Finalize();

    FREE(res_liberd);
    FREE(res_FO);

    return 0;
}
