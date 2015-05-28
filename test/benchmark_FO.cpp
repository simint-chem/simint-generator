#include <cstdio>

#include "eri/eri.h"
#include "boys/boys.h"
#include "test/common.hpp"

int main(int argc, char ** argv)
{
    // set up the function pointers
    Init_Test();

    if(argc != 7)
    {
        printf("Give me 6 arguments! I got %d\n", argc-1);
        return 1;
    }

    // files to read
    std::string basedir(argv[1]);

    // number of times to run
    const int ntest = atoi(argv[2]);

    std::array<int, 4> am = { atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]) };
    VecQuartet gshells(ReadQuartets(am, basedir, true));

    std::array<int, 4> nshell{gshells[0].size(), gshells[1].size(), 
                              gshells[2].size(), gshells[3].size()};


    const int nshell1234 = nshell[0] * nshell[1] * nshell[2] * nshell[3];
    const int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]); 

    /* Storage of integrals */
    double * res_ints = (double *)ALLOC(ncart * nshell1234 * sizeof(double));

    // set up the shell pairs here
    struct multishell_pair P = create_multishell_pair(nshell[0], gshells[0].data(),
                                                      nshell[1], gshells[1].data());
    struct multishell_pair Q = create_multishell_pair(nshell[2], gshells[2].data(),
                                                      nshell[3], gshells[3].data());

    Boys_Init();


    // Calculate the integral NTEST times 
    for(int n = 0; n < ntest; n++)
        Integral_FO(P, Q, res_ints);


    // compare at the end
    double * res_ref = (double *)ALLOC(ncart * nshell1234 * sizeof(double));
    ReadValeevIntegrals(basedir, am, res_ref);
    std::pair<double, double> err = CalcError(res_ints, res_ref, nshell1234 * ncart);
    free(res_ref);

    printf("\n");
    printf("Max abs error: %22.8e\n", err.first);
    printf("Max rel error: %22.8e\n", err.second);
    printf("\n");

    FreeQuartets(gshells);

    Boys_Finalize();

    free_multishell_pair(P);
    free_multishell_pair(Q);

    FREE(res_ints);

    return 0;
}
