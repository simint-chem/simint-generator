#include <cstdio>

#include "eri/eri.h"
#include "boys/boys.h"
#include "test/common.hpp"

int main(int argc, char ** argv)
{
    // set up the function pointers
    Init_Test();

    if(argc != 3)
    {
        printf("Give me 2 argument! I got %d\n", argc-1);
        return 1;
    }

    // files to read
    std::string basfile(argv[1]);

    // number of times to run
    const int ntest = atoi(argv[2]);

    // read in the shell info
    std::map<int, AlignedGaussianVec> shellmap = ReadBasis(basfile);

    // normalize
    for(const auto & it : shellmap)
        normalize_gaussian_shells(it.second.size(), it.second.data());

    // find the max dimensions
    int maxnprim = 0;
    int maxam = 0;
    int maxel = 0;
    for(auto & it : shellmap)
    {
        const int nca = NCART(it.first);
        const int nsh = it.second.size();
        const int n = nca * nsh;
        if(n > maxel)
            maxel = n;

        if(it.first > maxam)
            maxam = it.first;

        for(auto & it2 : it.second)
        {
            if(it2.nprim > maxnprim)
                maxnprim = it2.nprim;
        }
    }

    maxel = pow(maxel, 4);

    /* Storage of integrals */
    double * res_ints = (double *)ALLOC(maxel * sizeof(double));


    for(const auto & it_i : shellmap)
    for(const auto & it_j : shellmap)
    for(const auto & it_k : shellmap)
    for(const auto & it_l : shellmap)
    {
        std::array<int, 4> am{it_i.first, it_j.first, it_k.first, it_l.first};

        if(!ValidQuartet(am))
            continue;

        
        const int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);   



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
