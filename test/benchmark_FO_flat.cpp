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
    ShellMap shellmap = ReadBasis(basfile);

    // normalize
    for(auto & it : shellmap)
        normalize_gaussian_shells(it.second.size(), it.second.data());

    // find the max dimensions
    std::array<int, 3> maxparams = FindMapMaxParams(shellmap);
    const int maxam = maxparams[0];
    //const int maxnprim = maxparams[1];
    const int maxsize = maxparams[2];

    /* Storage of integrals */
    double * res_ints = (double *)ALLOC(maxsize * sizeof(double));

    // initialize stuff
    // nothing needs initializing!


    #ifdef BENCHMARK_VALIDATE
    double * res_ref = (double *)ALLOC(maxsize * sizeof(double));
    RefIntegralReader refint(basfile);
    #endif


    // loop ntest times
    for(int n = 0; n < ntest; n++)
    {
        #ifdef BENCHMARK_VALIDATE
        // move the reader back to the beginning of the file
        refint.Reset();
        printf("\n");
        printf("Run %d\n", n);
        #endif
        for(int i = 0; i <= maxam; i++)
        for(int j = 0; j <= maxam; j++)
        for(int k = 0; k <= maxam; k++)
        for(int l = 0; l <= maxam; l++)
        {
            if(!ValidQuartet(i, j, k, l))
                continue;

            const AlignedGaussianVec & it_i = shellmap[i];
            const AlignedGaussianVec & it_j = shellmap[j];
            const AlignedGaussianVec & it_k = shellmap[k];
            const AlignedGaussianVec & it_l = shellmap[l];


            // set up shell pairs
            struct multishell_pair_flat P = create_multishell_pair_flat(it_i.size(), it_i.data(),
                                                                        it_j.size(), it_j.data());
            struct multishell_pair_flat Q = create_multishell_pair_flat(it_k.size(), it_k.data(),
                                                                        it_l.size(), it_l.data());
            // actually calculate
            Integral_FO_flat(P, Q, res_ints);


            #ifdef BENCHMARK_VALIDATE
            const int ncart1234 = NCART(i) * NCART(j) * NCART(k) * NCART(l);
            const int nshell1 = shellmap[i].size();
            const int nshell2 = shellmap[j].size();
            const int nshell3 = shellmap[k].size();
            const int nshell4 = shellmap[l].size();
            const int nshell1234 = nshell1 * nshell2 * nshell3 * nshell4;
            const int arrlen = nshell1234 * ncart1234;
            refint.ReadNext(res_ref, arrlen);
            Chop(res_ints, arrlen);
            Chop(res_ref, arrlen);
            std::pair<double, double> err = CalcError(res_ints, res_ref, arrlen);
            printf("( %d %d | %d %d ) MaxAbsErr: %10.3e   MaxRelErr: %10.3e\n", i, j, k, l, err.first, err.second);
            #endif


            free_multishell_pair_flat(P);
            free_multishell_pair_flat(Q);
        }
    }

    FreeShellMap(shellmap);
    FREE(res_ints);


    #ifdef BENCHMARK_VALIDATE
    FREE(res_ref);
    #endif
    
    // Finalize stuff 
    // Nothing to finalize!

    return 0;
}
