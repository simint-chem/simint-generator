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

    // loop ntest times
    for(int n = 0; n < ntest; n++)
    {
        for(int i = 0; i <= maxam; i++)
        for(int j = 0; j <= maxam; j++)
        for(int k = 0; k <= maxam; k++)
        for(int l = 0; l <= maxam; l++)
        {
            std::array<int, 4> am{i, j, k, l};

            if(!ValidQuartet(am))
                continue;

            const AlignedGaussianVec & it_i = shellmap[i];
            const AlignedGaussianVec & it_j = shellmap[j];
            const AlignedGaussianVec & it_k = shellmap[k];
            const AlignedGaussianVec & it_l = shellmap[l];


            // set up shell pairs
            struct multishell_pair P = create_multishell_pair(it_i.size(), it_i.data(),
                                                              it_j.size(), it_j.data());
            struct multishell_pair Q = create_multishell_pair(it_k.size(), it_k.data(),
                                                              it_l.size(), it_l.data());
            // actually calculate
            Integral_FO(P, Q, res_ints);

            free_multishell_pair(P);
            free_multishell_pair(Q);
        }
    }

    FreeShellMap(shellmap);
    FREE(res_ints);

    return 0;
}
