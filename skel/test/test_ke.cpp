#include <cstdio>
#include <atomic>
#include <iostream>
#include <cmath>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "simint/simint.h"
#include "test/Common.hpp"
#include "test/ValeevRef.hpp"


#define SIMINT_SCREEN 0
#define SIMINT_SCREEN_TOL 0.0


int main(int argc, char ** argv)
{
    // set up the function pointers
    simint_init();

    // parse command line
    if(argc != 2)
    {
        printf("Give me 1 argument! I got %d\n", argc-1);
        return 1;
    }

    // basis functions file to read
    std::string basfile(argv[1]);

    // read in the shell info
    ShellMap shellmap = ReadBasis(basfile).first;

    // normalize the original
    for(auto & it : shellmap)
    {
        simint_normalize_shells(it.second.size(), it.second.data());
        for(const auto & it2 : it.second)
        {
            std::cout << it2.am << "\n";
            for(int i = 0; i < it2.nprim; i++)
                printf("%24.16e  %24.16e\n", it2.alpha[i], it2.coef[i]);
        }
        std::cout << "\n";
    }

    // find the max dimensions
    std::pair<int, int> maxparams = FindMaxParams(shellmap);
    const int maxam = maxparams.first;
    const int max_ncart = ( (maxam+1)*(maxam+2) )/2;
    const int maxsize = maxparams.second * max_ncart;

    // get the number of threads
    #ifdef _OPENMP
        const int nthread = omp_get_max_threads();
    #else
        const int nthread = 1;
    #endif

    /* Storage of integrals */
    double * all_res_simint = (double *)SIMINT_ALLOC(nthread * maxsize * sizeof(double));


    // loop over all AM quartets, choosing only valid ones
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    {
        const int ncart12 = NCART(i) * NCART(j);

        for(size_t a = 0; a < shellmap[i].size(); a++)
        for(size_t b = 0; b < shellmap[j].size(); b++)
        {
            #ifdef _OPENMP
                const int ithread = omp_get_thread_num();
            #else
                const int ithread = 0;
            #endif

            double * res_simint = all_res_simint + ithread * maxsize;

            ////////////////////////////
            // Calculate the integrals
            ////////////////////////////
            int simint_ret = simint_compute_ke(&shellmap[i][a], &shellmap[j][b], res_simint);

            // if the return is < 0, it didn't calculate anything
            // (everything was screened)
            if(simint_ret < 0)
                std::fill(res_simint, res_simint + ncart12, 0.0);

            #pragma omp critical
            {
                printf("AM: %d %d\n", i, j);
                printf("SHELL: %lu %lu\n", a, b);
                for(int n = 0; n < ncart12; n++) 
                    printf("  %4d   ->   %30.16e\n", n, res_simint[n]);
                printf("\n");
            }

        } // end threaded loop over a,b

    }

    simint_finalize();

    SIMINT_FREE(all_res_simint);

    return 0;
}
