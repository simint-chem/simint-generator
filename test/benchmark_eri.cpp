#include <cstdio>

#include <omp.h>

#include "eri/eri.h"
#include "boys/boys.h"
#include "test/common.hpp"


#define CLOCK(ticks, time) do {                                 \
    volatile unsigned int a, d;                              \
    __asm__ __volatile__("rdtsc" : "=a" (a), "=d" (d) : );   \
    (ticks) = ((unsigned long long) a)|(((unsigned long long)d)<<32); \
    (time) = (ticks) / 3700000000.;                              \
  } while(0)


int main(int argc, char ** argv)
{
    // set up the function pointers
    Init_Test();

    if(argc < 3 || argc > 4)
    {
        printf("Give me 2 or 3 arguments! I got %d\n", argc-1);
        return 1;
    }

    // files to read
    std::string basfile(argv[1]);

    // number of times to run
    const int ntest = atoi(argv[2]);

    // number of threads
    int nthread = 1;

    if(argc == 4)
        nthread = atoi(argv[3]);

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
    std::vector<double *> res_ints(nthread);
    for(int i = 0; i < nthread; i++)
        res_ints[i] = (double *)ALLOC(maxsize * sizeof(double));

    // initialize stuff
    // nothing needs initializing!


    #ifdef BENCHMARK_VALIDATE
    double * res_ref = (double *)ALLOC(maxsize * sizeof(double));
    RefIntegralReader refint(basfile);
    #endif

    // Timing header
    printf("%5s %13s %12s   %12s   %16s  %15s    %12s\n", "ID",
                           "Quartet", "NCont", "NPrim", "Ticks", "Clock", "Ticks/Prim");

    // loop ntest times omp threads
    #pragma omp parallel for num_threads(nthread)
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


            int ithread = omp_get_thread_num();


            const AlignedGaussianVec & it_i = shellmap[i];
            const AlignedGaussianVec & it_j = shellmap[j];
            const AlignedGaussianVec & it_k = shellmap[k];
            const AlignedGaussianVec & it_l = shellmap[l];


            // set up shell pairs
            struct multishell_pair P = create_multishell_pair(it_i.size(), it_i.data(),
                                                              it_j.size(), it_j.data());
            struct multishell_pair Q = create_multishell_pair(it_k.size(), it_k.data(),
                                                              it_l.size(), it_l.data());
            // for timing
            double wallclock0, wallclock1;
            unsigned long long ticks0, ticks1;

            // actually calculate
            CLOCK(ticks0, wallclock0); 
            Integral(P, Q, res_ints[ithread]);
            CLOCK(ticks1, wallclock1); 

            unsigned long ncont = (unsigned long)(P.nshell12) * (unsigned long)(Q.nshell12);
            unsigned long nprim = (unsigned long)(P.nprim) * (unsigned long)(Q.nprim);

            printf("[%3d] ( %d %d | %d %d ) %12lu   %12lu   %16lu  (%8.3f secs)    %12.3f\n",
                                                                          ithread,
                                                                          i, j, k, l,
                                                                          ncont, nprim,
                                                                          ticks1 - ticks0,
                                                                          wallclock1 - wallclock0,
                                                                          (double)(ticks1-ticks0)/(double)(nprim));


            #ifdef BENCHMARK_VALIDATE
            const int ncart1234 = NCART(i) * NCART(j) * NCART(k) * NCART(l);
            const int nshell1 = shellmap[i].size();
            const int nshell2 = shellmap[j].size();
            const int nshell3 = shellmap[k].size();
            const int nshell4 = shellmap[l].size();
            const int nshell1234 = nshell1 * nshell2 * nshell3 * nshell4;
            const int arrlen = nshell1234 * ncart1234;
            refint.ReadNext(res_ref, arrlen);
            Chop(res_ints[i], arrlen);
            Chop(res_ref, arrlen);
            std::pair<double, double> err = CalcError(res_ints[i], res_ref, arrlen);

            printf("( %d %d | %d %d ) MaxAbsErr: %10.3e   MaxRelErr: %10.3e\n", i, j, k, l, err.first, err.second);
            #endif


            free_multishell_pair(P);
            free_multishell_pair(Q);
        }
    }

    FreeShellMap(shellmap);
    for(int i = 0; i < nthread; i++)
        FREE(res_ints[i]);


    #ifdef BENCHMARK_VALIDATE
    FREE(res_ref);
    #endif
    
    // Finalize stuff 
    // Nothing to finalize!

    return 0;
}
