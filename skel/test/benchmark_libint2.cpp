#include <cstdio>
#include <memory>

#include <omp.h>

#include "vectorization/vectorization.h"
#include "eri/eri.h"
#include "boys/boys.h"
#include "test/common.hpp"
#include "test/Libint2.hpp"
#include "test/timer.h"

#ifdef BENCHMARK_VALIDATE
#include "test/valeev.hpp"
#endif

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

    #ifdef BENCHMARK_VALIDATE
    printf("Disabling threading since BENCHMARK_VALIDATE is set\n");
    nthread = 1;
    #endif

    if(argc == 4)
        nthread = atoi(argv[3]);

    // read in the shell info
    ShellMap shellmap = ReadBasis(basfile);

    // normalize
    for(auto & it : shellmap)
        normalize_gaussian_shells(it.second.size(), it.second.data());

    // find the max dimensions
    std::array<int, 3> maxparams = FindMapMaxParams(shellmap);
    const int maxam = (maxparams[0] > MAXAM ? MAXAM : maxparams[0]);
    const int maxnprim = maxparams[1];
    const int maxsize = maxparams[2];

    /* Storage of integrals */
    std::vector<double *> res_ints(nthread);
    for(int i = 0; i < nthread; i++)
        res_ints[i] = (double *)ALLOC(maxsize * sizeof(double));


    // initialize stuff
    LIBINT2_PREFIXED_NAME(libint2_static_init)();
    std::vector<std::unique_ptr<Libint2_ERI>> alleri(nthread);
    for(auto & it : alleri)
        it = std::unique_ptr<Libint2_ERI>(new Libint2_ERI(maxam, maxnprim));


    #ifdef BENCHMARK_VALIDATE
    Valeev_Init();
    double * res_ref = (double *)ALLOC(maxsize * sizeof(double));
    #endif

    // Timing header
    printf("%5s %13s %12s   %12s   %16s   %12s\n", "ID",
                           "Quartet", "NCont", "NPrim", "Ticks", "Ticks/Prim");

    // loop ntest times
    #pragma omp parallel for num_threads(nthread)
    for(int n = 0; n < ntest; n++)
    {
        int ithread = omp_get_thread_num();
        Libint2_ERI * eri = alleri[ithread].get();

        #ifdef BENCHMARK_VALIDATE
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

            const int nshell3 = shellmap[k].size();
            const int nshell4 = shellmap[l].size();

            gaussian_shell const * const C = &shellmap[k][0];
            gaussian_shell const * const D = &shellmap[l][0];
            struct multishell_pair Q = create_multishell_pair(nshell3, C, nshell4, D);

            TimerType time_total = 0;
            size_t nprim_total = 0;
            size_t nshell1234_total = 0;


            // do one shell pair at a time on the bra side
            for(size_t a = 0; a < shellmap[i].size(); a++)
            for(size_t b = 0; b < shellmap[j].size(); b++)
            {
                const int nshell1 = 1; 
                const int nshell2 = 1; 
                const size_t nshell1234 = nshell1 * nshell2 * nshell3 * nshell4;

                gaussian_shell const * const A = &shellmap[i][a];
                gaussian_shell const * const B = &shellmap[j][b];

                struct multishell_pair P = create_multishell_pair(nshell1, A, nshell2, B);


                // actually calculate
                time_total += eri->Integrals(P, Q, res_ints[ithread]);


                #ifdef BENCHMARK_VALIDATE
                const int ncart1234 = NCART(i) * NCART(j) * NCART(j) * NCART(l);
                const int arrlen = nshell1234 * ncart1234;
                ValeevIntegrals(A, nshell1,
                                B, nshell2,
                                C, nshell3,
                                D, nshell4,
                                res_ref, false);
                std::pair<double, double> err = CalcError(res_ints[ithread], res_ref, arrlen);
                printf("( %d %d | %d %d ) MaxAbsErr: %10.3e   MaxRelErr: %10.3e\n", i, j, k, l, err.first, err.second);
                #endif
                    
                nshell1234_total += nshell1234;

                // calculate the number of primitives
                for(int p = 0; p < nshell1; p++)
                for(int q = 0; q < nshell2; q++)
                for(int r = 0; r < nshell3; r++)
                for(int s = 0; s < nshell4; s++)
                    nprim_total += A[p].nprim * B[q].nprim * C[r].nprim * D[s].nprim;

                free_multishell_pair(P);
            }

            free_multishell_pair(Q);

            printf("[%3d] ( %d %d | %d %d ) %12lu   %12lu   %16llu   %12.3f\n",
                                                                          ithread,
                                                                          i, j, k, l,
                                                                          nshell1234_total, nprim_total,
                                                                          time_total,
                                                                          (double)(time_total)/(double)(nprim_total));
        }
    }

    FreeShellMap(shellmap);
    for(int i = 0; i < nthread; i++)
        FREE(res_ints[i]);


    #ifdef BENCHMARK_VALIDATE
    FREE(res_ref);
    #endif
   
    return 0;
}
