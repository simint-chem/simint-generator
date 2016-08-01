#include <cstdio>

#include "simint/vectorization/vectorization.h"
#include "simint/eri/eri.h"
#include "test/Simint.hpp"
#include "test/Common.hpp"
#include "test/Timer.h"

#include <omp.h>

#ifdef BENCHMARK_VALIDATE
#include "test/ValeevRef.hpp"
#endif


int main(int argc, char ** argv)
{
    // set up the function pointers
    Simint_Init();

    if(argc != 2)
    {
        printf("Give me 1 argument! I got %d\n", argc-1);
        return 1;
    }

    // files to read
    std::string basfile(argv[1]);

    // number of threads
    const int nthread = omp_get_max_threads();

    // read in the shell info
    ShellMap shellmap = ReadBasis(basfile);

    // normalize
    for(auto & it : shellmap)
        simint_normalize_shells(it.second.size(), it.second.data());


    // find the max dimensions
    std::array<int, 3> maxparams = FindMapMaxParams(shellmap);
    const int maxam = (maxparams[0] > SIMINT_ERI_MAXAM ? SIMINT_ERI_MAXAM : maxparams[0]);
    const int maxsize = maxparams[2];

    /* Storage of integrals */
    double * all_res_ints = (double *)ALLOC(nthread * maxsize * sizeof(double));

    /* contracted workspace */
    double * all_simint_work = (double *)ALLOC(nthread * SIMINT_ERI_MAX_WORKMEM);


    // initialize stuff
    // nothing needs initializing!

    #ifdef BENCHMARK_VALIDATE
    ValeevRef_Init();
    double * all_res_ref = (double *)ALLOC(nthread * maxsize * sizeof(double));
    #endif

    PrintTimingHeader();

    // running totals
    size_t ncont1234_total = 0;
    size_t nprim1234_total = 0;
    size_t nshell1234_total = 0;
    TimeContrib time_total;


    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    for(int k = 0; k <= maxam; k++)
    for(int l = 0; l <= maxam; l++)
    {
        if(!ValidQuartet({i,j,k,l}))
            continue;

        // total amount of time for this AM quartet
        TimeContrib time_am;

        const size_t nshell3 = shellmap[k].size();
        const size_t nshell4 = shellmap[l].size();

        simint_shell const * const C = &shellmap[k][0];
        simint_shell const * const D = &shellmap[l][0];

        // time creation of Q
        TimerType time_pair_34_0, time_pair_34_1 = 0;
        CLOCK(time_pair_34_0);
        struct simint_multi_shellpair Q = simint_create_multi_shellpair(nshell3, C, nshell4, D);
        CLOCK(time_pair_34_1);
        time_am.fill_shell_pair += time_pair_34_1 - time_pair_34_0;


        #ifdef BENCHMARK_VALIDATE
        std::pair<double, double> err{0.0, 0.0};
        #endif

        // running totals for this am
        size_t nprim1234_am = 0;
        size_t nshell1234_am = 0;
        size_t ncont1234_am = 0;

	    const auto & shellmap_i = shellmap[i];
	    const auto & shellmap_j = shellmap[j];

        const size_t i_size = shellmap_i.size();
        const size_t j_size = shellmap_j.size();
	    const size_t brasize = i_size * j_size;

        // do one shell pair at a time on the bra side
        #pragma omp parallel for schedule(dynamic)
        for(size_t ab = 0; ab < brasize; ab++)
        {
	        size_t a = ab / j_size;
    	    size_t b = ab % j_size;

            const size_t nshell1 = 1;
            const size_t nshell2 = 1;

            simint_shell const * const A = &shellmap_i[a];
            simint_shell const * const B = &shellmap_j[b];

            // time creation of P
            TimerType time_pair_12_0, time_pair_12_1;
            CLOCK(time_pair_12_0);
            struct simint_multi_shellpair P = simint_create_multi_shellpair(nshell1, A, nshell2, B);
            CLOCK(time_pair_12_1);

            int ithread = omp_get_thread_num();
            double * const res_ints = all_res_ints + ithread * maxsize;
            double * const simint_work = all_simint_work + ithread * SIMINT_ERI_MAX_WORKSIZE;

            // actually calculate
            TimerType my_time_am = Simint_Integral(P, Q, simint_work, res_ints);

            // acutal number of primitives and shells calculated
            // TODO - replace with return values from Integrals
            const size_t nprim1234 = P.nprim * Q.nprim;
            const size_t nshell1234 = P.nshell12 * Q.nshell12;
            const size_t ncart1234 = NCART(i) * NCART(j) * NCART(k) * NCART(l);
            const size_t ncont1234 = nshell1234 * ncart1234;

            #ifdef BENCHMARK_VALIDATE
            double * res_ref = all_res_ref + ithread * maxsize;

            ValeevRef_Integrals(A, nshell1,
                                B, nshell2,
                                C, nshell3,
                                D, nshell4,
                                res_ref, false);
            std::pair<double, double> err2 = CalcError(res_ints, res_ref, ncont1234);

            #pragma omp critical
            {
                err.first = std::max(err.first, err2.first);
                err.second = std::max(err.second, err2.second);
            }
            #endif

            // free this here, since we are done with it
            simint_free_multi_shellpair(P);


            #pragma omp critical
            {
                // add primitive and shell count to running totals for this am
                ncont1234_am += ncont1234;
                nprim1234_am += nprim1234;
                nshell1234_am += nshell1234;
                time_am.integrals += my_time_am;
                time_am.fill_shell_pair += time_pair_12_1 - time_pair_12_0;
            }
        }

        simint_free_multi_shellpair(Q);

		// add primitive and shell count to overall running totals
		ncont1234_total += ncont1234_am;
		nprim1234_total += nprim1234_am;
		nshell1234_total += nshell1234_am;
		time_total += time_am;

        PrintAMTimingInfo(i, j, k, l, nshell1234_am, nprim1234_am, time_am);


        #ifdef BENCHMARK_VALIDATE
        printf("( %d %d | %d %d ) MaxAbsErr: %10.3e   MaxRelErr: %10.3e\n", i, j, k, l, err.first, err.second);
        #endif
    }

    printf("\n");
    printf("Calculated %ld contracted integrals\n", static_cast<long>(ncont1234_total));
    printf("Calculated %ld contracted shells\n", static_cast<long>(nshell1234_total));
    printf("Calculated %ld primitive integrals\n", static_cast<long>(nprim1234_total));
    printf("Total ticks to calculate all:  %llu\n", time_total.TotalTime());
    printf("\n");

    FreeShellMap(shellmap);

    FREE(all_res_ints);
    FREE(all_simint_work);


    #ifdef BENCHMARK_VALIDATE
    FREE(all_res_ref);
    ValeevRef_Finalize();
    #endif
    
    // Finalize stuff 
    Simint_Finalize();

    return 0;
}
