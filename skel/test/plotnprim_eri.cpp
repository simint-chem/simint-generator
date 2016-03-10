#include <cstdio>
#include <sstream>
#include <cmath>

#include "simint/vectorization/vectorization.h"
#include "simint/eri/eri.h"
#include "test/Simint.hpp"
#include "test/Common.hpp"
#include "test/Timer.h"


#ifdef BENCHMARK_VALIDATE
#include "test/valeev.hpp"
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


    // read in the shell info
    ShellMap shellmap = ReadBasis(basfile);

    // normalize
    for(auto & it : shellmap)
        normalize_gaussian_shells(it.second.size(), it.second.data());


    // find the max dimensions
    std::array<int, 3> maxparams = FindMapMaxParams(shellmap);
    const int maxam = (maxparams[0] > SIMINT_ERI_MAXAM ? SIMINT_ERI_MAXAM : maxparams[0]);
    const int maxsize = maxparams[2];

    /* Storage of integrals */
    double * res_ints = (double *)ALLOC(maxsize * sizeof(double)); 

    /* contracted workspace */
    double * simint_work = (double *)ALLOC(SIMINT_ERI_MAX_WORKMEM);

    // initialize stuff
    // nothing needs initializing!


    #ifdef BENCHMARK_VALIDATE
    ValeevRef_Init();
    double * res_ref = (double *)ALLOC(maxsize * sizeof(double));
    #endif

    // Timing header
    printf("%13s %12s   %12s   %16s   %16s   %12s\n",
                           "Quartet", "NCont", "NPrim", "Ticks(Prep)", "Ticks(Ints)", "Ticks/Prim");

    // running totals
    size_t ncont1234_total = 0;
    size_t nprim1234_total = 0;
    size_t nshell1234_total = 0;
    std::pair<TimerType, TimerType> time_total{0,0};


    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    for(int k = 0; k <= maxam; k++)
    for(int l = 0; l <= maxam; l++)
    {
        if(!ValidQuartet({i,j,k,l}))
            continue;

        std::map<size_t, std::vector<size_t>> ticksperprim;

        // total amount of time for this AM quartet
        TimerType time_am = 0;


        #ifdef BENCHMARK_VALIDATE
        std::pair<double, double> err{0.0, 0.0};
        #endif

        // running totals for this am
        size_t nprim1234_am = 0;
        size_t nshell1234_am = 0;
        size_t ncont1234_am = 0;

        // do one shell pair at a time on the bra side
        for(size_t a = 0; a < shellmap[i].size(); a++)
        for(size_t b = 0; b < shellmap[j].size(); b++)
        for(size_t c = 0; c < shellmap[k].size(); c++)
        for(size_t d = 0; d < shellmap[l].size(); d++)
        {
            const size_t nshell1 = 1;
            const size_t nshell2 = 1;
            const size_t nshell3 = 1;
            const size_t nshell4 = 1;


            gaussian_shell const * const A = &shellmap[i][a];
            gaussian_shell const * const B = &shellmap[j][b];
            gaussian_shell const * const C = &shellmap[k][c];
            gaussian_shell const * const D = &shellmap[l][d];

            // time creation of P
            TimerType time_pair_12_0 = 0;
            TimerType time_pair_12_1 = 0;
            CLOCK(time_pair_12_0);
            struct multishell_pair P = create_multishell_pair(nshell1, A, nshell2, B);
            CLOCK(time_pair_12_1);
            time_am += time_pair_12_1 - time_pair_12_0;

            // time creation of Q
            TimerType time_pair_34_0 = 0;
            TimerType time_pair_34_1 = 0;
            CLOCK(time_pair_34_0);
            struct multishell_pair Q = create_multishell_pair(nshell3, C, nshell4, D);
            CLOCK(time_pair_34_1);
            time_am += time_pair_34_1 - time_pair_34_0;

            // actually calculate
            TimerType qtime = Simint_Integral(P, Q, simint_work, res_ints);
            time_am += qtime;

            // acutal number of primitives and shells calculated
            // TODO - replace with return values from Integrals
            const size_t nprim1234 = P.nprim * Q.nprim;
            const size_t nshell1234 = P.nshell12 * Q.nshell12;
            const size_t ncart1234 = NCART(i) * NCART(j) * NCART(k) * NCART(l);
            const size_t ncont1234 = nshell1234 * ncart1234;

            #ifdef BENCHMARK_VALIDATE
            ValeevRef_Integrals(A, nshell1,
                                B, nshell2,
                                C, nshell3,
                                D, nshell4,
                                res_ref, false);
            std::pair<double, double> err2 = CalcError(res_ints, res_ref, ncont1234);

            err.first = std::max(err.first, err2.first);
            err.second = std::max(err.second, err2.second);
            #endif


            free_multishell_pair(P);
            free_multishell_pair(Q);

            // add primitive and shell count to running totals for this am
            ncont1234_am += ncont1234;
            nprim1234_am += nprim1234;
            nshell1234_am += nshell1234;

            // add timings
            ticksperprim[nprim1234].push_back(qtime);
        }

        // add primitive and shell count to overall running totals
        ncont1234_total += ncont1234_am;
        nprim1234_total += nprim1234_am;
        nshell1234_total += nshell1234_am;
        time_total.first += 0;
        time_total.second += time_am;

        printf("( %d %d | %d %d ) %12lu   %12lu   %16llu   %16llu   %12.3f\n",
                                                                      i, j, k, l,
                                                                      nshell1234_am, nprim1234_am,
                                                                      0ull, time_am,
                                                                      (double)(time_am)/(double)(nprim1234_am));


        // write out the times
        std::stringstream fname;
        fname << "eri_nprim_" << i << "_" << j << "_" << k << "_" << l << ".out";
        std::ofstream nprimout(fname.str().c_str());
        for(const auto & it : ticksperprim)
        {
            double avg = 0;
            for(auto it2 : it.second)
                avg += static_cast<double>(it2);
            avg /= static_cast<double>(it.second.size());

            double stddev = 0;
            for(auto it2 : it.second)
                stddev += (it2 - avg)*(it2-avg);
            stddev = sqrt(stddev/(static_cast<double>(it.second.size()) - 1));
            nprimout << it.first << "    " << avg << "    " << stddev << "\n";
        }

        #ifdef BENCHMARK_VALIDATE
        printf("( %d %d | %d %d ) MaxAbsErr: %10.3e   MaxRelErr: %10.3e\n", i, j, k, l, err.first, err.second);
        #endif
    }

    printf("\n");
    printf("Calculated %ld contracted integrals\n", static_cast<long>(ncont1234_total));
    printf("Calculated %ld contracted shells\n", static_cast<long>(nshell1234_total));
    printf("Calculated %ld primitive integrals\n", static_cast<long>(nprim1234_total));
    printf("Total ticks to calculate all:  %llu\n", time_total.first + time_total.second);
    printf("\n");

    FreeShellMap(shellmap);
    FREE(res_ints);


    #ifdef BENCHMARK_VALIDATE
    FREE(res_ref);
    ValeevRef_Finalize();
    #endif
    
    // Finalize stuff 
    Simint_Finalize();

    return 0;
}
