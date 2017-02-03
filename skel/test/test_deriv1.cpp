#include <cstdio>
#include <atomic>
#include <cmath>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "simint/simint.h"
#include "test/Common.hpp"
#include "test/ValeevRef.hpp"


#define SIMINT_SCREEN 0
#define SIMINT_SCREEN_TOL 0.0


typedef std::array<int, 4> QAM;
typedef std::pair<double, double> ErrorPair;

// Prog -> { QAM -> { Absolute error, Relative Error } }
typedef std::map<QAM, ErrorPair> ErrorMap;

static void UpdateErrorMap(ErrorMap & m, QAM am, ErrorPair p)
{
    ErrorPair val = m.at(am);
    val.first = std::max(val.first, p.first);
    val.second = std::max(val.second, p.second);
    m[am] = val;
}


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
    ShellMap shellmap = ReadBasis(basfile);


    // normalize the original
    for(auto & it : shellmap)
        simint_normalize_shells(it.second.size(), it.second.data());


    // find the max dimensions
    std::pair<int, int> maxparams = FindMaxParams(shellmap);
    const int maxam = (maxparams.first > SIMINT_OSTEI_DERIV1_MAXAM ? SIMINT_OSTEI_DERIV1_MAXAM : maxparams.first);
    const int max_ncart = ( (maxam+1)*(maxam+2) )/2;
    const int maxsize = maxparams.second * maxparams.second * max_ncart * max_ncart * 12;

    // get the number of threads
    #ifdef _OPENMP
        const int nthread = omp_get_max_threads();
    #else
        const int nthread = 1;
    #endif

    /* contracted workspace */
    double * all_simint_work = (double *)SIMINT_ALLOC(nthread * SIMINT_OSTEI_MAX_WORKMEM);

    /* Storage of integrals */
    double * all_res_simint = (double *)SIMINT_ALLOC(nthread * maxsize * sizeof(double));
    double * all_res_valeev = (double *)SIMINT_ALLOC(nthread * maxsize * sizeof(double));


    // Map containing the errors
    // QAM -> (absolute error, relative error)
    ErrorMap errors;


    // initialize stuff
    ValeevRef_Init();

    // Set initial errors to zero
    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    for(int k = 0; k <= maxam; k++)
    for(int l = 0; l <= maxam; l++)
        errors[{{i,j,k,l}}] = {0,0};


    // Print the header for the final results table
    printf("\n");
    printf("%17s  %10s    %10s\n", "Quartet", "MaxErr", "MaxRelErr");

    // Number of contracted integrals calculated
    std::atomic<long> ncont(0);


    // loop over all AM quartets, choosing only valid ones
    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    for(int k = 0; k <= maxam; k++)
    for(int l = 0; l <= maxam; l++)
    {
        const int nshell3 = shellmap[k].size();
        const int nshell4 = shellmap[l].size();

        // create a multishell pair for all | k l ) quartets
        struct simint_multi_shellpair Q;
        simint_initialize_multi_shellpair(&Q);
        simint_create_multi_shellpair(nshell3, shellmap[k].data(),
                                      nshell4, shellmap[l].data(), &Q,
                                      SIMINT_SCREEN);


        // do bra one at a time
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(size_t a = 0; a < shellmap[i].size(); a++)
        for(size_t b = 0; b < shellmap[j].size(); b++)
        {
            #ifdef _OPENMP
                const int ithread = omp_get_thread_num();
            #else
                const int ithread = 0;
            #endif

            double * res_simint = all_res_simint + ithread * maxsize;
            double * simint_work = all_simint_work + ithread * SIMINT_OSTEI_MAX_WORKSIZE;
            double * res_valeev = all_res_valeev + ithread * maxsize;

            const int nshell1 = 1;
            const int nshell2 = 1;

            struct simint_multi_shellpair P;
            simint_initialize_multi_shellpair(&P);
            simint_create_multi_shellpair(nshell1, &shellmap[i][a],
                                          nshell2, &shellmap[j][b], &P,
                                          SIMINT_SCREEN);


            // acutal number of primitives and shells that
            // will be calculated
            const int ncart1234 = NCART(i) * NCART(j) * NCART(k) * NCART(l);
            const int nshell1234 = P.nshell12 * Q.nshell12;
            const int ncont1234 = nshell1234 * ncart1234;


            /////////////////////////////////
            // Calculate valeev
            // reference integrals
            /////////////////////////////////
            ValeevRef_Integrals(&shellmap[i][a], nshell1,
                                &shellmap[j][b], nshell2,
                                &shellmap[k][0], nshell3,
                                &shellmap[l][0], nshell4,
                                res_valeev, 1, false);



            ////////////////////////////
            // Calculate the integrals
            ////////////////////////////
            int simint_ret = simint_compute_eri_deriv(1, &P, &Q, SIMINT_SCREEN_TOL, simint_work, res_simint);

            // if the return is < 0, it didn't calculate anything
            // (everything was screened)
            if(simint_ret < 0)
                std::fill(res_simint, res_simint + nshell1234*ncart1234*12, 0.0);

            /////////////////////////////////
            // Update the error map
            /////////////////////////////////
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                UpdateErrorMap(errors, {{i, j, k, l}}, CalcError(res_simint, res_valeev, ncont1234*12));
            }

/*
            // For debugging
            int m = 0;
            for(int m1 = 0; m1 < nshell1; m1++)
            for(int m2 = 0; m2 < nshell2; m2++)
            for(int m3 = 0; m3 < nshell3; m3++)
            for(int m4 = 0; m4 < nshell4; m4++)
            for(int c1 = 0; c1 < NCART(i); c1++)
            for(int c2 = 0; c2 < NCART(j); c2++)
            for(int c3 = 0; c3 < NCART(k); c3++)
            for(int c4 = 0; c4 < NCART(l); c4++)
            for(int center = 0; center < 4; center++)
            for(int d = 0; d < 3; d++, m++)
            {
                double diff_simint = fabs(res_valeev[m] - res_simint[m]);
                double rdiff_simint = fabs(diff_simint / res_valeev[m]);


//                if( (diff_simint > 1e-14 && rdiff_simint > 1e-8) )
                {
                    printf(" [%d/%d] %d %d %d %d : %d %d %d %d : %d %d", ithread, nthread, int(a+m1), int(b+m2), m3, m4, c1, c2, c3, c4, center, d);

                    printf("   %25.16e  %25.16e  %25.16e", res_simint[m], res_valeev[m], rdiff_simint);
                    printf("\n");
                }
            }
*/

            simint_free_multi_shellpair(&P);

            ncont += ncont1234;


        } // end threaded loop over a,b


        simint_free_multi_shellpair(&Q);

        // print out the errors for this quartet
        printf("( %2d %2d | %2d %2d )  ", i, j, k, l);

        // should be the same order as the header, right?
        double abserr = errors.at({{i, j, k, l}}).first;
        double relerr = errors.at({{i, j, k, l}}).second;
        bool bad = (abserr > 1e-14 && relerr > 1e-8);
        printf("  %10.3e  %10.3e  %s\n", abserr, relerr, bad ? "***" : "");
    }

    printf("\n");
    printf("Calculated %ld contracted integrals\n", static_cast<long>(ncont));
    printf("\n");

    FreeShellMap(shellmap);
    ValeevRef_Finalize();
    simint_finalize();

    SIMINT_FREE(all_res_simint);
    SIMINT_FREE(all_simint_work);
    SIMINT_FREE(all_res_valeev);

    return 0;
}
