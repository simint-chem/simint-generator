#include <cstdio>
#include <random>

#include "constants.h"
#include "boys/boys.h"
#include "boys/boys_split.h"
#include "boys/boys_FO.h"

#ifdef TESTS_USE_LIBINT2
#include <boys.h> // this is from libint
#endif

#define MAXN 2
#define NTEST 500

#define MAXX 1e12

int main(int argc, char ** argv)
{
    int nval = atoi(argv[1]);
    int ncount = atoi(argv[2]);

    if(nval > MAXN)
    {
        printf("\nError - N value too high!\n");
        return 1;
    }

    Boys_Init(MAXX, MAXN+7);

    double * xval = (double *)ALLOC(ncount*sizeof(double));

    double * res_split = (double *)ALLOC(ncount*sizeof(double));
    double * res_FO = (double *)ALLOC(ncount*sizeof(double));

    #ifdef TESTS_USE_LIBINT2
    double * res_libint_cheby = (double *)ALLOC(ncount*sizeof(double));
    double * res_libint_taylor = (double *)ALLOC(ncount*sizeof(double));

    libint2::FmEval_Chebyshev3 libint_eval_cheby(MAXN);
    libint2::FmEval_Taylor<double, 7> libint_eval_taylor(MAXN, 1e-16);
    #endif

    // set up the RNG
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.0, MAXX); 

    for(int i = 0; i < NTEST; i++)
    {
        // random data
        for(int j = 0; j < ncount; j++)
            xval[j] = dist(gen);

        for(int j = 0; j < ncount; j++)
        {
            double tmp[MAXN+1];
            Boys_F_split(tmp, nval, xval[j]);
            res_split[j] = tmp[nval]; 
        }

        for(int j = 0; j < ncount; j++)
        {
            double tmp[MAXN+1];
            Boys_F_FO(tmp, nval, xval[j]);
            res_FO[j] = tmp[nval]; 
        }

        #ifdef TESTS_USE_LIBINT2
        for(int j = 0; j < ncount; j++)
        {
            double tmp[MAXN+1];
            libint_eval_cheby.eval(tmp, xval[j], nval); 
            res_libint_cheby[j] = tmp[nval]; 
        }

        for(int j = 0; j < ncount; j++)
        {
            double tmp[MAXN+1];
            libint_eval_taylor.eval(tmp, xval[j], nval); 
            res_libint_taylor[j] = tmp[nval]; 
        }
        #endif

        printf("%8.4e %8.4e", res_split[0], res_FO[0]);

        #ifdef TESTS_USE_LIBINT2
        printf(" %8.4e %8.4e", res_libint_cheby[0], res_libint_taylor[0]);
        #endif

        printf("\n");

        
    }

    Boys_Finalize();

    FREE(res_split);
    FREE(res_FO);

    #ifdef TESTS_USE_LIBINT2
    FREE(res_libint_cheby);
    FREE(res_libint_taylor);
    #endif

    return 0;
}
