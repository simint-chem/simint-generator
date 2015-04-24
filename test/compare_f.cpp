#include <cstdio>
#include <cmath>

#include "constants.h"
#include "boys/boys.h"
#include "boys/boys_split.h"
#include "boys/boys_FO.h"

#include "valeev.hpp"

#ifdef TESTS_USE_LIBINT2
#include <boys.h> // this is from libint
#endif

#define MAXN 2

int main(int argc, char ** argv)
{
    int nval = atoi(argv[1]);
    double maxx = atof(argv[2]);

    if(nval > MAXN)
    {
        printf("\nError - N value too high!\n");
        return 1;
    }

    Valeev_Init();
    Boys_Init(maxx, MAXN+7);

    double maxerr_taylor = 0.0;
    double maxerr_FO = 0.0;

    double maxrelerr_taylor = 0.0;
    double maxrelerr_FO = 0.0;

    double valeev[MAXN+1];
    double my_taylor[MAXN+1];
    double my_FO[MAXN+1];

    #ifdef TESTS_USE_LIBINT2
    double maxerr_libint_cheby = 0.0;
    double maxerr_libint_taylor = 0.0;

    double maxrelerr_libint_cheby = 0.0;
    double maxrelerr_libint_taylor = 0.0;

    double libint_cheby[MAXN+1];
    double libint_taylor[MAXN+1];
    libint2::FmEval_Chebyshev3 libint_eval_cheby(MAXN);
    libint2::FmEval_Taylor<double, 7> libint_eval_taylor(MAXN, 1e-16);
    #endif

    for(double x = 0.0; x < maxx; x += 0.00005)
    {
        Valeev_F(valeev, nval, x);
        Boys_F_taylor(my_taylor, nval, x);
        Boys_F_FO(my_FO, nval, x);

        double diff_taylor = fabs(my_taylor[nval] - valeev[nval]);
        double diff_FO     = fabs(my_FO[nval]     - valeev[nval]);
        double rel_taylor = diff_taylor / valeev[nval];
        double rel_FO = diff_FO / valeev[nval];

        if(diff_FO > maxerr_FO)
            maxerr_FO = diff_FO;
        if(diff_taylor > maxerr_taylor)
            maxerr_taylor = diff_taylor;

        if(rel_taylor > maxrelerr_taylor)
            maxrelerr_taylor = rel_taylor; 
        if(rel_FO > maxrelerr_FO)
            maxrelerr_FO = rel_FO; 

        #ifdef TESTS_USE_LIBINT2
        libint_eval_cheby.eval(libint_cheby, x, nval);
        libint_eval_taylor.eval(libint_taylor, x, nval);

        double diff_libint_cheby  = fabs(libint_cheby[nval]  - valeev[nval]);
        double diff_libint_taylor = fabs(libint_taylor[nval] - valeev[nval]);
        double rel_libint_taylor = diff_libint_taylor / valeev[nval];
        double rel_libint_cheby = diff_libint_cheby / valeev[nval];

        if(diff_libint_taylor > maxerr_libint_taylor)
            maxerr_libint_taylor = diff_libint_taylor;
        if(diff_libint_cheby > maxerr_libint_cheby)
            maxerr_libint_cheby = diff_libint_cheby;

        if(rel_libint_taylor > maxrelerr_libint_taylor)
            maxrelerr_libint_taylor = rel_libint_taylor;
        if(rel_libint_cheby > maxrelerr_libint_cheby)
            maxrelerr_libint_cheby = rel_libint_cheby;
        #endif


        //printf("%12.8f     %25.17e  %25.17e  %25.17e\n", x, valeev[nval], my_taylor[nval], my_FO[nval]);
    }

    printf("\n");
    printf("***Max abserr - taylor: %12.8e\n", maxerr_taylor);
    printf("***Max relerr - taylor: %12.7e\n", maxrelerr_taylor);
    printf("\n");
    printf("***Max abserr -     FO: %12.8e\n", maxerr_FO);
    printf("***Max relerr -     FO: %12.8e\n", maxrelerr_FO);
    printf("\n");

    #ifdef TESTS_USE_LIBINT2
    printf("***Max abserr -  libint cheby: %12.8e\n", maxerr_libint_cheby);
    printf("***Max relerr -  libint cheby: %12.8e\n", maxrelerr_libint_cheby);
    printf("\n");
    printf("***Max abserr - libint taylor: %12.8e\n", maxerr_libint_taylor);
    printf("***Max relerr - libint taylor: %12.7e\n", maxrelerr_libint_taylor);
    printf("\n");
    #endif

    printf("\n");

    Valeev_Finalize();
    Boys_Finalize();

    return 0;
}
