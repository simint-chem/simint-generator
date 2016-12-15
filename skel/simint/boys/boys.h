#pragma once

#include "simint/vectorization/vectorization.h"

#include "simint/boys/boys_taylor.h"
#include "simint/boys/boys_shortgrid.h"
#include "simint/boys/boys_long.h"

#ifdef __cplusplus
extern "C" {
#endif

static inline
void boys_F_split_single(double * restrict F,
                         double const * restrict x,
                         int n)
{
    for(int i = 0; i < SIMINT_SIMD_LEN; i++)
    {
        if(x[i] < BOYS_SHORTGRID_MAXX)
            F[i] = boys_F_taylor_single(x[i], n);
        else
            F[i] = boys_F_long_single(x[i], n); 
    }
}


static inline
void boys_F_split_small_n(double * restrict F,
                          double const * restrict x,
                          int n)
{
    // n is small - just do it all of them via
    // lookup or longfac (no recursion)

    for(int i = 0; i < SIMINT_SIMD_LEN; i++)
    {
        if(x[i] < BOYS_SHORTGRID_MAXX)
            boys_F_taylor(F + i, x[i], n);
        else
            boys_F_long(F + i, x[i], n); 
    }
}


static inline
void boys_F_split_large_n(SIMINT_DBLTYPE * restrict F,
                          SIMINT_DBLTYPE const * restrict x,
                          int n)
{
    // n is large - do only the highest, then recur down

    boys_F_split_single((double *)(F + n), (double *)(x), n);

    // factors for the recursion
    const SIMINT_DBLTYPE x2 = SIMINT_MUL(SIMINT_SET1(2.0), (*x));
    const SIMINT_DBLTYPE ex = SIMINT_EXP(SIMINT_NEG(*x));

    // now recur down
    for(int n2 = n-1; n2 >= 0; n2--)
    {
        const SIMINT_DBLTYPE den = SIMINT_SET1(1.0 / (2.0 * n2 + 1));

        //F[n2] = den * (x2 * F[(n2+1)] + ex);
        F[n2] = SIMINT_MUL(den, ( SIMINT_FMA(x2, F[(n2+1)], ex)));
    }
}

static inline
void boys_F_split(SIMINT_DBLTYPE * restrict F,
                  SIMINT_DBLTYPE const * restrict x,
                  int n)
{
    if(n > 1)
        boys_F_split_large_n(F, x, n);
    else
        boys_F_split_small_n((double *)F, (double *)x, n);
}



#ifdef __cplusplus
}
#endif

