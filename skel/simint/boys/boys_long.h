#pragma once

#include "simint/boys/boys_longfac.h"
#include "simint/vectorization/vectorization.h"

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline
void boys_F_long(double * restrict F, double x, int n)
{
    const double x1 = 1.0/x;
    double x2 = sqrt(x1);

     for(int i = 0; i <= n; i++)
     {
        *F = boys_longfac[i] * x2;
        x2 *= x1;
        F += SIMINT_SIMD_LEN;
     }
}

static inline
double boys_F_long_single(double x, int n)
{
    const double p = -(2*n+1);
    const double x2 = pow(x, p);
    return boys_longfac[n] * sqrt(x2);
}

static inline
void boys_F_long_vec(SIMINT_DBLTYPE * restrict F,
                     SIMINT_DBLTYPE x,
                     int n)
{
    const SIMINT_DBLTYPE x1 = SIMINT_DIV(SIMINT_DBLSET1(1.0), x);
    SIMINT_DBLTYPE x2 = SIMINT_SQRT(x1);

    for(int i = 0; i <= n; i++)
    {
        const SIMINT_DBLTYPE lfac = SIMINT_DBLSET1(boys_longfac[i]);
        F[i] = SIMINT_MUL(lfac, x2);
        x2 = SIMINT_MUL(x2, x1);
    }
}

static inline
SIMINT_DBLTYPE boys_F_long_single_vec(SIMINT_DBLTYPE x, int n)
{
    const SIMINT_DBLTYPE p = SIMINT_DBLSET1(-((double)n+0.5));
    const SIMINT_DBLTYPE x2 = SIMINT_POW(x, p);
    const SIMINT_DBLTYPE lfac = SIMINT_DBLSET1(boys_longfac[n]);
    return SIMINT_MUL(lfac, x2);
}


#ifdef __cplusplus
}
#endif

