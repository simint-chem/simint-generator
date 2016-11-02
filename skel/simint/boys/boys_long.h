#pragma once

#include "simint/boys/boys_longfac.h"
#include "simint/vectorization/vectorization.h"

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

// from boys_longfac.h
//extern double const boys_longfac[BOYS_LONGFAC_MAXN];

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


#ifdef __cplusplus
}
#endif

