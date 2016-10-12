#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "simint/boys/boys_shortgrid.h"
#include "simint/vectorization/vectorization.h"

extern double boys_shortgrid[BOYS_SHORTGRID_NPOINT][BOYS_SHORTGRID_MAXN+1];

static inline
void boys_F_taylor(double * restrict F, int n, double x)
{
    const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(x+BOYS_SHORTGRID_LOOKUPFAC2));
    const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
    const double dx = xi-x;   // -delta x

    double const * restrict gridpts = &(boys_shortgrid[lookup_idx][0]);

    for(int i = 0; i <= n; ++i)
    {
        const double f0xi = gridpts[i];
        const double f1xi = gridpts[i+1];
        const double f2xi = gridpts[i+2];
        const double f3xi = gridpts[i+3];
        const double f4xi = gridpts[i+4];
        const double f5xi = gridpts[i+5];
        const double f6xi = gridpts[i+6];
        const double f7xi = gridpts[i+7];

        F[0] = f0xi
               + dx * (                  f1xi
               + dx * ( (1.0/2.0   )   * f2xi
               + dx * ( (1.0/6.0   )   * f3xi
               + dx * ( (1.0/24.0  )   * f4xi
               + dx * ( (1.0/120.0 )   * f5xi
               + dx * ( (1.0/720.0 )   * f6xi
               + dx * ( (1.0/5040.0)   * f7xi
               )))))));

        F += SIMINT_SIMD_LEN;
    }
}


#ifdef __cplusplus
}
#endif

