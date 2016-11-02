#pragma once

#include "simint/boys/boys_shortgrid.h"
#include "simint/vectorization/vectorization.h"

#ifdef __cplusplus
extern "C" {
#endif

extern double boys_shortgrid[BOYS_SHORTGRID_NPOINT][BOYS_SHORTGRID_MAXN+1];

static inline
void boys_F_taylor(double * restrict F, double x, int n)
{

    const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(x+BOYS_SHORTGRID_LOOKUPFAC2));
    const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
    const double dx = xi-x;   // -delta x

    double const * restrict gridpts = &(boys_shortgrid[lookup_idx][0]);

    for(int i = 0; i <= n; ++i)
    {
        double fxi[8];
        for(int fi = 0; fi < 8; fi++)
            fxi[fi] = gridpts[i+fi];

        F[i*SIMINT_SIMD_LEN] = fxi[0]
               + dx * (                  fxi[1]
               + dx * ( (1.0/2.0   )   * fxi[2]
               + dx * ( (1.0/6.0   )   * fxi[3]
               + dx * ( (1.0/24.0  )   * fxi[4]
               + dx * ( (1.0/120.0 )   * fxi[5]
               + dx * ( (1.0/720.0 )   * fxi[6]
               + dx * ( (1.0/5040.0)   * fxi[7]
               )))))));
    }
}

static inline
double boys_F_taylor_single(double x, int n)
{
    double fxi[8];

    const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(x+BOYS_SHORTGRID_LOOKUPFAC2));
    const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
    const double dx = xi-x;   // -delta x

    double const * restrict gridpts = &(boys_shortgrid[lookup_idx][0]);

    for(int fi = 0; fi < 8; fi++)
        fxi[fi] = gridpts[n+fi];

    return fxi[0]
           + dx * (                  fxi[1]
           + dx * ( (1.0/2.0   )   * fxi[2]
           + dx * ( (1.0/6.0   )   * fxi[3]
           + dx * ( (1.0/24.0  )   * fxi[4]
           + dx * ( (1.0/120.0 )   * fxi[5]
           + dx * ( (1.0/720.0 )   * fxi[6]
           + dx * ( (1.0/5040.0)   * fxi[7]
           )))))));
}


#ifdef __cplusplus
}
#endif

