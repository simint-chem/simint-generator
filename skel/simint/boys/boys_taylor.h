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

    #ifndef SIMINT_SCALAR
    #pragma omp simd
    #endif
    for(int i = 0; i <= n; ++i)
    {
        double const * restrict gridpts2 = gridpts + i;

        F[i*SIMINT_SIMD_LEN] = gridpts2[0]
               + dx * (                  gridpts2[1]
               + dx * ( (1.0/2.0   )   * gridpts2[2]
               + dx * ( (1.0/6.0   )   * gridpts2[3]
               + dx * ( (1.0/24.0  )   * gridpts2[4]
               + dx * ( (1.0/120.0 )   * gridpts2[5]
               + dx * ( (1.0/720.0 )   * gridpts2[6]
               + dx * ( (1.0/5040.0)   * gridpts2[7]
               )))))));
    }
}

static inline
double boys_F_taylor_single(double x, int n)
{
    const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(x+BOYS_SHORTGRID_LOOKUPFAC2));
    const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
    const double dx = xi-x;   // -delta x

    double const * restrict gridpts = &(boys_shortgrid[lookup_idx][n]);

    return gridpts[0]
           + dx * (                  gridpts[1]
           + dx * ( (1.0/2.0   )   * gridpts[2]
           + dx * ( (1.0/6.0   )   * gridpts[3]
           + dx * ( (1.0/24.0  )   * gridpts[4]
           + dx * ( (1.0/120.0 )   * gridpts[5]
           + dx * ( (1.0/720.0 )   * gridpts[6]
           + dx * ( (1.0/5040.0)   * gridpts[7]
           )))))));
}

static inline
void boys_F_taylor_vec(SIMINT_DBLTYPE * restrict F, SIMINT_DBLTYPE x, int n)
{
    double * restrict Fd = (double *)F;
    double * restrict xd = (double *)&x;

    #ifndef SIMINT_SCALAR
      // workaround for GCC missing simdlen??
      #if defined __clang__ || defined __INTEL_COMPILER
        #pragma omp simd simdlen(SIMINT_SIMD_LEN)
      #else
        #pragma omp simd
      #endif
    #endif
    for(int v = 0; v < SIMINT_SIMD_LEN; v++)
    {
        const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(xd[v]+BOYS_SHORTGRID_LOOKUPFAC2));
        const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
        const double dx = xi-xd[v];   // -delta x

        double const * restrict gridpts = &(boys_shortgrid[lookup_idx][0]);

        for(int i = 0; i <= n; ++i)
        {
            double const * restrict gridpts2 = gridpts + i;

            Fd[i*SIMINT_SIMD_LEN + v] = gridpts2[0]
                   + dx * (                  gridpts2[1]
                   + dx * ( (1.0/2.0   )   * gridpts2[2]
                   + dx * ( (1.0/6.0   )   * gridpts2[3]
                   + dx * ( (1.0/24.0  )   * gridpts2[4]
                   + dx * ( (1.0/120.0 )   * gridpts2[5]
                   + dx * ( (1.0/720.0 )   * gridpts2[6]
                   + dx * ( (1.0/5040.0)   * gridpts2[7]
                   )))))));
        }
    }
}

static inline
SIMINT_DBLTYPE boys_F_taylor_single_vec(SIMINT_DBLTYPE x, int n)
{
    SIMINT_DBLTYPE ret;
    double * retd = (double *)&ret;
    double * restrict xd = (double *)&x;

    #ifndef SIMINT_SCALAR
      // workaround for GCC missing simdlen??
      #if defined __clang__ || defined __INTEL_COMPILER
        #pragma omp simd simdlen(SIMINT_SIMD_LEN)
      #else
        #pragma omp simd
      #endif
    #endif
    for(int v = 0; v < SIMINT_SIMD_LEN; v++)
    {
        const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(xd[v]+BOYS_SHORTGRID_LOOKUPFAC2));
        const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
        const double dx = xi-xd[v];   // -delta x

        double const * restrict gridpts = &(boys_shortgrid[lookup_idx][n]);

        retd[v] = gridpts[0]
               + dx * (                  gridpts[1]
               + dx * ( (1.0/2.0   )   * gridpts[2]
               + dx * ( (1.0/6.0   )   * gridpts[3]
               + dx * ( (1.0/24.0  )   * gridpts[4]
               + dx * ( (1.0/120.0 )   * gridpts[5]
               + dx * ( (1.0/720.0 )   * gridpts[6]
               + dx * ( (1.0/5040.0)   * gridpts[7]
               )))))));
    }

    return ret;
}


#ifdef __cplusplus
}
#endif

