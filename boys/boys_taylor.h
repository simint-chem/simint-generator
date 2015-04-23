#ifndef BOYS_TAYLOR_H
#define BOYS_TAYLOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "vectorization.h"
#include "boys/boys_taylorgrid.h"

extern double const * const restrict * const restrict boys_grid;

/////////////////////////
// Inline functions
/////////////////////////
inline double Boys_F0_taylor(double x)
{
    ASSUME_ALIGN(boys_grid);

    const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(x+BOYS_SHORTGRID_LOOKUPFAC2));
    const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
    const double dx = xi-x;   // -delta x

    double const * const restrict gridpts = &(boys_grid[lookup_idx][0]);

/*
    const double f0xi = gridpts[0];
    const double f1xi = gridpts[1];
    const double f2xi = gridpts[2];
    const double f3xi = gridpts[3];
    const double f4xi = gridpts[4];
    const double f5xi = gridpts[5];
    const double f6xi = gridpts[6];
    const double f7xi = gridpts[7];

    return f0xi
           + dx * (                  f1xi
           + dx * ( (1.0/2.0   )   * f2xi
           + dx * ( (1.0/6.0   )   * f3xi
           + dx * ( (1.0/24.0  )   * f4xi
           + dx * ( (1.0/120.0 )   * f5xi
           + dx * ( (1.0/720.0 )   * f6xi
           + dx * ( (1.0/5040.0)   * f7xi
           )))))));
*/
    const double f0xi = gridpts[0];
    const double f1xi = gridpts[1];
    const double f2xi = gridpts[2];
    const double f3xi = gridpts[3];
    const double f4xi = gridpts[4];
    const double f5xi = gridpts[5];
    const double f6xi = gridpts[6];

    return f0xi
           + dx * (                  f1xi
           + dx * ( (1.0/2.0   )   * f2xi
           + dx * ( (1.0/6.0   )   * f3xi
           + dx * ( (1.0/24.0  )   * f4xi
           + dx * ( (1.0/120.0 )   * f5xi
           + dx * ( (1.0/720.0 )   * f6xi
           ))))));
}

inline void Boys_F_taylor(double * const restrict F, int n, double x)
{
    ASSUME_ALIGN(boys_grid);

    const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(x+BOYS_SHORTGRID_LOOKUPFAC2));
    const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
    const double dx = xi-x;   // -delta x

    double const * const restrict gridpts = &(boys_grid[lookup_idx][0]);

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

        F[i] = f0xi
               + dx * (                  f1xi
               + dx * ( (1.0/2.0   )   * f2xi
               + dx * ( (1.0/6.0   )   * f3xi
               + dx * ( (1.0/24.0  )   * f4xi
               + dx * ( (1.0/120.0 )   * f5xi
               + dx * ( (1.0/720.0 )   * f6xi
               + dx * ( (1.0/5040.0)   * f7xi
               )))))));
    }
}


#ifdef __cplusplus
}
#endif

#endif
