#ifndef BOYS_LONG_H
#define BOYS_LONG_H

#include <math.h>
#include "boys/boys_longfac.h"

extern double const boys_longfac[BOYS_LONGFAC_MAXN];

inline double Boys_F0_long(double x)
{
    return boys_longfac[0] / sqrt(x);
}

inline void Boys_F_long(double * const restrict F, int n, double x)
{
    const double x1 = 1.0/x;
    double x2 = sqrt(x1);

    for(int i = 0; i <= n; i++)
    {
        F[i] = boys_longfac[i] * x2;
        x2 *= x1;
    }
}

#endif
