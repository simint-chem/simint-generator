#ifndef BOYS_LONG_H
#define BOYS_LONG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include "boys/boys_longfac.h"

extern double const boys_longfac[BOYS_LONGFAC_MAXN];

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


#ifdef __cplusplus
}
#endif

#endif
