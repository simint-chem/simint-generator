#ifndef BOYS_SPLIT_H
#define BOYS_SPLIT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "boys/boys_taylor.h"
#include "boys/boys_long.h"

/////////////////////////
// Inline functions
/////////////////////////
inline void Boys_F_split(double * const restrict F, int n, double x)
{
    if(x < BOYS_SHORTGRID_MAXX)
        Boys_F_taylor(F, n, x);
    else
        Boys_F_long(F, n, x); 
}

inline void Boys_F_split_simd(double * const restrict Farr, int n, double const * const restrict xvec)
{
    for(int i = 0; i < SIMINT_SIMD_ALIGN_DBL; i++)
    {
        if(xvec[i] < BOYS_SHORTGRID_MAXX)
            Boys_F_taylor_simd(Farr + i, n, xvec[i]);
        else
            Boys_F_long_simd(Farr + i, n, xvec[i]); 
    }
}

#ifdef __cplusplus
}
#endif

#endif
