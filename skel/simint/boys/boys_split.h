#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "simint/boys/boys_taylor.h"
#include "simint/boys/boys_long.h"

/////////////////////////
// Inline functions
/////////////////////////
static inline
void Boys_F_split(double * restrict F, int n, double x)
{
    if(x < BOYS_SHORTGRID_MAXX)
        Boys_F_taylor(F, n, x);
    else
        Boys_F_long(F, n, x); 
}

static inline
void Boys_F_split_simd(double * restrict Farr, int n, double const * restrict xvec)
{
    for(int i = 0; i < SIMINT_SIMD_LEN; i++)
    {
        if(xvec[i] < BOYS_SHORTGRID_MAXX)
            Boys_F_taylor_simd(Farr + i, n, xvec[i]);
        else
            Boys_F_long_simd(Farr + i, n, xvec[i]); 
    }
}


static inline
void Boys_F_split_single(double * restrict Farr, int n, double const * restrict xvec)
{
    for(int i = 0; i < SIMINT_SIMD_LEN; i++)
    {
        if(xvec[i] < BOYS_SHORTGRID_MAXX)
            Farr[i] = Boys_F_taylor_single(n, xvec[i]);
        else
            Farr[i] = Boys_F_long_single(n, xvec[i]); 
    }
}

#ifdef __cplusplus
}
#endif

