#ifndef BOYS_SPLIT_H
#define BOYS_SPLIT_H

#include "boys/boys_taylor.h"
#include "boys/boys_long.h"

/////////////////////////
// Inline functions
/////////////////////////
inline double Boys_F0_split(double x)
{
    ASSUME_ALIGN(boys_grid);

    if(x < BOYS_SHORTGRID_MAXX)
        return Boys_F0_taylor(x);
    else
        return Boys_F0_long(x); 
}

inline void Boys_F_split(double * const restrict F, int n, double x)
{
    if(x < BOYS_SHORTGRID_MAXX)
        Boys_F_taylor(F, n, x);
    else
        Boys_F_long(F, n, x); 
}


#endif
