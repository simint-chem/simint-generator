#include "simint/boys/boys.h"
#include "simint/boys/boys_taylor.h"
#include "simint/boys/boys_long.h"

void boys_F_split(double * restrict Farr, int n, double const * restrict xvec)
{
    for(int i = 0; i < SIMINT_SIMD_LEN; i++)
    {
        if(xvec[i] < BOYS_SHORTGRID_MAXX)
            boys_F_taylor(Farr + i, n, xvec[i]);
        else
            boys_F_long(Farr + i, n, xvec[i]); 
    }
}

