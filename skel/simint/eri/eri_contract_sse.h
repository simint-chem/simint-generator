#pragma once

#include "simint/vectorization/vectorization.h"

static inline
void contract(int ncart,
              int const * restrict offsets,
              __m128d const * restrict PRIM_INT,
              double * restrict PRIM_PTR)
{
    int n, np;

    for(n = 0; n < SIMINT_SIMD_LEN; ++n)
    {
        double const * restrict prim_int_tmp = (double *)PRIM_INT + n;
        double * restrict prim_ptr_tmp = PRIM_PTR + offsets[n]*ncart;

        for(np = 0; np < ncart; ++np)
        {
            prim_ptr_tmp[np] += *prim_int_tmp;
            prim_int_tmp += SIMINT_SIMD_LEN;
        }
    }
}



static inline
void contract_all(int ncart,
                  __m128d const * restrict PRIM_INT,
                  double * restrict PRIM_PTR)
{
    //TODO efficient algorithm for this
    int offsets[] = {0, 0};
    return contract(ncart, offsets, PRIM_INT, PRIM_PTR);
}
