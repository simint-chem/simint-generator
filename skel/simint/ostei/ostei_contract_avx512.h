#pragma once

#include "simint/vectorization/vectorization.h"

#ifdef __cplusplus
extern "C" {
#endif


static inline
void contract(int ncart,
              int const * restrict offsets,
              __m512d const * restrict PRIM_INT,
              double * restrict PRIM_PTR)
{
    for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
    {
        double const * restrict prim_int_tmp = (double *)PRIM_INT + n;
        double * restrict prim_ptr_tmp = PRIM_PTR + offsets[n]*ncart;

        for(int np = 0; np < ncart; ++np)
        {
            prim_ptr_tmp[np] += *prim_int_tmp;
            prim_int_tmp += SIMINT_SIMD_LEN;
        }
    }
}


static inline
void contract_all(int ncart,
                  __m512d const * restrict PRIM_INT,
                  double * restrict PRIM_PTR)
{
    #ifndef SIMINT_GCC

    for(int np = 0; np < ncart; np++)
        PRIM_PTR[np] += _mm512_reduce_add_pd(PRIM_INT[np]);

    #else

    // GCC is missing some instructions from the above block
    int offsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    contract(ncart, offsets, PRIM_INT, PRIM_PTR);

    #endif
}


static inline
void contract_fac(int ncart,
                  int const * restrict offsets,
                  __m512d const * restrict factor,
                  __m512d const * restrict PRIM_INT,
                  double * restrict PRIM_PTR)
{
    for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
    {
        double const * restrict prim_int_tmp = (double *)PRIM_INT + n;
        double * restrict prim_ptr_tmp = PRIM_PTR + offsets[n]*ncart;
        double factor_tmp = *((double *)factor + n);

        for(int np = 0; np < ncart; ++np)
        {
            prim_ptr_tmp[np] += factor_tmp * (*prim_int_tmp);
            prim_int_tmp += SIMINT_SIMD_LEN;
        }
    }
}


static inline
void contract_all_fac(int ncart,
                      __m512d const * restrict PRIM_INT,
                      __m512d const * restrict factor,
                      double * restrict PRIM_PTR)
{
    #ifndef SIMINT_GCC

    for(int np = 0; np < ncart; np++)
        PRIM_PTR[np] += _mm512_reduce_add_pd(_mm512_mul_pd((*factor), PRIM_INT[np]));

    #else

    // GCC is missing some instructions from the above block
    int offsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    contract_fac(ncart, offsets, factor, PRIM_INT, PRIM_PTR);

    #endif
}


#ifdef __cplusplus
}
#endif

