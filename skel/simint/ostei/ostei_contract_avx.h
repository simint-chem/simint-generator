#pragma once

#include "simint/vectorization/vectorization.h"

#ifdef __cplusplus
extern "C" {
#endif


static inline
void contract(int ncart,
              int const * restrict offsets,
              __m256d const * restrict PRIM_INT,
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
                  __m256d const * restrict PRIM_INT,
                  double * restrict PRIM_PTR)
{
    #ifndef SIMINT_GCC

    int n, n4;
    const int nbatch = ncart/4;

    for(n = 0, n4 = 0; n < nbatch; n++, n4 += 4)
    {
        __m256d t1 = _mm256_hadd_pd(PRIM_INT[n4],   PRIM_INT[n4+1]);
        __m256d t2 = _mm256_hadd_pd(PRIM_INT[n4+2], PRIM_INT[n4+3]);
        __m256d t3 = _mm256_set_m128d(  _mm256_extractf128_pd(t2, 0) + _mm256_extractf128_pd(t2, 1),
                                        _mm256_extractf128_pd(t1, 0) + _mm256_extractf128_pd(t1, 1));
        _mm256_storeu_pd(PRIM_PTR + n4, _mm256_loadu_pd(PRIM_PTR + n4) + t3);
    }

    const int left = ncart % 4;
    const int done = ncart - left;
    for(n = done; n < ncart; n++)
    {
        union double4 tmp = { PRIM_INT[n] };
        PRIM_PTR[n] += tmp.d[0] + tmp.d[1] + tmp.d[2] + tmp.d[3];
    }

    #else

    // GCC is missing some instructions from the above block
    int offsets[4] = {0, 0, 0, 0};
    contract(ncart, offsets, PRIM_INT, PRIM_PTR);

    #endif
}


static inline
void contract_fac(int ncart,
                  int const * restrict offsets,
                  __m256d const * restrict factor,
                  __m256d const * restrict PRIM_INT,
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
                      __m256d const * restrict PRIM_INT,
                      __m256d const * restrict factor,
                      double * restrict PRIM_PTR)
{
    #ifdef SIMINT_INTEL
    int n, n4;
    const int nbatch = ncart/4;

    for(n = 0, n4 = 0; n < nbatch; n++, n4 += 4)
    {
        __m256d t1 = _mm256_hadd_pd(PRIM_INT[n4],   PRIM_INT[n4+1]);
        __m256d t2 = _mm256_hadd_pd(PRIM_INT[n4+2], PRIM_INT[n4+3]);
        __m256d t3 = _mm256_set_m128d(  _mm256_extractf128_pd(t2, 0) + _mm256_extractf128_pd(t2, 1),
                                        _mm256_extractf128_pd(t1, 0) + _mm256_extractf128_pd(t1, 1));
        _mm256_storeu_pd(PRIM_PTR + n4, (_mm256_loadu_pd(PRIM_PTR + n4) + t3)*(*factor));
    }

    const int left = ncart % 4;
    const int done = ncart - left;
    for(n = done; n < ncart; n++)
    {
        union double4 tmp = (union double4)((*factor)*PRIM_INT[n]);
        PRIM_PTR[n] += tmp.d[0] + tmp.d[1] + tmp.d[2] + tmp.d[3];
    }

    #else

    // GCC is missing some instructions from the above block
    int offsets[4] = {0, 0, 0, 0};
    contract_fac(ncart, offsets, factor, PRIM_INT, PRIM_PTR);

    #endif
}


#ifdef __cplusplus
}
#endif

