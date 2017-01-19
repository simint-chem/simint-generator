#pragma once

#include <immintrin.h>
#include <stdint.h>
#include <math.h>

#include "simint/vectorization/intrinsics_avx.h"

#ifdef __cplusplus
extern "C" {
#endif

union double8
{
    __m512d v;
    double d[8];
};


// Missing GCC vectorized exp
static inline __m512d simint_exp_vec8(__m512d x)
{
    union double8 u = { x };
    union double8 res;
    for(int i = 0; i < 8; i++)
        res.d[i] = exp(u.d[i]);
    return res.v;
}

#if defined SIMINT_AVX512 || defined SIMINT_MICAVX512

    #define SIMINT_SIMD_LEN 8

    #define SIMINT_DBLTYPE         __m512d
    #define SIMINT_UNIONTYPE       union double8
    #define SIMINT_DBLLOAD(p,i)    _mm512_load_pd((p) + (i))
    #define SIMINT_DBLSET1(a)      _mm512_set1_pd((a))
    #define SIMINT_NEG(a)          (SIMINT_MUL((a), (SIMINT_DBLSET1(-1.0)))) 
    #define SIMINT_ADD(a,b)        _mm512_add_pd((a), (b))
    #define SIMINT_SUB(a,b)        _mm512_sub_pd((a), (b))
    #define SIMINT_MUL(a,b)        _mm512_mul_pd((a), (b))
    #define SIMINT_DIV(a,b)        _mm512_div_pd((a), (b))
    #define SIMINT_SQRT(a)         _mm512_sqrt_pd((a))
    #define SIMINT_FMADD(a,b,c)    _mm512_fmadd_pd((a), (b), (c))
    #define SIMINT_FMSUB(a,b,c)    _mm512_fmsub_pd((a), (b), (c))

    #ifdef SIMINT_INTEL
        #define SIMINT_EXP  _mm512_exp_pd
    #else
        #define SIMINT_EXP  simint_exp_vec8
    #endif



    ////////////////////////////////////////
    // Special functions
    ////////////////////////////////////////

    static inline
    void contract(int ncart,
                  int const * restrict offsets,
                  __m512d const * restrict PRIM_INT,
                  double * restrict PRIM_PTR)
    {
        for(int np = 0; np < ncart; ++np)
        {
            SIMINT_UNIONTYPE vtmp = { PRIM_INT[np] };

            for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
                PRIM_PTR[offsets[n]*ncart] += vtmp.d[n]; 
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


    static inline
    double vector_max(__m512d v)
    {
        return _mm512_reduce_max_pd(v);
    }

#endif // defined SIMINT_AVX512 || defined SIMINT_MICAVX512

#ifdef __cplusplus
}
#endif

