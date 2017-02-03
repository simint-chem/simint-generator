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


// Missing GCC vectorized exp and pow
static inline __m512d simint_exp_vec8(__m512d x)
{
    union double8 u = { x };
    union double8 res;
    for(int i = 0; i < 8; i++)
        res.d[i] = exp(u.d[i]);
    return res.v;
}

static inline __m512d simint_pow_vec8(__m512d a, __m512d p)
{
    union double8 ua = { a };
    union double8 up = { p };
    union double8 res;
    for(int i = 0; i < 8; i++)
        res.d[i] = pow(ua.d[i], up.d[i]);
    return res.v;
}

#if defined SIMINT_AVX512 || defined SIMINT_MICAVX512

    #define SIMINT_SIMD_LEN 8

    #define SIMINT_DBLTYPE         __m512d
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

    #if defined __INTEL_COMPILER 
        #define SIMINT_EXP(a)       _mm512_exp_pd((a))
        #define SIMINT_POW(a,p)     _mm512_pow_pd((a), (p))
    #else
        #define SIMINT_EXP(a)       simint_exp_vec8((a))
        #define SIMINT_POW(a,p)     simint_pow_vec8((a), (p))
    #endif



    ////////////////////////////////////////
    // Special functions
    ////////////////////////////////////////

    static inline
    void contract(int ncart,
                  int const * restrict offsets,
                  __m512d const * restrict src,
                  double * restrict dest)
    {
        for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
        {
            double const * restrict src_tmp = (double *)src + n;
            double * restrict dest_tmp = dest + offsets[n]*ncart;

            for(int np = 0; np < ncart; ++np)
            {
                dest_tmp[np] += *src_tmp;
                src_tmp += SIMINT_SIMD_LEN;
            }
        }
    }


    static inline
    void contract_all(int ncart,
                      __m512d const * restrict src,
                      double * restrict dest)
    {
        #if defined __clang__ || defined __INTEL_COMPILER

        for(int np = 0; np < ncart; np++)
            dest[np] += _mm512_reduce_add_pd(src[np]);

        #else

        int offsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        contract(ncart, offsets, src, dest);

        #endif
    }


    static inline
    void contract_fac(int ncart,
                      const __m512d factor,
                      int const * restrict offsets,
                      __m512d const * restrict src,
                      double * restrict dest)
    {
        for(int np = 0; np < ncart; ++np)
        {
            union double8 vtmp = { SIMINT_MUL(src[np], factor) };

            for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
                dest[offsets[n]*ncart+np] += vtmp.d[n]; 
        }
    }


    static inline
    void contract_all_fac(int ncart,
                          const __m512d factor,
                          __m512d const * restrict src,
                          double * restrict dest)
    {
        #if defined __clang__ || defined __INTEL_COMPILER

        for(int np = 0; np < ncart; np++)
            dest[np] += _mm512_reduce_add_pd(_mm512_mul_pd(factor, src[np]));

        #else

        int offsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        contract_fac(ncart, offsets, factor, src, dest);

        #endif
    }


    static inline
    double vector_min(__m512d v)
    {
        #if defined __clang__ || defined __INTEL_COMPILER
            return _mm512_reduce_min_pd(v);
        #else
            union double8 u = { v };
            double min = u.d[0];
            for(int i = 1; i < 8; i++)
                min = (u.d[i] < min ? u.d[i] : min);
            return min;
        #endif
    }


    static inline
    double vector_max(__m512d v)
    {
        #if defined __clang__ || defined __INTEL_COMPILER
            return _mm512_reduce_max_pd(v);
        #else
            union double8 u = { v };
            double max = u.d[0];
            for(int i = 1; i < 8; i++)
                max = (u.d[i] > max ? u.d[i] : max);
            return max;
        #endif
    }

    static inline
    __m512d mask_load(int nlane, double * memaddr)
    {
        union double8 u = { _mm512_load_pd(memaddr) };
        for(int n = nlane; n < SIMINT_SIMD_LEN; n++)
            u.d[n] = 0.0;
        return u.v;
    }

#endif // defined SIMINT_AVX512 || defined SIMINT_MICAVX512

#ifdef __cplusplus
}
#endif

