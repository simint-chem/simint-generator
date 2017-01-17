#pragma once

#include <immintrin.h>
#include <stdint.h>
#include <math.h>

#include "simint/vectorization/intrinsics_avx.h"

union double8
{
    __m512d d_512;
    __m256d d_256[2];
    __m128d d_128[4];
    double d[8];
};


// Missing GCC vectorized exp
static inline __m512d simint_exp_vec8(__m512d x)
{
    union double8 u = { x };
    union double8 res;
    for(int i = 0; i < 8; i++)
        res.d[i] = exp(u.d[i]);
    return res.d_512;
}

#if defined SIMINT_AVX512 || defined SIMINT_MICAVX512

    #define SIMINT_SIMD_LEN 8

    #define SIMINT_DBLTYPE         __m512d
    #define SIMINT_UNIONTYPE       union double8
    #define SIMINT_UNIONMEMBER(a)  (a.d_512)
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

#endif

