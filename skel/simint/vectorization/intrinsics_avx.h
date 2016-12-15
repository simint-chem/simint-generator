#pragma once

#include <immintrin.h>
#include <stdint.h>
#include <math.h>

#include "simint/vectorization/intrinsics_sse.h"

union double4
{
    __m256d d_256;
    __m128d d_128[2];
    double d[4];
};


// Missing GCC vectorized exp
static inline __m256d simint_exp_vec4(__m256d x)
{
    union double4 u = { x };
    union double4 res;
    for(int i = 0; i < 4; i++)
        res.d[i] = exp(u.d[i]);
    return res.d_256;
}


#if defined SIMINT_AVX

    #define SIMINT_SIMD_LEN 4

    #define SIMINT_DBLTYPE         __m256d
    #define SIMINT_UNIONTYPE       union double4
    #define SIMINT_UNIONMEMBER(a)  (a.d_256)
    #define SIMINT_DBLLOAD(p,i)    _mm256_load_pd((p) + (i))
    #define SIMINT_DBLSET1(a)      _mm256_set1_pd((a))
    #define SIMINT_NEG(a)          (-(a))
    #define SIMINT_ADD(a,b)        _mm256_add_pd((a), (b))
    #define SIMINT_SUB(a,b)        _mm256_sub_pd((a), (b))
    #define SIMINT_MUL(a,b)        _mm256_mul_pd((a), (b))
    #define SIMINT_DIV(a,b)        _mm256_div_pd((a), (b))
    #define SIMINT_SQRT(a)         _mm256_sqrt_pd((a))

    #ifdef SIMINT_FMA
      #define SIMINT_FMADD(a,b,c)  _mm256_fmadd_pd((a), (b), (c))
      #define SIMINT_FMSUB(a,b,c)  _mm256_fmsub_pd((a), (b), (c))
    #else
      #define SIMINT_FMADD(a,b,c)  SIMINT_ADD(SIMINT_MUL((a),(b)),(c))
      #define SIMINT_FMSUB(a,b,c)  SIMINT_SUB(SIMINT_MUL((a),(b)),(c))
    #endif

    #ifdef SIMINT_INTEL
        #define SIMINT_EXP  _mm256_exp_pd
    #else
        #define SIMINT_EXP  simint_exp_vec4
    #endif

#endif

