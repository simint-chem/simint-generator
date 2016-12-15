#pragma once

#include <immintrin.h>
#include <math.h>

union double2
{
    __m128d d_128;
    double d[2];
};


// Missing GCC vectorized exp
static inline __m128d simint_exp_vec2(__m128d x)
{
    union double2 u = { x };
    union double2 res;
    for(int i = 0; i < 2; i++)
        res.d[i] = exp(u.d[i]);
    return res.d_128;
}

#if defined SIMINT_SSE

    #define SIMINT_SIMD_LEN 2

    #define SIMINT_DBLTYPE         __m128d
    #define SIMINT_UNIONTYPE       union double2
    #define SIMINT_UNIONMEMBER(a)  (a.d_128)
    #define SIMINT_DBLLOAD(p,i)    _mm_load_pd((p) + (i))
    #define SIMINT_DBLSET1(a)      _mm_set1_pd((a))
    #define SIMINT_NEG(a)          (-(a))
    #define SIMINT_ADD(a,b)        _mm_add_pd((a), (b))
    #define SIMINT_SUB(a,b)        _mm_sub_pd((a), (b))
    #define SIMINT_MUL(a,b)        _mm_mul_pd((a), (b))
    #define SIMINT_DIV(a,b)        _mm_div_pd((a), (b))
    #define SIMINT_SQRT(a)         _mm_sqrt_pd((a))

    #ifdef SIMINT_FMA
      #define SIMINT_FMADD(a,b,c)  _mm_fmadd_pd((a), (b), (c))
      #define SIMINT_FMSUB(a,b,c)  _mm_fmsub_pd((a), (b), (c))
    #else
      #define SIMINT_FMADD(a,b,c)  SIMINT_ADD(SIMINT_MUL((a),(b)),(c))
      #define SIMINT_FMSUB(a,b,c)  SIMINT_SUB(SIMINT_MUL((a),(b)),(c))
    #endif

    #ifdef SIMINT_INTEL
        #define SIMINT_EXP  _mm_exp_pd
    #else
        #define SIMINT_EXP  simint_exp_vec2
    #endif

#endif

