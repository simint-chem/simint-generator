#pragma once

#include <immintrin.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

union double2
{
    __m128d v;
    double d[2];
};


// Missing GCC vectorized exp
static inline __m128d simint_exp_vec2(__m128d x)
{
    union double2 u = { x };
    union double2 res;
    for(int i = 0; i < 2; i++)
        res.d[i] = exp(u.d[i]);
    return res.v;
}

#if defined SIMINT_SSE

    #define SIMINT_SIMD_LEN 2

    #define SIMINT_DBLTYPE         __m128d
    #define SIMINT_UNIONTYPE       union double2
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



    ////////////////////////////////////////
    // Special functions
    ////////////////////////////////////////
    static inline
    void contract(int ncart,
                  int const * restrict offsets,
                  __m128d const * restrict PRIM_INT,
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
                      __m128d const * restrict PRIM_INT,
                      double * restrict PRIM_PTR)
    {
        //TODO efficient algorithm for this
        int offsets[2] = {0, 0};
        contract(ncart, offsets, PRIM_INT, PRIM_PTR);
    }


    static inline
    double vector_max(__m128d v)
    {
        SIMINT_UNIONTYPE m = { v };
        double max = 0.0;
        for(int n = 0; n < SIMINT_SIMD_LEN; n++)
            max = (m.d[n] > max ? m.d[n] : max);
        return max;
    }

#endif // defined SIMINT_SSE

#ifdef __cplusplus
}
#endif

