#pragma once

#include <immintrin.h>
#include <stdint.h>
#include <math.h>

#include "simint/vectorization/intrinsics_sse.h"


#ifdef __cplusplus
extern "C" {
#endif

union double4
{
    __m256d v;
    double d[4];
};


// Missing GCC vectorized exp
static inline __m256d simint_exp_vec4(__m256d x)
{
    union double4 u = { x };
    union double4 res;
    for(int i = 0; i < 4; i++)
        res.d[i] = exp(u.d[i]);
    return res.v;
}


#if defined SIMINT_AVX

    #define SIMINT_SIMD_LEN 4

    #define SIMINT_DBLTYPE         __m256d
    #define SIMINT_UNIONTYPE       union double4
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


    ////////////////////////////////////////
    // Special functions
    ////////////////////////////////////////

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
        for(int np = 0; np < ncart; ++np)
        {
            SIMINT_UNIONTYPE vtmp = { SIMINT_MUL(PRIM_INT[np], *factor) };

            for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
                PRIM_PTR[offsets[n]*ncart] += vtmp.d[n]; 
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
            __m256d t1 = _mm256_hadd_pd(SIMINT_MUL((*factor), PRIM_INT[n4]), SIMINT_MUL((*factor), PRIM_INT[n4+1]));
            __m256d t2 = _mm256_hadd_pd(SIMINT_MUL((*factor), PRIM_INT[n4+2]), SIMINT_MUL((*factor), PRIM_INT[n4+3]));
            __m256d t3 = _mm256_set_m128d(  _mm256_extractf128_pd(t2, 0) + _mm256_extractf128_pd(t2, 1),
                                            _mm256_extractf128_pd(t1, 0) + _mm256_extractf128_pd(t1, 1));
            _mm256_storeu_pd(PRIM_PTR + n4, (_mm256_loadu_pd(PRIM_PTR + n4) + t3));
        }

        const int left = ncart % 4;
        const int done = ncart - left;
        for(n = done; n < ncart; n++)
        {
            __m256d vtmp = SIMINT_MUL((*factor), PRIM_INT[n]);
            union double4 tmp = { vtmp };
            PRIM_PTR[n] += tmp.d[0] + tmp.d[1] + tmp.d[2] + tmp.d[3];
        }

        #else

        // GCC is missing some instructions from the above block
        int offsets[4] = {0, 0, 0, 0};
        contract_fac(ncart, offsets, factor, PRIM_INT, PRIM_PTR);

        #endif
    }


    static inline
    double vector_max(__m256d v)
    {
        SIMINT_UNIONTYPE m = { v };
        double max = 0.0;
        for(int n = 0; n < SIMINT_SIMD_LEN; n++)
            max = (m.d[n] > max ? m.d[n] : max);
        return max;
    }

#endif // defined SIMINT_AVX


#ifdef __cplusplus
}
#endif
