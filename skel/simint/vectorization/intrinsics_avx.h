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


// Missing GCC vectorized exp and pow
static inline __m256d simint_exp_vec4(__m256d x)
{
    union double4 u = { x };
    union double4 res;
    for(int i = 0; i < 4; i++)
        res.d[i] = exp(u.d[i]);
    return res.v;
}

static inline __m256d simint_pow_vec4(__m256d a, __m256d p)
{
    union double4 ua = { a };
    union double4 up = { p };
    union double4 res;
    for(int i = 0; i < 4; i++)
        res.d[i] = pow(ua.d[i], up.d[i]);
    return res.v;
}


#if defined SIMINT_AVX || defined SIMINT_AVXFMA

    #define SIMINT_SIMD_LEN 4

    #define SIMINT_DBLTYPE         __m256d
    #define SIMINT_DBLLOAD(p,i)    _mm256_load_pd((p) + (i))
    #define SIMINT_DBLSET1(a)      _mm256_set1_pd((a))
    #define SIMINT_NEG(a)          (-(a))
    #define SIMINT_ADD(a,b)        _mm256_add_pd((a), (b))
    #define SIMINT_SUB(a,b)        _mm256_sub_pd((a), (b))
    #define SIMINT_MUL(a,b)        _mm256_mul_pd((a), (b))
    #define SIMINT_DIV(a,b)        _mm256_div_pd((a), (b))
    #define SIMINT_SQRT(a)         _mm256_sqrt_pd((a))

    #ifdef SIMINT_AVXFMA
      #define SIMINT_FMADD(a,b,c)  _mm256_fmadd_pd((a), (b), (c))
      #define SIMINT_FMSUB(a,b,c)  _mm256_fmsub_pd((a), (b), (c))
    #else
      #define SIMINT_FMADD(a,b,c)  SIMINT_ADD(SIMINT_MUL((a),(b)),(c))
      #define SIMINT_FMSUB(a,b,c)  SIMINT_SUB(SIMINT_MUL((a),(b)),(c))
    #endif

    #if defined __INTEL_COMPILER 
        #define SIMINT_EXP(a)       _mm256_exp_pd((a))
        #define SIMINT_POW(a,p)     _mm256_pow_pd((a), (p))
    #else
        #define SIMINT_EXP(a)       simint_exp_vec4((a))
        #define SIMINT_POW(a,p)     simint_pow_vec4((a), (p))
    #endif


    ////////////////////////////////////////
    // Special functions
    ////////////////////////////////////////

    static inline
    void contract(int ncart,
                  int const * restrict offsets,
                  __m256d const * restrict src,
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
                      __m256d const * restrict src,
                      double * restrict dest)
    {
        #if defined __clang__ || defined __INTEL_COMPILER

        int n, n4;
        const int nbatch = ncart/4;

        for(n = 0, n4 = 0; n < nbatch; n++, n4 += 4)
        {
            __m256d t1 = _mm256_hadd_pd(src[n4],   src[n4+1]);
            __m256d t2 = _mm256_hadd_pd(src[n4+2], src[n4+3]);
            __m256d t3 = _mm256_set_m128d(  _mm256_extractf128_pd(t2, 0) + _mm256_extractf128_pd(t2, 1),
                                            _mm256_extractf128_pd(t1, 0) + _mm256_extractf128_pd(t1, 1));
            _mm256_storeu_pd(dest + n4, _mm256_loadu_pd(dest + n4) + t3);
        }

        const int left = ncart % 4;
        const int done = ncart - left;
        for(n = done; n < ncart; n++)
        {
            union double4 tmp = { src[n] };
            dest[n] += tmp.d[0] + tmp.d[1] + tmp.d[2] + tmp.d[3];
        }

        #else

        // GCC is missing some instructions from the above block
        int offsets[4] = {0, 0, 0, 0};
        contract(ncart, offsets, src, dest);

        #endif
    }


    static inline
    void contract_fac(int ncart,
                      const __m256d factor,
                      int const * restrict offsets,
                      __m256d const * restrict src,
                      double * restrict dest)
    {
        for(int np = 0; np < ncart; ++np)
        {
            union double4 vtmp = { SIMINT_MUL(src[np], factor) };

            for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
                dest[offsets[n]*ncart+np] += vtmp.d[n]; 
        }
    }


    static inline
    void contract_all_fac(int ncart,
                          const __m256d factor,
                          __m256d const * restrict src,
                          double * restrict dest)
    {
        #if defined __clang__ || defined __INTEL_COMPILER

        int n, n4;
        const int nbatch = ncart/4;

        for(n = 0, n4 = 0; n < nbatch; n++, n4 += 4)
        {
            __m256d t1 = _mm256_hadd_pd(SIMINT_MUL(factor, src[n4]), SIMINT_MUL(factor, src[n4+1]));
            __m256d t2 = _mm256_hadd_pd(SIMINT_MUL(factor, src[n4+2]), SIMINT_MUL(factor, src[n4+3]));
            __m256d t3 = _mm256_set_m128d(  _mm256_extractf128_pd(t2, 0) + _mm256_extractf128_pd(t2, 1),
                                            _mm256_extractf128_pd(t1, 0) + _mm256_extractf128_pd(t1, 1));
            _mm256_storeu_pd(dest + n4, (_mm256_loadu_pd(dest + n4) + t3));
        }

        const int left = ncart % 4;
        const int done = ncart - left;
        for(n = done; n < ncart; n++)
        {
            __m256d vtmp = SIMINT_MUL(factor, src[n]);
            union double4 tmp = { vtmp };
            dest[n] += tmp.d[0] + tmp.d[1] + tmp.d[2] + tmp.d[3];
        }

        #else

        // GCC is missing some instructions from the above block
        int offsets[4] = {0, 0, 0, 0};
        contract_fac(ncart, factor, offsets, src, dest);

        #endif
    }


    static inline
    double vector_min(__m256d v)
    {
        union double4 m = { v };
        double min = m.d[0];
        for(int n = 1; n < SIMINT_SIMD_LEN; n++)
            min = (m.d[n] < min ? m.d[n] : min);
        return min;
    }


    static inline
    double vector_max(__m256d v)
    {
        union double4 m = { v };
        double max = m.d[0];
        for(int n = 1; n < SIMINT_SIMD_LEN; n++)
            max = (m.d[n] > max ? m.d[n] : max);
        return max;
    }

    static inline
    __m256d mask_load(int nlane, double * memaddr)
    {
        union double4 u = { _mm256_load_pd(memaddr) };
        for(int n = nlane; n < SIMINT_SIMD_LEN; n++)
            u.d[n] = 0.0;
        return u.v;
    }

#endif // defined SIMINT_AVX


#ifdef __cplusplus
}
#endif
