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


// Missing GCC vectorized exp and pow
static inline __m128d simint_exp_vec2(__m128d x)
{
    union double2 u = { x };
    union double2 res;
    for(int i = 0; i < 2; i++)
        res.d[i] = exp(u.d[i]);
    return res.v;
}

static inline __m128d simint_pow_vec2(__m128d a, __m128d p)
{
    union double2 ua = { a };
    union double2 up = { p };
    union double2 res;
    for(int i = 0; i < 2; i++)
        res.d[i] = pow(ua.d[i], up.d[i]);
    return res.v;
}

#if defined SIMINT_SSE

    #define SIMINT_SIMD_LEN 2

    #define SIMINT_DBLTYPE         __m128d
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

    #if defined __INTEL_COMPILER 
        #define SIMINT_EXP(a)       _mm_exp_pd((a))
        #define SIMINT_POW(a,p)     _mm_pow_pd((a), (p))
    #else
        #define SIMINT_EXP(a)       simint_exp_vec2((a))
        #define SIMINT_POW(a,p)     simint_pow_vec2((a), (p))
    #endif



    ////////////////////////////////////////
    // Special functions
    ////////////////////////////////////////
    static inline
    void contract(int ncart,
                  int const * restrict offsets,
                  __m128d const * restrict src,
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
                      __m128d const * restrict src,
                      double * restrict dest)
    {
        //TODO efficient algorithm for this
        int offsets[2] = {0, 0};
        contract(ncart, offsets, src, dest);
    }

    static inline
    void contract_fac(int ncart,
                      __m128d factor,
                      int const * restrict offsets,
                      __m128d const * restrict src,
                      double * restrict dest)
    {
        for(int np = 0; np < ncart; ++np)
        {
            union double2 vtmp = { SIMINT_MUL(src[np], factor) };

            for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
                dest[offsets[n]*ncart+np] += vtmp.d[n]; 
        }
    }


    static inline
    void contract_all_fac(int ncart,
                          __m128d factor,
                          __m128d const * restrict src,
                          double * restrict dest)
    {
        //TODO efficient algorithm for this
        int offsets[2] = {0, 0};
        contract(ncart, offsets, src, dest);
    }


    static inline
    double vector_min(__m128d v)
    {
        union double2 m = { v };
        double min = m.d[0];
        for(int n = 1; n < SIMINT_SIMD_LEN; n++)
            min = (m.d[n] < min ? m.d[n] : min);
        return min;
    }


    static inline
    double vector_max(__m128d v)
    {
        union double2 m = { v };
        double max = 0.0;
        for(int n = 0; n < SIMINT_SIMD_LEN; n++)
            max = (m.d[n] > max ? m.d[n] : max);
        return max;
    }


    static inline
    __m128d mask_load(int nlane, double * memaddr)
    {
        union double2 u = { _mm_load_pd(memaddr) };
        for(int n = nlane; n < SIMINT_SIMD_LEN; n++)
            u.d[n] = 0.0;
        return u.v;
    }

#endif // defined SIMINT_SSE

#ifdef __cplusplus
}
#endif

