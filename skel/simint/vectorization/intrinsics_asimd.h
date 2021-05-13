#pragma once

#include <arm_neon.h>
#include <stdint.h>
#include <math.h>

#ifdef __cplusplus
#include "simint/cpp_restrict.hpp"
extern "C" {
#endif

union simint_double2
{
    float64x2_t v;
    double d[2];
};

// TODO: ARM compiler has its vectorized math library

static inline float64x2_t simint_exp_vec2(float64x2_t x)
{
    union simint_double2 u = { x };
    union simint_double2 res;
    for(int i = 0; i < 2; i++)
        res.d[i] = exp(u.d[i]);
    return res.v;
}

static inline float64x2_t simint_pow_vec2(float64x2_t a, float64x2_t p)
{
    union simint_double2 ua = { a };
    union simint_double2 up = { p };
    union simint_double2 res;
    for(int i = 0; i < 2; i++)
        res.d[i] = pow(ua.d[i], up.d[i]);
    return res.v;
}

#if defined SIMINT_ASIMD

    #define SIMINT_SIMD_LEN 2

    #define SIMINT_DBLTYPE         float64x2_t
    #define SIMINT_DBLLOAD(p,i)    vld1q_f64((p) + (i))
    #define SIMINT_DBLSET1(a)      vdupq_n_f64((a))
    #define SIMINT_NEG(a)          vnegq_f64((a))
    #define SIMINT_ADD(a,b)        vaddq_f64((a), (b))
    #define SIMINT_SUB(a,b)        vsubq_f64((a), (b))
    #define SIMINT_MUL(a,b)        vmulq_f64((a), (b))
    #define SIMINT_DIV(a,b)        vdivq_f64((a), (b))
    #define SIMINT_SQRT(a)         vsqrtq_f64((a))

    #define SIMINT_FMADD(a,b,c)    vfmaq_f64((c), (a), (b))
    #define SIMINT_FMSUB(a,b,c)    vnegq_f64(vfmsq_f64((c), (a), (b)))

    #define SIMINT_EXP(a)          simint_exp_vec2((a))
    #define SIMINT_POW(a,p)        simint_pow_vec2((a), (p))


    ////////////////////////////////////////
    // Special functions
    ////////////////////////////////////////

    static inline
    void contract(int ncart,
                  int const * restrict offsets,
                  float64x2_t const * restrict src,
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
                      float64x2_t const * restrict src,
                      double * restrict dest)
    {
        for (int np = 0; np < ncart; np++)
            dest[np] += vaddvq_f64(src[np]);
    }


    static inline
    void contract_fac(int ncart,
                      const float64x2_t factor,
                      int const * restrict offsets,
                      float64x2_t const * restrict src,
                      double * restrict dest)
    {
        for(int np = 0; np < ncart; ++np)
        {
            union simint_double2 vtmp = { SIMINT_MUL(src[np], factor) };

            for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
                dest[offsets[n]*ncart+np] += vtmp.d[n]; 
        }
    }


    static inline
    void contract_all_fac(int ncart,
                          const float64x2_t factor,
                          float64x2_t const * restrict src,
                          double * restrict dest)
    {
        for(int np = 0; np < ncart; np++)
            dest[np] += vaddvq_f64(SIMINT_MUL(factor, src[np]));
    }


    static inline
    double vector_min(float64x2_t v)
    {
        union simint_double2 m = { v };
        double min = m.d[0];
        for(int n = 1; n < SIMINT_SIMD_LEN; n++)
            min = (m.d[n] < min ? m.d[n] : min);
        return min;
    }


    static inline
    double vector_max(float64x2_t v)
    {
        union simint_double2 m = { v };
        double max = m.d[0];
        for(int n = 1; n < SIMINT_SIMD_LEN; n++)
            max = (m.d[n] > max ? m.d[n] : max);
        return max;
    }

    static inline
    float64x2_t mask_load(int nlane, double * memaddr)
    {
        union simint_double2 u = { vld1q_f64(memaddr) };
        for(int n = nlane; n < SIMINT_SIMD_LEN; n++)
            u.d[n] = 0.0;
        return u.v;
    }
    
    //#define SIMINT_PRIM_SCREEN_STAT
    static inline
    int count_prim_screen_survival(float64x2_t screen_val, const double screen_tol)
    {
        union simint_double2 u = {screen_val};
        int res = 0;
        for (int i = 0; i < SIMINT_SIMD_LEN; i++)
            if (u.d[i] >= screen_tol) res++;
        return res;
    }

#endif // defined SIMINT_ASIMD


#ifdef __cplusplus
}
#endif
