#pragma once

#include <arm_sve.h>
#include <stdint.h>
#include <math.h>

#ifdef __cplusplus
#include "simint/cpp_restrict.hpp"
extern "C" {
#endif

#define PTRUE32B    svptrue_b32()
#define PTRUE64B    svptrue_b64()

#if __ARM_FEATURE_SVE_BITS==0
#warning __ARM_FEATURE_SVE_BITS==0 detected, assuming SVE length == 512 bits
#undef   __ARM_FEATURE_SVE_BITS
#define  __ARM_FEATURE_SVE_BITS 512
#endif

#if __ARM_FEATURE_SVE_BITS==128
#define SIMD_LEN_S  4
#define SIMD_LEN_D  2
#define USE_SVE128
#endif 

#if __ARM_FEATURE_SVE_BITS==256
#define SIMD_LEN_S  8
#define SIMD_LEN_D  4
#define USE_SVE256
#endif 

#if __ARM_FEATURE_SVE_BITS==512
#define SIMD_LEN_S  16
#define SIMD_LEN_D  8
#define USE_SVE512
#endif 

#if __ARM_FEATURE_SVE_BITS==1024
#define SIMD_LEN_S  32
#define SIMD_LEN_D  16
#define USE_SVE1024
#endif 

#if __ARM_FEATURE_SVE_BITS==2048
#define SIMD_LEN_S  64
#define SIMD_LEN_D  32
#define USE_SVE2048
#endif 

typedef svfloat32_t vec_s __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
typedef svfloat64_t vec_d __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));

union simint_double_sve
{
    vec_d v;
    double d[SIMD_LEN_D];
};

// TODO: ARM compiler has its vectorized math library

static inline vec_d simint_exp_vec(vec_d x)
{
    union simint_double_sve u = { x };
    union simint_double_sve res;
    for(int i = 0; i < SIMD_LEN_D; i++)
        res.d[i] = exp(u.d[i]);
    return res.v;
}

static inline vec_d simint_pow_vec(vec_d a, vec_d p)
{
    union simint_double_sve ua = { a };
    union simint_double_sve up = { p };
    union simint_double_sve res;
    for(int i = 0; i < SIMD_LEN_D; i++)
        res.d[i] = pow(ua.d[i], up.d[i]);
    return res.v;
}

#if defined SIMINT_SVE

    #define SIMINT_SIMD_LEN SIMD_LEN_D

    #define SIMINT_DBLTYPE         vec_d
    #define SIMINT_DBLLOAD(p,i)    svld1_f64(PTRUE64B, (p) + (i))
    #define SIMINT_DBLSET1(a)      svdup_f64_z(PTRUE64B, (a))
    #define SIMINT_NEG(a)          svneg_f64_z(PTRUE64B, (a))
    #define SIMINT_ADD(a,b)        svadd_f64_z(PTRUE64B, (a), (b))
    #define SIMINT_SUB(a,b)        svsub_f64_z(PTRUE64B, (a), (b))
    #define SIMINT_MUL(a,b)        svmul_f64_z(PTRUE64B, (a), (b))
    #define SIMINT_DIV(a,b)        svdiv_f64_z(PTRUE64B, (a), (b))
    #define SIMINT_SQRT(a)         svsqrt_f64_z(PTRUE64B, (a))

    #define SIMINT_FMADD(a,b,c)    svmad_f64_z(PTRUE64B, (a), (b), (c))
    #define SIMINT_FMSUB(a,b,c)    svnmsb_f64_z(PTRUE64B, (a), (b), (c))

    #define SIMINT_EXP(a)          simint_exp_vec((a))
    #define SIMINT_POW(a,p)        simint_pow_vec((a), (p))


    ////////////////////////////////////////
    // Special functions
    ////////////////////////////////////////

    static inline
    void contract(int ncart,
                  int const * restrict offsets,
                  vec_d const * restrict src,
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
                      vec_d const * restrict src,
                      double * restrict dest)
    {
        for (int np = 0; np < ncart; np++)
            dest[np] += svaddv_f64(PTRUE64B, src[np]);
    }


    static inline
    void contract_fac(int ncart,
                      const vec_d factor,
                      int const * restrict offsets,
                      vec_d const * restrict src,
                      double * restrict dest)
    {
        for(int np = 0; np < ncart; ++np)
        {
            vec_d vtmp0 = SIMINT_MUL(src[np], factor);
            union simint_double_sve vtmp = { vtmp0 };

            for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
                dest[offsets[n]*ncart+np] += vtmp.d[n]; 
        }
    }


    static inline
    void contract_all_fac(int ncart,
                          const vec_d factor,
                          vec_d const * restrict src,
                          double * restrict dest)
    {
        for(int np = 0; np < ncart; np++)
            dest[np] += svaddv_f64(PTRUE64B, SIMINT_MUL(factor, src[np]));
    }


    static inline
    double vector_min(vec_d v)
    {
        return svminv_f64(PTRUE64B, v);
    }


    static inline
    double vector_max(vec_d v)
    {
        return svmaxv_f64(PTRUE64B, v);
    }

    static inline
    vec_d mask_load(int nlane, double * memaddr)
    {
        svbool_t pred = svwhilelt_b64(0, nlane);
        vec_d ret = svld1_f64(pred, memaddr);
        return ret;
    }
    
    //#define SIMINT_PRIM_SCREEN_STAT
    static inline
    int count_prim_screen_survival(vec_d screen_val, const double screen_tol)
    {
        union simint_double_sve u = {screen_val};
        int res = 0;
        for (int i = 0; i < SIMINT_SIMD_LEN; i++)
            if (u.d[i] >= screen_tol) res++;
        return res;
    }

#endif // defined SIMINT_SVE


#ifdef __cplusplus
}
#endif
