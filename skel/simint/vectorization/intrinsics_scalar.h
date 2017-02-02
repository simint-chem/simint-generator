#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#if defined SIMINT_SCALAR

    #define SIMINT_SIMD_LEN 1

    #define SIMINT_DBLTYPE         double
    #define SIMINT_DBLLOAD(p,i)    ((p)[(i)])
    #define SIMINT_DBLSET1(a)      ((a)) 
    #define SIMINT_NEG(a)          (-(a))
    #define SIMINT_ADD(a,b)        ((a)+(b))
    #define SIMINT_SUB(a,b)        ((a)-(b))
    #define SIMINT_MUL(a,b)        ((a)*(b))
    #define SIMINT_DIV(a,b)        ((a)/(b))
    #define SIMINT_SQRT(a)         sqrt((a))
    #define SIMINT_FMADD(a,b,c)    SIMINT_ADD(SIMINT_MUL((a),(b)),(c))
    #define SIMINT_FMSUB(a,b,c)    SIMINT_SUB(SIMINT_MUL((a),(b)),(c))
    #define SIMINT_EXP(a)          exp((a))
    #define SIMINT_POW(a, p)       pow((a), (p))


    ////////////////////////////////////////
    // Special functions
    ////////////////////////////////////////

    static inline
    void contract(int ncart,
                  int const * restrict offsets,
                  double const * restrict src,
                  double * restrict dest)
    {
        for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
        {
            double const * restrict src_tmp = src + n;
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
                      double const * restrict src,
                      double * restrict dest)
    {
        int offsets[1] = {0};
        contract(ncart, offsets, src, dest);
    }

    static inline
    void contract_fac(int ncart,
                      double factor,
                      int const * restrict offsets,
                      double const * restrict src,
                      double * restrict dest)
    {
        for(int n = 0; n < SIMINT_SIMD_LEN; ++n)
        {
            double const * restrict src_tmp = src + n;
            double * restrict dest_tmp = dest + offsets[n]*ncart;

            for(int np = 0; np < ncart; ++np)
            {
                dest_tmp[np] += factor*(*src_tmp);
                src_tmp += SIMINT_SIMD_LEN;
            }
        }
    }

    static inline
    void contract_all_fac(int ncart,
                          double factor,
                          double const * restrict src,
                          double * restrict dest)
    {
        int offsets[1] = {0};
        contract_fac(ncart, factor, offsets, src, dest);
    }

    static inline
    double vector_min(double v)
    {
        return v;
    }


    static inline
    double vector_max(double v)
    {
        return v;
    }


    static inline
    double mask_load(int nlane, double * memaddr)
    {
        if(nlane == 1)
            return *memaddr;
        else
            return 0.0;
    }



#endif // defined SIMINT_SCALAR

#ifdef __cplusplus
}
#endif
