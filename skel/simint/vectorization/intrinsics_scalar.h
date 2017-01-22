#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#if defined SIMINT_SCALAR

    #define SIMINT_SIMD_LEN 1

    #define SIMINT_DBLTYPE         double
    //#define SIMINT_UNIONTYPE       double
    #define SIMINT_DBLLOAD(p,i)    ((p)[(i)])
    #define SIMINT_DBLSET1(a)      ((a)) 
    #define SIMINT_NEG(a)          (-(a))
    #define SIMINT_ADD(a,b)        ((a)+(b))
    #define SIMINT_SUB(a,b)        ((a)-(b))
    #define SIMINT_MUL(a,b)        ((a)*(b))
    #define SIMINT_DIV(a,b)        ((a)/(b))
    #define SIMINT_SQRT(a)         sqrt((a))
    #define SIMINT_POW(a, p)       pow((a), (p))
    #define SIMINT_FMADD(a,b,c)    SIMINT_ADD(SIMINT_MUL((a),(b)),(c))
    #define SIMINT_FMSUB(a,b,c)    SIMINT_SUB(SIMINT_MUL((a),(b)),(c))
    #define SIMINT_EXP(a)          exp((a))


    ////////////////////////////////////////
    // Special functions
    ////////////////////////////////////////

    static inline
    void contract(int ncart,
                  double const * restrict PRIM_INT,
                  double * restrict PRIM_PTR)
    {
        for(int np = 0; np < ncart; ++np)
            PRIM_PTR[np] += PRIM_INT[np]; 
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


#endif // defined SIMINT_SCALAR

#ifdef __cplusplus
}
#endif
