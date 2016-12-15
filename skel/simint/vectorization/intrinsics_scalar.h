#pragma once

#if defined SIMINT_SCALAR

    #define SIMINT_SIMD_LEN 1

    #define SIMINT_DBLTYPE         double
    //#define SIMINT_UNIONTYPE       double
    //#define SIMINT_UNIONMEMBER(a)  (a)
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

#endif

