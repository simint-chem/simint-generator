#ifndef BOYS_H
#define BOYS_H

#include "vectorization.h"


#define BOYS_INTERP_ORDER 7
#define BOYS_MAX_DFAC 50


void Boys_Init(void);
void Boys_Finalize(void);
void Boys_F(double * F, int n, double x);



#ifdef SIMINT_SIMD
  #pragma omp declare simd simdlen(SIMD_LEN)
#endif
void Boys_F_short(double * F, int n, double x);



#ifdef SIMINT_SIMD
  #pragma omp declare simd simdlen(SIMD_LEN)
#endif
void Boys_F_long(double * F, int n, double x);



double Boys_Max(const int ncenter,
                const double * X, const double * Y, const double * Z,
                const int * n_prim_per_center, const double * alpha);


#endif
