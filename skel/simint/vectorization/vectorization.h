#pragma once

// Automatically generated. Includes SIMINT_AVX, SIMINT_SSE, etc
#include "simint/simint_config.h"

#if defined SIMINT_AVX
  #include "simint/vectorization/intrinsics_avx.h"
  #define SIMINT_SIMD_LEN 4
#elif defined SIMINT_SSE
  #include "simint/vectorization/intrinsics_sse.h"
  #define SIMINT_SIMD_LEN 2
#else
  #include "simint/vectorization/intrinsics_scalar.h"
  #define SIMINT_SIMD_LEN 1
#endif

#define SIMINT_SIMD_ALIGN_DBL (SIMINT_SIMD_LEN*8)
#define SIMINT_SIMD_ALIGN_INT (SIMINT_SIMD_LEN*sizeof(int))

// "Max" alignment. Can be used to allocate memory of mixed types
// (ie ints and doubles can be allocated this way)
#define SIMINT_SIMD_ALIGN SIMINT_SIMD_ALIGN_DBL


#define ASSUME_ALIGN_DBL(x)  __assume_aligned((x), SIMINT_SIMD_ALIGN_DBL)
#define ASSUME_ALIGN_INT(x)  __assume_aligned((x), SIMINT_SIMD_ALIGN_INT)
#define ASSUME __assume


// Aligned memory allocation
#define SIMINT_ALLOC(x) _mm_malloc((x), SIMINT_SIMD_ALIGN)
#define SIMINT_FREE(x) _mm_free((x))


// round up the number of elements to the nearest boundary
#define SIMINT_SIMD_ROUND(x) ((x + ((SIMINT_SIMD_LEN-1))) & (~(SIMINT_SIMD_LEN-1)))


// align an array
// ie double somearr[200] SIMINT_ALIGN_ARRAY
#define SIMINT_ALIGN_ARRAY_DBL __attribute__((aligned(SIMINT_SIMD_ALIGN_DBL)))
#define SIMINT_ALIGN_ARRAY_INT __attribute__((aligned(SIMINT_SIMD_ALIGN_INT)))


// Number of shells to use in a batch
#define SIMINT_NSHELL_SIMD (2*SIMINT_SIMD_LEN)


#if defined SIMINT_AVX
    #define SIMINT_DBLTYPE  __m256d
    #define SIMINT_SET1     _mm256_set1_pd

    #ifdef SIMINT_INTEL
        #define SIMINT_EXP  _mm256_exp_pd
    #else
        #define SIMINT_EXP  simint_exp_vec4
    #endif

#elif defined SIMINT_SSE
    #define SIMINT_DBLTYPE  __m128d
    #define SIMINT_SET1     _mm_set1_pd

    #ifdef SIMINT_INTEL
        #define SIMINT_EXP  _mm_exp_pd
    #else
        #define SIMINT_EXP  simint_exp_vec2
    #endif

#else
    #define SIMINT_DBLTYPE  double
    #define SIMINT_SET1
    #define SIMINT_EXP    exp

#endif
