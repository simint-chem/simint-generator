#pragma once

#include <stdlib.h>

// Defines SIMINT_AVX, etc
#include "simint/vectorization/vector_config.h"

#if defined SIMINT_AVX512 || defined SIMINT_MICAVX512
  #include "simint/vectorization/intrinsics_avx512.h"
#elif defined SIMINT_AVX || defined SIMINT_AVXFMA
  #include "simint/vectorization/intrinsics_avx.h"
#elif defined SIMINT_SSE
  #include "simint/vectorization/intrinsics_sse.h"
#elif defined SIMINT_SCALAR
  #include "simint/vectorization/intrinsics_scalar.h"
#else
  #error Vector type is not set
#endif

#define SIMINT_SIMD_ALIGN_DBL (SIMINT_SIMD_LEN*8)
#define SIMINT_SIMD_ALIGN_INT (SIMINT_SIMD_LEN*sizeof(int))

// "Max" alignment. Can be used to allocate memory of mixed types
// (ie ints and doubles can be allocated this way)
#define SIMINT_SIMD_ALIGN SIMINT_SIMD_ALIGN_DBL


#if defined __INTEL_COMPILER
    #define SIMINT_ASSUME_ALIGN_DBL(x)  __assume_aligned((x), SIMINT_SIMD_ALIGN_DBL)
    #define SIMINT_ASSUME_ALIGN_INT(x)  __assume_aligned((x), SIMINT_SIMD_ALIGN_INT)
#elif defined __clang__
    // TODO - find these
    #define SIMINT_ASSUME_ALIGN_DBL(x)
    #define SIMINT_ASSUME_ALIGN_INT(x)
#elif defined __GNUC__
    // TODO - find these
    #define SIMINT_ASSUME_ALIGN_DBL(x)
    #define SIMINT_ASSUME_ALIGN_INT(x)
#else
    #define SIMINT_ASSUME_ALIGN_DBL(x)
    #define SIMINT_ASSUME_ALIGN_INT(x)
#endif



// Aligned memory allocation
#ifdef SIMINT_SCALAR
  #define SIMINT_ALLOC(x) malloc((x))
  #define SIMINT_FREE(x) free((x))
#else
  #define SIMINT_ALLOC(x) _mm_malloc((x), SIMINT_SIMD_ALIGN)
  #define SIMINT_FREE(x) _mm_free((x))
#endif


// round up the number of elements to the nearest boundary
#define SIMINT_SIMD_ROUND(x) ((x + ((SIMINT_SIMD_LEN-1))) & (~(SIMINT_SIMD_LEN-1)))


// align an array
// ie double somearr[200] SIMINT_ALIGN_ARRAY
#define SIMINT_ALIGN_ARRAY_DBL __attribute__((aligned(SIMINT_SIMD_ALIGN_DBL)))
#define SIMINT_ALIGN_ARRAY_INT __attribute__((aligned(SIMINT_SIMD_ALIGN_INT)))


// Number of shells to use in a batch
#define SIMINT_NSHELL_SIMD (2*SIMINT_SIMD_LEN)

