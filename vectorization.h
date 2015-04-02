#ifndef SIMINT_VECTORIZATION_H
#define SIMINT_VECTORIZATION_H

#include <stdlib.h>

#if defined SIMINT_AVX

  #define SIMD_LEN 4
  #define SIMD_ALIGN 32

  #define SIMD_ALIGN_DBL (SIMD_ALIGN/8)

  #define ASSUME_ALIGN(x)  __assume_aligned((x), SIMD_ALIGN)
  #define ASSUME __assume

  // Aligned memory allocation
  #define ALLOC(x) _mm_malloc((x), SIMD_ALIGN)
  #define FREE(x) _mm_free((x))

  // round up to the nearest SIMD_ALIGN boundary
  #define SIMD_ROUND(x) ((x + ((SIMD_ALIGN-1))) & (~(SIMD_ALIGN-1)))
  #define SIMD_ROUND_DBL(x) ((x + ((SIMD_ALIGN_DBL-1))) & (~(SIMD_ALIGN_DBL-1)))

#elif defined SIMINT_MIC

  #define SIMD_LEN 8
  #define SIMD_ALIGN 64

  #define SIMD_ALIGN_DBL (SIMD_ALIGN/8)

  #define ASSUME_ALIGN(x)  __assume_aligned((x), SIMD_ALIGN)
  #define ASSUME __assume

  // Aligned memory allocation
  #define ALLOC(x) _mm_malloc((x), SIMD_ALIGN)
  #define FREE(x) _mm_free((x))

  // round up to the nearest SIMD_ALIGN boundary
  #define SIMD_ROUND(x) ((x + ((SIMD_ALIGN-1))) & (~(SIMD_ALIGN-1)))
  #define SIMD_ROUND_DBL(x) ((x + ((SIMD_ALIGN_DBL-1))) & (~(SIMD_ALIGN_DBL-1)))

#else

  #define SIMD_LEN 0
  #define SIMD_ALIGN 1
  #define SIMD_ALIGN_DBL 1

  // empty; do nothing
  #define ASSUME_ALIGN(x)
  #define ASSUME

  #define ALLOC(x) malloc((x))
  #define FREE(x) free((x))
  #define SIMD_ROUND(x) (x)
  #define SIMD_ROUND(x) (x)
  #define SIMD_ROUND_DBL(x) (x)

#endif


#endif
