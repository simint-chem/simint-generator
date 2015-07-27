#ifndef SIMINT_VECTORIZATION_H
#define SIMINT_VECTORIZATION_H

// SIMINT_SIMD_LEN should be defined on the command line
#ifndef SIMINT_SIMD_LEN
  #error "SIMINT_SIMD_LEN is not defined!"
#endif

#include <stdlib.h>

#define SIMINT_SIMD_ALIGN_DBL (SIMINT_SIMD_LEN*8)
#define SIMINT_SIMD_ALIGN_INT32 (SIMINT_SIMD_LEN*sizeof(int32_t))  // also used for uint32_t

// "Max" alignment. Can be used to allocate memory of mixed types
// (ie ints and doubles can be allocated this way)
#define SIMINT_SIMD_ALIGN SIMINT_SIMD_ALIGN_DBL


#define ASSUME_ALIGN_DBL(x)  __assume_aligned((x), SIMINT_SIMD_ALIGN_DBL)
#define ASSUME_ALIGN_INT32(x)  __assume_aligned((x), SIMINT_SIMD_ALIGN_INT32)
#define ASSUME __assume


// Aligned memory allocation
#define ALLOC(x) _mm_malloc((x), SIMINT_SIMD_ALIGN)
#define FREE(x) _mm_free((x))


// round up the number of elements to the nearest boundary
#define SIMINT_SIMD_ROUND(x) ((x + ((SIMINT_SIMD_LEN-1))) & (~(SIMINT_SIMD_LEN-1)))


// align an array
// ie double somearr[200] SIMINT_ALIGN_ARRAY
#define SIMINT_ALIGN_ARRAY_DBL __attribute__((aligned(SIMINT_SIMD_ALIGN_DBL)))
#define SIMINT_ALIGN_ARRAY_INT32 __attribute__((aligned(SIMINT_SIMD_ALIGN_INT32)))


// Number of shells to use in a batch
#define SIMINT_NSHELL_SIMD (2*SIMINT_SIMD_LEN)




#endif
