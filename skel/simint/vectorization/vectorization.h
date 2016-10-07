#pragma once

// Autogenerated
#include "simint/simint_config.h"

#define SIMINT_SIMD_ALIGN_DBL (SIMINT_SIMD_LEN*8)
#define SIMINT_SIMD_ALIGN_INT (SIMINT_SIMD_LEN*sizeof(int))

// "Max" alignment. Can be used to allocate memory of mixed types
// (ie ints and doubles can be allocated this way)
#define SIMINT_SIMD_ALIGN SIMINT_SIMD_ALIGN_DBL


#define ASSUME_ALIGN_DBL(x)  __assume_aligned((x), SIMINT_SIMD_ALIGN_DBL)
#define ASSUME_ALIGN_INT(x)  __assume_aligned((x), SIMINT_SIMD_ALIGN_INT)
#define ASSUME __assume


// Aligned memory allocation
#define ALLOC(x) _mm_malloc((x), SIMINT_SIMD_ALIGN)
#define FREE(x) _mm_free((x))


// round up the number of elements to the nearest boundary
#define SIMINT_SIMD_ROUND(x) ((x + ((SIMINT_SIMD_LEN-1))) & (~(SIMINT_SIMD_LEN-1)))


// align an array
// ie double somearr[200] SIMINT_ALIGN_ARRAY
#define SIMINT_ALIGN_ARRAY_DBL __attribute__((aligned(SIMINT_SIMD_ALIGN_DBL)))
#define SIMINT_ALIGN_ARRAY_INT __attribute__((aligned(SIMINT_SIMD_ALIGN_INT)))


// Number of shells to use in a batch
#define SIMINT_NSHELL_SIMD (2*SIMINT_SIMD_LEN)
