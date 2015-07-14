#ifndef SIMINT_VECTORIZATION_H
#define SIMINT_VECTORIZATION_H

// This file included below is auto generated from the 
// eri_generator program
// It defines SIMINT_SIMD_LEN
#include "vectorization/vectorization_generated.h"


#include <stdlib.h>


#define SIMINT_SIMD_ALIGN 8*SIMINT_SIMD_LEN
#define SIMINT_SIMD_ALIGN_DBL (SIMINT_SIMD_ALIGN/8)

#define ASSUME_ALIGN(x)  __assume_aligned((x), SIMINT_SIMD_ALIGN)
#define ASSUME __assume

// Aligned memory allocation
#define ALLOC(x) _mm_malloc((x), SIMINT_SIMD_ALIGN)
#define FREE(x) _mm_free((x))

// round up to the nearest SIMD_ALIGN boundary
#define SIMINT_SIMD_ROUND(x) ((x + ((SIMINT_SIMD_ALIGN-1))) & (~(SIMINT_SIMD_ALIGN-1)))
#define SIMINT_SIMD_ROUND_DBL(x) ((x + ((SIMINT_SIMD_ALIGN_DBL-1))) & (~(SIMINT_SIMD_ALIGN_DBL-1)))

// align an array
// ie double somearr[200] SIMINT_ALIGN_ARRAY
#define SIMINT_ALIGN_ARRAY __attribute__((aligned(SIMINT_SIMD_ALIGN)))


#define SIMINT_NSHELL_SIMD (2*SIMINT_SIMD_LEN)




#endif
