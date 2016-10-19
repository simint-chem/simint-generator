#pragma once

#include "simint/vectorization/vectorization.h"

#ifdef __cplusplus
extern "C" {
#endif


#define SIMINT_OSTEI_MAXAM 2
#define SIMINT_OSTEI_MAX_WORKSIZE ((SIMINT_SIMD_ROUND(SIMINT_NSHELL_SIMD * 961)))
#define SIMINT_OSTEI_MAX_WORKMEM (SIMINT_OSTEI_MAX_WORKSIZE * sizeof(double))



#ifdef __cplusplus
}
#endif

