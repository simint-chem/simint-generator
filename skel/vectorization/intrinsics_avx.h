#ifndef SIMINT_AVX_H
#define SIMINT_AVX_H

#include <immintrin.h>

#include "vectorization/intrinsics_sse.h"

union double4
{
    __m256d d_256;
    __m128d d_128[2];
    double d[4];
};

#endif
