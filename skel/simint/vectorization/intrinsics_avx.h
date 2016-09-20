#pragma once

#include <immintrin.h>
#include <stdint.h>

#include "simint/vectorization/intrinsics_sse.h"

union double4
{
    __m256d d_256;
    __m128d d_128[2];
    double d[4];
};

