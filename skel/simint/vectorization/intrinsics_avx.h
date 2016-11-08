#pragma once

#include <immintrin.h>
#include <stdint.h>
#include <math.h>

#include "simint/vectorization/intrinsics_sse.h"

union double4
{
    __m256d d_256;
    __m128d d_128[2];
    double d[4];
};


// Missing GCC vectorized exp
static inline __m256d simint_exp_vec4(__m256d x)
{
    union double4 u = { x };
    union double4 res;
    for(int i = 0; i < 4; i++)
        res.d[i] = exp(u.d[i]);
    return res.d_256;
}

