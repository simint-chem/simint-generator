#pragma once

#include <immintrin.h>
#include <stdint.h>
#include <math.h>

#include "simint/vectorization/intrinsics_avx.h"

union double8
{
    __m512d d_512;
    __m256d d_256[2];
    __m128d d_128[4];
    double d[8];
};


// Missing GCC vectorized exp
static inline __m512d simint_exp_vec8(__m512d x)
{
    union double8 u = { x };
    union double8 res;
    for(int i = 0; i < 8; i++)
        res.d[i] = exp(u.d[i]);
    return res.d_512;
}

