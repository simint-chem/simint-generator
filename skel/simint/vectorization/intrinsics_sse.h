#pragma once

#include <immintrin.h>
#include <math.h>

union double2
{
    __m128d d_128;
    double d[2];
};


// Missing GCC vectorized exp
static inline __m128d simint_exp_vec2(__m128d x)
{
    union double2 u = { x };
    union double2 res;
    for(int i = 0; i < 2; i++)
        res.d[i] = exp(u.d[i]);
    return res.d_128;
}

