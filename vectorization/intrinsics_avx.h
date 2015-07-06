#ifndef SIMINT_AVX_H
#define SIMINT_AVX_H

#include <immintrin.h>

union double4
{
    __m256d d;
    __m128d d_128[2];
    double v[4];
};

#endif
