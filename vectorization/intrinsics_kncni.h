#ifndef SIMINT_INTRINSICS_AVX512_H
#define SIMINT_INTRINSICS_AVX512_H

#include <immintrin.h>

union double8
{
    __m512d d;
    __m256d d_256[2];
    __m128d d_128[4];
    double v[8];
};

#endif
