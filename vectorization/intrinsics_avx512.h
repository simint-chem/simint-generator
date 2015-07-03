#ifndef SIMINT_INTRINSICS_AVX512_H
#define SIMINT_INTRINSICS_AVX512_H

#include <immintrin.h>

union double8
{
    __m512d d;
    double v[8];
};

#endif
