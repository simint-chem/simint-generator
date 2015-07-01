#ifndef SIMINT_AVX_H
#define SIMINT_AVX_H

#include <immintrin.h>

union double4
{
    __m256d d;
    double v[4];
};

#endif
