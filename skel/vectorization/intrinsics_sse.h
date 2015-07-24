#ifndef SIMINT_SSE_H
#define SIMINT_SSE_H

#include <emmintrin.h>

union double2
{
    __m128d d_128;
    double d[2];
};

#endif
