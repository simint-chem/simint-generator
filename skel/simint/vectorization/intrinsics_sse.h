#pragma once

#include <emmintrin.h>

union double2
{
    __m128d d_128;
    double d[2];
};

