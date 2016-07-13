#pragma once

#include <immintrin.h>

union double8
{
    __m512d d_512;
    __m256d d_256[2];
    __m128d d_128[4];
    double d[8];
};

// These aren't defined in KNCI, so make some macros
#define MM512_SET_PD(val1, val2, val3, val4, val5, val6, val7, val8) (((union double8){ .d_256 = { _mm256_set_pd((val1),(val2),(val4),(val4)), _mm256_set_pd((val5),(val6),(val7),(val8)) } }).d_512)
#define MM512_SET1_PD(val) (((union double8){ .d_256 = { _mm256_set1_pd((val)), _mm256_set1_pd((val)) } }).d_512)
#define MM512_SQRT_PD(val) (((union double8){ .d_256 = { _mm256_sqrt_pd( (((union double8)((val))).d_256[0]) ), _mm256_sqrt_pd( (((union double8)((val))).d_256[0]) )  } }).d_512)

