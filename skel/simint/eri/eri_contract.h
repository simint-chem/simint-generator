#pragma once

#include "simint/vectorization/vectorization.h"

#ifdef SIMINT_AVX
    #include "simint/eri/eri_contract_avx.h"
#elif defined SIMINT_SSE
    #include "simint/eri/eri_contract_sse.h"
#else
    #include "simint/eri/eri_contract_scalar.h"
#endif
