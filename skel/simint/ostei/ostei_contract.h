#pragma once

#include "simint/vectorization/vectorization.h"

#if defined SIMINT_AVX512 || defined SIMINT_MICAVX512
    #include "simint/ostei/ostei_contract_avx512.h"
#elif defined SIMINT_AVX
    #include "simint/ostei/ostei_contract_avx.h"
#elif defined SIMINT_SSE
    #include "simint/ostei/ostei_contract_sse.h"
#else
    #include "simint/ostei/ostei_contract_scalar.h"
#endif
