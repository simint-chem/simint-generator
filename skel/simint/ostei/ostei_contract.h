#pragma once

#include "simint/vectorization/vectorization.h"

#ifdef SIMINT_AVX
    #include "simint/ostei/ostei_contract_avx.h"
#elif defined SIMINT_SSE
    #include "simint/ostei/ostei_contract_sse.h"
#else
    #include "simint/ostei/ostei_contract_scalar.h"
#endif
