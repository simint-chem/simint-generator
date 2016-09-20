#pragma once

#include "simint/vectorization/vectorization.h"

#ifdef SIMINT_AVX
    #include "simint/eri/eri_contract_avx.h"
#else
    #include "simint/eri/eri_contract_scalar.h"
#endif
