#include "constants.h"
#include "test/Libint2.hpp"
#include "test/libint2/libint2.h"
#include "boys/boys_split.h"
#include "test/common.hpp"
#include "vectorization/vectorization.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

// Disable intel warnings
// 193 : zero used for undefined preprocessing identifier "LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18
//       Not much I can do about that 
#ifdef __INTEL_COMPILER
    #pragma warning(disable:193)
#endif



Libint2_ERI::Libint2_ERI(int maxam, size_t maxnprim)
{
    size_t size = maxnprim * maxnprim * maxnprim * maxnprim;

    erival_.resize(size);
    LIBINT2_PREFIXED_NAME(libint2_init_eri)(erival_.data(), maxam, 0); 
}


Libint2_ERI::~Libint2_ERI()
{
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(erival_.data());
}



TimerInfo Libint2_ERI::Integrals(struct multishell_pair P,
                                 struct multishell_pair Q,
                                 double * integrals)
{
    unsigned long long totaltime = libint2_integrals(erival_.data(), P, Q, integrals);
    return {totaltime, 0.0};
}

