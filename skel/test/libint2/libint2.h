#ifndef TEST_LIBINT_LIBINT2_H
#define TEST_LIBINT_LIBINT2_H

#include "shell/shell.h"


// Disable intel diagnostics for libint
// These happen in the libint header, so there's not much I can do about them
// 193 : zero used for undefined preprocessing identifier "LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18
// 1418: external function definition with no prior declaration
//   82: storage class is not first
// 2259: non-pointer conversion from "double" to "float" may lose significant bits  
#ifdef __INTEL_COMPILER
    #pragma warning(push)
    #pragma warning(disable:193)
    #pragma warning(disable:1418)
    #pragma warning(disable:82)
    #pragma warning(disable:2259)
#endif
#include <libint2.h>
#ifdef __INTEL_COMPILER
    #pragma warning(pop)
#endif

#ifdef __cplusplus
extern "C" {
#endif

unsigned long long libint2_integrals(Libint_eri_t * erival, struct multishell_pair P, struct multishell_pair Q, double * integrals);


#ifdef __cplusplus
}
#endif

#endif
