#ifndef ERI_H
#define ERI_H

#ifdef __cplusplus
extern "C" {
#endif

#include "eri/gen/eri_generated.h"


// Disable this warning:
// remark #2620: attribute appears more than once
// No idea what it means, and from what I can tell the
// syntax is correct
#ifdef __INTEL_COMPILER
    #pragma warning(push)
    #pragma warning(disable:2620)
#endif

#include "eri/gen/hrr_generated.h"
#include "eri/gen/vrr_generated.h"
#include "eri/gen/et_generated.h"

#ifdef __INTEL_COMPILER
    #pragma warning(pop)
#endif


typedef int (*simint_erifunc_sharedwork)(struct multishell_pair const, struct multishell_pair const, double * const restrict, double * const restrict);
typedef int (*simint_erifunc)(struct multishell_pair const, struct multishell_pair const, double * const restrict);

extern simint_erifunc_sharedwork  simint_compute_eri_sharedwork[SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1];
extern simint_erifunc             simint_compute_eri[SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1];


#ifdef __cplusplus
}
#endif

#endif
