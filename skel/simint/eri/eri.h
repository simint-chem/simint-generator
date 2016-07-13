#ifndef SIMINT__ERI_H
#define SIMINT__ERI_H

#ifdef __cplusplus
extern "C" {
#endif

#include "simint/eri/gen/eri_generated.h"


// Disable this warning:
// remark #2620: attribute appears more than once
// No idea what it means, and from what I can tell the
// syntax is correct
#ifdef __INTEL_COMPILER
    #pragma warning(push)
    #pragma warning(disable:2620)
#endif

#include "simint/eri/gen/hrr_generated.h"
#include "simint/eri/gen/vrr_generated.h"
#include "simint/eri/gen/et_generated.h"

#ifdef __INTEL_COMPILER
    #pragma warning(pop)
#endif


// Typedefs for eri functions
typedef int (*simint_erifunc_sharedwork)(struct multishell_pair const, struct multishell_pair const, double * const restrict, double * const restrict);
typedef int (*simint_erifunc)(struct multishell_pair const, struct multishell_pair const, double * const restrict);


// Stores pointers to the eri functions
extern simint_erifunc_sharedwork  simint_erifunc_sharedwork_array[SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1];
extern simint_erifunc             simint_erifunc_array[SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1];




// Convenience
inline int simint_compute_eri(struct multishell_pair const P,
                              struct multishell_pair const Q,
                              double * const restrict integrals)
{
    return simint_erifunc_array[P.am1][P.am2][Q.am1][Q.am2](P, Q, integrals);
}

inline int simint_compute_eri_sharedwork(struct multishell_pair const P,
                                         struct multishell_pair const Q,
                                         double * const restrict work,
                                         double * const restrict integrals)
{
    return simint_erifunc_sharedwork_array[P.am1][P.am2][Q.am1][Q.am2](P, Q, work, integrals);
}


#ifdef __cplusplus
}
#endif

#endif
