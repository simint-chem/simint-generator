#pragma once

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


//! A pointer to a function that calculates ERI utilizing a shared workspace
typedef int (*simint_erifunc_sharedwork)(struct simint_multi_shellpair const,
                                         struct simint_multi_shellpair const,
                                         double * restrict,
                                         double * restrict);

//! A pointer to a function that calculates ERI that allocates its own workspace
typedef int (*simint_erifunc)(struct simint_multi_shellpair const,
                              struct simint_multi_shellpair const,
                              double * restrict);



/*! \brief Compute an eri given shell pair information
 *
 * \param [in] P The shell pairs for the bra side of the integral 
 * \param [in] Q The shell pairs for the ket side of the integral
 * \param [inout] integrals Storage for the final integrals. Since size information
 *                          is not passed, you are expected to ensure that this buffer
 *                          is large enough
 */
int simint_compute_eri(struct simint_multi_shellpair const * P,
                       struct simint_multi_shellpair const * Q,
                       double * restrict integrals);


/*! \brief Compute an eri given shell pair information
 *
 * \param [in] P The shell pairs for the bra side of the integral 
 * \param [in] Q The shell pairs for the ket side of the integral
 * \param [in] work Workspace to use in calculating the integrals
 * \param [inout] integrals Storage for the final integrals. Since size information
 *                          is not passed, you are expected to ensure that this buffer
 *                          is large enough
 */
int simint_compute_eri_sharedwork(struct simint_multi_shellpair const * P,
                                  struct simint_multi_shellpair const * Q,
                                  double * restrict work,
                                  double * restrict integrals);


#ifdef __cplusplus
}
#endif

