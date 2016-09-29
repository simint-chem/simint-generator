#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "simint/eri/gen/eri_generated.h"


//! A pointer to a function that calculates ERI utilizing a shared workspace
typedef int (*simint_erifunc_sharedwork)(struct simint_multi_shellpair const,
                                         struct simint_multi_shellpair const,
                                         double,
                                         double * restrict,
                                         double * restrict);

//! A pointer to a function that calculates ERI that allocates its own workspace
typedef int (*simint_erifunc)(struct simint_multi_shellpair const,
                              struct simint_multi_shellpair const,
                              double,
                              double * restrict);



/*! \brief Compute an eri given shell pair information
 *
 * \param [in] P The shell pairs for the bra side of the integral 
 * \param [in] Q The shell pairs for the ket side of the integral
 * \param [in] screen_tol Tolerance for screening (set to zero to disable)
 * \param [inout] integrals Storage for the final integrals. Since size information
 *                          is not passed, you are expected to ensure that this buffer
 *                          is large enough
 */
int simint_compute_eri(struct simint_multi_shellpair const * P,
                       struct simint_multi_shellpair const * Q,
                       double screen_tol,
                       double * restrict integrals);


/*! \brief Compute an eri given shell pair information
 *
 * \param [in] P The shell pairs for the bra side of the integral 
 * \param [in] Q The shell pairs for the ket side of the integral
 * \param [in] screen_tol Tolerance for screening (set to zero to disable)
 * \param [in] work Workspace to use in calculating the integrals
 * \param [inout] integrals Storage for the final integrals. Since size information
 *                          is not passed, you are expected to ensure that this buffer
 *                          is large enough
 */
int simint_compute_eri_sharedwork(struct simint_multi_shellpair const * P,
                                  struct simint_multi_shellpair const * Q,
                                  double screen_tol,
                                  double * restrict work,
                                  double * restrict integrals);


#ifdef __cplusplus
}
#endif

