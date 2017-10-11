#pragma once

#include "simint/shell/shell.h"
#include "simint/ostei/ostei_config.h"

#ifdef __cplusplus
#include "simint/cpp_restrict.hpp"
extern "C" {
#endif


/*! \brief Compute an ostei given shell pair information
 *
 * \param [in] P The shell pairs for the bra side of the integral 
 * \param [in] Q The shell pairs for the ket side of the integral
 * \param [in] screen_tol Tolerance for screening (set to zero to disable)
 * \param [in] work Workspace to use in calculating the integrals
 * \param [inout] integrals Storage for the final integrals. Since size information
 *                          is not passed, you are expected to ensure that this buffer
 *                          is large enough
 */
int simint_compute_eri(struct simint_multi_shellpair const * P,
                       struct simint_multi_shellpair const * Q,
                       double screen_tol,
                       double * restrict work,
                       double * restrict integrals);

/*! \brief Compute an ostei given shell pair information
 *
 * \param [in] deriv Order of the derivative to compute
 * \param [in] P The shell pairs for the bra side of the integral 
 * \param [in] Q The shell pairs for the ket side of the integral
 * \param [in] screen_tol Tolerance for screening (set to zero to disable)
 * \param [in] work Workspace to use in calculating the integrals
 * \param [inout] integrals Storage for the final integrals. Since size information
 *                          is not passed, you are expected to ensure that this buffer
 *                          is large enough
 */
int simint_compute_eri_deriv(int deriv,
                             struct simint_multi_shellpair const * P,
                             struct simint_multi_shellpair const * Q,
                             double screen_tol,
                             double * restrict work,
                             double * restrict integrals);


/*! \brief Get the required size of the workspace required (number of elements)
 *
 * \param [in] derorder Order of the derivative (0 = no derivative, 1 = first derivative)
 * \param [in] maxam Maximum angular momentum to be used in an ERI calculation
 * \return Minimum size of the workspace required (as number of double-precision elements)
 */
size_t simint_eri_worksize(int derorder, int maxam);


/*! \brief Get the required size of the workspace required (in bytes)
 *
 * \param [in] derorder Order of the derivative (0 = no derivative, 1 = first derivative)
 * \param [in] maxam Maximum angular momentum to be used in an ERI calculation
 * \return Minimum size of the workspace required (in bytes)
 */
size_t simint_eri_workmem(int derorder, int maxam);


#ifdef __cplusplus
}
#endif
