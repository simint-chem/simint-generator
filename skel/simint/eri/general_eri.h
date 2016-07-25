#pragma once

#include "simint/shell/shell.h"
#include "simint/vectorization/vectorization.h"

#ifdef __cplusplus
extern "C" {
#endif
 
/*! \brief Compute an eri via a general, unoptimized algorithm
 *
 * \param [in] P The shell pairs for the bra side of the integral 
 * \param [in] Q The shell pairs for the ket side of the integral
 * \param [in] work Workspace to use in calculating the integrals
 * \param [inout] integrals Storage for the final integrals. Since size information
 *                          is not passed, you are expected to ensure that this buffer
 *                          is large enough
 */
int eri_sharedwork_X_s_s_s(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__X_s_s_s);

int eri_sharedwork_s_X_s_s(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__s_X_s_s);

int eri_sharedwork_s_s_X_s(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__s_s_X_s);

int eri_sharedwork_s_s_s_X(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__s_s_s_X);

int eri_sharedwork_X_s_X_s(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__X_s_X_s);

int eri_sharedwork_s_X_s_X(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__s_X_s_X);

int eri_sharedwork_s_X_X_s(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__s_X_X_s);

int eri_sharedwork_X_s_s_X(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__X_s_s_X);

int eri_X_s_s_s(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__X_s_s_s);

int eri_s_X_s_s(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__s_X_s_s);

int eri_s_s_X_s(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__s_s_X_s);

int eri_s_s_s_X(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__s_s_s_X);

int eri_X_s_X_s(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__X_s_X_s);

int eri_s_X_s_X(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__s_X_s_X);

int eri_s_X_X_s(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__s_X_X_s);

int eri_X_s_s_X(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__X_s_s_X);



/*! \brief Compute X_s_s_s via a single vrr step
 *
 * \todo writeme
 */
void general_vrr1(int i, int N,
                  __m256d one_over_2p, __m256d a_over_p,
                  __m256d PA[3], __m256d PQ[3],
                  __m256d * theta1, __m256d * theta2,
                  __m256d * output);


/*! \brief Compute X_s_X_s via a single vrr step
 *
 * \todo writeme
 */
void general_vrr2(int i, int k, int N,
                  __m256d one_over_2q, __m256d one_over_2pq,
                  __m256d a_over_q,
                  __m256d QC[3], __m256d PQ[3],
                  __m256d * theta1, __m256d * theta2, __m256d * theta3,
                  __m256d * output);

#ifdef __cplusplus
}
#endif
