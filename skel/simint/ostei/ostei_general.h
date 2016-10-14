#pragma once

#include "simint/shell/shell.h"
#include "simint/vectorization/vectorization.h"

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Compute a vrr step for incrementing a single center,
 *         all other centers being zero
 *
 * \todo writeme
 */
void ostei_general_vrr1(int i, int num_n,
                        __m256d one_over_2p, __m256d a_over_p,
                        const __m256d aop_PQ[3],
                        const __m256d PA[3],
                        __m256d const * restrict theta1,
                        __m256d const * restrict theta2,
                        __m256d * restrict output);


/*! \brief Compute X_s_s_s via a single vrr step
 *
 * \todo writeme
 */
static inline
void ostei_general_vrr1_i(int i, int num_n,
                          __m256d one_over_2p, __m256d a_over_p,
                          const __m256d aop_PQ[3],
                          const __m256d PA[3],
                        __m256d const * restrict theta1,
                        __m256d const * restrict theta2,
                        __m256d * restrict output)
{
    ostei_general_vrr1(i, num_n, one_over_2p, a_over_p, aop_PQ, PA,
                 theta1, theta2, output);
}

/*! \brief Compute s_X_s_s via a single vrr step
 *
 * \todo writeme
 */
static inline
void ostei_general_vrr1_j(int j, int num_n,
                          __m256d one_over_2p, __m256d a_over_p,
                          const __m256d aop_PQ[3],
                          const __m256d PB[3],
                        __m256d const * restrict theta1,
                        __m256d const * restrict theta2,
                        __m256d * restrict output)
{
    ostei_general_vrr1(j, num_n, one_over_2p, a_over_p, aop_PQ, PB,
                 theta1, theta2, output);
}

/*! \brief Compute s_s_X_s via a single vrr step
 *
 * \todo writeme
 */
static inline
void ostei_general_vrr1_k(int k, int num_n,
                          __m256d one_over_2q, __m256d a_over_q,
                          const __m256d aoq_PQ[3],
                          const __m256d QC[3],
                          __m256d const * restrict theta1,
                          __m256d const * restrict theta2,
                          __m256d * restrict output)
{
    ostei_general_vrr1(k, num_n, one_over_2q, a_over_q, aoq_PQ, QC,
                 theta1, theta2, output);
}

/*! \brief Compute s_s_s_X via a single vrr step
 *
 * \todo writeme
 */
static inline
void ostei_general_vrr1_l(int l, int num_n,
                          __m256d one_over_2q, __m256d a_over_q,
                          const __m256d aoq_PQ[3],
                          const __m256d QD[3],
                          __m256d const * restrict theta1,
                          __m256d const * restrict theta2,
                          __m256d * restrict output)
{
    ostei_general_vrr1(l, num_n, one_over_2q, a_over_q, aoq_PQ, QD,
                 theta1, theta2, output);
}


/*! \brief Compute X_X_X_X via a single vrr step on i
 *
 * \todo writeme
 */
void ostei_general_vrr_i(int i, int j, int k, int l, int num_n,
                         __m256d one_over_2p, __m256d a_over_p,
                         __m256d one_over_2pq,
                         const __m256d aop_PQ[3],
                         const __m256d PA[3],
                         __m256d const * restrict theta1,
                         __m256d const * restrict theta2,
                         __m256d const * restrict theta3,
                         __m256d const * restrict theta4,
                         __m256d const * restrict theta5,
                         __m256d * restrict output);

/*! \brief Compute X_X_X_X via a single vrr step on j
 *
 * \todo writeme
 */
void ostei_general_vrr_j(int i, int j, int k, int l, int num_n,
                         __m256d one_over_2p, __m256d a_over_p,
                         __m256d one_over_2pq,
                         const __m256d aop_PQ[3],
                         const __m256d PB[3],
                         __m256d const * restrict theta1,
                         __m256d const * restrict theta2,
                         __m256d const * restrict theta3,
                         __m256d const * restrict theta4,
                         __m256d const * restrict theta5,
                         __m256d * restrict output);

/*! \brief Compute X_X_X_X via a single vrr step on k
 *
 * \todo writeme
 */
void ostei_general_vrr_k(int i, int j, int k, int l, int num_n,
                         __m256d one_over_2q, __m256d a_over_q,
                         __m256d one_over_2pq,
                         const __m256d aoq_PQ[3],
                         const __m256d QC[3],
                         __m256d const * restrict theta1,
                         __m256d const * restrict theta2,
                         __m256d const * restrict theta3,
                         __m256d const * restrict theta4,
                         __m256d const * restrict theta5,
                         __m256d * restrict output);

/*! \brief Compute X_X_X_X via a single vrr step on l
 *
 * \todo writeme
 */
void ostei_general_vrr_l(int i, int j, int k, int l, int num_n,
                         __m256d one_over_2q, __m256d a_over_q,
                         __m256d one_over_2pq,
                         const __m256d aoq_PQ[3],
                         const __m256d QD[3],
                         __m256d const * restrict theta1,
                         __m256d const * restrict theta2,
                         __m256d const * restrict theta3,
                         __m256d const * restrict theta4,
                         __m256d const * restrict theta5,
                         __m256d * restrict output);

#ifdef __cplusplus
}
#endif
