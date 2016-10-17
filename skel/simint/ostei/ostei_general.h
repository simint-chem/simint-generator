#pragma once

#include "simint/shell/shell.h"
#include "simint/vectorization/vectorization.h"

#ifdef SIMINT_AVX
  #define DBLTYPE __m256d
#elif defined SIMINT_SSE
  #define DBLTYPE __m128d
#elif defined SIMINT_SCALAR
  #define DBLTYPE double
#endif

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Compute a vrr step for incrementing a single center,
 *         all other centers being zero
 *
 * \todo writeme
 */
void ostei_general_vrr1(int i, int num_n,
                        DBLTYPE one_over_2p, DBLTYPE a_over_p,
                        const DBLTYPE aop_PQ[3],
                        const DBLTYPE PA[3],
                        DBLTYPE const * restrict theta1,
                        DBLTYPE const * restrict theta2,
                        DBLTYPE * restrict output);


/*! \brief Compute X_s_s_s via a single vrr step
 *
 * \todo writeme
 */
static inline
void ostei_general_vrr1_i(int i, int num_n,
                          DBLTYPE one_over_2p, DBLTYPE a_over_p,
                          const DBLTYPE aop_PQ[3],
                          const DBLTYPE PA[3],
                        DBLTYPE const * restrict theta1,
                        DBLTYPE const * restrict theta2,
                        DBLTYPE * restrict output)
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
                          DBLTYPE one_over_2p, DBLTYPE a_over_p,
                          const DBLTYPE aop_PQ[3],
                          const DBLTYPE PB[3],
                        DBLTYPE const * restrict theta1,
                        DBLTYPE const * restrict theta2,
                        DBLTYPE * restrict output)
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
                          DBLTYPE one_over_2q, DBLTYPE a_over_q,
                          const DBLTYPE aoq_PQ[3],
                          const DBLTYPE QC[3],
                          DBLTYPE const * restrict theta1,
                          DBLTYPE const * restrict theta2,
                          DBLTYPE * restrict output)
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
                          DBLTYPE one_over_2q, DBLTYPE a_over_q,
                          const DBLTYPE aoq_PQ[3],
                          const DBLTYPE QD[3],
                          DBLTYPE const * restrict theta1,
                          DBLTYPE const * restrict theta2,
                          DBLTYPE * restrict output)
{
    ostei_general_vrr1(l, num_n, one_over_2q, a_over_q, aoq_PQ, QD,
                 theta1, theta2, output);
}


/*! \brief Compute X_X_X_X via a single vrr step on i
 *
 * \todo writeme
 */
void ostei_general_vrr_i(int i, int j, int k, int l, int num_n,
                         DBLTYPE one_over_2p, DBLTYPE a_over_p,
                         DBLTYPE one_over_2pq,
                         const DBLTYPE aop_PQ[3],
                         const DBLTYPE PA[3],
                         DBLTYPE const * restrict theta1,
                         DBLTYPE const * restrict theta2,
                         DBLTYPE const * restrict theta3,
                         DBLTYPE const * restrict theta4,
                         DBLTYPE const * restrict theta5,
                         DBLTYPE * restrict output);

/*! \brief Compute X_X_X_X via a single vrr step on j
 *
 * \todo writeme
 */
void ostei_general_vrr_j(int i, int j, int k, int l, int num_n,
                         DBLTYPE one_over_2p, DBLTYPE a_over_p,
                         DBLTYPE one_over_2pq,
                         const DBLTYPE aop_PQ[3],
                         const DBLTYPE PB[3],
                         DBLTYPE const * restrict theta1,
                         DBLTYPE const * restrict theta2,
                         DBLTYPE const * restrict theta3,
                         DBLTYPE const * restrict theta4,
                         DBLTYPE const * restrict theta5,
                         DBLTYPE * restrict output);

/*! \brief Compute X_X_X_X via a single vrr step on k
 *
 * \todo writeme
 */
void ostei_general_vrr_k(int i, int j, int k, int l, int num_n,
                         DBLTYPE one_over_2q, DBLTYPE a_over_q,
                         DBLTYPE one_over_2pq,
                         const DBLTYPE aoq_PQ[3],
                         const DBLTYPE QC[3],
                         DBLTYPE const * restrict theta1,
                         DBLTYPE const * restrict theta2,
                         DBLTYPE const * restrict theta3,
                         DBLTYPE const * restrict theta4,
                         DBLTYPE const * restrict theta5,
                         DBLTYPE * restrict output);

/*! \brief Compute X_X_X_X via a single vrr step on l
 *
 * \todo writeme
 */
void ostei_general_vrr_l(int i, int j, int k, int l, int num_n,
                         DBLTYPE one_over_2q, DBLTYPE a_over_q,
                         DBLTYPE one_over_2pq,
                         const DBLTYPE aoq_PQ[3],
                         const DBLTYPE QD[3],
                         DBLTYPE const * restrict theta1,
                         DBLTYPE const * restrict theta2,
                         DBLTYPE const * restrict theta3,
                         DBLTYPE const * restrict theta4,
                         DBLTYPE const * restrict theta5,
                         DBLTYPE * restrict output);


/*! \brief Compute X_X_X_X via a single hrr step on i
 *
 * \todo writeme
 */
void ostei_general_hrr_i(int i, int j, int k, int l,
                         const double AB[3],
                         const double * theta1,
                         const double * theta2,
                         double * output);

/*! \brief Compute X_X_X_X via a single hrr step on j
 *
 * \todo writeme
 */
void ostei_general_hrr_j(int i, int j, int k, int l,
                         const double AB[3],
                         const double * theta1,
                         const double * theta2,
                         double * output);

/*! \brief Compute X_X_X_X via a single hrr step on k
 *
 * \todo writeme
 */
void ostei_general_hrr_k(int i, int j, int k, int l,
                         const double CD[3],
                         const double * theta1,
                         const double * theta2,
                         double * output);

/*! \brief Compute X_X_X_X via a single hrr step on l
 *
 * \todo writeme
 */
void ostei_general_hrr_l(int i, int j, int k, int l,
                         const double CD[3],
                         const double * theta1,
                         const double * theta2,
                         double * output);

#ifdef __cplusplus
}
#endif
