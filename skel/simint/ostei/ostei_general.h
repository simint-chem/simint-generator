#pragma once

#include "simint/ostei/ostei.h"
#include "simint/vectorization/vectorization.h"


#ifdef __cplusplus
extern "C" {
#endif

//////////////////////////////////////////////////////
// General VRR
//////////////////////////////////////////////////////


/*! \brief Compute a vrr step for incrementing a single center,
 *         all other centers being zero
 *
 * \todo writeme
 */
void ostei_general_vrr1(int i, int num_n,
                        SIMINT_DBLTYPE one_over_2p, SIMINT_DBLTYPE a_over_p,
                        const SIMINT_DBLTYPE aop_PQ[3],
                        const SIMINT_DBLTYPE PA[3],
                        SIMINT_DBLTYPE const * restrict theta1,
                        SIMINT_DBLTYPE const * restrict theta2,
                        SIMINT_DBLTYPE * restrict output);


/*! \brief Compute X_s_s_s via a single vrr step
 *
 * \todo writeme
 */
static inline
void ostei_general_vrr1_I(int i, int num_n,
                          SIMINT_DBLTYPE one_over_2p, SIMINT_DBLTYPE a_over_p,
                          const SIMINT_DBLTYPE aop_PQ[3],
                          const SIMINT_DBLTYPE PA[3],
                          SIMINT_DBLTYPE const * restrict theta1,
                          SIMINT_DBLTYPE const * restrict theta2,
                          SIMINT_DBLTYPE * restrict output)
{
    ostei_general_vrr1(i, num_n, one_over_2p, a_over_p, aop_PQ, PA,
                 theta1, theta2, output);
}

/*! \brief Compute s_X_s_s via a single vrr step
 *
 * \todo writeme
 */
static inline
void ostei_general_vrr1_J(int j, int num_n,
                          SIMINT_DBLTYPE one_over_2p, SIMINT_DBLTYPE a_over_p,
                          const SIMINT_DBLTYPE aop_PQ[3],
                          const SIMINT_DBLTYPE PB[3],
                          SIMINT_DBLTYPE const * restrict theta1,
                          SIMINT_DBLTYPE const * restrict theta2,
                          SIMINT_DBLTYPE * restrict output)
{
    ostei_general_vrr1(j, num_n, one_over_2p, a_over_p, aop_PQ, PB,
                 theta1, theta2, output);
}

/*! \brief Compute s_s_X_s via a single vrr step
 *
 * \todo writeme
 */
static inline
void ostei_general_vrr1_K(int k, int num_n,
                          SIMINT_DBLTYPE one_over_2q, SIMINT_DBLTYPE a_over_q,
                          const SIMINT_DBLTYPE aoq_PQ[3],
                          const SIMINT_DBLTYPE QC[3],
                          SIMINT_DBLTYPE const * restrict theta1,
                          SIMINT_DBLTYPE const * restrict theta2,
                          SIMINT_DBLTYPE * restrict output)
{
    ostei_general_vrr1(k, num_n, one_over_2q, a_over_q, aoq_PQ, QC,
                 theta1, theta2, output);
}

/*! \brief Compute s_s_s_X via a single vrr step
 *
 * \todo writeme
 */
static inline
void ostei_general_vrr1_L(int l, int num_n,
                          SIMINT_DBLTYPE one_over_2q, SIMINT_DBLTYPE a_over_q,
                          const SIMINT_DBLTYPE aoq_PQ[3],
                          const SIMINT_DBLTYPE QD[3],
                          SIMINT_DBLTYPE const * restrict theta1,
                          SIMINT_DBLTYPE const * restrict theta2,
                          SIMINT_DBLTYPE * restrict output)
{
    ostei_general_vrr1(l, num_n, one_over_2q, a_over_q, aoq_PQ, QD,
                 theta1, theta2, output);
}


/*! \brief Compute X_X_X_X via a single vrr step on i
 *
 * \todo writeme
 */
void ostei_general_vrr_I(int i, int j, int k, int l, int num_n,
                         SIMINT_DBLTYPE one_over_2p, SIMINT_DBLTYPE a_over_p,
                         SIMINT_DBLTYPE one_over_2pq,
                         const SIMINT_DBLTYPE aop_PQ[3],
                         const SIMINT_DBLTYPE PA[3],
                         SIMINT_DBLTYPE const * restrict theta1,
                         SIMINT_DBLTYPE const * restrict theta2,
                         SIMINT_DBLTYPE const * restrict theta3,
                         SIMINT_DBLTYPE const * restrict theta4,
                         SIMINT_DBLTYPE const * restrict theta5,
                         SIMINT_DBLTYPE * restrict output);

/*! \brief Compute X_X_X_X via a single vrr step on j
 *
 * \todo writeme
 */
void ostei_general_vrr_J(int i, int j, int k, int l, int num_n,
                         SIMINT_DBLTYPE one_over_2p, SIMINT_DBLTYPE a_over_p,
                         SIMINT_DBLTYPE one_over_2pq,
                         const SIMINT_DBLTYPE aop_PQ[3],
                         const SIMINT_DBLTYPE PB[3],
                         SIMINT_DBLTYPE const * restrict theta1,
                         SIMINT_DBLTYPE const * restrict theta2,
                         SIMINT_DBLTYPE const * restrict theta3,
                         SIMINT_DBLTYPE const * restrict theta4,
                         SIMINT_DBLTYPE const * restrict theta5,
                         SIMINT_DBLTYPE * restrict output);

/*! \brief Compute X_X_X_X via a single vrr step on k
 *
 * \todo writeme
 */
void ostei_general_vrr_K(int i, int j, int k, int l, int num_n,
                         SIMINT_DBLTYPE one_over_2q, SIMINT_DBLTYPE a_over_q,
                         SIMINT_DBLTYPE one_over_2pq,
                         const SIMINT_DBLTYPE aoq_PQ[3],
                         const SIMINT_DBLTYPE QC[3],
                         SIMINT_DBLTYPE const * restrict theta1,
                         SIMINT_DBLTYPE const * restrict theta2,
                         SIMINT_DBLTYPE const * restrict theta3,
                         SIMINT_DBLTYPE const * restrict theta4,
                         SIMINT_DBLTYPE const * restrict theta5,
                         SIMINT_DBLTYPE * restrict output);

/*! \brief Compute X_X_X_X via a single vrr step on l
 *
 * \todo writeme
 */
void ostei_general_vrr_L(int i, int j, int k, int l, int num_n,
                         SIMINT_DBLTYPE one_over_2q, SIMINT_DBLTYPE a_over_q,
                         SIMINT_DBLTYPE one_over_2pq,
                         const SIMINT_DBLTYPE aoq_PQ[3],
                         const SIMINT_DBLTYPE QD[3],
                         SIMINT_DBLTYPE const * restrict theta1,
                         SIMINT_DBLTYPE const * restrict theta2,
                         SIMINT_DBLTYPE const * restrict theta3,
                         SIMINT_DBLTYPE const * restrict theta4,
                         SIMINT_DBLTYPE const * restrict theta5,
                         SIMINT_DBLTYPE * restrict output);


//////////////////////////////////////////////////////
// General HRR
//////////////////////////////////////////////////////

/*! \brief Compute X_X_X_X via a single hrr step on i
 *
 * \todo writeme
 */
void ostei_general_hrr_I(int i, int j, int k, int l,
                         const double AB[3],
                         const double * theta1,
                         const double * theta2,
                         double * output);

/*! \brief Compute X_X_X_X via a single hrr step on j
 *
 * \todo writeme
 */
void ostei_general_hrr_J(int i, int j, int k, int l,
                         const double AB[3],
                         const double * theta1,
                         const double * theta2,
                         double * output);

/*! \brief Compute X_X_X_X via a single hrr step on k
 *
 * \todo writeme
 */
void ostei_general_hrr_K(int i, int j, int k, int l,
                         const double CD[3],
                         const double * theta1,
                         const double * theta2,
                         double * output);

/*! \brief Compute X_X_X_X via a single hrr step on l
 *
 * \todo writeme
 */
void ostei_general_hrr_L(int i, int j, int k, int l,
                         const double CD[3],
                         const double * theta1,
                         const double * theta2,
                         double * output);


#ifdef __cplusplus
}
#endif
