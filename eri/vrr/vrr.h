#ifndef VRR__H
#define VRR__H

//////////////////////////////////////////////
// VRR functions
//////////////////////////////////////////////

#include "vectorization.h"




// VRR to obtain AUX_INT__p_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__p_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p,
                        double * const restrict AUX_INT__p_s_s_s,
                        double const * const restrict AUX_INT__s_s_s_s);



// VRR to obtain AUX_INT__d_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__d_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__d_s_s_s,
                        double const * const restrict AUX_INT__p_s_s_s,
                        double const * const restrict AUX_INT__s_s_s_s);



// VRR to obtain AUX_INT__f_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__f_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__f_s_s_s,
                        double const * const restrict AUX_INT__d_s_s_s,
                        double const * const restrict AUX_INT__p_s_s_s);



// VRR to obtain AUX_INT__g_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__g_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__g_s_s_s,
                        double const * const restrict AUX_INT__f_s_s_s,
                        double const * const restrict AUX_INT__d_s_s_s);



// VRR to obtain AUX_INT__h_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__h_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__h_s_s_s,
                        double const * const restrict AUX_INT__g_s_s_s,
                        double const * const restrict AUX_INT__f_s_s_s);



// VRR to obtain AUX_INT__i_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__i_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__i_s_s_s,
                        double const * const restrict AUX_INT__h_s_s_s,
                        double const * const restrict AUX_INT__g_s_s_s);



// VRR to obtain AUX_INT__j_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__j_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__j_s_s_s,
                        double const * const restrict AUX_INT__i_s_s_s,
                        double const * const restrict AUX_INT__h_s_s_s);



// VRR to obtain AUX_INT__k_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__k_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__k_s_s_s,
                        double const * const restrict AUX_INT__j_s_s_s,
                        double const * const restrict AUX_INT__i_s_s_s);



// VRR to obtain AUX_INT__l_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__l_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__l_s_s_s,
                        double const * const restrict AUX_INT__k_s_s_s,
                        double const * const restrict AUX_INT__j_s_s_s);



// VRR to obtain AUX_INT__m_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__m_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__m_s_s_s,
                        double const * const restrict AUX_INT__l_s_s_s,
                        double const * const restrict AUX_INT__k_s_s_s);



// VRR to obtain AUX_INT__n_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__n_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__n_s_s_s,
                        double const * const restrict AUX_INT__m_s_s_s,
                        double const * const restrict AUX_INT__l_s_s_s);



// VRR to obtain AUX_INT__o_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__o_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__o_s_s_s,
                        double const * const restrict AUX_INT__n_s_s_s,
                        double const * const restrict AUX_INT__m_s_s_s);



// VRR to obtain AUX_INT__q_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__q_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__q_s_s_s,
                        double const * const restrict AUX_INT__o_s_s_s,
                        double const * const restrict AUX_INT__n_s_s_s);



// VRR to obtain AUX_INT__r_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__r_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__r_s_s_s,
                        double const * const restrict AUX_INT__q_s_s_s,
                        double const * const restrict AUX_INT__o_s_s_s);



// VRR to obtain AUX_INT__t_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__t_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__t_s_s_s,
                        double const * const restrict AUX_INT__r_s_s_s,
                        double const * const restrict AUX_INT__q_s_s_s);



// VRR to obtain AUX_INT__u_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__u_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__u_s_s_s,
                        double const * const restrict AUX_INT__t_s_s_s,
                        double const * const restrict AUX_INT__r_s_s_s);



// VRR to obtain AUX_INT__v_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__v_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__v_s_s_s,
                        double const * const restrict AUX_INT__u_s_s_s,
                        double const * const restrict AUX_INT__t_s_s_s);



// VRR to obtain AUX_INT__w_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__w_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__w_s_s_s,
                        double const * const restrict AUX_INT__v_s_s_s,
                        double const * const restrict AUX_INT__u_s_s_s);



// VRR to obtain AUX_INT__x_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__x_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__x_s_s_s,
                        double const * const restrict AUX_INT__w_s_s_s,
                        double const * const restrict AUX_INT__v_s_s_s);



// VRR to obtain AUX_INT__y_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN)
void VRR_AUX_INT__y_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__y_s_s_s,
                        double const * const restrict AUX_INT__x_s_s_s,
                        double const * const restrict AUX_INT__w_s_s_s);
#endif
