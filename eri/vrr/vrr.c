//////////////////////////////////////////////
// VRR functions
//////////////////////////////////////////////

#include "vectorization.h"




// VRR to obtain AUX_INT__p_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_p(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p,
           double * const restrict AUX_INT__p_s_s_s,
           double const * const restrict AUX_INT__s_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__p_s_s_s[num_m * 3];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //p_1_0_0 : STEP: x
                        AUX_INT__p_s_s_s[m * 3 + 0] = P_PA_x * AUX_INT__s_s_s_s[m * 1 + 0] - aop_PQ_x * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //p_0_1_0 : STEP: y
                        AUX_INT__p_s_s_s[m * 3 + 1] = P_PA_y * AUX_INT__s_s_s_s[m * 1 + 0] - aop_PQ_y * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //p_0_0_1 : STEP: z
                        AUX_INT__p_s_s_s[m * 3 + 2] = P_PA_z * AUX_INT__s_s_s_s[m * 1 + 0] - aop_PQ_z * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                    }
}



// VRR to obtain AUX_INT__d_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_d(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__d_s_s_s,
           double const * const restrict AUX_INT__p_s_s_s,
           double const * const restrict AUX_INT__s_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__d_s_s_s[num_m * 6];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //d_2_0_0 : STEP: x
                        AUX_INT__d_s_s_s[m * 6 + 0] = P_PA_x * AUX_INT__p_s_s_s[m * 3 + 0] - aop_PQ_x * AUX_INT__p_s_s_s[(m+1) * 3 + 0]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                        //d_1_1_0 : STEP: y
                        AUX_INT__d_s_s_s[m * 6 + 1] = P_PA_y * AUX_INT__p_s_s_s[m * 3 + 0] - aop_PQ_y * AUX_INT__p_s_s_s[(m+1) * 3 + 0];

                        //d_1_0_1 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 2] = P_PA_z * AUX_INT__p_s_s_s[m * 3 + 0] - aop_PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 0];

                        //d_0_2_0 : STEP: y
                        AUX_INT__d_s_s_s[m * 6 + 3] = P_PA_y * AUX_INT__p_s_s_s[m * 3 + 1] - aop_PQ_y * AUX_INT__p_s_s_s[(m+1) * 3 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                        //d_0_1_1 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 4] = P_PA_z * AUX_INT__p_s_s_s[m * 3 + 1] - aop_PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 1];

                        //d_0_0_2 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 5] = P_PA_z * AUX_INT__p_s_s_s[m * 3 + 2] - aop_PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                    }
}



// VRR to obtain AUX_INT__f_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_f(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__f_s_s_s,
           double const * const restrict AUX_INT__d_s_s_s,
           double const * const restrict AUX_INT__p_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__f_s_s_s[num_m * 10];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //f_3_0_0 : STEP: x
                        AUX_INT__f_s_s_s[m * 10 + 0] = P_PA_x * AUX_INT__d_s_s_s[m * 6 + 0] - aop_PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 0]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  0] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 0] );

                        //f_2_1_0 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 1] = P_PA_y * AUX_INT__d_s_s_s[m * 6 + 0] - aop_PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 0];

                        //f_2_0_1 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 2] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 0] - aop_PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 0];

                        //f_1_2_0 : STEP: x
                        AUX_INT__f_s_s_s[m * 10 + 3] = P_PA_x * AUX_INT__d_s_s_s[m * 6 + 3] - aop_PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 3];

                        //f_1_1_1 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 4] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 1] - aop_PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 1];

                        //f_1_0_2 : STEP: x
                        AUX_INT__f_s_s_s[m * 10 + 5] = P_PA_x * AUX_INT__d_s_s_s[m * 6 + 5] - aop_PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 5];

                        //f_0_3_0 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 6] = P_PA_y * AUX_INT__d_s_s_s[m * 6 + 3] - aop_PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  1] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 1] );

                        //f_0_2_1 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 7] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 3] - aop_PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 3];

                        //f_0_1_2 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 8] = P_PA_y * AUX_INT__d_s_s_s[m * 6 + 5] - aop_PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 5];

                        //f_0_0_3 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 9] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 5] - aop_PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  2] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 2] );

                    }
}



// VRR to obtain AUX_INT__g_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_g(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__g_s_s_s,
           double const * const restrict AUX_INT__f_s_s_s,
           double const * const restrict AUX_INT__d_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__g_s_s_s[num_m * 15];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //g_4_0_0 : STEP: x
                        AUX_INT__g_s_s_s[m * 15 + 0] = P_PA_x * AUX_INT__f_s_s_s[m * 10 + 0] - aop_PQ_x * AUX_INT__f_s_s_s[(m+1) * 10 + 0]
                                      + 3 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  0] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 0] );

                        //g_3_1_0 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 1] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 0] - aop_PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 0];

                        //g_3_0_1 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 2] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 0] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 0];

                        //g_2_2_0 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 3] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 1] - aop_PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  0] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 0] );

                        //g_2_1_1 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 4] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 1] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 1];

                        //g_2_0_2 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 5] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 2] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  0] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 0] );

                        //g_1_3_0 : STEP: x
                        AUX_INT__g_s_s_s[m * 15 + 6] = P_PA_x * AUX_INT__f_s_s_s[m * 10 + 6] - aop_PQ_x * AUX_INT__f_s_s_s[(m+1) * 10 + 6];

                        //g_1_2_1 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 7] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 3] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 3];

                        //g_1_1_2 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 8] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 5] - aop_PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 5];

                        //g_1_0_3 : STEP: x
                        AUX_INT__g_s_s_s[m * 15 + 9] = P_PA_x * AUX_INT__f_s_s_s[m * 10 + 9] - aop_PQ_x * AUX_INT__f_s_s_s[(m+1) * 10 + 9];

                        //g_0_4_0 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 10] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 6] - aop_PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  3] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 3] );

                        //g_0_3_1 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 11] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 6] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 6];

                        //g_0_2_2 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 12] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 7] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  3] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 3] );

                        //g_0_1_3 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 13] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 9] - aop_PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 9];

                        //g_0_0_4 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 14] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 9] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  5] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 5] );

                    }
}



// VRR to obtain AUX_INT__h_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_h(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__h_s_s_s,
           double const * const restrict AUX_INT__g_s_s_s,
           double const * const restrict AUX_INT__f_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__h_s_s_s[num_m * 21];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //h_5_0_0 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 0] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 0] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 0]
                                      + 4 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  0] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 0] );

                        //h_4_1_0 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 1] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 0] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 0];

                        //h_4_0_1 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 2] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 0] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 0];

                        //h_3_2_0 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 3] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 1] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  0] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 0] );

                        //h_3_1_1 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 4] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 1] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 1];

                        //h_3_0_2 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 5] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 2] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  0] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 0] );

                        //h_2_3_0 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 6] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 6] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 6]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  6] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 6] );

                        //h_2_2_1 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 7] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 3] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 3];

                        //h_2_1_2 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 8] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 5] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 5];

                        //h_2_0_3 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 9] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 9] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 9]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  9] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 9] );

                        //h_1_4_0 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 10] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 10] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 10];

                        //h_1_3_1 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 11] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 6] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 6];

                        //h_1_2_2 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 12] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 12] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 12];

                        //h_1_1_3 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 13] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 9] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 9];

                        //h_1_0_4 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 14] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 14] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 14];

                        //h_0_5_0 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 15] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 10] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  6] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 6] );

                        //h_0_4_1 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 16] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 10] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 10];

                        //h_0_3_2 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 17] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 11] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  6] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 6] );

                        //h_0_2_3 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 18] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 13] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  9] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 9] );

                        //h_0_1_4 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 19] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 14] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 14];

                        //h_0_0_5 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 20] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 14] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  9] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 9] );

                    }
}



// VRR to obtain AUX_INT__i_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_i(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__i_s_s_s,
           double const * const restrict AUX_INT__h_s_s_s,
           double const * const restrict AUX_INT__g_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__i_s_s_s[num_m * 28];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //i_6_0_0 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 0] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 0] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 0]
                                      + 5 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  0] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 0] );

                        //i_5_1_0 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 1] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 0] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 0];

                        //i_5_0_1 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 2] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 0] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 0];

                        //i_4_2_0 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 3] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 1] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  0] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 0] );

                        //i_4_1_1 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 4] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 1] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 1];

                        //i_4_0_2 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 5] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 2] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  0] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 0] );

                        //i_3_3_0 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 6] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 3] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  1] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 1] );

                        //i_3_2_1 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 7] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 3] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 3];

                        //i_3_1_2 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 8] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 5] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 5];

                        //i_3_0_3 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 9] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 5] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  2] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 2] );

                        //i_2_4_0 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 10] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 10] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 10]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  10] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 10] );

                        //i_2_3_1 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 11] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 6] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 6];

                        //i_2_2_2 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 12] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 7] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  3] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 3] );

                        //i_2_1_3 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 13] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 9] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 9];

                        //i_2_0_4 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 14] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 14] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 14]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  14] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 14] );

                        //i_1_5_0 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 15] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 15] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 15];

                        //i_1_4_1 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 16] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 10] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 10];

                        //i_1_3_2 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 17] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 17] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 17];

                        //i_1_2_3 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 18] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 18] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 18];

                        //i_1_1_4 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 19] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 14] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 14];

                        //i_1_0_5 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 20] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 20] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 20];

                        //i_0_6_0 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 21] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 15] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  10] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 10] );

                        //i_0_5_1 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 22] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 15] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 15];

                        //i_0_4_2 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 23] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 16] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  10] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 10] );

                        //i_0_3_3 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 24] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 17] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  11] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 11] );

                        //i_0_2_4 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 25] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 19] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  14] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 14] );

                        //i_0_1_5 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 26] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 20] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 20];

                        //i_0_0_6 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 27] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 20] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  14] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 14] );

                    }
}



// VRR to obtain AUX_INT__j_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_j(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__j_s_s_s,
           double const * const restrict AUX_INT__i_s_s_s,
           double const * const restrict AUX_INT__h_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__j_s_s_s[num_m * 36];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //j_7_0_0 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 0] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 0] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 0]
                                      + 6 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  0] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 0] );

                        //j_6_1_0 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 1] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 0] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 0];

                        //j_6_0_1 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 2] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 0] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 0];

                        //j_5_2_0 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 3] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 1] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  0] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 0] );

                        //j_5_1_1 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 4] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 1] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 1];

                        //j_5_0_2 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 5] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 2] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  0] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 0] );

                        //j_4_3_0 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 6] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 3] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  1] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 1] );

                        //j_4_2_1 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 7] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 3] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 3];

                        //j_4_1_2 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 8] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 5] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 5];

                        //j_4_0_3 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 9] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 5] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  2] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 2] );

                        //j_3_4_0 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 10] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 10] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 10]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  10] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 10] );

                        //j_3_3_1 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 11] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 6] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 6];

                        //j_3_2_2 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 12] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 7] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  3] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 3] );

                        //j_3_1_3 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 13] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 9] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 9];

                        //j_3_0_4 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 14] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 14] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 14]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  14] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 14] );

                        //j_2_5_0 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 15] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 15] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 15]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  15] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 15] );

                        //j_2_4_1 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 16] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 10] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 10];

                        //j_2_3_2 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 17] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 11] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  6] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 6] );

                        //j_2_2_3 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 18] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 13] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  9] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 9] );

                        //j_2_1_4 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 19] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 14] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 14];

                        //j_2_0_5 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 20] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 20] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 20]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  20] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 20] );

                        //j_1_6_0 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 21] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 21] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 21];

                        //j_1_5_1 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 22] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 15] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 15];

                        //j_1_4_2 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 23] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 23] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 23];

                        //j_1_3_3 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 24] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 24] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 24];

                        //j_1_2_4 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 25] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 25] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 25];

                        //j_1_1_5 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 26] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 20] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 20];

                        //j_1_0_6 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 27] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 27] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 27];

                        //j_0_7_0 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 28] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 21] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  15] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 15] );

                        //j_0_6_1 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 29] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 21] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 21];

                        //j_0_5_2 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 30] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 22] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  15] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 15] );

                        //j_0_4_3 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 31] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 23] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  16] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 16] );

                        //j_0_3_4 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 32] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 25] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  19] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 19] );

                        //j_0_2_5 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 33] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 26] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  20] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 20] );

                        //j_0_1_6 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 34] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 27] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 27];

                        //j_0_0_7 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 35] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 27] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  20] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 20] );

                    }
}



// VRR to obtain AUX_INT__k_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_k(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__k_s_s_s,
           double const * const restrict AUX_INT__j_s_s_s,
           double const * const restrict AUX_INT__i_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__k_s_s_s[num_m * 45];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //k_8_0_0 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 0] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 0] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 0]
                                      + 7 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  0] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 0] );

                        //k_7_1_0 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 1] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 0] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 0];

                        //k_7_0_1 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 2] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 0] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 0];

                        //k_6_2_0 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 3] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 1] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  0] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 0] );

                        //k_6_1_1 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 4] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 1] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 1];

                        //k_6_0_2 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 5] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 2] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  0] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 0] );

                        //k_5_3_0 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 6] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 3] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  1] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 1] );

                        //k_5_2_1 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 7] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 3] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 3];

                        //k_5_1_2 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 8] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 5] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 5];

                        //k_5_0_3 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 9] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 5] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  2] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 2] );

                        //k_4_4_0 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 10] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 6] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  3] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 3] );

                        //k_4_3_1 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 11] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 6] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 6];

                        //k_4_2_2 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 12] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 7] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  3] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 3] );

                        //k_4_1_3 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 13] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 9] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 9];

                        //k_4_0_4 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 14] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 9] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  5] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 5] );

                        //k_3_5_0 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 15] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 15] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 15]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  15] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 15] );

                        //k_3_4_1 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 16] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 10] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 10];

                        //k_3_3_2 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 17] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 11] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  6] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 6] );

                        //k_3_2_3 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 18] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 13] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  9] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 9] );

                        //k_3_1_4 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 19] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 14] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 14];

                        //k_3_0_5 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 20] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 20] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 20]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  20] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 20] );

                        //k_2_6_0 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 21] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 21] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 21]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  21] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 21] );

                        //k_2_5_1 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 22] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 15] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 15];

                        //k_2_4_2 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 23] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 16] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  10] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 10] );

                        //k_2_3_3 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 24] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 24] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 24]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  24] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 24] );

                        //k_2_2_4 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 25] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 19] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  14] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 14] );

                        //k_2_1_5 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 26] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 20] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 20];

                        //k_2_0_6 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 27] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 27] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 27]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  27] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 27] );

                        //k_1_7_0 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 28] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 28] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 28];

                        //k_1_6_1 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 29] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 21] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 21];

                        //k_1_5_2 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 30] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 30] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 30];

                        //k_1_4_3 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 31] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 31] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 31];

                        //k_1_3_4 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 32] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 32] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 32];

                        //k_1_2_5 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 33] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 33] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 33];

                        //k_1_1_6 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 34] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 27] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 27];

                        //k_1_0_7 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 35] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 35] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 35];

                        //k_0_8_0 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 36] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 28] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  21] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 21] );

                        //k_0_7_1 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 37] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 28] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 28];

                        //k_0_6_2 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 38] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 29] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  21] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 21] );

                        //k_0_5_3 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 39] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 30] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  22] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 22] );

                        //k_0_4_4 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 40] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 31] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  23] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 23] );

                        //k_0_3_5 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 41] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 33] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  26] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 26] );

                        //k_0_2_6 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 42] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 34] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  27] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 27] );

                        //k_0_1_7 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 43] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 35] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 35];

                        //k_0_0_8 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 44] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 35] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  27] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 27] );

                    }
}



// VRR to obtain AUX_INT__l_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_l(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__l_s_s_s,
           double const * const restrict AUX_INT__k_s_s_s,
           double const * const restrict AUX_INT__j_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__l_s_s_s[num_m * 55];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //l_9_0_0 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 0] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 0] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 0]
                                      + 8 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  0] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 0] );

                        //l_8_1_0 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 1] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 0] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 0];

                        //l_8_0_1 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 2] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 0] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 0];

                        //l_7_2_0 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 3] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 1] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  0] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 0] );

                        //l_7_1_1 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 4] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 1] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 1];

                        //l_7_0_2 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 5] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 2] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  0] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 0] );

                        //l_6_3_0 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 6] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 3] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  1] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 1] );

                        //l_6_2_1 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 7] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 3] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 3];

                        //l_6_1_2 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 8] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 5] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 5];

                        //l_6_0_3 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 9] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 5] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  2] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 2] );

                        //l_5_4_0 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 10] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 6] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  3] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 3] );

                        //l_5_3_1 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 11] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 6] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 6];

                        //l_5_2_2 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 12] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 7] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  3] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 3] );

                        //l_5_1_3 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 13] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 9] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 9];

                        //l_5_0_4 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 14] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 9] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  5] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 5] );

                        //l_4_5_0 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 15] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 15] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 15]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  15] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 15] );

                        //l_4_4_1 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 16] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 10] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 10];

                        //l_4_3_2 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 17] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 11] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  6] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 6] );

                        //l_4_2_3 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 18] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 13] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  9] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 9] );

                        //l_4_1_4 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 19] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 14] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 14];

                        //l_4_0_5 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 20] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 20] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 20]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  20] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 20] );

                        //l_3_6_0 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 21] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 21] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 21]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  21] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 21] );

                        //l_3_5_1 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 22] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 15] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 15];

                        //l_3_4_2 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 23] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 16] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  10] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 10] );

                        //l_3_3_3 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 24] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 17] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  11] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 11] );

                        //l_3_2_4 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 25] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 19] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  14] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 14] );

                        //l_3_1_5 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 26] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 20] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 20];

                        //l_3_0_6 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 27] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 27] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 27]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  27] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 27] );

                        //l_2_7_0 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 28] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 28] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 28]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  28] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 28] );

                        //l_2_6_1 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 29] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 21] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 21];

                        //l_2_5_2 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 30] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 22] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  15] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 15] );

                        //l_2_4_3 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 31] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 31] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 31]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  31] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 31] );

                        //l_2_3_4 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 32] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 32] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 32]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  32] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 32] );

                        //l_2_2_5 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 33] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 26] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  20] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 20] );

                        //l_2_1_6 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 34] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 27] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 27];

                        //l_2_0_7 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 35] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 35] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 35]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  35] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 35] );

                        //l_1_8_0 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 36] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 36] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 36];

                        //l_1_7_1 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 37] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 28] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 28];

                        //l_1_6_2 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 38] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 38] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 38];

                        //l_1_5_3 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 39] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 39] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 39];

                        //l_1_4_4 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 40] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 40] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 40];

                        //l_1_3_5 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 41] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 41] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 41];

                        //l_1_2_6 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 42] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 42] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 42];

                        //l_1_1_7 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 43] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 35] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 35];

                        //l_1_0_8 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 44] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 44] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 44];

                        //l_0_9_0 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 45] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 36] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 36]
                                      + 8 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  28] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 28] );

                        //l_0_8_1 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 46] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 36] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 36];

                        //l_0_7_2 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 47] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 37] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  28] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 28] );

                        //l_0_6_3 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 48] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 38] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  29] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 29] );

                        //l_0_5_4 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 49] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 39] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  30] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 30] );

                        //l_0_4_5 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 50] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 41] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  33] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 33] );

                        //l_0_3_6 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 51] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 42] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  34] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 34] );

                        //l_0_2_7 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 52] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 43] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  35] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 35] );

                        //l_0_1_8 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 53] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 44] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 44];

                        //l_0_0_9 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 54] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 44] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 44]
                                      + 8 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  35] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 35] );

                    }
}



// VRR to obtain AUX_INT__m_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_m(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__m_s_s_s,
           double const * const restrict AUX_INT__l_s_s_s,
           double const * const restrict AUX_INT__k_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__m_s_s_s[num_m * 66];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //m_10_0_0 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 0] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 0] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 0]
                                      + 9 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  0] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 0] );

                        //m_9_1_0 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 1] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 0] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 0];

                        //m_9_0_1 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 2] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 0] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 0];

                        //m_8_2_0 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 3] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 1] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  0] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 0] );

                        //m_8_1_1 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 4] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 1] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 1];

                        //m_8_0_2 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 5] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 2] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  0] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 0] );

                        //m_7_3_0 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 6] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 3] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  1] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 1] );

                        //m_7_2_1 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 7] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 3] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 3];

                        //m_7_1_2 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 8] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 5] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 5];

                        //m_7_0_3 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 9] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 5] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  2] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 2] );

                        //m_6_4_0 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 10] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 6] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  3] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 3] );

                        //m_6_3_1 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 11] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 6] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 6];

                        //m_6_2_2 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 12] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 7] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  3] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 3] );

                        //m_6_1_3 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 13] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 9] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 9];

                        //m_6_0_4 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 14] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 9] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  5] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 5] );

                        //m_5_5_0 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 15] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 10] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  6] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 6] );

                        //m_5_4_1 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 16] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 10] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 10];

                        //m_5_3_2 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 17] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 11] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  6] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 6] );

                        //m_5_2_3 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 18] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 13] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  9] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 9] );

                        //m_5_1_4 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 19] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 14] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 14];

                        //m_5_0_5 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 20] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 14] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  9] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 9] );

                        //m_4_6_0 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 21] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 21] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 21]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  21] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 21] );

                        //m_4_5_1 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 22] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 15] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 15];

                        //m_4_4_2 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 23] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 16] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  10] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 10] );

                        //m_4_3_3 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 24] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 17] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  11] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 11] );

                        //m_4_2_4 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 25] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 19] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  14] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 14] );

                        //m_4_1_5 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 26] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 20] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 20];

                        //m_4_0_6 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 27] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 27] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 27]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  27] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 27] );

                        //m_3_7_0 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 28] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 28] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 28]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  28] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 28] );

                        //m_3_6_1 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 29] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 21] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 21];

                        //m_3_5_2 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 30] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 22] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  15] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 15] );

                        //m_3_4_3 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 31] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 23] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  16] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 16] );

                        //m_3_3_4 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 32] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 25] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  19] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 19] );

                        //m_3_2_5 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 33] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 26] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  20] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 20] );

                        //m_3_1_6 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 34] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 27] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 27];

                        //m_3_0_7 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 35] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 35] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 35]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  35] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 35] );

                        //m_2_8_0 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 36] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 36] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 36]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  36] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 36] );

                        //m_2_7_1 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 37] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 28] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 28];

                        //m_2_6_2 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 38] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 29] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  21] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 21] );

                        //m_2_5_3 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 39] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 39] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 39]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  39] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 39] );

                        //m_2_4_4 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 40] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 40] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 40]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  40] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 40] );

                        //m_2_3_5 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 41] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 41] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 41]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  41] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 41] );

                        //m_2_2_6 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 42] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 34] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  27] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 27] );

                        //m_2_1_7 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 43] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 35] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 35];

                        //m_2_0_8 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 44] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 44] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 44]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  44] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 44] );

                        //m_1_9_0 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 45] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 45] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 45];

                        //m_1_8_1 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 46] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 36] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 36];

                        //m_1_7_2 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 47] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 47] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 47];

                        //m_1_6_3 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 48] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 48] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 48];

                        //m_1_5_4 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 49] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 49] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 49];

                        //m_1_4_5 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 50] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 50] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 50];

                        //m_1_3_6 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 51] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 51] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 51];

                        //m_1_2_7 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 52] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 52] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 52];

                        //m_1_1_8 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 53] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 44] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 44];

                        //m_1_0_9 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 54] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 54] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 54];

                        //m_0_10_0 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 55] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 45] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 45]
                                      + 9 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  36] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 36] );

                        //m_0_9_1 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 56] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 45] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 45];

                        //m_0_8_2 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 57] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 46] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  36] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 36] );

                        //m_0_7_3 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 58] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 47] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  37] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 37] );

                        //m_0_6_4 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 59] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 48] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  38] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 38] );

                        //m_0_5_5 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 60] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 49] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  39] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 39] );

                        //m_0_4_6 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 61] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 51] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  42] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 42] );

                        //m_0_3_7 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 62] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 52] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  43] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 43] );

                        //m_0_2_8 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 63] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 53] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  44] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 44] );

                        //m_0_1_9 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 64] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 54] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 54];

                        //m_0_0_10 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 65] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 54] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 54]
                                      + 9 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  44] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 44] );

                    }
}



// VRR to obtain AUX_INT__n_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_n(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__n_s_s_s,
           double const * const restrict AUX_INT__m_s_s_s,
           double const * const restrict AUX_INT__l_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__n_s_s_s[num_m * 78];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //n_11_0_0 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 0] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 0] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 0]
                                      + 10 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  0] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 0] );

                        //n_10_1_0 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 1] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 0] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 0];

                        //n_10_0_1 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 2] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 0] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 0];

                        //n_9_2_0 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 3] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 1] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  0] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 0] );

                        //n_9_1_1 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 4] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 1] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 1];

                        //n_9_0_2 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 5] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 2] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  0] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 0] );

                        //n_8_3_0 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 6] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 3] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  1] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 1] );

                        //n_8_2_1 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 7] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 3] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 3];

                        //n_8_1_2 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 8] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 5] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 5];

                        //n_8_0_3 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 9] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 5] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  2] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 2] );

                        //n_7_4_0 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 10] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 6] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  3] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 3] );

                        //n_7_3_1 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 11] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 6] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 6];

                        //n_7_2_2 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 12] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 7] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  3] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 3] );

                        //n_7_1_3 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 13] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 9] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 9];

                        //n_7_0_4 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 14] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 9] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  5] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 5] );

                        //n_6_5_0 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 15] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 10] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  6] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 6] );

                        //n_6_4_1 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 16] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 10] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 10];

                        //n_6_3_2 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 17] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 11] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  6] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 6] );

                        //n_6_2_3 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 18] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 13] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  9] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 9] );

                        //n_6_1_4 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 19] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 14] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 14];

                        //n_6_0_5 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 20] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 14] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  9] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 9] );

                        //n_5_6_0 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 21] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 21] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 21]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  21] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 21] );

                        //n_5_5_1 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 22] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 15] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 15];

                        //n_5_4_2 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 23] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 16] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  10] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 10] );

                        //n_5_3_3 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 24] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 17] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  11] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 11] );

                        //n_5_2_4 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 25] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 19] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  14] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 14] );

                        //n_5_1_5 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 26] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 20] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 20];

                        //n_5_0_6 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 27] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 27] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 27]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  27] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 27] );

                        //n_4_7_0 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 28] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 28] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 28]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  28] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 28] );

                        //n_4_6_1 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 29] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 21] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 21];

                        //n_4_5_2 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 30] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 22] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  15] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 15] );

                        //n_4_4_3 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 31] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 23] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  16] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 16] );

                        //n_4_3_4 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 32] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 25] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  19] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 19] );

                        //n_4_2_5 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 33] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 26] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  20] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 20] );

                        //n_4_1_6 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 34] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 27] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 27];

                        //n_4_0_7 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 35] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 35] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 35]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  35] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 35] );

                        //n_3_8_0 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 36] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 36] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 36]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  36] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 36] );

                        //n_3_7_1 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 37] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 28] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 28];

                        //n_3_6_2 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 38] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 29] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  21] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 21] );

                        //n_3_5_3 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 39] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 30] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  22] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 22] );

                        //n_3_4_4 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 40] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 40] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 40]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  40] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 40] );

                        //n_3_3_5 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 41] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 33] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  26] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 26] );

                        //n_3_2_6 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 42] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 34] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  27] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 27] );

                        //n_3_1_7 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 43] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 35] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 35];

                        //n_3_0_8 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 44] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 44] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 44]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  44] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 44] );

                        //n_2_9_0 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 45] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 45] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 45]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  45] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 45] );

                        //n_2_8_1 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 46] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 36] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 36];

                        //n_2_7_2 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 47] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 37] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  28] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 28] );

                        //n_2_6_3 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 48] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 48] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 48]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  48] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 48] );

                        //n_2_5_4 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 49] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 49] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 49]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  49] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 49] );

                        //n_2_4_5 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 50] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 50] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 50]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  50] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 50] );

                        //n_2_3_6 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 51] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 51] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 51]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  51] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 51] );

                        //n_2_2_7 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 52] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 43] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  35] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 35] );

                        //n_2_1_8 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 53] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 44] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 44];

                        //n_2_0_9 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 54] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 54] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 54]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  54] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 54] );

                        //n_1_10_0 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 55] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 55] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 55];

                        //n_1_9_1 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 56] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 45] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 45];

                        //n_1_8_2 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 57] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 57] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 57];

                        //n_1_7_3 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 58] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 58] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 58];

                        //n_1_6_4 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 59] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 59] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 59];

                        //n_1_5_5 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 60] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 60] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 60];

                        //n_1_4_6 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 61] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 61] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 61];

                        //n_1_3_7 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 62] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 62] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 62];

                        //n_1_2_8 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 63] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 63] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 63];

                        //n_1_1_9 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 64] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 54] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 54];

                        //n_1_0_10 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 65] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 65] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 65];

                        //n_0_11_0 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 66] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 55] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 55]
                                      + 10 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  45] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 45] );

                        //n_0_10_1 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 67] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 55] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 55];

                        //n_0_9_2 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 68] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 56] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  45] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 45] );

                        //n_0_8_3 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 69] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 57] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  46] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 46] );

                        //n_0_7_4 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 70] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 58] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  47] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 47] );

                        //n_0_6_5 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 71] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 59] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 59]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  48] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 48] );

                        //n_0_5_6 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 72] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 61] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 61]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  51] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 51] );

                        //n_0_4_7 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 73] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 62] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  52] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 52] );

                        //n_0_3_8 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 74] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 63] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  53] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 53] );

                        //n_0_2_9 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 75] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 64] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  54] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 54] );

                        //n_0_1_10 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 76] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 65] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 65];

                        //n_0_0_11 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 77] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 65] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 65]
                                      + 10 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  54] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 54] );

                    }
}



// VRR to obtain AUX_INT__o_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_o(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__o_s_s_s,
           double const * const restrict AUX_INT__n_s_s_s,
           double const * const restrict AUX_INT__m_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__o_s_s_s[num_m * 91];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //o_12_0_0 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 0] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 0] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 0]
                                      + 11 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  0] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 0] );

                        //o_11_1_0 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 1] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 0] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 0];

                        //o_11_0_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 2] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 0] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 0];

                        //o_10_2_0 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 3] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 1] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  0] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 0] );

                        //o_10_1_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 4] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 1] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 1];

                        //o_10_0_2 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 5] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 2] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  0] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 0] );

                        //o_9_3_0 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 6] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 3] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  1] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 1] );

                        //o_9_2_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 7] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 3] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 3];

                        //o_9_1_2 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 8] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 5] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 5];

                        //o_9_0_3 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 9] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 5] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  2] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 2] );

                        //o_8_4_0 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 10] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 6] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  3] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 3] );

                        //o_8_3_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 11] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 6] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 6];

                        //o_8_2_2 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 12] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 7] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  3] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 3] );

                        //o_8_1_3 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 13] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 9] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 9];

                        //o_8_0_4 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 14] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 9] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  5] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 5] );

                        //o_7_5_0 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 15] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 10] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  6] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 6] );

                        //o_7_4_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 16] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 10] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 10];

                        //o_7_3_2 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 17] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 11] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  6] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 6] );

                        //o_7_2_3 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 18] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 13] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  9] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 9] );

                        //o_7_1_4 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 19] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 14] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 14];

                        //o_7_0_5 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 20] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 14] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  9] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 9] );

                        //o_6_6_0 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 21] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 15] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  10] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 10] );

                        //o_6_5_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 22] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 15] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 15];

                        //o_6_4_2 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 23] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 16] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  10] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 10] );

                        //o_6_3_3 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 24] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 17] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  11] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 11] );

                        //o_6_2_4 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 25] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 19] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  14] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 14] );

                        //o_6_1_5 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 26] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 20] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 20];

                        //o_6_0_6 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 27] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 20] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  14] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 14] );

                        //o_5_7_0 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 28] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 28] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 28]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  28] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 28] );

                        //o_5_6_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 29] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 21] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 21];

                        //o_5_5_2 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 30] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 22] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  15] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 15] );

                        //o_5_4_3 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 31] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 23] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  16] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 16] );

                        //o_5_3_4 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 32] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 25] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  19] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 19] );

                        //o_5_2_5 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 33] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 26] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  20] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 20] );

                        //o_5_1_6 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 34] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 27] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 27];

                        //o_5_0_7 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 35] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 35] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 35]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  35] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 35] );

                        //o_4_8_0 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 36] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 36] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 36]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  36] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 36] );

                        //o_4_7_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 37] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 28] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 28];

                        //o_4_6_2 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 38] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 29] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  21] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 21] );

                        //o_4_5_3 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 39] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 30] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  22] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 22] );

                        //o_4_4_4 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 40] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 31] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  23] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 23] );

                        //o_4_3_5 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 41] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 33] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  26] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 26] );

                        //o_4_2_6 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 42] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 34] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  27] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 27] );

                        //o_4_1_7 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 43] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 35] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 35];

                        //o_4_0_8 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 44] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 44] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 44]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  44] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 44] );

                        //o_3_9_0 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 45] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 45] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 45]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  45] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 45] );

                        //o_3_8_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 46] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 36] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 36];

                        //o_3_7_2 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 47] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 37] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  28] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 28] );

                        //o_3_6_3 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 48] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 38] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  29] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 29] );

                        //o_3_5_4 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 49] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 49] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 49]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  49] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 49] );

                        //o_3_4_5 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 50] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 50] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 50]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  50] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 50] );

                        //o_3_3_6 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 51] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 42] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  34] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 34] );

                        //o_3_2_7 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 52] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 43] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  35] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 35] );

                        //o_3_1_8 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 53] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 44] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 44];

                        //o_3_0_9 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 54] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 54] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 54]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  54] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 54] );

                        //o_2_10_0 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 55] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 55] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 55]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  55] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 55] );

                        //o_2_9_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 56] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 45] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 45];

                        //o_2_8_2 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 57] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 46] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  36] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 36] );

                        //o_2_7_3 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 58] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 58] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 58]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  58] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 58] );

                        //o_2_6_4 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 59] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 59] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 59]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  59] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 59] );

                        //o_2_5_5 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 60] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 60] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 60]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  60] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 60] );

                        //o_2_4_6 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 61] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 61] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 61]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  61] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 61] );

                        //o_2_3_7 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 62] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 62] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 62]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  62] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 62] );

                        //o_2_2_8 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 63] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 53] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  44] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 44] );

                        //o_2_1_9 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 64] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 54] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 54];

                        //o_2_0_10 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 65] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 65] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 65]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  65] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 65] );

                        //o_1_11_0 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 66] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 66] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 66];

                        //o_1_10_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 67] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 55] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 55];

                        //o_1_9_2 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 68] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 68] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 68];

                        //o_1_8_3 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 69] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 69] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 69];

                        //o_1_7_4 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 70] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 70] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 70];

                        //o_1_6_5 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 71] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 71] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 71];

                        //o_1_5_6 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 72] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 72] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 72];

                        //o_1_4_7 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 73] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 73] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 73];

                        //o_1_3_8 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 74] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 74] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 74];

                        //o_1_2_9 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 75] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 75] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 75];

                        //o_1_1_10 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 76] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 65] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 65];

                        //o_1_0_11 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 77] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 77] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 77];

                        //o_0_12_0 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 78] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 66] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 66]
                                      + 11 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  55] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 55] );

                        //o_0_11_1 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 79] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 66] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 66];

                        //o_0_10_2 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 80] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 67] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  55] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 55] );

                        //o_0_9_3 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 81] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 68] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  56] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 56] );

                        //o_0_8_4 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 82] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 69] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 69]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  57] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 57] );

                        //o_0_7_5 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 83] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 70] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 70]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  58] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 58] );

                        //o_0_6_6 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 84] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 71] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 71]
                                      + 5 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  59] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 59] );

                        //o_0_5_7 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 85] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 73] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 73]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  62] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 62] );

                        //o_0_4_8 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 86] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 74] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 74]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  63] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 63] );

                        //o_0_3_9 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 87] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 75] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  64] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 64] );

                        //o_0_2_10 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 88] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 76] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  65] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 65] );

                        //o_0_1_11 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 89] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 77] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 77];

                        //o_0_0_12 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 90] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 77] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 77]
                                      + 11 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  65] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 65] );

                    }
}



// VRR to obtain AUX_INT__q_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_q(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__q_s_s_s,
           double const * const restrict AUX_INT__o_s_s_s,
           double const * const restrict AUX_INT__n_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__q_s_s_s[num_m * 105];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //q_13_0_0 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 0] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 0] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 0]
                                      + 12 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  0] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 0] );

                        //q_12_1_0 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 1] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 0] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 0];

                        //q_12_0_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 2] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 0] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 0];

                        //q_11_2_0 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 3] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 1] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  0] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 0] );

                        //q_11_1_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 4] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 1] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 1];

                        //q_11_0_2 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 5] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 2] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  0] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 0] );

                        //q_10_3_0 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 6] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 3] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  1] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 1] );

                        //q_10_2_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 7] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 3] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 3];

                        //q_10_1_2 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 8] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 5] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 5];

                        //q_10_0_3 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 9] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 5] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  2] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 2] );

                        //q_9_4_0 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 10] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 6] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  3] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 3] );

                        //q_9_3_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 11] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 6] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 6];

                        //q_9_2_2 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 12] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 7] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  3] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 3] );

                        //q_9_1_3 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 13] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 9] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 9];

                        //q_9_0_4 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 14] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 9] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  5] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 5] );

                        //q_8_5_0 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 15] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 10] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  6] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 6] );

                        //q_8_4_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 16] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 10] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 10];

                        //q_8_3_2 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 17] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 11] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  6] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 6] );

                        //q_8_2_3 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 18] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 13] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  9] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 9] );

                        //q_8_1_4 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 19] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 14] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 14];

                        //q_8_0_5 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 20] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 14] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  9] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 9] );

                        //q_7_6_0 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 21] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 15] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  10] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 10] );

                        //q_7_5_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 22] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 15] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 15];

                        //q_7_4_2 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 23] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 16] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  10] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 10] );

                        //q_7_3_3 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 24] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 17] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  11] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 11] );

                        //q_7_2_4 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 25] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 19] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  14] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 14] );

                        //q_7_1_5 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 26] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 20] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 20];

                        //q_7_0_6 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 27] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 20] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  14] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 14] );

                        //q_6_7_0 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 28] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 28] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 28]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  28] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 28] );

                        //q_6_6_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 29] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 21] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 21];

                        //q_6_5_2 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 30] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 22] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  15] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 15] );

                        //q_6_4_3 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 31] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 23] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  16] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 16] );

                        //q_6_3_4 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 32] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 25] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  19] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 19] );

                        //q_6_2_5 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 33] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 26] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  20] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 20] );

                        //q_6_1_6 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 34] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 27] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 27];

                        //q_6_0_7 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 35] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 35] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 35]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  35] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 35] );

                        //q_5_8_0 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 36] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 36] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 36]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  36] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 36] );

                        //q_5_7_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 37] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 28] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 28];

                        //q_5_6_2 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 38] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 29] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  21] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 21] );

                        //q_5_5_3 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 39] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 30] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  22] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 22] );

                        //q_5_4_4 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 40] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 31] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  23] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 23] );

                        //q_5_3_5 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 41] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 33] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  26] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 26] );

                        //q_5_2_6 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 42] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 34] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  27] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 27] );

                        //q_5_1_7 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 43] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 35] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 35];

                        //q_5_0_8 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 44] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 44] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 44]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  44] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 44] );

                        //q_4_9_0 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 45] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 45] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 45]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  45] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 45] );

                        //q_4_8_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 46] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 36] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 36];

                        //q_4_7_2 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 47] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 37] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  28] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 28] );

                        //q_4_6_3 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 48] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 38] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  29] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 29] );

                        //q_4_5_4 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 49] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 39] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  30] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 30] );

                        //q_4_4_5 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 50] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 41] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  33] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 33] );

                        //q_4_3_6 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 51] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 42] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  34] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 34] );

                        //q_4_2_7 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 52] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 43] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  35] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 35] );

                        //q_4_1_8 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 53] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 44] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 44];

                        //q_4_0_9 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 54] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 54] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 54]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  54] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 54] );

                        //q_3_10_0 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 55] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 55] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 55]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  55] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 55] );

                        //q_3_9_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 56] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 45] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 45];

                        //q_3_8_2 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 57] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 46] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  36] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 36] );

                        //q_3_7_3 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 58] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 47] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  37] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 37] );

                        //q_3_6_4 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 59] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 59] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 59]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  59] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 59] );

                        //q_3_5_5 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 60] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 60] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 60]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  60] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 60] );

                        //q_3_4_6 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 61] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 61] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 61]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  61] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 61] );

                        //q_3_3_7 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 62] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 52] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  43] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 43] );

                        //q_3_2_8 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 63] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 53] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  44] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 44] );

                        //q_3_1_9 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 64] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 54] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 54];

                        //q_3_0_10 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 65] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 65] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 65]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  65] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 65] );

                        //q_2_11_0 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 66] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 66] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 66]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  66] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 66] );

                        //q_2_10_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 67] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 55] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 55];

                        //q_2_9_2 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 68] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 56] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  45] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 45] );

                        //q_2_8_3 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 69] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 69] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 69]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  69] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 69] );

                        //q_2_7_4 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 70] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 70] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 70]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  70] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 70] );

                        //q_2_6_5 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 71] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 71] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 71]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  71] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 71] );

                        //q_2_5_6 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 72] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 72] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 72]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  72] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 72] );

                        //q_2_4_7 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 73] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 73] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 73]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  73] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 73] );

                        //q_2_3_8 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 74] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 74] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 74]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  74] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 74] );

                        //q_2_2_9 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 75] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 64] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  54] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 54] );

                        //q_2_1_10 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 76] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 65] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 65];

                        //q_2_0_11 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 77] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 77] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 77]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  77] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 77] );

                        //q_1_12_0 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 78] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 78] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 78];

                        //q_1_11_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 79] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 66] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 66];

                        //q_1_10_2 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 80] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 80] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 80];

                        //q_1_9_3 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 81] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 81] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 81];

                        //q_1_8_4 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 82] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 82] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 82];

                        //q_1_7_5 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 83] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 83] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 83];

                        //q_1_6_6 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 84] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 84] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 84];

                        //q_1_5_7 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 85] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 85] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 85];

                        //q_1_4_8 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 86] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 86] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 86];

                        //q_1_3_9 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 87] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 87] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 87];

                        //q_1_2_10 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 88] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 88] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 88];

                        //q_1_1_11 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 89] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 77] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 77];

                        //q_1_0_12 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 90] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 90] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 90];

                        //q_0_13_0 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 91] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 78] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 78]
                                      + 12 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  66] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 66] );

                        //q_0_12_1 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 92] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 78] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 78];

                        //q_0_11_2 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 93] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 79] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  66] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 66] );

                        //q_0_10_3 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 94] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 80] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 80]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  67] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 67] );

                        //q_0_9_4 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 95] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 81] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 81]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  68] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 68] );

                        //q_0_8_5 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 96] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 82] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 82]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  69] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 69] );

                        //q_0_7_6 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 97] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 83] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 83]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  70] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 70] );

                        //q_0_6_7 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 98] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 85] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 85]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  73] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 73] );

                        //q_0_5_8 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 99] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 86] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 86]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  74] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 74] );

                        //q_0_4_9 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 100] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 87] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 87]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  75] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 75] );

                        //q_0_3_10 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 101] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 88] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 88]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  76] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 76] );

                        //q_0_2_11 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 102] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 89] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  77] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 77] );

                        //q_0_1_12 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 103] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 90] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 90];

                        //q_0_0_13 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 104] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 90] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 90]
                                      + 12 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  77] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 77] );

                    }
}



// VRR to obtain AUX_INT__r_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_r(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__r_s_s_s,
           double const * const restrict AUX_INT__q_s_s_s,
           double const * const restrict AUX_INT__o_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__r_s_s_s[num_m * 120];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //r_14_0_0 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 0] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 0] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 0]
                                      + 13 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  0] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 0] );

                        //r_13_1_0 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 1] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 0] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 0];

                        //r_13_0_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 2] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 0] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 0];

                        //r_12_2_0 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 3] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 1] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  0] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 0] );

                        //r_12_1_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 4] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 1] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 1];

                        //r_12_0_2 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 5] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 2] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  0] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 0] );

                        //r_11_3_0 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 6] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 3] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  1] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 1] );

                        //r_11_2_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 7] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 3] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 3];

                        //r_11_1_2 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 8] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 5] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 5];

                        //r_11_0_3 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 9] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 5] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  2] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 2] );

                        //r_10_4_0 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 10] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 6] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  3] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 3] );

                        //r_10_3_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 11] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 6] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 6];

                        //r_10_2_2 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 12] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 7] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  3] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 3] );

                        //r_10_1_3 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 13] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 9] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 9];

                        //r_10_0_4 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 14] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 9] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  5] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 5] );

                        //r_9_5_0 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 15] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 10] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  6] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 6] );

                        //r_9_4_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 16] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 10] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 10];

                        //r_9_3_2 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 17] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 11] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  6] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 6] );

                        //r_9_2_3 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 18] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 13] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  9] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 9] );

                        //r_9_1_4 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 19] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 14] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 14];

                        //r_9_0_5 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 20] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 14] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  9] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 9] );

                        //r_8_6_0 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 21] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 15] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  10] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 10] );

                        //r_8_5_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 22] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 15] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 15];

                        //r_8_4_2 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 23] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 16] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  10] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 10] );

                        //r_8_3_3 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 24] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 17] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  11] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 11] );

                        //r_8_2_4 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 25] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 19] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  14] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 14] );

                        //r_8_1_5 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 26] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 20] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 20];

                        //r_8_0_6 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 27] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 20] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  14] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 14] );

                        //r_7_7_0 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 28] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 21] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  15] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 15] );

                        //r_7_6_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 29] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 21] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 21];

                        //r_7_5_2 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 30] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 22] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  15] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 15] );

                        //r_7_4_3 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 31] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 23] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  16] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 16] );

                        //r_7_3_4 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 32] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 25] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  19] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 19] );

                        //r_7_2_5 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 33] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 26] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  20] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 20] );

                        //r_7_1_6 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 34] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 27] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 27];

                        //r_7_0_7 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 35] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 27] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  20] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 20] );

                        //r_6_8_0 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 36] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 36] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 36]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  36] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 36] );

                        //r_6_7_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 37] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 28] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 28];

                        //r_6_6_2 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 38] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 29] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  21] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 21] );

                        //r_6_5_3 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 39] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 30] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  22] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 22] );

                        //r_6_4_4 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 40] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 31] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  23] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 23] );

                        //r_6_3_5 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 41] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 33] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  26] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 26] );

                        //r_6_2_6 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 42] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 34] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  27] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 27] );

                        //r_6_1_7 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 43] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 35] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 35];

                        //r_6_0_8 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 44] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 44] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 44]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  44] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 44] );

                        //r_5_9_0 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 45] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 45] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 45]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  45] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 45] );

                        //r_5_8_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 46] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 36] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 36];

                        //r_5_7_2 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 47] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 37] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  28] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 28] );

                        //r_5_6_3 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 48] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 38] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  29] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 29] );

                        //r_5_5_4 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 49] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 39] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  30] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 30] );

                        //r_5_4_5 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 50] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 41] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  33] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 33] );

                        //r_5_3_6 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 51] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 42] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  34] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 34] );

                        //r_5_2_7 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 52] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 43] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  35] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 35] );

                        //r_5_1_8 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 53] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 44] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 44];

                        //r_5_0_9 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 54] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 54] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 54]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  54] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 54] );

                        //r_4_10_0 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 55] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 55] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 55]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  55] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 55] );

                        //r_4_9_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 56] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 45] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 45];

                        //r_4_8_2 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 57] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 46] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  36] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 36] );

                        //r_4_7_3 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 58] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 47] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  37] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 37] );

                        //r_4_6_4 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 59] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 48] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  38] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 38] );

                        //r_4_5_5 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 60] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 60] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 60]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  60] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 60] );

                        //r_4_4_6 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 61] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 51] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  42] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 42] );

                        //r_4_3_7 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 62] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 52] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  43] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 43] );

                        //r_4_2_8 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 63] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 53] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  44] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 44] );

                        //r_4_1_9 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 64] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 54] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 54];

                        //r_4_0_10 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 65] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 65] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 65]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  65] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 65] );

                        //r_3_11_0 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 66] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 66] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 66]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  66] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 66] );

                        //r_3_10_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 67] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 55] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 55];

                        //r_3_9_2 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 68] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 56] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  45] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 45] );

                        //r_3_8_3 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 69] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 57] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  46] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 46] );

                        //r_3_7_4 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 70] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 70] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 70]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  70] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 70] );

                        //r_3_6_5 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 71] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 71] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 71]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  71] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 71] );

                        //r_3_5_6 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 72] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 72] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 72]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  72] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 72] );

                        //r_3_4_7 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 73] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 73] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 73]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  73] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 73] );

                        //r_3_3_8 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 74] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 63] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  53] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 53] );

                        //r_3_2_9 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 75] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 64] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  54] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 54] );

                        //r_3_1_10 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 76] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 65] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 65];

                        //r_3_0_11 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 77] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 77] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 77]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  77] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 77] );

                        //r_2_12_0 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 78] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 78] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 78]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  78] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 78] );

                        //r_2_11_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 79] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 66] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 66];

                        //r_2_10_2 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 80] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 67] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  55] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 55] );

                        //r_2_9_3 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 81] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 81] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 81]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  81] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 81] );

                        //r_2_8_4 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 82] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 82] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 82]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  82] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 82] );

                        //r_2_7_5 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 83] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 83] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 83]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  83] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 83] );

                        //r_2_6_6 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 84] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 84] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 84]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  84] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 84] );

                        //r_2_5_7 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 85] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 85] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 85]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  85] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 85] );

                        //r_2_4_8 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 86] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 86] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 86]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  86] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 86] );

                        //r_2_3_9 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 87] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 87] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 87]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  87] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 87] );

                        //r_2_2_10 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 88] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 76] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  65] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 65] );

                        //r_2_1_11 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 89] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 77] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 77];

                        //r_2_0_12 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 90] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 90] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 90]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  90] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 90] );

                        //r_1_13_0 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 91] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 91] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 91];

                        //r_1_12_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 92] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 78] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 78];

                        //r_1_11_2 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 93] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 93] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 93];

                        //r_1_10_3 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 94] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 94] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 94];

                        //r_1_9_4 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 95] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 95] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 95];

                        //r_1_8_5 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 96] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 96] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 96];

                        //r_1_7_6 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 97] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 97] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 97];

                        //r_1_6_7 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 98] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 98] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 98];

                        //r_1_5_8 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 99] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 99] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 99];

                        //r_1_4_9 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 100] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 100] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 100];

                        //r_1_3_10 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 101] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 101] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 101];

                        //r_1_2_11 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 102] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 102] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 102];

                        //r_1_1_12 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 103] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 90] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 90];

                        //r_1_0_13 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 104] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 104] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 104];

                        //r_0_14_0 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 105] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 91] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 91]
                                      + 13 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  78] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 78] );

                        //r_0_13_1 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 106] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 91] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 91];

                        //r_0_12_2 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 107] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 92] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 92]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  78] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 78] );

                        //r_0_11_3 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 108] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 93] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 93]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  79] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 79] );

                        //r_0_10_4 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 109] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 94] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 94]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  80] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 80] );

                        //r_0_9_5 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 110] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 95] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 95]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  81] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 81] );

                        //r_0_8_6 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 111] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 96] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 96]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  82] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 82] );

                        //r_0_7_7 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 112] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 97] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 97]
                                      + 6 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  83] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 83] );

                        //r_0_6_8 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 113] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 99] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 99]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  86] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 86] );

                        //r_0_5_9 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 114] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 100] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 100]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  87] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 87] );

                        //r_0_4_10 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 115] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 101] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 101]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  88] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 88] );

                        //r_0_3_11 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 116] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 102] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 102]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  89] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 89] );

                        //r_0_2_12 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 117] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 103] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 103]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  90] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 90] );

                        //r_0_1_13 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 118] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 104] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 104];

                        //r_0_0_14 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 119] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 104] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 104]
                                      + 13 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  90] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 90] );

                    }
}



// VRR to obtain AUX_INT__t_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_t(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__t_s_s_s,
           double const * const restrict AUX_INT__r_s_s_s,
           double const * const restrict AUX_INT__q_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__t_s_s_s[num_m * 136];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //t_15_0_0 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 0] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 0] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 0]
                                      + 14 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  0] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 0] );

                        //t_14_1_0 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 1] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 0] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 0];

                        //t_14_0_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 2] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 0] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 0];

                        //t_13_2_0 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 3] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 1] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  0] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 0] );

                        //t_13_1_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 4] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 1] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 1];

                        //t_13_0_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 5] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 2] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  0] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 0] );

                        //t_12_3_0 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 6] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 3] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  1] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 1] );

                        //t_12_2_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 7] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 3] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 3];

                        //t_12_1_2 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 8] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 5] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 5];

                        //t_12_0_3 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 9] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 5] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  2] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 2] );

                        //t_11_4_0 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 10] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 6] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  3] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 3] );

                        //t_11_3_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 11] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 6] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 6];

                        //t_11_2_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 12] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 7] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  3] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 3] );

                        //t_11_1_3 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 13] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 9] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 9];

                        //t_11_0_4 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 14] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 9] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  5] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 5] );

                        //t_10_5_0 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 15] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 10] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  6] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 6] );

                        //t_10_4_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 16] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 10] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 10];

                        //t_10_3_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 17] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 11] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  6] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 6] );

                        //t_10_2_3 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 18] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 13] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  9] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 9] );

                        //t_10_1_4 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 19] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 14] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 14];

                        //t_10_0_5 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 20] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 14] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  9] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 9] );

                        //t_9_6_0 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 21] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 15] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  10] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 10] );

                        //t_9_5_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 22] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 15] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 15];

                        //t_9_4_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 23] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 16] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  10] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 10] );

                        //t_9_3_3 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 24] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 17] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  11] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 11] );

                        //t_9_2_4 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 25] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 19] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  14] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 14] );

                        //t_9_1_5 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 26] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 20] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 20];

                        //t_9_0_6 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 27] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 20] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  14] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 14] );

                        //t_8_7_0 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 28] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 21] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  15] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 15] );

                        //t_8_6_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 29] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 21] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 21];

                        //t_8_5_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 30] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 22] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  15] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 15] );

                        //t_8_4_3 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 31] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 23] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  16] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 16] );

                        //t_8_3_4 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 32] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 25] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  19] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 19] );

                        //t_8_2_5 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 33] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 26] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  20] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 20] );

                        //t_8_1_6 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 34] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 27] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 27];

                        //t_8_0_7 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 35] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 27] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  20] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 20] );

                        //t_7_8_0 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 36] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 36] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 36]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  36] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 36] );

                        //t_7_7_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 37] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 28] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 28];

                        //t_7_6_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 38] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 29] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  21] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 21] );

                        //t_7_5_3 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 39] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 30] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  22] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 22] );

                        //t_7_4_4 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 40] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 31] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  23] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 23] );

                        //t_7_3_5 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 41] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 33] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  26] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 26] );

                        //t_7_2_6 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 42] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 34] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  27] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 27] );

                        //t_7_1_7 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 43] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 35] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 35];

                        //t_7_0_8 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 44] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 44] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 44]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  44] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 44] );

                        //t_6_9_0 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 45] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 45] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 45]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  45] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 45] );

                        //t_6_8_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 46] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 36] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 36];

                        //t_6_7_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 47] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 37] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  28] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 28] );

                        //t_6_6_3 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 48] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 38] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  29] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 29] );

                        //t_6_5_4 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 49] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 39] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  30] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 30] );

                        //t_6_4_5 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 50] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 41] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  33] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 33] );

                        //t_6_3_6 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 51] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 42] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  34] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 34] );

                        //t_6_2_7 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 52] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 43] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  35] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 35] );

                        //t_6_1_8 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 53] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 44] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 44];

                        //t_6_0_9 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 54] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 54] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 54]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  54] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 54] );

                        //t_5_10_0 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 55] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 55] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 55]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  55] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 55] );

                        //t_5_9_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 56] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 45] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 45];

                        //t_5_8_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 57] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 46] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  36] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 36] );

                        //t_5_7_3 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 58] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 47] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  37] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 37] );

                        //t_5_6_4 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 59] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 48] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  38] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 38] );

                        //t_5_5_5 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 60] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 49] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  39] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 39] );

                        //t_5_4_6 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 61] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 51] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  42] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 42] );

                        //t_5_3_7 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 62] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 52] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  43] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 43] );

                        //t_5_2_8 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 63] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 53] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  44] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 44] );

                        //t_5_1_9 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 64] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 54] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 54];

                        //t_5_0_10 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 65] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 65] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 65]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  65] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 65] );

                        //t_4_11_0 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 66] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 66] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 66]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  66] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 66] );

                        //t_4_10_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 67] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 55] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 55];

                        //t_4_9_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 68] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 56] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  45] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 45] );

                        //t_4_8_3 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 69] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 57] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  46] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 46] );

                        //t_4_7_4 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 70] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 58] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  47] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 47] );

                        //t_4_6_5 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 71] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 71] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 71]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  71] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 71] );

                        //t_4_5_6 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 72] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 72] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 72]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  72] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 72] );

                        //t_4_4_7 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 73] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 62] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  52] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 52] );

                        //t_4_3_8 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 74] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 63] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  53] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 53] );

                        //t_4_2_9 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 75] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 64] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  54] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 54] );

                        //t_4_1_10 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 76] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 65] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 65];

                        //t_4_0_11 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 77] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 77] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 77]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  77] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 77] );

                        //t_3_12_0 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 78] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 78] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 78]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  78] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 78] );

                        //t_3_11_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 79] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 66] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 66];

                        //t_3_10_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 80] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 67] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  55] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 55] );

                        //t_3_9_3 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 81] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 68] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  56] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 56] );

                        //t_3_8_4 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 82] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 82] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 82]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  82] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 82] );

                        //t_3_7_5 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 83] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 83] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 83]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  83] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 83] );

                        //t_3_6_6 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 84] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 84] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 84]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  84] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 84] );

                        //t_3_5_7 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 85] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 85] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 85]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  85] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 85] );

                        //t_3_4_8 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 86] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 86] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 86]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  86] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 86] );

                        //t_3_3_9 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 87] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 75] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  64] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 64] );

                        //t_3_2_10 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 88] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 76] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  65] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 65] );

                        //t_3_1_11 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 89] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 77] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 77];

                        //t_3_0_12 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 90] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 90] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 90]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  90] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 90] );

                        //t_2_13_0 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 91] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 91] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 91]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  91] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 91] );

                        //t_2_12_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 92] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 78] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 78];

                        //t_2_11_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 93] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 79] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  66] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 66] );

                        //t_2_10_3 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 94] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 94] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 94]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  94] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 94] );

                        //t_2_9_4 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 95] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 95] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 95]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  95] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 95] );

                        //t_2_8_5 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 96] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 96] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 96]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  96] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 96] );

                        //t_2_7_6 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 97] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 97] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 97]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  97] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 97] );

                        //t_2_6_7 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 98] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 98] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 98]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  98] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 98] );

                        //t_2_5_8 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 99] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 99] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 99]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  99] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 99] );

                        //t_2_4_9 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 100] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 100] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 100]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  100] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 100] );

                        //t_2_3_10 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 101] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 101] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 101]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  101] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 101] );

                        //t_2_2_11 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 102] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 89] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  77] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 77] );

                        //t_2_1_12 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 103] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 90] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 90];

                        //t_2_0_13 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 104] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 104] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 104]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  104] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 104] );

                        //t_1_14_0 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 105] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 105] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 105];

                        //t_1_13_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 106] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 91] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 91];

                        //t_1_12_2 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 107] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 107] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 107];

                        //t_1_11_3 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 108] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 108] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 108];

                        //t_1_10_4 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 109] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 109] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 109];

                        //t_1_9_5 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 110] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 110] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 110];

                        //t_1_8_6 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 111] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 111] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 111];

                        //t_1_7_7 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 112] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 112] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 112];

                        //t_1_6_8 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 113] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 113] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 113];

                        //t_1_5_9 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 114] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 114] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 114];

                        //t_1_4_10 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 115] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 115] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 115];

                        //t_1_3_11 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 116] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 116] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 116];

                        //t_1_2_12 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 117] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 117] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 117];

                        //t_1_1_13 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 118] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 104] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 104];

                        //t_1_0_14 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 119] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 119] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 119];

                        //t_0_15_0 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 120] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 105] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 105]
                                      + 14 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  91] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 91] );

                        //t_0_14_1 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 121] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 105] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 105];

                        //t_0_13_2 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 122] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 106] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 106]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  91] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 91] );

                        //t_0_12_3 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 123] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 107] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 107]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  92] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 92] );

                        //t_0_11_4 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 124] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 108] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 108]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  93] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 93] );

                        //t_0_10_5 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 125] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 109] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 109]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  94] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 94] );

                        //t_0_9_6 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 126] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 110] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 110]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  95] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 95] );

                        //t_0_8_7 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 127] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 111] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 111]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  96] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 96] );

                        //t_0_7_8 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 128] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 113] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 113]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  99] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 99] );

                        //t_0_6_9 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 129] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 114] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 114]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  100] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 100] );

                        //t_0_5_10 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 130] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 115] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 115]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  101] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 101] );

                        //t_0_4_11 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 131] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 116] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 116]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  102] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 102] );

                        //t_0_3_12 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 132] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 117] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 117]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  103] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 103] );

                        //t_0_2_13 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 133] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 118] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 118]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  104] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 104] );

                        //t_0_1_14 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 134] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 119] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 119];

                        //t_0_0_15 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 135] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 119] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 119]
                                      + 14 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  104] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 104] );

                    }
}



// VRR to obtain AUX_INT__u_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_u(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__u_s_s_s,
           double const * const restrict AUX_INT__t_s_s_s,
           double const * const restrict AUX_INT__r_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__u_s_s_s[num_m * 153];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //u_16_0_0 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 0] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 0] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 0]
                                      + 15 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  0] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 0] );

                        //u_15_1_0 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 1] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 0] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 0];

                        //u_15_0_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 2] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 0] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 0];

                        //u_14_2_0 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 3] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 1] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  0] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 0] );

                        //u_14_1_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 4] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 1] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 1];

                        //u_14_0_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 5] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 2] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  0] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 0] );

                        //u_13_3_0 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 6] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 3] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  1] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 1] );

                        //u_13_2_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 7] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 3] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 3];

                        //u_13_1_2 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 8] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 5] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 5];

                        //u_13_0_3 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 9] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 5] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  2] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 2] );

                        //u_12_4_0 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 10] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 6] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  3] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 3] );

                        //u_12_3_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 11] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 6] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 6];

                        //u_12_2_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 12] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 7] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  3] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 3] );

                        //u_12_1_3 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 13] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 9] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 9];

                        //u_12_0_4 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 14] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 9] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  5] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 5] );

                        //u_11_5_0 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 15] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 10] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  6] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 6] );

                        //u_11_4_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 16] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 10] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 10];

                        //u_11_3_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 17] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 11] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  6] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 6] );

                        //u_11_2_3 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 18] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 13] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  9] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 9] );

                        //u_11_1_4 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 19] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 14] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 14];

                        //u_11_0_5 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 20] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 14] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  9] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 9] );

                        //u_10_6_0 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 21] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 15] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  10] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 10] );

                        //u_10_5_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 22] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 15] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 15];

                        //u_10_4_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 23] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 16] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  10] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 10] );

                        //u_10_3_3 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 24] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 17] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  11] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 11] );

                        //u_10_2_4 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 25] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 19] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  14] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 14] );

                        //u_10_1_5 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 26] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 20] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 20];

                        //u_10_0_6 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 27] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 20] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  14] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 14] );

                        //u_9_7_0 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 28] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 21] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  15] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 15] );

                        //u_9_6_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 29] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 21] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 21];

                        //u_9_5_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 30] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 22] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  15] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 15] );

                        //u_9_4_3 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 31] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 23] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  16] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 16] );

                        //u_9_3_4 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 32] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 25] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  19] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 19] );

                        //u_9_2_5 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 33] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 26] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  20] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 20] );

                        //u_9_1_6 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 34] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 27] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 27];

                        //u_9_0_7 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 35] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 27] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  20] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 20] );

                        //u_8_8_0 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 36] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 28] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  21] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 21] );

                        //u_8_7_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 37] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 28] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 28];

                        //u_8_6_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 38] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 29] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  21] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 21] );

                        //u_8_5_3 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 39] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 30] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  22] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 22] );

                        //u_8_4_4 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 40] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 31] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  23] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 23] );

                        //u_8_3_5 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 41] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 33] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  26] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 26] );

                        //u_8_2_6 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 42] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 34] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  27] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 27] );

                        //u_8_1_7 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 43] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 35] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 35];

                        //u_8_0_8 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 44] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 35] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  27] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 27] );

                        //u_7_9_0 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 45] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 45] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 45]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  45] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 45] );

                        //u_7_8_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 46] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 36] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 36];

                        //u_7_7_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 47] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 37] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  28] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 28] );

                        //u_7_6_3 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 48] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 38] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  29] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 29] );

                        //u_7_5_4 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 49] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 39] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  30] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 30] );

                        //u_7_4_5 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 50] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 41] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  33] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 33] );

                        //u_7_3_6 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 51] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 42] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  34] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 34] );

                        //u_7_2_7 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 52] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 43] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  35] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 35] );

                        //u_7_1_8 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 53] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 44] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 44];

                        //u_7_0_9 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 54] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 54] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 54]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  54] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 54] );

                        //u_6_10_0 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 55] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 55] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 55]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  55] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 55] );

                        //u_6_9_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 56] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 45] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 45];

                        //u_6_8_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 57] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 46] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  36] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 36] );

                        //u_6_7_3 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 58] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 47] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  37] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 37] );

                        //u_6_6_4 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 59] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 48] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  38] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 38] );

                        //u_6_5_5 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 60] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 49] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  39] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 39] );

                        //u_6_4_6 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 61] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 51] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  42] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 42] );

                        //u_6_3_7 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 62] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 52] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  43] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 43] );

                        //u_6_2_8 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 63] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 53] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  44] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 44] );

                        //u_6_1_9 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 64] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 54] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 54];

                        //u_6_0_10 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 65] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 65] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 65]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  65] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 65] );

                        //u_5_11_0 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 66] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 66] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 66]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  66] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 66] );

                        //u_5_10_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 67] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 55] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 55];

                        //u_5_9_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 68] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 56] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  45] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 45] );

                        //u_5_8_3 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 69] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 57] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  46] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 46] );

                        //u_5_7_4 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 70] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 58] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  47] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 47] );

                        //u_5_6_5 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 71] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 59] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 59]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  48] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 48] );

                        //u_5_5_6 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 72] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 61] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 61]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  51] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 51] );

                        //u_5_4_7 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 73] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 62] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  52] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 52] );

                        //u_5_3_8 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 74] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 63] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  53] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 53] );

                        //u_5_2_9 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 75] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 64] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  54] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 54] );

                        //u_5_1_10 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 76] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 65] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 65];

                        //u_5_0_11 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 77] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 77] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 77]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  77] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 77] );

                        //u_4_12_0 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 78] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 78] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 78]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  78] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 78] );

                        //u_4_11_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 79] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 66] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 66];

                        //u_4_10_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 80] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 67] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  55] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 55] );

                        //u_4_9_3 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 81] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 68] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  56] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 56] );

                        //u_4_8_4 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 82] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 69] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 69]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  57] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 57] );

                        //u_4_7_5 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 83] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 83] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 83]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  83] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 83] );

                        //u_4_6_6 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 84] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 84] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 84]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  84] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 84] );

                        //u_4_5_7 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 85] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 85] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 85]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  85] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 85] );

                        //u_4_4_8 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 86] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 74] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 74]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  63] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 63] );

                        //u_4_3_9 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 87] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 75] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  64] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 64] );

                        //u_4_2_10 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 88] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 76] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  65] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 65] );

                        //u_4_1_11 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 89] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 77] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 77];

                        //u_4_0_12 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 90] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 90] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 90]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  90] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 90] );

                        //u_3_13_0 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 91] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 91] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 91]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  91] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 91] );

                        //u_3_12_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 92] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 78] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 78];

                        //u_3_11_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 93] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 79] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  66] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 66] );

                        //u_3_10_3 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 94] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 80] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 80]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  67] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 67] );

                        //u_3_9_4 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 95] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 95] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 95]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  95] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 95] );

                        //u_3_8_5 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 96] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 96] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 96]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  96] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 96] );

                        //u_3_7_6 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 97] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 97] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 97]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  97] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 97] );

                        //u_3_6_7 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 98] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 98] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 98]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  98] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 98] );

                        //u_3_5_8 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 99] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 99] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 99]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  99] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 99] );

                        //u_3_4_9 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 100] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 100] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 100]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  100] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 100] );

                        //u_3_3_10 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 101] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 88] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 88]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  76] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 76] );

                        //u_3_2_11 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 102] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 89] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  77] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 77] );

                        //u_3_1_12 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 103] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 90] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 90];

                        //u_3_0_13 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 104] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 104] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 104]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  104] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 104] );

                        //u_2_14_0 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 105] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 105] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 105]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  105] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 105] );

                        //u_2_13_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 106] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 91] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 91];

                        //u_2_12_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 107] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 92] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 92]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  78] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 78] );

                        //u_2_11_3 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 108] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 108] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 108]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  108] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 108] );

                        //u_2_10_4 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 109] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 109] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 109]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  109] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 109] );

                        //u_2_9_5 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 110] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 110] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 110]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  110] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 110] );

                        //u_2_8_6 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 111] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 111] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 111]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  111] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 111] );

                        //u_2_7_7 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 112] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 112] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 112]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  112] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 112] );

                        //u_2_6_8 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 113] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 113] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 113]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  113] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 113] );

                        //u_2_5_9 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 114] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 114] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 114]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  114] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 114] );

                        //u_2_4_10 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 115] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 115] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 115]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  115] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 115] );

                        //u_2_3_11 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 116] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 116] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 116]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  116] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 116] );

                        //u_2_2_12 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 117] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 103] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 103]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  90] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 90] );

                        //u_2_1_13 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 118] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 104] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 104];

                        //u_2_0_14 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 119] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 119] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 119]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  119] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 119] );

                        //u_1_15_0 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 120] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 120] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 120];

                        //u_1_14_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 121] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 105] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 105];

                        //u_1_13_2 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 122] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 122] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 122];

                        //u_1_12_3 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 123] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 123] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 123];

                        //u_1_11_4 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 124] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 124] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 124];

                        //u_1_10_5 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 125] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 125] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 125];

                        //u_1_9_6 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 126] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 126] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 126];

                        //u_1_8_7 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 127] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 127] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 127];

                        //u_1_7_8 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 128] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 128] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 128];

                        //u_1_6_9 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 129] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 129] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 129];

                        //u_1_5_10 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 130] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 130] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 130];

                        //u_1_4_11 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 131] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 131] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 131];

                        //u_1_3_12 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 132] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 132] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 132];

                        //u_1_2_13 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 133] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 133] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 133];

                        //u_1_1_14 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 134] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 119] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 119];

                        //u_1_0_15 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 135] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 135] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 135];

                        //u_0_16_0 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 136] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 120] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 120]
                                      + 15 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  105] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 105] );

                        //u_0_15_1 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 137] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 120] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 120];

                        //u_0_14_2 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 138] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 121] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 121]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  105] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 105] );

                        //u_0_13_3 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 139] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 122] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 122]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  106] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 106] );

                        //u_0_12_4 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 140] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 123] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 123]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  107] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 107] );

                        //u_0_11_5 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 141] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 124] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 124]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  108] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 108] );

                        //u_0_10_6 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 142] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 125] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 125]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  109] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 109] );

                        //u_0_9_7 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 143] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 126] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 126]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  110] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 110] );

                        //u_0_8_8 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 144] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 127] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 127]
                                      + 7 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  111] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 111] );

                        //u_0_7_9 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 145] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 129] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 129]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  114] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 114] );

                        //u_0_6_10 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 146] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 130] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 130]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  115] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 115] );

                        //u_0_5_11 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 147] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 131] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 131]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  116] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 116] );

                        //u_0_4_12 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 148] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 132] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 132]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  117] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 117] );

                        //u_0_3_13 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 149] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 133] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 133]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  118] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 118] );

                        //u_0_2_14 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 150] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 134] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 134]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  119] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 119] );

                        //u_0_1_15 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 151] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 135] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 135];

                        //u_0_0_16 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 152] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 135] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 135]
                                      + 15 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  119] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 119] );

                    }
}



// VRR to obtain AUX_INT__v_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_v(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__v_s_s_s,
           double const * const restrict AUX_INT__u_s_s_s,
           double const * const restrict AUX_INT__t_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__v_s_s_s[num_m * 171];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //v_17_0_0 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 0] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 0] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 0]
                                      + 16 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  0] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 0] );

                        //v_16_1_0 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 1] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 0] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 0];

                        //v_16_0_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 2] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 0] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 0];

                        //v_15_2_0 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 3] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 1] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  0] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 0] );

                        //v_15_1_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 4] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 1] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 1];

                        //v_15_0_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 5] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 2] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  0] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 0] );

                        //v_14_3_0 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 6] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 3] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  1] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 1] );

                        //v_14_2_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 7] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 3] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 3];

                        //v_14_1_2 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 8] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 5] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 5];

                        //v_14_0_3 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 9] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 5] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  2] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 2] );

                        //v_13_4_0 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 10] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 6] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  3] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 3] );

                        //v_13_3_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 11] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 6] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 6];

                        //v_13_2_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 12] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 7] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  3] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 3] );

                        //v_13_1_3 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 13] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 9] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 9];

                        //v_13_0_4 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 14] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 9] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  5] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 5] );

                        //v_12_5_0 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 15] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 10] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  6] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 6] );

                        //v_12_4_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 16] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 10] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 10];

                        //v_12_3_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 17] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 11] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  6] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 6] );

                        //v_12_2_3 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 18] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 13] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  9] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 9] );

                        //v_12_1_4 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 19] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 14] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 14];

                        //v_12_0_5 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 20] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 14] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  9] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 9] );

                        //v_11_6_0 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 21] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 15] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  10] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 10] );

                        //v_11_5_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 22] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 15] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 15];

                        //v_11_4_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 23] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 16] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  10] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 10] );

                        //v_11_3_3 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 24] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 17] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  11] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 11] );

                        //v_11_2_4 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 25] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 19] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  14] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 14] );

                        //v_11_1_5 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 26] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 20] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 20];

                        //v_11_0_6 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 27] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 20] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  14] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 14] );

                        //v_10_7_0 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 28] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 21] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  15] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 15] );

                        //v_10_6_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 29] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 21] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 21];

                        //v_10_5_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 30] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 22] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  15] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 15] );

                        //v_10_4_3 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 31] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 23] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  16] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 16] );

                        //v_10_3_4 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 32] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 25] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  19] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 19] );

                        //v_10_2_5 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 33] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 26] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  20] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 20] );

                        //v_10_1_6 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 34] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 27] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 27];

                        //v_10_0_7 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 35] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 27] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  20] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 20] );

                        //v_9_8_0 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 36] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 28] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  21] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 21] );

                        //v_9_7_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 37] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 28] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 28];

                        //v_9_6_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 38] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 29] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  21] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 21] );

                        //v_9_5_3 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 39] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 30] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  22] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 22] );

                        //v_9_4_4 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 40] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 31] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  23] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 23] );

                        //v_9_3_5 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 41] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 33] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  26] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 26] );

                        //v_9_2_6 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 42] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 34] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  27] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 27] );

                        //v_9_1_7 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 43] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 35] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 35];

                        //v_9_0_8 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 44] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 35] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  27] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 27] );

                        //v_8_9_0 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 45] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 45] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 45]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  45] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 45] );

                        //v_8_8_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 46] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 36] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 36];

                        //v_8_7_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 47] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 37] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  28] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 28] );

                        //v_8_6_3 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 48] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 38] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  29] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 29] );

                        //v_8_5_4 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 49] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 39] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  30] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 30] );

                        //v_8_4_5 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 50] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 41] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  33] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 33] );

                        //v_8_3_6 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 51] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 42] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  34] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 34] );

                        //v_8_2_7 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 52] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 43] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  35] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 35] );

                        //v_8_1_8 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 53] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 44] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 44];

                        //v_8_0_9 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 54] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 54] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 54]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  54] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 54] );

                        //v_7_10_0 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 55] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 55] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 55]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  55] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 55] );

                        //v_7_9_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 56] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 45] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 45];

                        //v_7_8_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 57] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 46] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  36] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 36] );

                        //v_7_7_3 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 58] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 47] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  37] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 37] );

                        //v_7_6_4 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 59] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 48] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  38] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 38] );

                        //v_7_5_5 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 60] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 49] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  39] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 39] );

                        //v_7_4_6 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 61] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 51] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  42] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 42] );

                        //v_7_3_7 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 62] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 52] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  43] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 43] );

                        //v_7_2_8 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 63] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 53] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  44] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 44] );

                        //v_7_1_9 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 64] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 54] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 54];

                        //v_7_0_10 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 65] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 65] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 65]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  65] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 65] );

                        //v_6_11_0 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 66] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 66] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 66]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  66] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 66] );

                        //v_6_10_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 67] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 55] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 55];

                        //v_6_9_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 68] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 56] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  45] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 45] );

                        //v_6_8_3 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 69] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 57] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  46] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 46] );

                        //v_6_7_4 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 70] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 58] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  47] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 47] );

                        //v_6_6_5 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 71] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 59] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 59]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  48] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 48] );

                        //v_6_5_6 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 72] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 61] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 61]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  51] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 51] );

                        //v_6_4_7 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 73] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 62] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  52] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 52] );

                        //v_6_3_8 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 74] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 63] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  53] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 53] );

                        //v_6_2_9 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 75] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 64] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  54] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 54] );

                        //v_6_1_10 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 76] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 65] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 65];

                        //v_6_0_11 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 77] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 77] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 77]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  77] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 77] );

                        //v_5_12_0 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 78] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 78] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 78]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  78] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 78] );

                        //v_5_11_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 79] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 66] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 66];

                        //v_5_10_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 80] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 67] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  55] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 55] );

                        //v_5_9_3 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 81] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 68] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  56] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 56] );

                        //v_5_8_4 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 82] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 69] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 69]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  57] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 57] );

                        //v_5_7_5 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 83] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 70] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 70]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  58] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 58] );

                        //v_5_6_6 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 84] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 84] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 84]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  84] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 84] );

                        //v_5_5_7 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 85] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 73] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 73]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  62] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 62] );

                        //v_5_4_8 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 86] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 74] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 74]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  63] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 63] );

                        //v_5_3_9 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 87] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 75] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  64] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 64] );

                        //v_5_2_10 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 88] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 76] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  65] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 65] );

                        //v_5_1_11 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 89] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 77] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 77];

                        //v_5_0_12 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 90] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 90] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 90]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  90] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 90] );

                        //v_4_13_0 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 91] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 91] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 91]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  91] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 91] );

                        //v_4_12_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 92] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 78] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 78];

                        //v_4_11_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 93] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 79] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  66] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 66] );

                        //v_4_10_3 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 94] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 80] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 80]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  67] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 67] );

                        //v_4_9_4 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 95] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 81] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 81]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  68] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 68] );

                        //v_4_8_5 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 96] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 96] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 96]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  96] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 96] );

                        //v_4_7_6 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 97] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 97] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 97]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  97] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 97] );

                        //v_4_6_7 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 98] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 98] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 98]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  98] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 98] );

                        //v_4_5_8 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 99] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 99] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 99]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  99] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 99] );

                        //v_4_4_9 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 100] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 87] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 87]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  75] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 75] );

                        //v_4_3_10 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 101] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 88] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 88]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  76] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 76] );

                        //v_4_2_11 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 102] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 89] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  77] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 77] );

                        //v_4_1_12 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 103] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 90] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 90];

                        //v_4_0_13 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 104] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 104] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 104]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  104] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 104] );

                        //v_3_14_0 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 105] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 105] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 105]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  105] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 105] );

                        //v_3_13_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 106] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 91] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 91];

                        //v_3_12_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 107] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 92] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 92]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  78] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 78] );

                        //v_3_11_3 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 108] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 93] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 93]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  79] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 79] );

                        //v_3_10_4 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 109] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 109] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 109]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  109] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 109] );

                        //v_3_9_5 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 110] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 110] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 110]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  110] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 110] );

                        //v_3_8_6 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 111] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 111] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 111]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  111] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 111] );

                        //v_3_7_7 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 112] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 112] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 112]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  112] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 112] );

                        //v_3_6_8 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 113] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 113] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 113]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  113] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 113] );

                        //v_3_5_9 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 114] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 114] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 114]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  114] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 114] );

                        //v_3_4_10 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 115] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 115] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 115]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  115] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 115] );

                        //v_3_3_11 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 116] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 102] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 102]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  89] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 89] );

                        //v_3_2_12 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 117] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 103] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 103]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  90] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 90] );

                        //v_3_1_13 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 118] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 104] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 104];

                        //v_3_0_14 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 119] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 119] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 119]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  119] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 119] );

                        //v_2_15_0 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 120] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 120] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 120]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  120] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 120] );

                        //v_2_14_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 121] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 105] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 105];

                        //v_2_13_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 122] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 106] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 106]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  91] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 91] );

                        //v_2_12_3 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 123] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 123] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 123]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  123] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 123] );

                        //v_2_11_4 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 124] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 124] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 124]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  124] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 124] );

                        //v_2_10_5 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 125] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 125] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 125]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  125] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 125] );

                        //v_2_9_6 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 126] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 126] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 126]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  126] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 126] );

                        //v_2_8_7 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 127] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 127] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 127]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  127] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 127] );

                        //v_2_7_8 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 128] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 128] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 128]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  128] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 128] );

                        //v_2_6_9 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 129] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 129] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 129]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  129] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 129] );

                        //v_2_5_10 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 130] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 130] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 130]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  130] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 130] );

                        //v_2_4_11 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 131] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 131] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 131]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  131] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 131] );

                        //v_2_3_12 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 132] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 132] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 132]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  132] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 132] );

                        //v_2_2_13 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 133] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 118] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 118]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  104] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 104] );

                        //v_2_1_14 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 134] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 119] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 119];

                        //v_2_0_15 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 135] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 135] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 135]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  135] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 135] );

                        //v_1_16_0 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 136] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 136] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 136];

                        //v_1_15_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 137] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 120] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 120];

                        //v_1_14_2 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 138] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 138] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 138];

                        //v_1_13_3 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 139] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 139] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 139];

                        //v_1_12_4 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 140] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 140] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 140];

                        //v_1_11_5 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 141] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 141] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 141];

                        //v_1_10_6 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 142] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 142] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 142];

                        //v_1_9_7 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 143] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 143] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 143];

                        //v_1_8_8 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 144] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 144] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 144];

                        //v_1_7_9 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 145] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 145] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 145];

                        //v_1_6_10 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 146] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 146] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 146];

                        //v_1_5_11 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 147] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 147] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 147];

                        //v_1_4_12 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 148] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 148] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 148];

                        //v_1_3_13 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 149] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 149] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 149];

                        //v_1_2_14 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 150] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 150] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 150];

                        //v_1_1_15 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 151] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 135] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 135];

                        //v_1_0_16 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 152] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 152] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 152];

                        //v_0_17_0 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 153] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 136] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 136]
                                      + 16 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  120] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 120] );

                        //v_0_16_1 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 154] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 136] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 136];

                        //v_0_15_2 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 155] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 137] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 137]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  120] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 120] );

                        //v_0_14_3 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 156] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 138] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 138]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  121] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 121] );

                        //v_0_13_4 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 157] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 139] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 139]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  122] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 122] );

                        //v_0_12_5 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 158] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 140] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 140]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  123] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 123] );

                        //v_0_11_6 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 159] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 141] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 141]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  124] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 124] );

                        //v_0_10_7 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 160] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 142] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 142]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  125] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 125] );

                        //v_0_9_8 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 161] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 143] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 143]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  126] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 126] );

                        //v_0_8_9 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 162] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 145] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 145]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  129] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 129] );

                        //v_0_7_10 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 163] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 146] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 146]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  130] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 130] );

                        //v_0_6_11 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 164] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 147] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 147]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  131] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 131] );

                        //v_0_5_12 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 165] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 148] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 148]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  132] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 132] );

                        //v_0_4_13 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 166] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 149] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 149]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  133] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 133] );

                        //v_0_3_14 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 167] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 150] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 150]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  134] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 134] );

                        //v_0_2_15 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 168] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 151] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 151]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  135] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 135] );

                        //v_0_1_16 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 169] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 152] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 152];

                        //v_0_0_17 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 170] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 152] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 152]
                                      + 16 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  135] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 135] );

                    }
}



// VRR to obtain AUX_INT__w_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_w(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__w_s_s_s,
           double const * const restrict AUX_INT__v_s_s_s,
           double const * const restrict AUX_INT__u_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__w_s_s_s[num_m * 190];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //w_18_0_0 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 0] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 0] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 0]
                                      + 17 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  0] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 0] );

                        //w_17_1_0 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 1] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 0] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 0];

                        //w_17_0_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 2] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 0] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 0];

                        //w_16_2_0 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 3] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 1] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  0] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 0] );

                        //w_16_1_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 4] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 1] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 1];

                        //w_16_0_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 5] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 2] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  0] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 0] );

                        //w_15_3_0 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 6] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 3] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  1] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 1] );

                        //w_15_2_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 7] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 3] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 3];

                        //w_15_1_2 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 8] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 5] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 5];

                        //w_15_0_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 9] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 5] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  2] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 2] );

                        //w_14_4_0 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 10] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 6] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  3] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 3] );

                        //w_14_3_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 11] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 6] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 6];

                        //w_14_2_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 12] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 7] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  3] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 3] );

                        //w_14_1_3 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 13] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 9] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 9];

                        //w_14_0_4 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 14] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 9] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  5] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 5] );

                        //w_13_5_0 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 15] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 10] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  6] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 6] );

                        //w_13_4_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 16] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 10] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 10];

                        //w_13_3_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 17] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 11] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  6] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 6] );

                        //w_13_2_3 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 18] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 13] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  9] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 9] );

                        //w_13_1_4 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 19] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 14] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 14];

                        //w_13_0_5 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 20] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 14] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  9] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 9] );

                        //w_12_6_0 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 21] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 15] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  10] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 10] );

                        //w_12_5_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 22] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 15] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 15];

                        //w_12_4_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 23] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 16] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  10] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 10] );

                        //w_12_3_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 24] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 17] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  11] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 11] );

                        //w_12_2_4 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 25] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 19] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  14] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 14] );

                        //w_12_1_5 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 26] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 20] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 20];

                        //w_12_0_6 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 27] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 20] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  14] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 14] );

                        //w_11_7_0 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 28] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 21] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  15] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 15] );

                        //w_11_6_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 29] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 21] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 21];

                        //w_11_5_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 30] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 22] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  15] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 15] );

                        //w_11_4_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 31] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 23] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  16] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 16] );

                        //w_11_3_4 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 32] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 25] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  19] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 19] );

                        //w_11_2_5 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 33] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 26] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  20] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 20] );

                        //w_11_1_6 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 34] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 27] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 27];

                        //w_11_0_7 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 35] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 27] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  20] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 20] );

                        //w_10_8_0 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 36] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 28] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  21] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 21] );

                        //w_10_7_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 37] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 28] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 28];

                        //w_10_6_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 38] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 29] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  21] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 21] );

                        //w_10_5_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 39] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 30] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  22] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 22] );

                        //w_10_4_4 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 40] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 31] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  23] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 23] );

                        //w_10_3_5 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 41] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 33] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  26] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 26] );

                        //w_10_2_6 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 42] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 34] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  27] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 27] );

                        //w_10_1_7 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 43] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 35] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 35];

                        //w_10_0_8 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 44] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 35] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  27] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 27] );

                        //w_9_9_0 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 45] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 36] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 36]
                                      + 8 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  28] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 28] );

                        //w_9_8_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 46] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 36] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 36];

                        //w_9_7_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 47] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 37] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  28] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 28] );

                        //w_9_6_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 48] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 38] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  29] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 29] );

                        //w_9_5_4 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 49] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 39] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  30] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 30] );

                        //w_9_4_5 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 50] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 41] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  33] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 33] );

                        //w_9_3_6 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 51] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 42] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  34] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 34] );

                        //w_9_2_7 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 52] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 43] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  35] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 35] );

                        //w_9_1_8 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 53] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 44] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 44];

                        //w_9_0_9 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 54] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 44] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 44]
                                      + 8 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  35] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 35] );

                        //w_8_10_0 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 55] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 55] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 55]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  55] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 55] );

                        //w_8_9_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 56] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 45] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 45];

                        //w_8_8_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 57] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 46] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  36] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 36] );

                        //w_8_7_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 58] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 47] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  37] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 37] );

                        //w_8_6_4 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 59] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 48] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  38] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 38] );

                        //w_8_5_5 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 60] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 49] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  39] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 39] );

                        //w_8_4_6 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 61] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 51] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  42] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 42] );

                        //w_8_3_7 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 62] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 52] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  43] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 43] );

                        //w_8_2_8 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 63] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 53] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  44] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 44] );

                        //w_8_1_9 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 64] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 54] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 54];

                        //w_8_0_10 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 65] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 65] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 65]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  65] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 65] );

                        //w_7_11_0 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 66] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 66] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 66]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  66] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 66] );

                        //w_7_10_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 67] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 55] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 55];

                        //w_7_9_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 68] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 56] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  45] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 45] );

                        //w_7_8_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 69] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 57] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  46] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 46] );

                        //w_7_7_4 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 70] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 58] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  47] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 47] );

                        //w_7_6_5 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 71] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 59] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 59]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  48] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 48] );

                        //w_7_5_6 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 72] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 61] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 61]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  51] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 51] );

                        //w_7_4_7 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 73] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 62] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  52] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 52] );

                        //w_7_3_8 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 74] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 63] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  53] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 53] );

                        //w_7_2_9 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 75] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 64] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  54] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 54] );

                        //w_7_1_10 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 76] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 65] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 65];

                        //w_7_0_11 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 77] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 77] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 77]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  77] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 77] );

                        //w_6_12_0 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 78] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 78] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 78]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  78] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 78] );

                        //w_6_11_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 79] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 66] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 66];

                        //w_6_10_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 80] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 67] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  55] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 55] );

                        //w_6_9_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 81] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 68] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  56] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 56] );

                        //w_6_8_4 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 82] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 69] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 69]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  57] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 57] );

                        //w_6_7_5 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 83] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 70] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 70]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  58] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 58] );

                        //w_6_6_6 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 84] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 71] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 71]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  59] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 59] );

                        //w_6_5_7 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 85] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 73] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 73]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  62] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 62] );

                        //w_6_4_8 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 86] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 74] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 74]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  63] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 63] );

                        //w_6_3_9 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 87] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 75] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  64] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 64] );

                        //w_6_2_10 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 88] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 76] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  65] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 65] );

                        //w_6_1_11 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 89] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 77] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 77];

                        //w_6_0_12 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 90] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 90] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 90]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  90] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 90] );

                        //w_5_13_0 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 91] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 91] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 91]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  91] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 91] );

                        //w_5_12_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 92] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 78] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 78];

                        //w_5_11_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 93] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 79] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  66] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 66] );

                        //w_5_10_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 94] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 80] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 80]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  67] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 67] );

                        //w_5_9_4 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 95] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 81] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 81]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  68] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 68] );

                        //w_5_8_5 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 96] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 82] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 82]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  69] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 69] );

                        //w_5_7_6 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 97] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 97] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 97]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  97] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 97] );

                        //w_5_6_7 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 98] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 98] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 98]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  98] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 98] );

                        //w_5_5_8 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 99] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 86] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 86]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  74] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 74] );

                        //w_5_4_9 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 100] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 87] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 87]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  75] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 75] );

                        //w_5_3_10 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 101] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 88] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 88]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  76] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 76] );

                        //w_5_2_11 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 102] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 89] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  77] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 77] );

                        //w_5_1_12 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 103] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 90] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 90];

                        //w_5_0_13 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 104] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 104] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 104]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  104] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 104] );

                        //w_4_14_0 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 105] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 105] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 105]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  105] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 105] );

                        //w_4_13_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 106] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 91] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 91];

                        //w_4_12_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 107] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 92] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 92]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  78] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 78] );

                        //w_4_11_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 108] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 93] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 93]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  79] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 79] );

                        //w_4_10_4 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 109] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 94] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 94]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  80] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 80] );

                        //w_4_9_5 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 110] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 110] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 110]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  110] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 110] );

                        //w_4_8_6 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 111] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 111] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 111]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  111] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 111] );

                        //w_4_7_7 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 112] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 112] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 112]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  112] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 112] );

                        //w_4_6_8 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 113] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 113] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 113]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  113] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 113] );

                        //w_4_5_9 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 114] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 114] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 114]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  114] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 114] );

                        //w_4_4_10 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 115] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 101] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 101]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  88] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 88] );

                        //w_4_3_11 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 116] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 102] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 102]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  89] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 89] );

                        //w_4_2_12 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 117] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 103] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 103]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  90] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 90] );

                        //w_4_1_13 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 118] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 104] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 104];

                        //w_4_0_14 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 119] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 119] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 119]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  119] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 119] );

                        //w_3_15_0 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 120] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 120] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 120]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  120] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 120] );

                        //w_3_14_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 121] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 105] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 105];

                        //w_3_13_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 122] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 106] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 106]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  91] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 91] );

                        //w_3_12_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 123] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 107] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 107]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  92] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 92] );

                        //w_3_11_4 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 124] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 124] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 124]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  124] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 124] );

                        //w_3_10_5 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 125] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 125] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 125]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  125] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 125] );

                        //w_3_9_6 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 126] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 126] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 126]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  126] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 126] );

                        //w_3_8_7 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 127] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 127] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 127]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  127] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 127] );

                        //w_3_7_8 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 128] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 128] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 128]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  128] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 128] );

                        //w_3_6_9 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 129] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 129] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 129]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  129] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 129] );

                        //w_3_5_10 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 130] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 130] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 130]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  130] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 130] );

                        //w_3_4_11 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 131] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 131] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 131]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  131] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 131] );

                        //w_3_3_12 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 132] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 117] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 117]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  103] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 103] );

                        //w_3_2_13 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 133] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 118] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 118]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  104] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 104] );

                        //w_3_1_14 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 134] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 119] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 119];

                        //w_3_0_15 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 135] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 135] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 135]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  135] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 135] );

                        //w_2_16_0 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 136] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 136] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 136]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  136] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 136] );

                        //w_2_15_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 137] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 120] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 120];

                        //w_2_14_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 138] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 121] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 121]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  105] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 105] );

                        //w_2_13_3 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 139] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 139] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 139]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  139] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 139] );

                        //w_2_12_4 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 140] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 140] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 140]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  140] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 140] );

                        //w_2_11_5 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 141] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 141] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 141]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  141] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 141] );

                        //w_2_10_6 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 142] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 142] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 142]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  142] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 142] );

                        //w_2_9_7 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 143] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 143] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 143]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  143] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 143] );

                        //w_2_8_8 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 144] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 144] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 144]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  144] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 144] );

                        //w_2_7_9 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 145] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 145] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 145]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  145] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 145] );

                        //w_2_6_10 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 146] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 146] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 146]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  146] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 146] );

                        //w_2_5_11 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 147] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 147] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 147]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  147] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 147] );

                        //w_2_4_12 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 148] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 148] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 148]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  148] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 148] );

                        //w_2_3_13 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 149] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 149] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 149]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  149] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 149] );

                        //w_2_2_14 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 150] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 134] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 134]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  119] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 119] );

                        //w_2_1_15 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 151] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 135] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 135];

                        //w_2_0_16 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 152] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 152] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 152]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  152] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 152] );

                        //w_1_17_0 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 153] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 153] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 153];

                        //w_1_16_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 154] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 136] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 136];

                        //w_1_15_2 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 155] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 155] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 155];

                        //w_1_14_3 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 156] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 156] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 156];

                        //w_1_13_4 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 157] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 157] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 157];

                        //w_1_12_5 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 158] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 158] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 158];

                        //w_1_11_6 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 159] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 159] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 159];

                        //w_1_10_7 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 160] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 160] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 160];

                        //w_1_9_8 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 161] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 161] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 161];

                        //w_1_8_9 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 162] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 162] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 162];

                        //w_1_7_10 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 163] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 163] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 163];

                        //w_1_6_11 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 164] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 164] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 164];

                        //w_1_5_12 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 165] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 165] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 165];

                        //w_1_4_13 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 166] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 166] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 166];

                        //w_1_3_14 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 167] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 167] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 167];

                        //w_1_2_15 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 168] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 168] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 168];

                        //w_1_1_16 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 169] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 152] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 152];

                        //w_1_0_17 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 170] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 170] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 170];

                        //w_0_18_0 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 171] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 153] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 153]
                                      + 17 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  136] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 136] );

                        //w_0_17_1 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 172] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 153] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 153];

                        //w_0_16_2 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 173] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 154] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 154]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  136] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 136] );

                        //w_0_15_3 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 174] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 155] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 155]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  137] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 137] );

                        //w_0_14_4 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 175] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 156] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 156]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  138] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 138] );

                        //w_0_13_5 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 176] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 157] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 157]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  139] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 139] );

                        //w_0_12_6 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 177] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 158] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 158]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  140] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 140] );

                        //w_0_11_7 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 178] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 159] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 159]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  141] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 141] );

                        //w_0_10_8 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 179] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 160] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 160]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  142] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 142] );

                        //w_0_9_9 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 180] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 161] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 161]
                                      + 8 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  143] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 143] );

                        //w_0_8_10 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 181] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 163] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 163]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  146] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 146] );

                        //w_0_7_11 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 182] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 164] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 164]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  147] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 147] );

                        //w_0_6_12 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 183] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 165] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 165]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  148] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 148] );

                        //w_0_5_13 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 184] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 166] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 166]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  149] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 149] );

                        //w_0_4_14 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 185] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 167] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 167]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  150] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 150] );

                        //w_0_3_15 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 186] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 168] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 168]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  151] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 151] );

                        //w_0_2_16 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 187] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 169] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 169]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  152] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 152] );

                        //w_0_1_17 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 188] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 170] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 170];

                        //w_0_0_18 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 189] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 170] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 170]
                                      + 17 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  152] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 152] );

                    }
}



// VRR to obtain AUX_INT__x_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_x(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__x_s_s_s,
           double const * const restrict AUX_INT__w_s_s_s,
           double const * const restrict AUX_INT__v_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__x_s_s_s[num_m * 210];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //x_19_0_0 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 0] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 0] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 0]
                                      + 18 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  0] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 0] );

                        //x_18_1_0 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 1] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 0] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 0];

                        //x_18_0_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 2] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 0] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 0];

                        //x_17_2_0 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 3] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 1] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  0] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 0] );

                        //x_17_1_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 4] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 1] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 1];

                        //x_17_0_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 5] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 2] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  0] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 0] );

                        //x_16_3_0 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 6] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 3] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  1] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 1] );

                        //x_16_2_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 7] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 3] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 3];

                        //x_16_1_2 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 8] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 5] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 5];

                        //x_16_0_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 9] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 5] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  2] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 2] );

                        //x_15_4_0 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 10] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 6] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  3] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 3] );

                        //x_15_3_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 11] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 6] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 6];

                        //x_15_2_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 12] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 7] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  3] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 3] );

                        //x_15_1_3 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 13] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 9] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 9];

                        //x_15_0_4 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 14] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 9] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  5] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 5] );

                        //x_14_5_0 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 15] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 10] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  6] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 6] );

                        //x_14_4_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 16] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 10] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 10];

                        //x_14_3_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 17] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 11] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  6] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 6] );

                        //x_14_2_3 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 18] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 13] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  9] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 9] );

                        //x_14_1_4 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 19] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 14] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 14];

                        //x_14_0_5 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 20] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 14] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  9] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 9] );

                        //x_13_6_0 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 21] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 15] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  10] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 10] );

                        //x_13_5_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 22] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 15] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 15];

                        //x_13_4_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 23] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 16] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  10] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 10] );

                        //x_13_3_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 24] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 17] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  11] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 11] );

                        //x_13_2_4 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 25] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 19] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  14] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 14] );

                        //x_13_1_5 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 26] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 20] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 20];

                        //x_13_0_6 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 27] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 20] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  14] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 14] );

                        //x_12_7_0 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 28] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 21] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  15] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 15] );

                        //x_12_6_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 29] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 21] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 21];

                        //x_12_5_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 30] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 22] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  15] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 15] );

                        //x_12_4_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 31] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 23] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  16] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 16] );

                        //x_12_3_4 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 32] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 25] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  19] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 19] );

                        //x_12_2_5 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 33] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 26] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  20] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 20] );

                        //x_12_1_6 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 34] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 27] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 27];

                        //x_12_0_7 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 35] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 27] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  20] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 20] );

                        //x_11_8_0 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 36] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 28] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  21] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 21] );

                        //x_11_7_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 37] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 28] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 28];

                        //x_11_6_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 38] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 29] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  21] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 21] );

                        //x_11_5_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 39] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 30] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  22] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 22] );

                        //x_11_4_4 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 40] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 31] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  23] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 23] );

                        //x_11_3_5 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 41] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 33] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  26] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 26] );

                        //x_11_2_6 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 42] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 34] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  27] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 27] );

                        //x_11_1_7 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 43] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 35] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 35];

                        //x_11_0_8 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 44] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 35] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  27] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 27] );

                        //x_10_9_0 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 45] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 36] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 36]
                                      + 8 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  28] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 28] );

                        //x_10_8_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 46] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 36] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 36];

                        //x_10_7_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 47] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 37] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  28] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 28] );

                        //x_10_6_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 48] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 38] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  29] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 29] );

                        //x_10_5_4 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 49] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 39] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  30] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 30] );

                        //x_10_4_5 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 50] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 41] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  33] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 33] );

                        //x_10_3_6 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 51] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 42] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  34] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 34] );

                        //x_10_2_7 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 52] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 43] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  35] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 35] );

                        //x_10_1_8 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 53] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 44] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 44];

                        //x_10_0_9 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 54] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 44] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 44]
                                      + 8 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  35] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 35] );

                        //x_9_10_0 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 55] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 55] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 55]
                                      + 8 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  55] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 55] );

                        //x_9_9_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 56] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 45] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 45];

                        //x_9_8_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 57] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 46] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  36] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 36] );

                        //x_9_7_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 58] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 47] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  37] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 37] );

                        //x_9_6_4 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 59] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 48] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  38] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 38] );

                        //x_9_5_5 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 60] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 49] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  39] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 39] );

                        //x_9_4_6 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 61] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 51] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  42] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 42] );

                        //x_9_3_7 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 62] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 52] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  43] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 43] );

                        //x_9_2_8 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 63] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 53] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  44] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 44] );

                        //x_9_1_9 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 64] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 54] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 54];

                        //x_9_0_10 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 65] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 65] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 65]
                                      + 8 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  65] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 65] );

                        //x_8_11_0 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 66] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 66] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 66]
                                      + 7 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  66] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 66] );

                        //x_8_10_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 67] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 55] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 55];

                        //x_8_9_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 68] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 56] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  45] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 45] );

                        //x_8_8_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 69] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 57] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  46] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 46] );

                        //x_8_7_4 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 70] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 58] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  47] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 47] );

                        //x_8_6_5 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 71] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 59] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 59]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  48] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 48] );

                        //x_8_5_6 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 72] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 61] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 61]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  51] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 51] );

                        //x_8_4_7 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 73] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 62] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  52] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 52] );

                        //x_8_3_8 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 74] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 63] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  53] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 53] );

                        //x_8_2_9 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 75] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 64] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  54] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 54] );

                        //x_8_1_10 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 76] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 65] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 65];

                        //x_8_0_11 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 77] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 77] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 77]
                                      + 7 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  77] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 77] );

                        //x_7_12_0 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 78] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 78] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 78]
                                      + 6 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  78] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 78] );

                        //x_7_11_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 79] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 66] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 66];

                        //x_7_10_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 80] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 67] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  55] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 55] );

                        //x_7_9_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 81] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 68] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  56] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 56] );

                        //x_7_8_4 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 82] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 69] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 69]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  57] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 57] );

                        //x_7_7_5 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 83] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 70] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 70]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  58] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 58] );

                        //x_7_6_6 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 84] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 71] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 71]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  59] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 59] );

                        //x_7_5_7 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 85] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 73] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 73]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  62] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 62] );

                        //x_7_4_8 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 86] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 74] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 74]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  63] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 63] );

                        //x_7_3_9 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 87] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 75] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  64] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 64] );

                        //x_7_2_10 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 88] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 76] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  65] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 65] );

                        //x_7_1_11 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 89] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 77] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 77];

                        //x_7_0_12 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 90] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 90] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 90]
                                      + 6 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  90] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 90] );

                        //x_6_13_0 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 91] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 91] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 91]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  91] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 91] );

                        //x_6_12_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 92] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 78] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 78];

                        //x_6_11_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 93] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 79] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  66] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 66] );

                        //x_6_10_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 94] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 80] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 80]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  67] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 67] );

                        //x_6_9_4 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 95] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 81] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 81]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  68] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 68] );

                        //x_6_8_5 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 96] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 82] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 82]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  69] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 69] );

                        //x_6_7_6 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 97] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 83] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 83]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  70] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 70] );

                        //x_6_6_7 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 98] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 85] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 85]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  73] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 73] );

                        //x_6_5_8 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 99] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 86] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 86]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  74] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 74] );

                        //x_6_4_9 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 100] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 87] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 87]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  75] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 75] );

                        //x_6_3_10 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 101] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 88] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 88]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  76] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 76] );

                        //x_6_2_11 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 102] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 89] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  77] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 77] );

                        //x_6_1_12 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 103] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 90] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 90];

                        //x_6_0_13 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 104] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 104] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 104]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  104] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 104] );

                        //x_5_14_0 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 105] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 105] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 105]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  105] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 105] );

                        //x_5_13_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 106] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 91] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 91];

                        //x_5_12_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 107] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 92] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 92]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  78] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 78] );

                        //x_5_11_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 108] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 93] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 93]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  79] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 79] );

                        //x_5_10_4 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 109] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 94] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 94]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  80] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 80] );

                        //x_5_9_5 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 110] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 95] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 95]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  81] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 81] );

                        //x_5_8_6 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 111] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 111] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 111]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  111] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 111] );

                        //x_5_7_7 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 112] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 112] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 112]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  112] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 112] );

                        //x_5_6_8 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 113] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 113] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 113]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  113] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 113] );

                        //x_5_5_9 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 114] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 100] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 100]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  87] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 87] );

                        //x_5_4_10 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 115] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 101] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 101]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  88] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 88] );

                        //x_5_3_11 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 116] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 102] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 102]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  89] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 89] );

                        //x_5_2_12 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 117] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 103] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 103]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  90] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 90] );

                        //x_5_1_13 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 118] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 104] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 104];

                        //x_5_0_14 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 119] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 119] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 119]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  119] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 119] );

                        //x_4_15_0 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 120] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 120] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 120]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  120] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 120] );

                        //x_4_14_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 121] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 105] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 105];

                        //x_4_13_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 122] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 106] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 106]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  91] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 91] );

                        //x_4_12_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 123] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 107] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 107]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  92] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 92] );

                        //x_4_11_4 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 124] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 108] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 108]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  93] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 93] );

                        //x_4_10_5 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 125] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 125] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 125]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  125] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 125] );

                        //x_4_9_6 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 126] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 126] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 126]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  126] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 126] );

                        //x_4_8_7 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 127] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 127] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 127]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  127] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 127] );

                        //x_4_7_8 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 128] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 128] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 128]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  128] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 128] );

                        //x_4_6_9 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 129] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 129] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 129]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  129] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 129] );

                        //x_4_5_10 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 130] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 130] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 130]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  130] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 130] );

                        //x_4_4_11 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 131] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 116] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 116]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  102] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 102] );

                        //x_4_3_12 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 132] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 117] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 117]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  103] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 103] );

                        //x_4_2_13 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 133] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 118] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 118]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  104] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 104] );

                        //x_4_1_14 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 134] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 119] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 119];

                        //x_4_0_15 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 135] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 135] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 135]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  135] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 135] );

                        //x_3_16_0 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 136] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 136] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 136]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  136] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 136] );

                        //x_3_15_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 137] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 120] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 120];

                        //x_3_14_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 138] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 121] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 121]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  105] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 105] );

                        //x_3_13_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 139] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 122] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 122]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  106] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 106] );

                        //x_3_12_4 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 140] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 140] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 140]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  140] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 140] );

                        //x_3_11_5 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 141] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 141] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 141]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  141] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 141] );

                        //x_3_10_6 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 142] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 142] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 142]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  142] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 142] );

                        //x_3_9_7 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 143] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 143] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 143]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  143] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 143] );

                        //x_3_8_8 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 144] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 144] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 144]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  144] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 144] );

                        //x_3_7_9 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 145] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 145] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 145]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  145] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 145] );

                        //x_3_6_10 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 146] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 146] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 146]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  146] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 146] );

                        //x_3_5_11 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 147] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 147] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 147]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  147] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 147] );

                        //x_3_4_12 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 148] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 148] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 148]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  148] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 148] );

                        //x_3_3_13 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 149] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 133] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 133]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  118] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 118] );

                        //x_3_2_14 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 150] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 134] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 134]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  119] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 119] );

                        //x_3_1_15 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 151] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 135] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 135];

                        //x_3_0_16 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 152] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 152] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 152]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  152] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 152] );

                        //x_2_17_0 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 153] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 153] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 153]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  153] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 153] );

                        //x_2_16_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 154] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 136] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 136];

                        //x_2_15_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 155] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 137] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 137]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  120] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 120] );

                        //x_2_14_3 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 156] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 156] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 156]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  156] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 156] );

                        //x_2_13_4 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 157] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 157] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 157]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  157] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 157] );

                        //x_2_12_5 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 158] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 158] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 158]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  158] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 158] );

                        //x_2_11_6 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 159] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 159] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 159]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  159] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 159] );

                        //x_2_10_7 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 160] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 160] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 160]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  160] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 160] );

                        //x_2_9_8 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 161] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 161] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 161]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  161] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 161] );

                        //x_2_8_9 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 162] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 162] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 162]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  162] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 162] );

                        //x_2_7_10 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 163] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 163] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 163]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  163] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 163] );

                        //x_2_6_11 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 164] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 164] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 164]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  164] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 164] );

                        //x_2_5_12 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 165] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 165] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 165]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  165] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 165] );

                        //x_2_4_13 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 166] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 166] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 166]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  166] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 166] );

                        //x_2_3_14 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 167] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 167] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 167]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  167] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 167] );

                        //x_2_2_15 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 168] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 151] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 151]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  135] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 135] );

                        //x_2_1_16 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 169] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 152] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 152];

                        //x_2_0_17 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 170] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 170] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 170]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  170] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 170] );

                        //x_1_18_0 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 171] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 171] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 171];

                        //x_1_17_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 172] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 153] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 153];

                        //x_1_16_2 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 173] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 173] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 173];

                        //x_1_15_3 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 174] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 174] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 174];

                        //x_1_14_4 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 175] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 175] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 175];

                        //x_1_13_5 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 176] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 176] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 176];

                        //x_1_12_6 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 177] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 177] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 177];

                        //x_1_11_7 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 178] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 178] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 178];

                        //x_1_10_8 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 179] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 179] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 179];

                        //x_1_9_9 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 180] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 180] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 180];

                        //x_1_8_10 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 181] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 181] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 181];

                        //x_1_7_11 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 182] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 182] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 182];

                        //x_1_6_12 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 183] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 183] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 183];

                        //x_1_5_13 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 184] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 184] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 184];

                        //x_1_4_14 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 185] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 185] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 185];

                        //x_1_3_15 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 186] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 186] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 186];

                        //x_1_2_16 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 187] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 187] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 187];

                        //x_1_1_17 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 188] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 170] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 170];

                        //x_1_0_18 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 189] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 189] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 189];

                        //x_0_19_0 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 190] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 171] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 171]
                                      + 18 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  153] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 153] );

                        //x_0_18_1 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 191] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 171] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 171];

                        //x_0_17_2 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 192] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 172] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 172]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  153] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 153] );

                        //x_0_16_3 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 193] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 173] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 173]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  154] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 154] );

                        //x_0_15_4 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 194] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 174] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 174]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  155] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 155] );

                        //x_0_14_5 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 195] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 175] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 175]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  156] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 156] );

                        //x_0_13_6 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 196] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 176] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 176]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  157] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 157] );

                        //x_0_12_7 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 197] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 177] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 177]
                                      + 6 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  158] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 158] );

                        //x_0_11_8 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 198] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 178] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 178]
                                      + 7 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  159] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 159] );

                        //x_0_10_9 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 199] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 179] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 179]
                                      + 8 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  160] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 160] );

                        //x_0_9_10 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 200] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 181] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 181]
                                      + 8 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  163] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 163] );

                        //x_0_8_11 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 201] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 182] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 182]
                                      + 7 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  164] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 164] );

                        //x_0_7_12 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 202] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 183] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 183]
                                      + 6 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  165] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 165] );

                        //x_0_6_13 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 203] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 184] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 184]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  166] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 166] );

                        //x_0_5_14 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 204] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 185] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 185]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  167] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 167] );

                        //x_0_4_15 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 205] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 186] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 186]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  168] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 168] );

                        //x_0_3_16 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 206] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 187] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 187]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  169] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 169] );

                        //x_0_2_17 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 207] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 188] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 188]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  170] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 170] );

                        //x_0_1_18 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 208] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 189] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 189];

                        //x_0_0_19 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 209] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 189] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 189]
                                      + 18 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  170] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 170] );

                    }
}



// VRR to obtain AUX_INT__y_s_s_s
#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_m)
void VRR_y(const int num_m,
           const double P_PA_x, const double P_PA_y, const double P_PA_z,
           const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
           const double a_over_p, const double one_over_2p,
           double * const restrict AUX_INT__y_s_s_s,
           double const * const restrict AUX_INT__x_s_s_s,
           double const * const restrict AUX_INT__w_s_s_s)
{
    int m = 0;
                    // Forming AUX_INT__y_s_s_s[num_m * 231];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //y_20_0_0 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 0] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 0] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 0]
                                      + 19 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  0] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 0] );

                        //y_19_1_0 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 1] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 0] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 0];

                        //y_19_0_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 2] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 0] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 0];

                        //y_18_2_0 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 3] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 1] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  0] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 0] );

                        //y_18_1_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 4] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 1] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 1];

                        //y_18_0_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 5] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 2] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  0] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 0] );

                        //y_17_3_0 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 6] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 3] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  1] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 1] );

                        //y_17_2_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 7] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 3] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 3];

                        //y_17_1_2 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 8] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 5] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 5];

                        //y_17_0_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 9] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 5] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  2] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 2] );

                        //y_16_4_0 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 10] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 6] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  3] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 3] );

                        //y_16_3_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 11] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 6] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 6];

                        //y_16_2_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 12] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 7] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  3] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 3] );

                        //y_16_1_3 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 13] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 9] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 9];

                        //y_16_0_4 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 14] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 9] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  5] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 5] );

                        //y_15_5_0 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 15] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 10] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  6] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 6] );

                        //y_15_4_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 16] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 10] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 10];

                        //y_15_3_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 17] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 11] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  6] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 6] );

                        //y_15_2_3 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 18] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 13] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  9] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 9] );

                        //y_15_1_4 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 19] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 14] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 14];

                        //y_15_0_5 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 20] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 14] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  9] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 9] );

                        //y_14_6_0 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 21] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 15] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  10] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 10] );

                        //y_14_5_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 22] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 15] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 15];

                        //y_14_4_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 23] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 16] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  10] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 10] );

                        //y_14_3_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 24] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 17] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  11] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 11] );

                        //y_14_2_4 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 25] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 19] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  14] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 14] );

                        //y_14_1_5 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 26] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 20] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 20];

                        //y_14_0_6 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 27] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 20] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  14] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 14] );

                        //y_13_7_0 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 28] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 21] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  15] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 15] );

                        //y_13_6_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 29] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 21] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 21];

                        //y_13_5_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 30] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 22] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  15] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 15] );

                        //y_13_4_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 31] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 23] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  16] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 16] );

                        //y_13_3_4 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 32] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 25] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  19] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 19] );

                        //y_13_2_5 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 33] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 26] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  20] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 20] );

                        //y_13_1_6 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 34] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 27] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 27];

                        //y_13_0_7 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 35] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 27] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  20] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 20] );

                        //y_12_8_0 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 36] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 28] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  21] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 21] );

                        //y_12_7_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 37] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 28] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 28];

                        //y_12_6_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 38] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 29] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  21] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 21] );

                        //y_12_5_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 39] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 30] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  22] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 22] );

                        //y_12_4_4 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 40] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 31] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  23] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 23] );

                        //y_12_3_5 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 41] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 33] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  26] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 26] );

                        //y_12_2_6 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 42] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 34] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  27] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 27] );

                        //y_12_1_7 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 43] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 35] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 35];

                        //y_12_0_8 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 44] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 35] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  27] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 27] );

                        //y_11_9_0 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 45] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 36] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 36]
                                      + 8 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  28] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 28] );

                        //y_11_8_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 46] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 36] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 36];

                        //y_11_7_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 47] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 37] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  28] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 28] );

                        //y_11_6_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 48] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 38] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  29] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 29] );

                        //y_11_5_4 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 49] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 39] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  30] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 30] );

                        //y_11_4_5 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 50] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 41] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  33] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 33] );

                        //y_11_3_6 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 51] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 42] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  34] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 34] );

                        //y_11_2_7 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 52] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 43] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  35] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 35] );

                        //y_11_1_8 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 53] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 44] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 44];

                        //y_11_0_9 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 54] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 44] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 44]
                                      + 8 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  35] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 35] );

                        //y_10_10_0 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 55] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 45] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 45]
                                      + 9 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  36] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 36] );

                        //y_10_9_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 56] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 45] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 45];

                        //y_10_8_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 57] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 46] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  36] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 36] );

                        //y_10_7_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 58] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 47] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  37] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 37] );

                        //y_10_6_4 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 59] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 48] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  38] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 38] );

                        //y_10_5_5 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 60] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 49] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  39] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 39] );

                        //y_10_4_6 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 61] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 51] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  42] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 42] );

                        //y_10_3_7 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 62] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 52] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  43] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 43] );

                        //y_10_2_8 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 63] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 53] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  44] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 44] );

                        //y_10_1_9 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 64] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 54] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 54];

                        //y_10_0_10 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 65] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 54] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 54]
                                      + 9 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  44] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 44] );

                        //y_9_11_0 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 66] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 66] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 66]
                                      + 8 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  66] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 66] );

                        //y_9_10_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 67] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 55] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 55];

                        //y_9_9_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 68] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 56] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  45] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 45] );

                        //y_9_8_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 69] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 57] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  46] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 46] );

                        //y_9_7_4 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 70] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 58] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  47] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 47] );

                        //y_9_6_5 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 71] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 59] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 59]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  48] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 48] );

                        //y_9_5_6 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 72] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 61] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 61]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  51] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 51] );

                        //y_9_4_7 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 73] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 62] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  52] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 52] );

                        //y_9_3_8 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 74] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 63] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  53] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 53] );

                        //y_9_2_9 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 75] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 64] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  54] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 54] );

                        //y_9_1_10 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 76] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 65] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 65];

                        //y_9_0_11 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 77] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 77] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 77]
                                      + 8 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  77] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 77] );

                        //y_8_12_0 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 78] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 78] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 78]
                                      + 7 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  78] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 78] );

                        //y_8_11_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 79] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 66] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 66];

                        //y_8_10_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 80] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 67] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  55] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 55] );

                        //y_8_9_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 81] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 68] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  56] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 56] );

                        //y_8_8_4 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 82] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 69] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 69]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  57] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 57] );

                        //y_8_7_5 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 83] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 70] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 70]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  58] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 58] );

                        //y_8_6_6 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 84] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 71] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 71]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  59] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 59] );

                        //y_8_5_7 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 85] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 73] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 73]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  62] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 62] );

                        //y_8_4_8 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 86] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 74] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 74]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  63] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 63] );

                        //y_8_3_9 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 87] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 75] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  64] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 64] );

                        //y_8_2_10 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 88] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 76] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  65] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 65] );

                        //y_8_1_11 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 89] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 77] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 77];

                        //y_8_0_12 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 90] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 90] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 90]
                                      + 7 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  90] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 90] );

                        //y_7_13_0 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 91] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 91] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 91]
                                      + 6 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  91] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 91] );

                        //y_7_12_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 92] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 78] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 78];

                        //y_7_11_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 93] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 79] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  66] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 66] );

                        //y_7_10_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 94] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 80] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 80]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  67] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 67] );

                        //y_7_9_4 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 95] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 81] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 81]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  68] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 68] );

                        //y_7_8_5 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 96] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 82] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 82]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  69] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 69] );

                        //y_7_7_6 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 97] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 83] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 83]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  70] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 70] );

                        //y_7_6_7 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 98] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 85] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 85]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  73] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 73] );

                        //y_7_5_8 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 99] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 86] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 86]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  74] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 74] );

                        //y_7_4_9 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 100] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 87] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 87]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  75] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 75] );

                        //y_7_3_10 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 101] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 88] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 88]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  76] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 76] );

                        //y_7_2_11 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 102] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 89] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  77] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 77] );

                        //y_7_1_12 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 103] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 90] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 90];

                        //y_7_0_13 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 104] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 104] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 104]
                                      + 6 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  104] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 104] );

                        //y_6_14_0 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 105] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 105] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 105]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  105] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 105] );

                        //y_6_13_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 106] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 91] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 91];

                        //y_6_12_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 107] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 92] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 92]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  78] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 78] );

                        //y_6_11_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 108] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 93] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 93]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  79] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 79] );

                        //y_6_10_4 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 109] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 94] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 94]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  80] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 80] );

                        //y_6_9_5 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 110] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 95] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 95]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  81] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 81] );

                        //y_6_8_6 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 111] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 96] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 96]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  82] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 82] );

                        //y_6_7_7 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 112] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 112] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 112]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  112] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 112] );

                        //y_6_6_8 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 113] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 99] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 99]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  86] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 86] );

                        //y_6_5_9 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 114] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 100] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 100]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  87] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 87] );

                        //y_6_4_10 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 115] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 101] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 101]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  88] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 88] );

                        //y_6_3_11 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 116] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 102] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 102]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  89] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 89] );

                        //y_6_2_12 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 117] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 103] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 103]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  90] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 90] );

                        //y_6_1_13 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 118] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 104] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 104];

                        //y_6_0_14 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 119] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 119] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 119]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  119] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 119] );

                        //y_5_15_0 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 120] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 120] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 120]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  120] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 120] );

                        //y_5_14_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 121] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 105] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 105];

                        //y_5_13_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 122] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 106] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 106]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  91] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 91] );

                        //y_5_12_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 123] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 107] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 107]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  92] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 92] );

                        //y_5_11_4 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 124] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 108] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 108]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  93] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 93] );

                        //y_5_10_5 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 125] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 109] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 109]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  94] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 94] );

                        //y_5_9_6 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 126] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 126] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 126]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  126] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 126] );

                        //y_5_8_7 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 127] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 127] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 127]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  127] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 127] );

                        //y_5_7_8 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 128] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 128] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 128]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  128] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 128] );

                        //y_5_6_9 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 129] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 129] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 129]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  129] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 129] );

                        //y_5_5_10 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 130] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 115] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 115]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  101] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 101] );

                        //y_5_4_11 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 131] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 116] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 116]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  102] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 102] );

                        //y_5_3_12 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 132] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 117] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 117]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  103] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 103] );

                        //y_5_2_13 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 133] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 118] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 118]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  104] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 104] );

                        //y_5_1_14 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 134] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 119] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 119];

                        //y_5_0_15 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 135] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 135] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 135]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  135] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 135] );

                        //y_4_16_0 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 136] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 136] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 136]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  136] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 136] );

                        //y_4_15_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 137] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 120] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 120];

                        //y_4_14_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 138] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 121] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 121]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  105] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 105] );

                        //y_4_13_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 139] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 122] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 122]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  106] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 106] );

                        //y_4_12_4 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 140] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 123] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 123]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  107] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 107] );

                        //y_4_11_5 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 141] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 141] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 141]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  141] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 141] );

                        //y_4_10_6 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 142] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 142] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 142]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  142] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 142] );

                        //y_4_9_7 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 143] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 143] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 143]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  143] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 143] );

                        //y_4_8_8 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 144] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 144] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 144]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  144] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 144] );

                        //y_4_7_9 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 145] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 145] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 145]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  145] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 145] );

                        //y_4_6_10 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 146] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 146] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 146]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  146] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 146] );

                        //y_4_5_11 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 147] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 147] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 147]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  147] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 147] );

                        //y_4_4_12 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 148] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 132] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 132]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  117] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 117] );

                        //y_4_3_13 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 149] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 133] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 133]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  118] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 118] );

                        //y_4_2_14 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 150] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 134] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 134]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  119] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 119] );

                        //y_4_1_15 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 151] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 135] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 135];

                        //y_4_0_16 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 152] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 152] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 152]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  152] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 152] );

                        //y_3_17_0 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 153] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 153] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 153]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  153] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 153] );

                        //y_3_16_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 154] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 136] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 136];

                        //y_3_15_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 155] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 137] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 137]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  120] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 120] );

                        //y_3_14_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 156] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 138] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 138]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  121] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 121] );

                        //y_3_13_4 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 157] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 157] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 157]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  157] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 157] );

                        //y_3_12_5 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 158] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 158] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 158]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  158] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 158] );

                        //y_3_11_6 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 159] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 159] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 159]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  159] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 159] );

                        //y_3_10_7 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 160] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 160] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 160]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  160] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 160] );

                        //y_3_9_8 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 161] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 161] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 161]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  161] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 161] );

                        //y_3_8_9 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 162] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 162] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 162]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  162] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 162] );

                        //y_3_7_10 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 163] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 163] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 163]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  163] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 163] );

                        //y_3_6_11 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 164] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 164] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 164]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  164] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 164] );

                        //y_3_5_12 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 165] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 165] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 165]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  165] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 165] );

                        //y_3_4_13 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 166] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 166] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 166]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  166] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 166] );

                        //y_3_3_14 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 167] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 150] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 150]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  134] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 134] );

                        //y_3_2_15 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 168] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 151] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 151]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  135] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 135] );

                        //y_3_1_16 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 169] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 152] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 152];

                        //y_3_0_17 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 170] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 170] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 170]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  170] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 170] );

                        //y_2_18_0 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 171] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 171] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 171]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  171] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 171] );

                        //y_2_17_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 172] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 153] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 153];

                        //y_2_16_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 173] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 154] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 154]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  136] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 136] );

                        //y_2_15_3 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 174] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 174] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 174]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  174] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 174] );

                        //y_2_14_4 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 175] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 175] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 175]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  175] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 175] );

                        //y_2_13_5 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 176] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 176] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 176]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  176] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 176] );

                        //y_2_12_6 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 177] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 177] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 177]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  177] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 177] );

                        //y_2_11_7 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 178] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 178] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 178]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  178] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 178] );

                        //y_2_10_8 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 179] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 179] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 179]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  179] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 179] );

                        //y_2_9_9 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 180] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 180] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 180]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  180] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 180] );

                        //y_2_8_10 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 181] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 181] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 181]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  181] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 181] );

                        //y_2_7_11 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 182] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 182] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 182]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  182] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 182] );

                        //y_2_6_12 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 183] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 183] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 183]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  183] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 183] );

                        //y_2_5_13 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 184] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 184] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 184]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  184] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 184] );

                        //y_2_4_14 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 185] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 185] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 185]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  185] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 185] );

                        //y_2_3_15 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 186] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 186] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 186]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  186] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 186] );

                        //y_2_2_16 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 187] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 169] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 169]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  152] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 152] );

                        //y_2_1_17 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 188] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 170] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 170];

                        //y_2_0_18 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 189] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 189] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 189]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  189] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 189] );

                        //y_1_19_0 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 190] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 190] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 190];

                        //y_1_18_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 191] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 171] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 171];

                        //y_1_17_2 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 192] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 192] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 192];

                        //y_1_16_3 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 193] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 193] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 193];

                        //y_1_15_4 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 194] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 194] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 194];

                        //y_1_14_5 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 195] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 195] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 195];

                        //y_1_13_6 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 196] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 196] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 196];

                        //y_1_12_7 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 197] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 197] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 197];

                        //y_1_11_8 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 198] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 198] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 198];

                        //y_1_10_9 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 199] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 199] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 199];

                        //y_1_9_10 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 200] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 200] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 200];

                        //y_1_8_11 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 201] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 201] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 201];

                        //y_1_7_12 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 202] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 202] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 202];

                        //y_1_6_13 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 203] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 203] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 203];

                        //y_1_5_14 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 204] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 204] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 204];

                        //y_1_4_15 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 205] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 205] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 205];

                        //y_1_3_16 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 206] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 206] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 206];

                        //y_1_2_17 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 207] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 207] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 207];

                        //y_1_1_18 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 208] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 189] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 189];

                        //y_1_0_19 : STEP: x
                        AUX_INT__y_s_s_s[m * 231 + 209] = P_PA_x * AUX_INT__x_s_s_s[m * 210 + 209] - aop_PQ_x * AUX_INT__x_s_s_s[(m+1) * 210 + 209];

                        //y_0_20_0 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 210] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 190] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 190]
                                      + 19 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  171] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 171] );

                        //y_0_19_1 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 211] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 190] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 190];

                        //y_0_18_2 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 212] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 191] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 191]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  171] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 171] );

                        //y_0_17_3 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 213] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 192] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 192]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  172] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 172] );

                        //y_0_16_4 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 214] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 193] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 193]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  173] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 173] );

                        //y_0_15_5 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 215] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 194] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 194]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  174] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 174] );

                        //y_0_14_6 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 216] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 195] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 195]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  175] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 175] );

                        //y_0_13_7 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 217] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 196] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 196]
                                      + 6 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  176] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 176] );

                        //y_0_12_8 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 218] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 197] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 197]
                                      + 7 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  177] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 177] );

                        //y_0_11_9 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 219] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 198] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 198]
                                      + 8 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  178] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 178] );

                        //y_0_10_10 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 220] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 199] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 199]
                                      + 9 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  179] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 179] );

                        //y_0_9_11 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 221] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 201] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 201]
                                      + 8 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  182] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 182] );

                        //y_0_8_12 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 222] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 202] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 202]
                                      + 7 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  183] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 183] );

                        //y_0_7_13 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 223] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 203] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 203]
                                      + 6 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  184] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 184] );

                        //y_0_6_14 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 224] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 204] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 204]
                                      + 5 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  185] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 185] );

                        //y_0_5_15 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 225] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 205] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 205]
                                      + 4 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  186] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 186] );

                        //y_0_4_16 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 226] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 206] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 206]
                                      + 3 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  187] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 187] );

                        //y_0_3_17 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 227] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 207] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 207]
                                      + 2 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  188] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 188] );

                        //y_0_2_18 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 228] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 208] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 208]
                                      + 1 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  189] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 189] );

                        //y_0_1_19 : STEP: y
                        AUX_INT__y_s_s_s[m * 231 + 229] = P_PA_y * AUX_INT__x_s_s_s[m * 210 + 209] - aop_PQ_y * AUX_INT__x_s_s_s[(m+1) * 210 + 209];

                        //y_0_0_20 : STEP: z
                        AUX_INT__y_s_s_s[m * 231 + 230] = P_PA_z * AUX_INT__x_s_s_s[m * 210 + 209] - aop_PQ_z * AUX_INT__x_s_s_s[(m+1) * 210 + 209]
                                      + 19 * one_over_2p * ( AUX_INT__w_s_s_s[m * 190 +  189] - a_over_p * AUX_INT__w_s_s_s[(m+1) * 190 + 189] );

                    }
}
