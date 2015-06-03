//////////////////////////////////////////////
// VRR functions
//////////////////////////////////////////////




// VRR to obtain AUX_INT__p_s_s_s
void VRR_AUX_INT__p_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__p_s_s_s,
                        double const * const restrict AUX_INT__s_s_s_s)
{
                    // Forming AUX_INT__p_s_s_s[num_m * 3];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //P_100 : STEP: x
                        AUX_INT__p_s_s_s[m * 3 + 0] = P_PA_x * AUX_INT__s_s_s_s[m * 1 + 0] - aop_PQ_x * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_010 : STEP: y
                        AUX_INT__p_s_s_s[m * 3 + 1] = P_PA_y * AUX_INT__s_s_s_s[m * 1 + 0] - aop_PQ_y * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_001 : STEP: z
                        AUX_INT__p_s_s_s[m * 3 + 2] = P_PA_z * AUX_INT__s_s_s_s[m * 1 + 0] - aop_PQ_z * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                    }
}



// VRR to obtain AUX_INT__d_s_s_s
void VRR_AUX_INT__d_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__d_s_s_s,
                        double const * const restrict AUX_INT__p_s_s_s,
                        double const * const restrict AUX_INT__s_s_s_s)
{
                    // Forming AUX_INT__d_s_s_s[num_m * 6];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //D_200 : STEP: x
                        AUX_INT__d_s_s_s[m * 6 + 0] = P_PA_x * AUX_INT__p_s_s_s[m * 3 + 0] - aop_PQ_x * AUX_INT__p_s_s_s[(m+1) * 3 + 0]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                        //D_110 : STEP: y
                        AUX_INT__d_s_s_s[m * 6 + 1] = P_PA_y * AUX_INT__p_s_s_s[m * 3 + 0] - aop_PQ_y * AUX_INT__p_s_s_s[(m+1) * 3 + 0];

                        //D_101 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 2] = P_PA_z * AUX_INT__p_s_s_s[m * 3 + 0] - aop_PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 0];

                        //D_020 : STEP: y
                        AUX_INT__d_s_s_s[m * 6 + 3] = P_PA_y * AUX_INT__p_s_s_s[m * 3 + 1] - aop_PQ_y * AUX_INT__p_s_s_s[(m+1) * 3 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                        //D_011 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 4] = P_PA_z * AUX_INT__p_s_s_s[m * 3 + 1] - aop_PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 1];

                        //D_002 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 5] = P_PA_z * AUX_INT__p_s_s_s[m * 3 + 2] - aop_PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                    }
}



// VRR to obtain AUX_INT__f_s_s_s
void VRR_AUX_INT__f_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__f_s_s_s,
                        double const * const restrict AUX_INT__d_s_s_s,
                        double const * const restrict AUX_INT__p_s_s_s)
{
                    // Forming AUX_INT__f_s_s_s[num_m * 10];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //F_300 : STEP: x
                        AUX_INT__f_s_s_s[m * 10 + 0] = P_PA_x * AUX_INT__d_s_s_s[m * 6 + 0] - aop_PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 0]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  0] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 0] );

                        //F_210 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 1] = P_PA_y * AUX_INT__d_s_s_s[m * 6 + 0] - aop_PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 0];

                        //F_201 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 2] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 0] - aop_PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 0];

                        //F_120 : STEP: x
                        AUX_INT__f_s_s_s[m * 10 + 3] = P_PA_x * AUX_INT__d_s_s_s[m * 6 + 3] - aop_PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 3];

                        //F_111 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 4] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 1] - aop_PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 1];

                        //F_102 : STEP: x
                        AUX_INT__f_s_s_s[m * 10 + 5] = P_PA_x * AUX_INT__d_s_s_s[m * 6 + 5] - aop_PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 5];

                        //F_030 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 6] = P_PA_y * AUX_INT__d_s_s_s[m * 6 + 3] - aop_PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  1] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 1] );

                        //F_021 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 7] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 3] - aop_PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 3];

                        //F_012 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 8] = P_PA_y * AUX_INT__d_s_s_s[m * 6 + 5] - aop_PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 5];

                        //F_003 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 9] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 5] - aop_PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  2] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 2] );

                    }
}



// VRR to obtain AUX_INT__g_s_s_s
void VRR_AUX_INT__g_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__g_s_s_s,
                        double const * const restrict AUX_INT__f_s_s_s,
                        double const * const restrict AUX_INT__d_s_s_s)
{
                    // Forming AUX_INT__g_s_s_s[num_m * 15];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //G_400 : STEP: x
                        AUX_INT__g_s_s_s[m * 15 + 0] = P_PA_x * AUX_INT__f_s_s_s[m * 10 + 0] - aop_PQ_x * AUX_INT__f_s_s_s[(m+1) * 10 + 0]
                                      + 3 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  0] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 0] );

                        //G_310 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 1] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 0] - aop_PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 0];

                        //G_301 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 2] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 0] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 0];

                        //G_220 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 3] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 1] - aop_PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  0] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 0] );

                        //G_211 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 4] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 1] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 1];

                        //G_202 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 5] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 2] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  0] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 0] );

                        //G_130 : STEP: x
                        AUX_INT__g_s_s_s[m * 15 + 6] = P_PA_x * AUX_INT__f_s_s_s[m * 10 + 6] - aop_PQ_x * AUX_INT__f_s_s_s[(m+1) * 10 + 6];

                        //G_121 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 7] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 3] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 3];

                        //G_112 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 8] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 5] - aop_PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 5];

                        //G_103 : STEP: x
                        AUX_INT__g_s_s_s[m * 15 + 9] = P_PA_x * AUX_INT__f_s_s_s[m * 10 + 9] - aop_PQ_x * AUX_INT__f_s_s_s[(m+1) * 10 + 9];

                        //G_040 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 10] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 6] - aop_PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  3] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 3] );

                        //G_031 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 11] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 6] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 6];

                        //G_022 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 12] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 7] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  3] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 3] );

                        //G_013 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 13] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 9] - aop_PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 9];

                        //G_004 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 14] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 9] - aop_PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  5] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 5] );

                    }
}



// VRR to obtain AUX_INT__h_s_s_s
void VRR_AUX_INT__h_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__h_s_s_s,
                        double const * const restrict AUX_INT__g_s_s_s,
                        double const * const restrict AUX_INT__f_s_s_s)
{
                    // Forming AUX_INT__h_s_s_s[num_m * 21];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //H_500 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 0] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 0] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 0]
                                      + 4 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  0] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 0] );

                        //H_410 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 1] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 0] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 0];

                        //H_401 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 2] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 0] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 0];

                        //H_320 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 3] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 1] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  0] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 0] );

                        //H_311 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 4] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 1] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 1];

                        //H_302 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 5] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 2] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  0] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 0] );

                        //H_230 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 6] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 6] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 6]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  6] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 6] );

                        //H_221 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 7] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 3] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 3];

                        //H_212 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 8] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 5] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 5];

                        //H_203 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 9] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 9] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 9]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  9] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 9] );

                        //H_140 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 10] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 10] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 10];

                        //H_131 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 11] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 6] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 6];

                        //H_122 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 12] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 12] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 12];

                        //H_113 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 13] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 9] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 9];

                        //H_104 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 14] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 14] - aop_PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 14];

                        //H_050 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 15] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 10] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  6] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 6] );

                        //H_041 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 16] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 10] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 10];

                        //H_032 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 17] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 11] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  6] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 6] );

                        //H_023 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 18] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 13] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  9] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 9] );

                        //H_014 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 19] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 14] - aop_PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 14];

                        //H_005 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 20] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 14] - aop_PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  9] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 9] );

                    }
}



// VRR to obtain AUX_INT__i_s_s_s
void VRR_AUX_INT__i_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__i_s_s_s,
                        double const * const restrict AUX_INT__h_s_s_s,
                        double const * const restrict AUX_INT__g_s_s_s)
{
                    // Forming AUX_INT__i_s_s_s[num_m * 28];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //I_600 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 0] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 0] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 0]
                                      + 5 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  0] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 0] );

                        //I_510 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 1] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 0] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 0];

                        //I_501 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 2] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 0] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 0];

                        //I_420 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 3] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 1] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  0] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 0] );

                        //I_411 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 4] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 1] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 1];

                        //I_402 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 5] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 2] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  0] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 0] );

                        //I_330 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 6] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 3] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  1] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 1] );

                        //I_321 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 7] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 3] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 3];

                        //I_312 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 8] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 5] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 5];

                        //I_303 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 9] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 5] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  2] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 2] );

                        //I_240 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 10] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 10] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 10]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  10] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 10] );

                        //I_231 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 11] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 6] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 6];

                        //I_222 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 12] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 7] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  3] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 3] );

                        //I_213 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 13] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 9] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 9];

                        //I_204 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 14] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 14] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 14]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  14] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 14] );

                        //I_150 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 15] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 15] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 15];

                        //I_141 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 16] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 10] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 10];

                        //I_132 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 17] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 17] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 17];

                        //I_123 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 18] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 18] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 18];

                        //I_114 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 19] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 14] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 14];

                        //I_105 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 20] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 20] - aop_PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 20];

                        //I_060 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 21] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 15] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  10] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 10] );

                        //I_051 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 22] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 15] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 15];

                        //I_042 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 23] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 16] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  10] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 10] );

                        //I_033 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 24] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 17] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  11] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 11] );

                        //I_024 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 25] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 19] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  14] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 14] );

                        //I_015 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 26] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 20] - aop_PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 20];

                        //I_006 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 27] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 20] - aop_PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  14] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 14] );

                    }
}



// VRR to obtain AUX_INT__j_s_s_s
void VRR_AUX_INT__j_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__j_s_s_s,
                        double const * const restrict AUX_INT__i_s_s_s,
                        double const * const restrict AUX_INT__h_s_s_s)
{
                    // Forming AUX_INT__j_s_s_s[num_m * 36];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //J_700 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 0] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 0] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 0]
                                      + 6 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  0] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 0] );

                        //J_610 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 1] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 0] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 0];

                        //J_601 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 2] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 0] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 0];

                        //J_520 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 3] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 1] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  0] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 0] );

                        //J_511 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 4] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 1] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 1];

                        //J_502 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 5] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 2] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  0] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 0] );

                        //J_430 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 6] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 3] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  1] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 1] );

                        //J_421 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 7] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 3] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 3];

                        //J_412 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 8] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 5] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 5];

                        //J_403 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 9] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 5] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  2] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 2] );

                        //J_340 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 10] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 10] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 10]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  10] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 10] );

                        //J_331 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 11] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 6] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 6];

                        //J_322 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 12] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 7] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  3] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 3] );

                        //J_313 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 13] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 9] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 9];

                        //J_304 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 14] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 14] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 14]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  14] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 14] );

                        //J_250 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 15] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 15] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 15]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  15] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 15] );

                        //J_241 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 16] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 10] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 10];

                        //J_232 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 17] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 11] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  6] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 6] );

                        //J_223 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 18] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 13] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  9] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 9] );

                        //J_214 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 19] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 14] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 14];

                        //J_205 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 20] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 20] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 20]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  20] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 20] );

                        //J_160 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 21] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 21] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 21];

                        //J_151 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 22] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 15] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 15];

                        //J_142 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 23] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 23] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 23];

                        //J_133 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 24] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 24] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 24];

                        //J_124 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 25] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 25] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 25];

                        //J_115 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 26] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 20] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 20];

                        //J_106 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 27] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 27] - aop_PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 27];

                        //J_070 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 28] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 21] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  15] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 15] );

                        //J_061 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 29] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 21] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 21];

                        //J_052 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 30] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 22] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  15] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 15] );

                        //J_043 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 31] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 23] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  16] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 16] );

                        //J_034 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 32] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 25] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  19] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 19] );

                        //J_025 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 33] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 26] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  20] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 20] );

                        //J_016 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 34] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 27] - aop_PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 27];

                        //J_007 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 35] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 27] - aop_PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  20] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 20] );

                    }
}



// VRR to obtain AUX_INT__k_s_s_s
void VRR_AUX_INT__k_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__k_s_s_s,
                        double const * const restrict AUX_INT__j_s_s_s,
                        double const * const restrict AUX_INT__i_s_s_s)
{
                    // Forming AUX_INT__k_s_s_s[num_m * 45];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //K_800 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 0] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 0] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 0]
                                      + 7 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  0] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 0] );

                        //K_710 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 1] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 0] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 0];

                        //K_701 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 2] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 0] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 0];

                        //K_620 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 3] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 1] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  0] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 0] );

                        //K_611 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 4] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 1] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 1];

                        //K_602 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 5] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 2] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  0] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 0] );

                        //K_530 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 6] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 3] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  1] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 1] );

                        //K_521 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 7] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 3] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 3];

                        //K_512 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 8] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 5] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 5];

                        //K_503 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 9] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 5] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  2] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 2] );

                        //K_440 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 10] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 6] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  3] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 3] );

                        //K_431 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 11] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 6] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 6];

                        //K_422 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 12] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 7] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  3] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 3] );

                        //K_413 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 13] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 9] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 9];

                        //K_404 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 14] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 9] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  5] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 5] );

                        //K_350 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 15] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 15] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 15]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  15] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 15] );

                        //K_341 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 16] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 10] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 10];

                        //K_332 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 17] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 11] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  6] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 6] );

                        //K_323 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 18] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 13] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  9] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 9] );

                        //K_314 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 19] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 14] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 14];

                        //K_305 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 20] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 20] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 20]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  20] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 20] );

                        //K_260 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 21] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 21] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 21]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  21] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 21] );

                        //K_251 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 22] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 15] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 15];

                        //K_242 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 23] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 16] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  10] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 10] );

                        //K_233 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 24] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 24] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 24]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  24] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 24] );

                        //K_224 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 25] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 19] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  14] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 14] );

                        //K_215 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 26] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 20] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 20];

                        //K_206 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 27] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 27] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 27]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  27] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 27] );

                        //K_170 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 28] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 28] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 28];

                        //K_161 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 29] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 21] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 21];

                        //K_152 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 30] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 30] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 30];

                        //K_143 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 31] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 31] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 31];

                        //K_134 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 32] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 32] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 32];

                        //K_125 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 33] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 33] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 33];

                        //K_116 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 34] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 27] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 27];

                        //K_107 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 35] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 35] - aop_PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 35];

                        //K_080 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 36] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 28] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  21] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 21] );

                        //K_071 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 37] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 28] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 28];

                        //K_062 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 38] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 29] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  21] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 21] );

                        //K_053 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 39] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 30] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  22] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 22] );

                        //K_044 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 40] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 31] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  23] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 23] );

                        //K_035 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 41] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 33] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  26] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 26] );

                        //K_026 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 42] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 34] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  27] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 27] );

                        //K_017 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 43] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 35] - aop_PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 35];

                        //K_008 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 44] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 35] - aop_PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  27] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 27] );

                    }
}



// VRR to obtain AUX_INT__l_s_s_s
void VRR_AUX_INT__l_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__l_s_s_s,
                        double const * const restrict AUX_INT__k_s_s_s,
                        double const * const restrict AUX_INT__j_s_s_s)
{
                    // Forming AUX_INT__l_s_s_s[num_m * 55];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //L_900 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 0] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 0] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 0]
                                      + 8 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  0] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 0] );

                        //L_810 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 1] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 0] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 0];

                        //L_801 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 2] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 0] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 0];

                        //L_720 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 3] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 1] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  0] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 0] );

                        //L_711 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 4] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 1] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 1];

                        //L_702 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 5] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 2] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  0] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 0] );

                        //L_630 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 6] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 3] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  1] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 1] );

                        //L_621 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 7] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 3] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 3];

                        //L_612 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 8] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 5] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 5];

                        //L_603 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 9] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 5] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  2] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 2] );

                        //L_540 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 10] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 6] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  3] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 3] );

                        //L_531 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 11] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 6] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 6];

                        //L_522 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 12] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 7] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  3] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 3] );

                        //L_513 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 13] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 9] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 9];

                        //L_504 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 14] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 9] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  5] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 5] );

                        //L_450 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 15] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 15] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 15]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  15] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 15] );

                        //L_441 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 16] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 10] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 10];

                        //L_432 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 17] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 11] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  6] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 6] );

                        //L_423 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 18] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 13] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  9] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 9] );

                        //L_414 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 19] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 14] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 14];

                        //L_405 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 20] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 20] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 20]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  20] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 20] );

                        //L_360 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 21] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 21] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 21]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  21] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 21] );

                        //L_351 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 22] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 15] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 15];

                        //L_342 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 23] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 16] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  10] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 10] );

                        //L_333 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 24] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 17] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  11] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 11] );

                        //L_324 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 25] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 19] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  14] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 14] );

                        //L_315 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 26] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 20] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 20];

                        //L_306 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 27] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 27] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 27]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  27] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 27] );

                        //L_270 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 28] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 28] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 28]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  28] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 28] );

                        //L_261 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 29] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 21] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 21];

                        //L_252 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 30] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 22] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  15] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 15] );

                        //L_243 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 31] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 31] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 31]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  31] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 31] );

                        //L_234 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 32] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 32] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 32]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  32] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 32] );

                        //L_225 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 33] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 26] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  20] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 20] );

                        //L_216 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 34] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 27] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 27];

                        //L_207 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 35] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 35] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 35]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  35] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 35] );

                        //L_180 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 36] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 36] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 36];

                        //L_171 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 37] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 28] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 28];

                        //L_162 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 38] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 38] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 38];

                        //L_153 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 39] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 39] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 39];

                        //L_144 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 40] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 40] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 40];

                        //L_135 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 41] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 41] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 41];

                        //L_126 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 42] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 42] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 42];

                        //L_117 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 43] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 35] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 35];

                        //L_108 : STEP: x
                        AUX_INT__l_s_s_s[m * 55 + 44] = P_PA_x * AUX_INT__k_s_s_s[m * 45 + 44] - aop_PQ_x * AUX_INT__k_s_s_s[(m+1) * 45 + 44];

                        //L_090 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 45] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 36] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 36]
                                      + 8 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  28] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 28] );

                        //L_081 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 46] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 36] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 36];

                        //L_072 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 47] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 37] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  28] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 28] );

                        //L_063 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 48] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 38] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  29] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 29] );

                        //L_054 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 49] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 39] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  30] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 30] );

                        //L_045 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 50] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 41] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  33] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 33] );

                        //L_036 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 51] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 42] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  34] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 34] );

                        //L_027 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 52] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 43] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  35] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 35] );

                        //L_018 : STEP: y
                        AUX_INT__l_s_s_s[m * 55 + 53] = P_PA_y * AUX_INT__k_s_s_s[m * 45 + 44] - aop_PQ_y * AUX_INT__k_s_s_s[(m+1) * 45 + 44];

                        //L_009 : STEP: z
                        AUX_INT__l_s_s_s[m * 55 + 54] = P_PA_z * AUX_INT__k_s_s_s[m * 45 + 44] - aop_PQ_z * AUX_INT__k_s_s_s[(m+1) * 45 + 44]
                                      + 8 * one_over_2p * ( AUX_INT__j_s_s_s[m * 36 +  35] - a_over_p * AUX_INT__j_s_s_s[(m+1) * 36 + 35] );

                    }
}



// VRR to obtain AUX_INT__m_s_s_s
void VRR_AUX_INT__m_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__m_s_s_s,
                        double const * const restrict AUX_INT__l_s_s_s,
                        double const * const restrict AUX_INT__k_s_s_s)
{
                    // Forming AUX_INT__m_s_s_s[num_m * 66];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //M_1000 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 0] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 0] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 0]
                                      + 9 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  0] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 0] );

                        //M_910 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 1] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 0] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 0];

                        //M_901 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 2] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 0] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 0];

                        //M_820 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 3] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 1] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  0] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 0] );

                        //M_811 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 4] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 1] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 1];

                        //M_802 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 5] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 2] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  0] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 0] );

                        //M_730 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 6] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 3] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  1] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 1] );

                        //M_721 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 7] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 3] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 3];

                        //M_712 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 8] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 5] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 5];

                        //M_703 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 9] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 5] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  2] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 2] );

                        //M_640 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 10] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 6] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  3] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 3] );

                        //M_631 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 11] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 6] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 6];

                        //M_622 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 12] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 7] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  3] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 3] );

                        //M_613 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 13] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 9] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 9];

                        //M_604 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 14] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 9] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  5] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 5] );

                        //M_550 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 15] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 10] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  6] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 6] );

                        //M_541 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 16] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 10] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 10];

                        //M_532 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 17] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 11] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  6] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 6] );

                        //M_523 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 18] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 13] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  9] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 9] );

                        //M_514 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 19] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 14] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 14];

                        //M_505 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 20] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 14] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  9] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 9] );

                        //M_460 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 21] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 21] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 21]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  21] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 21] );

                        //M_451 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 22] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 15] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 15];

                        //M_442 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 23] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 16] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  10] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 10] );

                        //M_433 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 24] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 17] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  11] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 11] );

                        //M_424 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 25] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 19] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  14] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 14] );

                        //M_415 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 26] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 20] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 20];

                        //M_406 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 27] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 27] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 27]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  27] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 27] );

                        //M_370 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 28] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 28] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 28]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  28] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 28] );

                        //M_361 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 29] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 21] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 21];

                        //M_352 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 30] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 22] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  15] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 15] );

                        //M_343 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 31] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 23] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  16] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 16] );

                        //M_334 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 32] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 25] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  19] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 19] );

                        //M_325 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 33] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 26] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  20] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 20] );

                        //M_316 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 34] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 27] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 27];

                        //M_307 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 35] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 35] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 35]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  35] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 35] );

                        //M_280 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 36] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 36] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 36]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  36] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 36] );

                        //M_271 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 37] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 28] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 28];

                        //M_262 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 38] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 29] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  21] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 21] );

                        //M_253 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 39] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 39] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 39]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  39] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 39] );

                        //M_244 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 40] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 40] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 40]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  40] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 40] );

                        //M_235 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 41] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 41] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 41]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  41] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 41] );

                        //M_226 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 42] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 34] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  27] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 27] );

                        //M_217 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 43] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 35] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 35];

                        //M_208 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 44] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 44] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 44]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  44] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 44] );

                        //M_190 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 45] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 45] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 45];

                        //M_181 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 46] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 36] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 36];

                        //M_172 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 47] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 47] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 47];

                        //M_163 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 48] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 48] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 48];

                        //M_154 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 49] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 49] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 49];

                        //M_145 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 50] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 50] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 50];

                        //M_136 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 51] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 51] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 51];

                        //M_127 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 52] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 52] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 52];

                        //M_118 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 53] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 44] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 44];

                        //M_109 : STEP: x
                        AUX_INT__m_s_s_s[m * 66 + 54] = P_PA_x * AUX_INT__l_s_s_s[m * 55 + 54] - aop_PQ_x * AUX_INT__l_s_s_s[(m+1) * 55 + 54];

                        //M_0100 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 55] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 45] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 45]
                                      + 9 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  36] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 36] );

                        //M_091 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 56] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 45] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 45];

                        //M_082 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 57] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 46] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  36] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 36] );

                        //M_073 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 58] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 47] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  37] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 37] );

                        //M_064 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 59] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 48] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  38] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 38] );

                        //M_055 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 60] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 49] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  39] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 39] );

                        //M_046 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 61] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 51] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  42] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 42] );

                        //M_037 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 62] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 52] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  43] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 43] );

                        //M_028 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 63] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 53] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  44] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 44] );

                        //M_019 : STEP: y
                        AUX_INT__m_s_s_s[m * 66 + 64] = P_PA_y * AUX_INT__l_s_s_s[m * 55 + 54] - aop_PQ_y * AUX_INT__l_s_s_s[(m+1) * 55 + 54];

                        //M_0010 : STEP: z
                        AUX_INT__m_s_s_s[m * 66 + 65] = P_PA_z * AUX_INT__l_s_s_s[m * 55 + 54] - aop_PQ_z * AUX_INT__l_s_s_s[(m+1) * 55 + 54]
                                      + 9 * one_over_2p * ( AUX_INT__k_s_s_s[m * 45 +  44] - a_over_p * AUX_INT__k_s_s_s[(m+1) * 45 + 44] );

                    }
}



// VRR to obtain AUX_INT__n_s_s_s
void VRR_AUX_INT__n_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__n_s_s_s,
                        double const * const restrict AUX_INT__m_s_s_s,
                        double const * const restrict AUX_INT__l_s_s_s)
{
                    // Forming AUX_INT__n_s_s_s[num_m * 78];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //N_1100 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 0] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 0] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 0]
                                      + 10 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  0] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 0] );

                        //N_1010 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 1] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 0] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 0];

                        //N_1001 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 2] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 0] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 0];

                        //N_920 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 3] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 1] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  0] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 0] );

                        //N_911 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 4] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 1] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 1];

                        //N_902 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 5] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 2] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  0] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 0] );

                        //N_830 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 6] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 3] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  1] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 1] );

                        //N_821 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 7] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 3] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 3];

                        //N_812 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 8] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 5] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 5];

                        //N_803 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 9] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 5] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  2] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 2] );

                        //N_740 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 10] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 6] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  3] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 3] );

                        //N_731 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 11] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 6] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 6];

                        //N_722 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 12] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 7] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  3] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 3] );

                        //N_713 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 13] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 9] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 9];

                        //N_704 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 14] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 9] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  5] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 5] );

                        //N_650 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 15] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 10] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  6] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 6] );

                        //N_641 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 16] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 10] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 10];

                        //N_632 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 17] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 11] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  6] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 6] );

                        //N_623 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 18] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 13] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  9] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 9] );

                        //N_614 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 19] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 14] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 14];

                        //N_605 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 20] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 14] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  9] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 9] );

                        //N_560 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 21] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 21] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 21]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  21] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 21] );

                        //N_551 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 22] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 15] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 15];

                        //N_542 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 23] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 16] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  10] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 10] );

                        //N_533 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 24] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 17] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  11] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 11] );

                        //N_524 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 25] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 19] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  14] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 14] );

                        //N_515 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 26] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 20] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 20];

                        //N_506 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 27] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 27] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 27]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  27] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 27] );

                        //N_470 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 28] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 28] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 28]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  28] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 28] );

                        //N_461 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 29] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 21] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 21];

                        //N_452 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 30] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 22] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  15] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 15] );

                        //N_443 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 31] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 23] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  16] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 16] );

                        //N_434 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 32] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 25] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  19] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 19] );

                        //N_425 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 33] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 26] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  20] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 20] );

                        //N_416 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 34] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 27] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 27];

                        //N_407 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 35] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 35] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 35]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  35] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 35] );

                        //N_380 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 36] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 36] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 36]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  36] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 36] );

                        //N_371 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 37] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 28] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 28];

                        //N_362 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 38] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 29] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  21] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 21] );

                        //N_353 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 39] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 30] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  22] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 22] );

                        //N_344 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 40] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 40] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 40]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  40] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 40] );

                        //N_335 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 41] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 33] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  26] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 26] );

                        //N_326 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 42] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 34] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  27] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 27] );

                        //N_317 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 43] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 35] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 35];

                        //N_308 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 44] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 44] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 44]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  44] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 44] );

                        //N_290 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 45] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 45] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 45]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  45] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 45] );

                        //N_281 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 46] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 36] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 36];

                        //N_272 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 47] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 37] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  28] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 28] );

                        //N_263 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 48] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 48] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 48]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  48] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 48] );

                        //N_254 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 49] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 49] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 49]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  49] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 49] );

                        //N_245 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 50] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 50] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 50]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  50] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 50] );

                        //N_236 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 51] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 51] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 51]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  51] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 51] );

                        //N_227 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 52] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 43] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  35] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 35] );

                        //N_218 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 53] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 44] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 44];

                        //N_209 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 54] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 54] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 54]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  54] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 54] );

                        //N_1100 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 55] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 55] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 55];

                        //N_191 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 56] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 45] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 45];

                        //N_182 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 57] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 57] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 57];

                        //N_173 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 58] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 58] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 58];

                        //N_164 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 59] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 59] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 59];

                        //N_155 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 60] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 60] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 60];

                        //N_146 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 61] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 61] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 61];

                        //N_137 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 62] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 62] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 62];

                        //N_128 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 63] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 63] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 63];

                        //N_119 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 64] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 54] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 54];

                        //N_1010 : STEP: x
                        AUX_INT__n_s_s_s[m * 78 + 65] = P_PA_x * AUX_INT__m_s_s_s[m * 66 + 65] - aop_PQ_x * AUX_INT__m_s_s_s[(m+1) * 66 + 65];

                        //N_0110 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 66] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 55] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 55]
                                      + 10 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  45] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 45] );

                        //N_0101 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 67] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 55] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 55];

                        //N_092 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 68] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 56] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  45] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 45] );

                        //N_083 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 69] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 57] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  46] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 46] );

                        //N_074 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 70] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 58] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  47] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 47] );

                        //N_065 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 71] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 59] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 59]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  48] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 48] );

                        //N_056 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 72] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 61] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 61]
                                      + 4 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  51] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 51] );

                        //N_047 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 73] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 62] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  52] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 52] );

                        //N_038 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 74] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 63] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  53] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 53] );

                        //N_029 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 75] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 64] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  54] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 54] );

                        //N_0110 : STEP: y
                        AUX_INT__n_s_s_s[m * 78 + 76] = P_PA_y * AUX_INT__m_s_s_s[m * 66 + 65] - aop_PQ_y * AUX_INT__m_s_s_s[(m+1) * 66 + 65];

                        //N_0011 : STEP: z
                        AUX_INT__n_s_s_s[m * 78 + 77] = P_PA_z * AUX_INT__m_s_s_s[m * 66 + 65] - aop_PQ_z * AUX_INT__m_s_s_s[(m+1) * 66 + 65]
                                      + 10 * one_over_2p * ( AUX_INT__l_s_s_s[m * 55 +  54] - a_over_p * AUX_INT__l_s_s_s[(m+1) * 55 + 54] );

                    }
}



// VRR to obtain AUX_INT__o_s_s_s
void VRR_AUX_INT__o_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__o_s_s_s,
                        double const * const restrict AUX_INT__n_s_s_s,
                        double const * const restrict AUX_INT__m_s_s_s)
{
                    // Forming AUX_INT__o_s_s_s[num_m * 91];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //O_1200 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 0] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 0] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 0]
                                      + 11 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  0] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 0] );

                        //O_1110 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 1] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 0] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 0];

                        //O_1101 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 2] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 0] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 0];

                        //O_1020 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 3] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 1] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  0] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 0] );

                        //O_1011 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 4] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 1] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 1];

                        //O_1002 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 5] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 2] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  0] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 0] );

                        //O_930 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 6] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 3] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  1] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 1] );

                        //O_921 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 7] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 3] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 3];

                        //O_912 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 8] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 5] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 5];

                        //O_903 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 9] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 5] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  2] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 2] );

                        //O_840 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 10] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 6] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  3] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 3] );

                        //O_831 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 11] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 6] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 6];

                        //O_822 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 12] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 7] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  3] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 3] );

                        //O_813 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 13] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 9] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 9];

                        //O_804 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 14] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 9] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  5] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 5] );

                        //O_750 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 15] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 10] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  6] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 6] );

                        //O_741 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 16] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 10] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 10];

                        //O_732 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 17] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 11] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  6] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 6] );

                        //O_723 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 18] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 13] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  9] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 9] );

                        //O_714 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 19] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 14] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 14];

                        //O_705 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 20] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 14] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  9] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 9] );

                        //O_660 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 21] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 15] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  10] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 10] );

                        //O_651 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 22] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 15] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 15];

                        //O_642 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 23] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 16] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  10] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 10] );

                        //O_633 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 24] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 17] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  11] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 11] );

                        //O_624 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 25] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 19] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  14] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 14] );

                        //O_615 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 26] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 20] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 20];

                        //O_606 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 27] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 20] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  14] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 14] );

                        //O_570 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 28] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 28] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 28]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  28] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 28] );

                        //O_561 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 29] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 21] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 21];

                        //O_552 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 30] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 22] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  15] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 15] );

                        //O_543 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 31] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 23] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  16] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 16] );

                        //O_534 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 32] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 25] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  19] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 19] );

                        //O_525 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 33] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 26] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  20] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 20] );

                        //O_516 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 34] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 27] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 27];

                        //O_507 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 35] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 35] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 35]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  35] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 35] );

                        //O_480 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 36] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 36] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 36]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  36] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 36] );

                        //O_471 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 37] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 28] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 28];

                        //O_462 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 38] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 29] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  21] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 21] );

                        //O_453 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 39] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 30] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  22] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 22] );

                        //O_444 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 40] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 31] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  23] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 23] );

                        //O_435 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 41] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 33] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  26] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 26] );

                        //O_426 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 42] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 34] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  27] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 27] );

                        //O_417 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 43] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 35] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 35];

                        //O_408 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 44] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 44] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 44]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  44] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 44] );

                        //O_390 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 45] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 45] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 45]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  45] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 45] );

                        //O_381 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 46] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 36] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 36];

                        //O_372 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 47] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 37] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  28] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 28] );

                        //O_363 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 48] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 38] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  29] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 29] );

                        //O_354 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 49] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 49] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 49]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  49] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 49] );

                        //O_345 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 50] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 50] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 50]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  50] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 50] );

                        //O_336 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 51] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 42] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  34] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 34] );

                        //O_327 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 52] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 43] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  35] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 35] );

                        //O_318 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 53] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 44] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 44];

                        //O_309 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 54] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 54] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 54]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  54] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 54] );

                        //O_2100 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 55] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 55] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 55]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  55] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 55] );

                        //O_291 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 56] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 45] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 45];

                        //O_282 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 57] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 46] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  36] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 36] );

                        //O_273 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 58] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 58] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 58]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  58] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 58] );

                        //O_264 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 59] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 59] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 59]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  59] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 59] );

                        //O_255 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 60] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 60] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 60]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  60] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 60] );

                        //O_246 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 61] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 61] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 61]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  61] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 61] );

                        //O_237 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 62] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 62] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 62]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  62] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 62] );

                        //O_228 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 63] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 53] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  44] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 44] );

                        //O_219 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 64] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 54] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 54];

                        //O_2010 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 65] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 65] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 65]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  65] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 65] );

                        //O_1110 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 66] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 66] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 66];

                        //O_1101 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 67] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 55] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 55];

                        //O_192 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 68] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 68] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 68];

                        //O_183 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 69] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 69] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 69];

                        //O_174 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 70] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 70] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 70];

                        //O_165 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 71] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 71] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 71];

                        //O_156 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 72] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 72] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 72];

                        //O_147 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 73] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 73] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 73];

                        //O_138 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 74] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 74] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 74];

                        //O_129 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 75] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 75] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 75];

                        //O_1110 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 76] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 65] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 65];

                        //O_1011 : STEP: x
                        AUX_INT__o_s_s_s[m * 91 + 77] = P_PA_x * AUX_INT__n_s_s_s[m * 78 + 77] - aop_PQ_x * AUX_INT__n_s_s_s[(m+1) * 78 + 77];

                        //O_0120 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 78] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 66] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 66]
                                      + 11 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  55] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 55] );

                        //O_0111 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 79] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 66] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 66];

                        //O_0102 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 80] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 67] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  55] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 55] );

                        //O_093 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 81] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 68] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  56] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 56] );

                        //O_084 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 82] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 69] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 69]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  57] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 57] );

                        //O_075 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 83] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 70] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 70]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  58] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 58] );

                        //O_066 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 84] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 71] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 71]
                                      + 5 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  59] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 59] );

                        //O_057 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 85] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 73] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 73]
                                      + 4 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  62] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 62] );

                        //O_048 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 86] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 74] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 74]
                                      + 3 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  63] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 63] );

                        //O_039 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 87] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 75] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  64] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 64] );

                        //O_0210 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 88] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 76] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  65] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 65] );

                        //O_0111 : STEP: y
                        AUX_INT__o_s_s_s[m * 91 + 89] = P_PA_y * AUX_INT__n_s_s_s[m * 78 + 77] - aop_PQ_y * AUX_INT__n_s_s_s[(m+1) * 78 + 77];

                        //O_0012 : STEP: z
                        AUX_INT__o_s_s_s[m * 91 + 90] = P_PA_z * AUX_INT__n_s_s_s[m * 78 + 77] - aop_PQ_z * AUX_INT__n_s_s_s[(m+1) * 78 + 77]
                                      + 11 * one_over_2p * ( AUX_INT__m_s_s_s[m * 66 +  65] - a_over_p * AUX_INT__m_s_s_s[(m+1) * 66 + 65] );

                    }
}



// VRR to obtain AUX_INT__q_s_s_s
void VRR_AUX_INT__q_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__q_s_s_s,
                        double const * const restrict AUX_INT__o_s_s_s,
                        double const * const restrict AUX_INT__n_s_s_s)
{
                    // Forming AUX_INT__q_s_s_s[num_m * 105];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //Q_1300 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 0] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 0] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 0]
                                      + 12 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  0] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 0] );

                        //Q_1210 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 1] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 0] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 0];

                        //Q_1201 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 2] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 0] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 0];

                        //Q_1120 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 3] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 1] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  0] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 0] );

                        //Q_1111 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 4] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 1] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 1];

                        //Q_1102 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 5] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 2] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  0] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 0] );

                        //Q_1030 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 6] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 3] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  1] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 1] );

                        //Q_1021 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 7] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 3] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 3];

                        //Q_1012 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 8] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 5] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 5];

                        //Q_1003 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 9] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 5] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  2] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 2] );

                        //Q_940 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 10] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 6] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  3] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 3] );

                        //Q_931 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 11] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 6] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 6];

                        //Q_922 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 12] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 7] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  3] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 3] );

                        //Q_913 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 13] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 9] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 9];

                        //Q_904 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 14] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 9] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  5] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 5] );

                        //Q_850 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 15] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 10] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  6] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 6] );

                        //Q_841 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 16] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 10] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 10];

                        //Q_832 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 17] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 11] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  6] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 6] );

                        //Q_823 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 18] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 13] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  9] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 9] );

                        //Q_814 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 19] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 14] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 14];

                        //Q_805 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 20] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 14] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  9] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 9] );

                        //Q_760 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 21] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 15] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  10] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 10] );

                        //Q_751 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 22] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 15] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 15];

                        //Q_742 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 23] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 16] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  10] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 10] );

                        //Q_733 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 24] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 17] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  11] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 11] );

                        //Q_724 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 25] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 19] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  14] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 14] );

                        //Q_715 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 26] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 20] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 20];

                        //Q_706 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 27] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 20] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  14] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 14] );

                        //Q_670 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 28] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 28] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 28]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  28] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 28] );

                        //Q_661 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 29] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 21] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 21];

                        //Q_652 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 30] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 22] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  15] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 15] );

                        //Q_643 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 31] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 23] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  16] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 16] );

                        //Q_634 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 32] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 25] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  19] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 19] );

                        //Q_625 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 33] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 26] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  20] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 20] );

                        //Q_616 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 34] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 27] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 27];

                        //Q_607 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 35] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 35] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 35]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  35] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 35] );

                        //Q_580 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 36] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 36] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 36]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  36] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 36] );

                        //Q_571 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 37] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 28] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 28];

                        //Q_562 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 38] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 29] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  21] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 21] );

                        //Q_553 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 39] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 30] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  22] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 22] );

                        //Q_544 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 40] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 31] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  23] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 23] );

                        //Q_535 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 41] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 33] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  26] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 26] );

                        //Q_526 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 42] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 34] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  27] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 27] );

                        //Q_517 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 43] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 35] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 35];

                        //Q_508 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 44] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 44] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 44]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  44] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 44] );

                        //Q_490 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 45] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 45] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 45]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  45] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 45] );

                        //Q_481 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 46] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 36] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 36];

                        //Q_472 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 47] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 37] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  28] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 28] );

                        //Q_463 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 48] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 38] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  29] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 29] );

                        //Q_454 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 49] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 39] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  30] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 30] );

                        //Q_445 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 50] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 41] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  33] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 33] );

                        //Q_436 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 51] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 42] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  34] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 34] );

                        //Q_427 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 52] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 43] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  35] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 35] );

                        //Q_418 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 53] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 44] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 44];

                        //Q_409 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 54] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 54] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 54]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  54] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 54] );

                        //Q_3100 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 55] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 55] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 55]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  55] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 55] );

                        //Q_391 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 56] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 45] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 45];

                        //Q_382 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 57] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 46] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  36] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 36] );

                        //Q_373 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 58] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 47] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  37] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 37] );

                        //Q_364 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 59] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 59] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 59]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  59] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 59] );

                        //Q_355 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 60] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 60] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 60]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  60] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 60] );

                        //Q_346 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 61] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 61] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 61]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  61] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 61] );

                        //Q_337 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 62] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 52] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  43] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 43] );

                        //Q_328 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 63] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 53] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  44] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 44] );

                        //Q_319 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 64] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 54] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 54];

                        //Q_3010 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 65] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 65] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 65]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  65] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 65] );

                        //Q_2110 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 66] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 66] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 66]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  66] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 66] );

                        //Q_2101 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 67] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 55] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 55];

                        //Q_292 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 68] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 56] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  45] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 45] );

                        //Q_283 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 69] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 69] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 69]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  69] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 69] );

                        //Q_274 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 70] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 70] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 70]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  70] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 70] );

                        //Q_265 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 71] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 71] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 71]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  71] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 71] );

                        //Q_256 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 72] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 72] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 72]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  72] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 72] );

                        //Q_247 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 73] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 73] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 73]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  73] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 73] );

                        //Q_238 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 74] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 74] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 74]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  74] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 74] );

                        //Q_229 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 75] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 64] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  54] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 54] );

                        //Q_2110 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 76] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 65] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 65];

                        //Q_2011 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 77] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 77] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 77]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  77] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 77] );

                        //Q_1120 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 78] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 78] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 78];

                        //Q_1111 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 79] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 66] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 66];

                        //Q_1102 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 80] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 80] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 80];

                        //Q_193 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 81] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 81] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 81];

                        //Q_184 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 82] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 82] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 82];

                        //Q_175 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 83] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 83] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 83];

                        //Q_166 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 84] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 84] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 84];

                        //Q_157 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 85] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 85] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 85];

                        //Q_148 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 86] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 86] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 86];

                        //Q_139 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 87] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 87] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 87];

                        //Q_1210 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 88] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 88] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 88];

                        //Q_1111 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 89] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 77] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 77];

                        //Q_1012 : STEP: x
                        AUX_INT__q_s_s_s[m * 105 + 90] = P_PA_x * AUX_INT__o_s_s_s[m * 91 + 90] - aop_PQ_x * AUX_INT__o_s_s_s[(m+1) * 91 + 90];

                        //Q_0130 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 91] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 78] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 78]
                                      + 12 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  66] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 66] );

                        //Q_0121 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 92] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 78] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 78];

                        //Q_0112 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 93] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 79] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  66] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 66] );

                        //Q_0103 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 94] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 80] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 80]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  67] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 67] );

                        //Q_094 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 95] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 81] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 81]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  68] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 68] );

                        //Q_085 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 96] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 82] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 82]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  69] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 69] );

                        //Q_076 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 97] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 83] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 83]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  70] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 70] );

                        //Q_067 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 98] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 85] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 85]
                                      + 5 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  73] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 73] );

                        //Q_058 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 99] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 86] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 86]
                                      + 4 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  74] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 74] );

                        //Q_049 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 100] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 87] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 87]
                                      + 3 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  75] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 75] );

                        //Q_0310 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 101] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 88] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 88]
                                      + 2 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  76] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 76] );

                        //Q_0211 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 102] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 89] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  77] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 77] );

                        //Q_0112 : STEP: y
                        AUX_INT__q_s_s_s[m * 105 + 103] = P_PA_y * AUX_INT__o_s_s_s[m * 91 + 90] - aop_PQ_y * AUX_INT__o_s_s_s[(m+1) * 91 + 90];

                        //Q_0013 : STEP: z
                        AUX_INT__q_s_s_s[m * 105 + 104] = P_PA_z * AUX_INT__o_s_s_s[m * 91 + 90] - aop_PQ_z * AUX_INT__o_s_s_s[(m+1) * 91 + 90]
                                      + 12 * one_over_2p * ( AUX_INT__n_s_s_s[m * 78 +  77] - a_over_p * AUX_INT__n_s_s_s[(m+1) * 78 + 77] );

                    }
}



// VRR to obtain AUX_INT__r_s_s_s
void VRR_AUX_INT__r_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__r_s_s_s,
                        double const * const restrict AUX_INT__q_s_s_s,
                        double const * const restrict AUX_INT__o_s_s_s)
{
                    // Forming AUX_INT__r_s_s_s[num_m * 120];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //R_1400 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 0] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 0] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 0]
                                      + 13 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  0] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 0] );

                        //R_1310 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 1] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 0] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 0];

                        //R_1301 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 2] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 0] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 0];

                        //R_1220 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 3] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 1] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  0] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 0] );

                        //R_1211 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 4] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 1] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 1];

                        //R_1202 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 5] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 2] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  0] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 0] );

                        //R_1130 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 6] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 3] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  1] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 1] );

                        //R_1121 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 7] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 3] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 3];

                        //R_1112 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 8] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 5] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 5];

                        //R_1103 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 9] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 5] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  2] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 2] );

                        //R_1040 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 10] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 6] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  3] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 3] );

                        //R_1031 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 11] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 6] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 6];

                        //R_1022 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 12] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 7] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  3] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 3] );

                        //R_1013 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 13] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 9] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 9];

                        //R_1004 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 14] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 9] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  5] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 5] );

                        //R_950 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 15] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 10] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  6] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 6] );

                        //R_941 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 16] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 10] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 10];

                        //R_932 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 17] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 11] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  6] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 6] );

                        //R_923 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 18] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 13] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  9] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 9] );

                        //R_914 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 19] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 14] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 14];

                        //R_905 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 20] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 14] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  9] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 9] );

                        //R_860 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 21] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 15] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  10] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 10] );

                        //R_851 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 22] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 15] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 15];

                        //R_842 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 23] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 16] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  10] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 10] );

                        //R_833 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 24] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 17] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  11] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 11] );

                        //R_824 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 25] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 19] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  14] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 14] );

                        //R_815 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 26] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 20] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 20];

                        //R_806 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 27] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 20] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  14] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 14] );

                        //R_770 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 28] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 21] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  15] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 15] );

                        //R_761 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 29] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 21] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 21];

                        //R_752 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 30] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 22] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  15] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 15] );

                        //R_743 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 31] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 23] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  16] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 16] );

                        //R_734 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 32] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 25] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  19] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 19] );

                        //R_725 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 33] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 26] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  20] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 20] );

                        //R_716 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 34] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 27] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 27];

                        //R_707 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 35] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 27] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  20] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 20] );

                        //R_680 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 36] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 36] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 36]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  36] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 36] );

                        //R_671 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 37] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 28] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 28];

                        //R_662 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 38] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 29] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  21] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 21] );

                        //R_653 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 39] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 30] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  22] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 22] );

                        //R_644 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 40] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 31] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  23] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 23] );

                        //R_635 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 41] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 33] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  26] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 26] );

                        //R_626 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 42] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 34] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  27] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 27] );

                        //R_617 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 43] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 35] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 35];

                        //R_608 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 44] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 44] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 44]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  44] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 44] );

                        //R_590 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 45] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 45] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 45]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  45] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 45] );

                        //R_581 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 46] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 36] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 36];

                        //R_572 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 47] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 37] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  28] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 28] );

                        //R_563 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 48] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 38] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  29] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 29] );

                        //R_554 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 49] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 39] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  30] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 30] );

                        //R_545 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 50] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 41] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  33] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 33] );

                        //R_536 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 51] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 42] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  34] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 34] );

                        //R_527 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 52] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 43] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  35] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 35] );

                        //R_518 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 53] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 44] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 44];

                        //R_509 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 54] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 54] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 54]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  54] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 54] );

                        //R_4100 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 55] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 55] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 55]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  55] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 55] );

                        //R_491 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 56] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 45] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 45];

                        //R_482 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 57] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 46] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  36] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 36] );

                        //R_473 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 58] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 47] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  37] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 37] );

                        //R_464 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 59] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 48] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  38] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 38] );

                        //R_455 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 60] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 60] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 60]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  60] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 60] );

                        //R_446 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 61] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 51] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  42] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 42] );

                        //R_437 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 62] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 52] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  43] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 43] );

                        //R_428 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 63] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 53] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  44] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 44] );

                        //R_419 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 64] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 54] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 54];

                        //R_4010 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 65] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 65] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 65]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  65] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 65] );

                        //R_3110 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 66] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 66] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 66]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  66] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 66] );

                        //R_3101 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 67] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 55] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 55];

                        //R_392 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 68] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 56] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  45] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 45] );

                        //R_383 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 69] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 57] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  46] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 46] );

                        //R_374 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 70] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 70] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 70]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  70] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 70] );

                        //R_365 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 71] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 71] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 71]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  71] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 71] );

                        //R_356 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 72] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 72] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 72]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  72] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 72] );

                        //R_347 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 73] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 73] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 73]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  73] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 73] );

                        //R_338 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 74] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 63] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  53] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 53] );

                        //R_329 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 75] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 64] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  54] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 54] );

                        //R_3110 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 76] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 65] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 65];

                        //R_3011 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 77] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 77] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 77]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  77] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 77] );

                        //R_2120 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 78] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 78] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 78]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  78] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 78] );

                        //R_2111 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 79] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 66] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 66];

                        //R_2102 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 80] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 67] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  55] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 55] );

                        //R_293 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 81] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 81] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 81]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  81] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 81] );

                        //R_284 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 82] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 82] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 82]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  82] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 82] );

                        //R_275 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 83] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 83] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 83]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  83] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 83] );

                        //R_266 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 84] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 84] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 84]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  84] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 84] );

                        //R_257 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 85] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 85] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 85]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  85] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 85] );

                        //R_248 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 86] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 86] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 86]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  86] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 86] );

                        //R_239 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 87] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 87] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 87]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  87] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 87] );

                        //R_2210 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 88] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 76] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  65] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 65] );

                        //R_2111 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 89] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 77] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 77];

                        //R_2012 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 90] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 90] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 90]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  90] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 90] );

                        //R_1130 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 91] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 91] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 91];

                        //R_1121 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 92] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 78] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 78];

                        //R_1112 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 93] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 93] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 93];

                        //R_1103 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 94] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 94] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 94];

                        //R_194 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 95] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 95] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 95];

                        //R_185 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 96] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 96] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 96];

                        //R_176 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 97] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 97] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 97];

                        //R_167 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 98] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 98] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 98];

                        //R_158 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 99] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 99] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 99];

                        //R_149 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 100] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 100] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 100];

                        //R_1310 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 101] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 101] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 101];

                        //R_1211 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 102] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 102] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 102];

                        //R_1112 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 103] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 90] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 90];

                        //R_1013 : STEP: x
                        AUX_INT__r_s_s_s[m * 120 + 104] = P_PA_x * AUX_INT__q_s_s_s[m * 105 + 104] - aop_PQ_x * AUX_INT__q_s_s_s[(m+1) * 105 + 104];

                        //R_0140 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 105] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 91] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 91]
                                      + 13 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  78] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 78] );

                        //R_0131 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 106] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 91] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 91];

                        //R_0122 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 107] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 92] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 92]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  78] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 78] );

                        //R_0113 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 108] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 93] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 93]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  79] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 79] );

                        //R_0104 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 109] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 94] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 94]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  80] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 80] );

                        //R_095 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 110] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 95] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 95]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  81] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 81] );

                        //R_086 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 111] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 96] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 96]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  82] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 82] );

                        //R_077 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 112] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 97] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 97]
                                      + 6 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  83] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 83] );

                        //R_068 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 113] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 99] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 99]
                                      + 5 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  86] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 86] );

                        //R_059 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 114] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 100] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 100]
                                      + 4 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  87] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 87] );

                        //R_0410 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 115] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 101] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 101]
                                      + 3 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  88] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 88] );

                        //R_0311 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 116] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 102] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 102]
                                      + 2 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  89] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 89] );

                        //R_0212 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 117] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 103] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 103]
                                      + 1 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  90] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 90] );

                        //R_0113 : STEP: y
                        AUX_INT__r_s_s_s[m * 120 + 118] = P_PA_y * AUX_INT__q_s_s_s[m * 105 + 104] - aop_PQ_y * AUX_INT__q_s_s_s[(m+1) * 105 + 104];

                        //R_0014 : STEP: z
                        AUX_INT__r_s_s_s[m * 120 + 119] = P_PA_z * AUX_INT__q_s_s_s[m * 105 + 104] - aop_PQ_z * AUX_INT__q_s_s_s[(m+1) * 105 + 104]
                                      + 13 * one_over_2p * ( AUX_INT__o_s_s_s[m * 91 +  90] - a_over_p * AUX_INT__o_s_s_s[(m+1) * 91 + 90] );

                    }
}



// VRR to obtain AUX_INT__t_s_s_s
void VRR_AUX_INT__t_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__t_s_s_s,
                        double const * const restrict AUX_INT__r_s_s_s,
                        double const * const restrict AUX_INT__q_s_s_s)
{
                    // Forming AUX_INT__t_s_s_s[num_m * 136];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //T_1500 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 0] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 0] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 0]
                                      + 14 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  0] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 0] );

                        //T_1410 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 1] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 0] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 0];

                        //T_1401 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 2] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 0] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 0];

                        //T_1320 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 3] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 1] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  0] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 0] );

                        //T_1311 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 4] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 1] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 1];

                        //T_1302 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 5] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 2] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  0] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 0] );

                        //T_1230 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 6] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 3] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  1] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 1] );

                        //T_1221 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 7] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 3] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 3];

                        //T_1212 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 8] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 5] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 5];

                        //T_1203 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 9] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 5] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  2] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 2] );

                        //T_1140 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 10] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 6] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  3] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 3] );

                        //T_1131 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 11] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 6] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 6];

                        //T_1122 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 12] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 7] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  3] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 3] );

                        //T_1113 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 13] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 9] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 9];

                        //T_1104 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 14] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 9] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  5] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 5] );

                        //T_1050 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 15] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 10] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  6] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 6] );

                        //T_1041 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 16] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 10] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 10];

                        //T_1032 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 17] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 11] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  6] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 6] );

                        //T_1023 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 18] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 13] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  9] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 9] );

                        //T_1014 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 19] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 14] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 14];

                        //T_1005 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 20] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 14] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  9] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 9] );

                        //T_960 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 21] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 15] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  10] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 10] );

                        //T_951 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 22] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 15] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 15];

                        //T_942 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 23] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 16] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  10] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 10] );

                        //T_933 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 24] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 17] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  11] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 11] );

                        //T_924 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 25] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 19] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  14] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 14] );

                        //T_915 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 26] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 20] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 20];

                        //T_906 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 27] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 20] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  14] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 14] );

                        //T_870 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 28] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 21] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  15] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 15] );

                        //T_861 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 29] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 21] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 21];

                        //T_852 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 30] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 22] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  15] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 15] );

                        //T_843 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 31] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 23] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  16] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 16] );

                        //T_834 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 32] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 25] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  19] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 19] );

                        //T_825 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 33] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 26] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  20] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 20] );

                        //T_816 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 34] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 27] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 27];

                        //T_807 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 35] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 27] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  20] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 20] );

                        //T_780 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 36] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 36] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 36]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  36] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 36] );

                        //T_771 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 37] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 28] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 28];

                        //T_762 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 38] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 29] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  21] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 21] );

                        //T_753 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 39] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 30] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  22] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 22] );

                        //T_744 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 40] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 31] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  23] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 23] );

                        //T_735 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 41] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 33] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  26] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 26] );

                        //T_726 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 42] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 34] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  27] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 27] );

                        //T_717 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 43] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 35] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 35];

                        //T_708 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 44] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 44] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 44]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  44] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 44] );

                        //T_690 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 45] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 45] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 45]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  45] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 45] );

                        //T_681 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 46] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 36] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 36];

                        //T_672 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 47] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 37] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  28] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 28] );

                        //T_663 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 48] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 38] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  29] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 29] );

                        //T_654 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 49] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 39] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  30] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 30] );

                        //T_645 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 50] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 41] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  33] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 33] );

                        //T_636 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 51] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 42] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  34] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 34] );

                        //T_627 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 52] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 43] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  35] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 35] );

                        //T_618 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 53] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 44] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 44];

                        //T_609 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 54] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 54] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 54]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  54] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 54] );

                        //T_5100 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 55] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 55] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 55]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  55] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 55] );

                        //T_591 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 56] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 45] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 45];

                        //T_582 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 57] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 46] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  36] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 36] );

                        //T_573 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 58] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 47] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  37] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 37] );

                        //T_564 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 59] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 48] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  38] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 38] );

                        //T_555 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 60] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 49] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  39] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 39] );

                        //T_546 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 61] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 51] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  42] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 42] );

                        //T_537 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 62] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 52] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  43] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 43] );

                        //T_528 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 63] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 53] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  44] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 44] );

                        //T_519 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 64] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 54] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 54];

                        //T_5010 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 65] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 65] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 65]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  65] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 65] );

                        //T_4110 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 66] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 66] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 66]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  66] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 66] );

                        //T_4101 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 67] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 55] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 55];

                        //T_492 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 68] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 56] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  45] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 45] );

                        //T_483 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 69] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 57] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  46] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 46] );

                        //T_474 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 70] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 58] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  47] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 47] );

                        //T_465 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 71] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 71] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 71]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  71] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 71] );

                        //T_456 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 72] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 72] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 72]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  72] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 72] );

                        //T_447 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 73] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 62] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  52] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 52] );

                        //T_438 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 74] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 63] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  53] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 53] );

                        //T_429 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 75] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 64] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  54] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 54] );

                        //T_4110 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 76] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 65] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 65];

                        //T_4011 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 77] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 77] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 77]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  77] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 77] );

                        //T_3120 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 78] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 78] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 78]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  78] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 78] );

                        //T_3111 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 79] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 66] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 66];

                        //T_3102 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 80] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 67] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  55] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 55] );

                        //T_393 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 81] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 68] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  56] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 56] );

                        //T_384 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 82] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 82] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 82]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  82] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 82] );

                        //T_375 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 83] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 83] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 83]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  83] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 83] );

                        //T_366 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 84] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 84] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 84]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  84] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 84] );

                        //T_357 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 85] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 85] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 85]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  85] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 85] );

                        //T_348 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 86] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 86] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 86]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  86] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 86] );

                        //T_339 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 87] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 75] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  64] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 64] );

                        //T_3210 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 88] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 76] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  65] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 65] );

                        //T_3111 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 89] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 77] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 77];

                        //T_3012 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 90] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 90] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 90]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  90] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 90] );

                        //T_2130 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 91] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 91] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 91]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  91] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 91] );

                        //T_2121 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 92] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 78] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 78];

                        //T_2112 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 93] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 79] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  66] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 66] );

                        //T_2103 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 94] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 94] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 94]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  94] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 94] );

                        //T_294 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 95] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 95] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 95]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  95] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 95] );

                        //T_285 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 96] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 96] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 96]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  96] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 96] );

                        //T_276 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 97] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 97] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 97]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  97] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 97] );

                        //T_267 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 98] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 98] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 98]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  98] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 98] );

                        //T_258 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 99] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 99] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 99]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  99] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 99] );

                        //T_249 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 100] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 100] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 100]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  100] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 100] );

                        //T_2310 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 101] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 101] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 101]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  101] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 101] );

                        //T_2211 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 102] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 89] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  77] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 77] );

                        //T_2112 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 103] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 90] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 90];

                        //T_2013 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 104] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 104] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 104]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  104] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 104] );

                        //T_1140 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 105] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 105] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 105];

                        //T_1131 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 106] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 91] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 91];

                        //T_1122 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 107] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 107] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 107];

                        //T_1113 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 108] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 108] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 108];

                        //T_1104 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 109] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 109] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 109];

                        //T_195 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 110] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 110] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 110];

                        //T_186 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 111] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 111] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 111];

                        //T_177 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 112] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 112] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 112];

                        //T_168 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 113] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 113] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 113];

                        //T_159 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 114] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 114] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 114];

                        //T_1410 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 115] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 115] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 115];

                        //T_1311 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 116] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 116] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 116];

                        //T_1212 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 117] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 117] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 117];

                        //T_1113 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 118] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 104] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 104];

                        //T_1014 : STEP: x
                        AUX_INT__t_s_s_s[m * 136 + 119] = P_PA_x * AUX_INT__r_s_s_s[m * 120 + 119] - aop_PQ_x * AUX_INT__r_s_s_s[(m+1) * 120 + 119];

                        //T_0150 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 120] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 105] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 105]
                                      + 14 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  91] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 91] );

                        //T_0141 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 121] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 105] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 105];

                        //T_0132 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 122] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 106] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 106]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  91] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 91] );

                        //T_0123 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 123] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 107] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 107]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  92] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 92] );

                        //T_0114 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 124] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 108] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 108]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  93] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 93] );

                        //T_0105 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 125] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 109] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 109]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  94] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 94] );

                        //T_096 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 126] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 110] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 110]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  95] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 95] );

                        //T_087 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 127] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 111] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 111]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  96] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 96] );

                        //T_078 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 128] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 113] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 113]
                                      + 6 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  99] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 99] );

                        //T_069 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 129] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 114] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 114]
                                      + 5 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  100] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 100] );

                        //T_0510 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 130] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 115] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 115]
                                      + 4 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  101] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 101] );

                        //T_0411 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 131] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 116] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 116]
                                      + 3 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  102] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 102] );

                        //T_0312 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 132] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 117] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 117]
                                      + 2 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  103] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 103] );

                        //T_0213 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 133] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 118] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 118]
                                      + 1 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  104] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 104] );

                        //T_0114 : STEP: y
                        AUX_INT__t_s_s_s[m * 136 + 134] = P_PA_y * AUX_INT__r_s_s_s[m * 120 + 119] - aop_PQ_y * AUX_INT__r_s_s_s[(m+1) * 120 + 119];

                        //T_0015 : STEP: z
                        AUX_INT__t_s_s_s[m * 136 + 135] = P_PA_z * AUX_INT__r_s_s_s[m * 120 + 119] - aop_PQ_z * AUX_INT__r_s_s_s[(m+1) * 120 + 119]
                                      + 14 * one_over_2p * ( AUX_INT__q_s_s_s[m * 105 +  104] - a_over_p * AUX_INT__q_s_s_s[(m+1) * 105 + 104] );

                    }
}



// VRR to obtain AUX_INT__u_s_s_s
void VRR_AUX_INT__u_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__u_s_s_s,
                        double const * const restrict AUX_INT__t_s_s_s,
                        double const * const restrict AUX_INT__r_s_s_s)
{
                    // Forming AUX_INT__u_s_s_s[num_m * 153];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //U_1600 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 0] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 0] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 0]
                                      + 15 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  0] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 0] );

                        //U_1510 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 1] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 0] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 0];

                        //U_1501 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 2] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 0] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 0];

                        //U_1420 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 3] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 1] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  0] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 0] );

                        //U_1411 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 4] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 1] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 1];

                        //U_1402 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 5] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 2] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  0] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 0] );

                        //U_1330 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 6] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 3] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  1] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 1] );

                        //U_1321 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 7] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 3] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 3];

                        //U_1312 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 8] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 5] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 5];

                        //U_1303 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 9] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 5] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  2] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 2] );

                        //U_1240 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 10] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 6] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  3] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 3] );

                        //U_1231 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 11] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 6] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 6];

                        //U_1222 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 12] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 7] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  3] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 3] );

                        //U_1213 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 13] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 9] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 9];

                        //U_1204 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 14] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 9] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  5] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 5] );

                        //U_1150 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 15] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 10] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  6] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 6] );

                        //U_1141 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 16] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 10] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 10];

                        //U_1132 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 17] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 11] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  6] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 6] );

                        //U_1123 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 18] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 13] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  9] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 9] );

                        //U_1114 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 19] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 14] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 14];

                        //U_1105 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 20] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 14] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  9] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 9] );

                        //U_1060 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 21] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 15] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  10] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 10] );

                        //U_1051 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 22] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 15] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 15];

                        //U_1042 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 23] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 16] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  10] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 10] );

                        //U_1033 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 24] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 17] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  11] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 11] );

                        //U_1024 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 25] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 19] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  14] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 14] );

                        //U_1015 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 26] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 20] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 20];

                        //U_1006 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 27] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 20] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  14] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 14] );

                        //U_970 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 28] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 21] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  15] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 15] );

                        //U_961 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 29] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 21] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 21];

                        //U_952 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 30] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 22] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  15] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 15] );

                        //U_943 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 31] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 23] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  16] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 16] );

                        //U_934 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 32] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 25] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  19] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 19] );

                        //U_925 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 33] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 26] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  20] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 20] );

                        //U_916 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 34] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 27] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 27];

                        //U_907 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 35] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 27] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  20] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 20] );

                        //U_880 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 36] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 28] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  21] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 21] );

                        //U_871 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 37] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 28] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 28];

                        //U_862 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 38] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 29] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  21] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 21] );

                        //U_853 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 39] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 30] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  22] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 22] );

                        //U_844 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 40] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 31] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  23] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 23] );

                        //U_835 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 41] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 33] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  26] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 26] );

                        //U_826 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 42] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 34] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  27] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 27] );

                        //U_817 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 43] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 35] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 35];

                        //U_808 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 44] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 35] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  27] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 27] );

                        //U_790 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 45] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 45] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 45]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  45] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 45] );

                        //U_781 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 46] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 36] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 36];

                        //U_772 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 47] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 37] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  28] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 28] );

                        //U_763 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 48] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 38] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  29] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 29] );

                        //U_754 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 49] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 39] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  30] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 30] );

                        //U_745 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 50] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 41] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  33] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 33] );

                        //U_736 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 51] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 42] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  34] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 34] );

                        //U_727 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 52] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 43] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  35] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 35] );

                        //U_718 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 53] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 44] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 44];

                        //U_709 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 54] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 54] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 54]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  54] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 54] );

                        //U_6100 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 55] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 55] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 55]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  55] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 55] );

                        //U_691 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 56] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 45] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 45];

                        //U_682 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 57] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 46] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  36] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 36] );

                        //U_673 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 58] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 47] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  37] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 37] );

                        //U_664 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 59] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 48] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  38] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 38] );

                        //U_655 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 60] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 49] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  39] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 39] );

                        //U_646 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 61] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 51] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  42] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 42] );

                        //U_637 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 62] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 52] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  43] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 43] );

                        //U_628 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 63] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 53] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  44] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 44] );

                        //U_619 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 64] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 54] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 54];

                        //U_6010 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 65] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 65] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 65]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  65] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 65] );

                        //U_5110 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 66] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 66] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 66]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  66] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 66] );

                        //U_5101 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 67] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 55] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 55];

                        //U_592 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 68] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 56] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  45] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 45] );

                        //U_583 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 69] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 57] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  46] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 46] );

                        //U_574 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 70] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 58] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  47] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 47] );

                        //U_565 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 71] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 59] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 59]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  48] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 48] );

                        //U_556 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 72] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 61] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 61]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  51] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 51] );

                        //U_547 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 73] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 62] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  52] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 52] );

                        //U_538 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 74] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 63] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  53] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 53] );

                        //U_529 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 75] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 64] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  54] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 54] );

                        //U_5110 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 76] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 65] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 65];

                        //U_5011 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 77] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 77] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 77]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  77] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 77] );

                        //U_4120 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 78] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 78] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 78]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  78] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 78] );

                        //U_4111 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 79] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 66] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 66];

                        //U_4102 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 80] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 67] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  55] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 55] );

                        //U_493 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 81] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 68] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  56] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 56] );

                        //U_484 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 82] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 69] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 69]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  57] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 57] );

                        //U_475 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 83] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 83] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 83]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  83] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 83] );

                        //U_466 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 84] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 84] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 84]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  84] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 84] );

                        //U_457 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 85] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 85] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 85]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  85] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 85] );

                        //U_448 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 86] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 74] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 74]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  63] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 63] );

                        //U_439 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 87] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 75] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  64] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 64] );

                        //U_4210 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 88] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 76] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  65] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 65] );

                        //U_4111 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 89] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 77] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 77];

                        //U_4012 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 90] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 90] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 90]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  90] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 90] );

                        //U_3130 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 91] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 91] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 91]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  91] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 91] );

                        //U_3121 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 92] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 78] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 78];

                        //U_3112 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 93] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 79] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  66] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 66] );

                        //U_3103 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 94] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 80] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 80]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  67] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 67] );

                        //U_394 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 95] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 95] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 95]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  95] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 95] );

                        //U_385 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 96] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 96] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 96]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  96] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 96] );

                        //U_376 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 97] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 97] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 97]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  97] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 97] );

                        //U_367 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 98] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 98] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 98]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  98] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 98] );

                        //U_358 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 99] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 99] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 99]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  99] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 99] );

                        //U_349 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 100] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 100] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 100]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  100] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 100] );

                        //U_3310 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 101] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 88] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 88]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  76] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 76] );

                        //U_3211 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 102] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 89] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  77] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 77] );

                        //U_3112 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 103] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 90] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 90];

                        //U_3013 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 104] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 104] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 104]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  104] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 104] );

                        //U_2140 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 105] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 105] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 105]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  105] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 105] );

                        //U_2131 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 106] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 91] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 91];

                        //U_2122 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 107] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 92] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 92]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  78] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 78] );

                        //U_2113 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 108] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 108] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 108]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  108] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 108] );

                        //U_2104 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 109] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 109] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 109]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  109] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 109] );

                        //U_295 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 110] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 110] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 110]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  110] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 110] );

                        //U_286 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 111] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 111] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 111]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  111] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 111] );

                        //U_277 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 112] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 112] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 112]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  112] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 112] );

                        //U_268 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 113] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 113] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 113]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  113] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 113] );

                        //U_259 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 114] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 114] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 114]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  114] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 114] );

                        //U_2410 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 115] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 115] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 115]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  115] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 115] );

                        //U_2311 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 116] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 116] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 116]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  116] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 116] );

                        //U_2212 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 117] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 103] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 103]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  90] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 90] );

                        //U_2113 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 118] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 104] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 104];

                        //U_2014 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 119] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 119] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 119]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  119] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 119] );

                        //U_1150 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 120] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 120] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 120];

                        //U_1141 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 121] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 105] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 105];

                        //U_1132 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 122] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 122] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 122];

                        //U_1123 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 123] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 123] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 123];

                        //U_1114 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 124] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 124] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 124];

                        //U_1105 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 125] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 125] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 125];

                        //U_196 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 126] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 126] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 126];

                        //U_187 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 127] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 127] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 127];

                        //U_178 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 128] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 128] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 128];

                        //U_169 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 129] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 129] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 129];

                        //U_1510 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 130] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 130] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 130];

                        //U_1411 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 131] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 131] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 131];

                        //U_1312 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 132] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 132] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 132];

                        //U_1213 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 133] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 133] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 133];

                        //U_1114 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 134] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 119] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 119];

                        //U_1015 : STEP: x
                        AUX_INT__u_s_s_s[m * 153 + 135] = P_PA_x * AUX_INT__t_s_s_s[m * 136 + 135] - aop_PQ_x * AUX_INT__t_s_s_s[(m+1) * 136 + 135];

                        //U_0160 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 136] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 120] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 120]
                                      + 15 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  105] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 105] );

                        //U_0151 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 137] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 120] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 120];

                        //U_0142 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 138] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 121] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 121]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  105] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 105] );

                        //U_0133 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 139] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 122] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 122]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  106] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 106] );

                        //U_0124 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 140] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 123] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 123]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  107] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 107] );

                        //U_0115 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 141] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 124] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 124]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  108] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 108] );

                        //U_0106 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 142] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 125] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 125]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  109] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 109] );

                        //U_097 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 143] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 126] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 126]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  110] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 110] );

                        //U_088 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 144] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 127] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 127]
                                      + 7 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  111] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 111] );

                        //U_079 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 145] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 129] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 129]
                                      + 6 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  114] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 114] );

                        //U_0610 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 146] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 130] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 130]
                                      + 5 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  115] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 115] );

                        //U_0511 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 147] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 131] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 131]
                                      + 4 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  116] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 116] );

                        //U_0412 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 148] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 132] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 132]
                                      + 3 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  117] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 117] );

                        //U_0313 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 149] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 133] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 133]
                                      + 2 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  118] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 118] );

                        //U_0214 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 150] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 134] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 134]
                                      + 1 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  119] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 119] );

                        //U_0115 : STEP: y
                        AUX_INT__u_s_s_s[m * 153 + 151] = P_PA_y * AUX_INT__t_s_s_s[m * 136 + 135] - aop_PQ_y * AUX_INT__t_s_s_s[(m+1) * 136 + 135];

                        //U_0016 : STEP: z
                        AUX_INT__u_s_s_s[m * 153 + 152] = P_PA_z * AUX_INT__t_s_s_s[m * 136 + 135] - aop_PQ_z * AUX_INT__t_s_s_s[(m+1) * 136 + 135]
                                      + 15 * one_over_2p * ( AUX_INT__r_s_s_s[m * 120 +  119] - a_over_p * AUX_INT__r_s_s_s[(m+1) * 120 + 119] );

                    }
}



// VRR to obtain AUX_INT__v_s_s_s
void VRR_AUX_INT__v_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__v_s_s_s,
                        double const * const restrict AUX_INT__u_s_s_s,
                        double const * const restrict AUX_INT__t_s_s_s)
{
                    // Forming AUX_INT__v_s_s_s[num_m * 171];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //V_1700 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 0] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 0] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 0]
                                      + 16 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  0] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 0] );

                        //V_1610 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 1] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 0] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 0];

                        //V_1601 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 2] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 0] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 0];

                        //V_1520 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 3] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 1] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  0] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 0] );

                        //V_1511 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 4] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 1] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 1];

                        //V_1502 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 5] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 2] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  0] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 0] );

                        //V_1430 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 6] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 3] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  1] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 1] );

                        //V_1421 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 7] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 3] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 3];

                        //V_1412 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 8] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 5] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 5];

                        //V_1403 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 9] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 5] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  2] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 2] );

                        //V_1340 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 10] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 6] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  3] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 3] );

                        //V_1331 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 11] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 6] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 6];

                        //V_1322 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 12] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 7] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  3] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 3] );

                        //V_1313 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 13] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 9] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 9];

                        //V_1304 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 14] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 9] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  5] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 5] );

                        //V_1250 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 15] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 10] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  6] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 6] );

                        //V_1241 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 16] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 10] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 10];

                        //V_1232 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 17] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 11] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  6] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 6] );

                        //V_1223 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 18] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 13] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  9] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 9] );

                        //V_1214 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 19] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 14] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 14];

                        //V_1205 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 20] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 14] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  9] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 9] );

                        //V_1160 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 21] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 15] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  10] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 10] );

                        //V_1151 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 22] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 15] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 15];

                        //V_1142 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 23] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 16] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  10] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 10] );

                        //V_1133 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 24] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 17] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  11] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 11] );

                        //V_1124 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 25] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 19] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  14] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 14] );

                        //V_1115 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 26] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 20] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 20];

                        //V_1106 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 27] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 20] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  14] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 14] );

                        //V_1070 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 28] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 21] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  15] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 15] );

                        //V_1061 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 29] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 21] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 21];

                        //V_1052 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 30] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 22] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  15] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 15] );

                        //V_1043 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 31] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 23] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  16] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 16] );

                        //V_1034 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 32] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 25] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  19] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 19] );

                        //V_1025 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 33] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 26] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  20] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 20] );

                        //V_1016 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 34] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 27] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 27];

                        //V_1007 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 35] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 27] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  20] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 20] );

                        //V_980 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 36] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 28] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  21] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 21] );

                        //V_971 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 37] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 28] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 28];

                        //V_962 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 38] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 29] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  21] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 21] );

                        //V_953 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 39] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 30] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  22] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 22] );

                        //V_944 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 40] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 31] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  23] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 23] );

                        //V_935 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 41] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 33] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  26] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 26] );

                        //V_926 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 42] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 34] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  27] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 27] );

                        //V_917 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 43] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 35] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 35];

                        //V_908 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 44] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 35] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  27] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 27] );

                        //V_890 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 45] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 45] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 45]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  45] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 45] );

                        //V_881 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 46] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 36] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 36];

                        //V_872 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 47] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 37] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  28] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 28] );

                        //V_863 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 48] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 38] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  29] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 29] );

                        //V_854 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 49] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 39] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  30] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 30] );

                        //V_845 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 50] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 41] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  33] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 33] );

                        //V_836 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 51] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 42] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  34] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 34] );

                        //V_827 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 52] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 43] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  35] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 35] );

                        //V_818 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 53] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 44] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 44];

                        //V_809 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 54] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 54] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 54]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  54] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 54] );

                        //V_7100 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 55] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 55] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 55]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  55] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 55] );

                        //V_791 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 56] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 45] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 45];

                        //V_782 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 57] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 46] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  36] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 36] );

                        //V_773 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 58] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 47] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  37] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 37] );

                        //V_764 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 59] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 48] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  38] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 38] );

                        //V_755 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 60] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 49] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  39] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 39] );

                        //V_746 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 61] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 51] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  42] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 42] );

                        //V_737 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 62] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 52] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  43] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 43] );

                        //V_728 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 63] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 53] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  44] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 44] );

                        //V_719 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 64] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 54] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 54];

                        //V_7010 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 65] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 65] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 65]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  65] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 65] );

                        //V_6110 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 66] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 66] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 66]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  66] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 66] );

                        //V_6101 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 67] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 55] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 55];

                        //V_692 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 68] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 56] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  45] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 45] );

                        //V_683 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 69] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 57] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  46] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 46] );

                        //V_674 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 70] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 58] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  47] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 47] );

                        //V_665 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 71] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 59] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 59]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  48] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 48] );

                        //V_656 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 72] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 61] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 61]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  51] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 51] );

                        //V_647 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 73] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 62] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  52] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 52] );

                        //V_638 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 74] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 63] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  53] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 53] );

                        //V_629 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 75] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 64] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  54] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 54] );

                        //V_6110 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 76] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 65] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 65];

                        //V_6011 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 77] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 77] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 77]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  77] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 77] );

                        //V_5120 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 78] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 78] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 78]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  78] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 78] );

                        //V_5111 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 79] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 66] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 66];

                        //V_5102 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 80] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 67] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  55] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 55] );

                        //V_593 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 81] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 68] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  56] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 56] );

                        //V_584 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 82] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 69] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 69]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  57] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 57] );

                        //V_575 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 83] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 70] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 70]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  58] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 58] );

                        //V_566 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 84] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 84] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 84]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  84] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 84] );

                        //V_557 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 85] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 73] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 73]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  62] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 62] );

                        //V_548 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 86] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 74] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 74]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  63] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 63] );

                        //V_539 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 87] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 75] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  64] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 64] );

                        //V_5210 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 88] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 76] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  65] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 65] );

                        //V_5111 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 89] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 77] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 77];

                        //V_5012 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 90] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 90] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 90]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  90] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 90] );

                        //V_4130 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 91] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 91] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 91]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  91] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 91] );

                        //V_4121 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 92] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 78] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 78];

                        //V_4112 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 93] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 79] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  66] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 66] );

                        //V_4103 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 94] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 80] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 80]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  67] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 67] );

                        //V_494 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 95] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 81] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 81]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  68] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 68] );

                        //V_485 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 96] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 96] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 96]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  96] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 96] );

                        //V_476 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 97] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 97] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 97]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  97] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 97] );

                        //V_467 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 98] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 98] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 98]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  98] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 98] );

                        //V_458 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 99] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 99] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 99]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  99] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 99] );

                        //V_449 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 100] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 87] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 87]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  75] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 75] );

                        //V_4310 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 101] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 88] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 88]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  76] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 76] );

                        //V_4211 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 102] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 89] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  77] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 77] );

                        //V_4112 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 103] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 90] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 90];

                        //V_4013 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 104] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 104] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 104]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  104] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 104] );

                        //V_3140 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 105] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 105] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 105]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  105] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 105] );

                        //V_3131 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 106] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 91] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 91];

                        //V_3122 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 107] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 92] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 92]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  78] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 78] );

                        //V_3113 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 108] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 93] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 93]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  79] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 79] );

                        //V_3104 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 109] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 109] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 109]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  109] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 109] );

                        //V_395 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 110] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 110] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 110]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  110] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 110] );

                        //V_386 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 111] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 111] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 111]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  111] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 111] );

                        //V_377 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 112] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 112] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 112]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  112] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 112] );

                        //V_368 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 113] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 113] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 113]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  113] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 113] );

                        //V_359 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 114] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 114] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 114]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  114] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 114] );

                        //V_3410 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 115] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 115] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 115]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  115] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 115] );

                        //V_3311 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 116] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 102] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 102]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  89] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 89] );

                        //V_3212 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 117] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 103] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 103]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  90] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 90] );

                        //V_3113 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 118] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 104] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 104];

                        //V_3014 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 119] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 119] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 119]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  119] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 119] );

                        //V_2150 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 120] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 120] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 120]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  120] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 120] );

                        //V_2141 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 121] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 105] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 105];

                        //V_2132 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 122] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 106] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 106]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  91] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 91] );

                        //V_2123 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 123] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 123] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 123]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  123] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 123] );

                        //V_2114 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 124] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 124] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 124]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  124] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 124] );

                        //V_2105 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 125] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 125] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 125]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  125] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 125] );

                        //V_296 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 126] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 126] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 126]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  126] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 126] );

                        //V_287 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 127] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 127] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 127]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  127] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 127] );

                        //V_278 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 128] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 128] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 128]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  128] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 128] );

                        //V_269 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 129] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 129] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 129]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  129] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 129] );

                        //V_2510 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 130] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 130] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 130]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  130] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 130] );

                        //V_2411 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 131] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 131] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 131]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  131] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 131] );

                        //V_2312 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 132] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 132] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 132]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  132] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 132] );

                        //V_2213 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 133] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 118] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 118]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  104] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 104] );

                        //V_2114 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 134] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 119] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 119];

                        //V_2015 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 135] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 135] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 135]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  135] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 135] );

                        //V_1160 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 136] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 136] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 136];

                        //V_1151 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 137] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 120] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 120];

                        //V_1142 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 138] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 138] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 138];

                        //V_1133 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 139] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 139] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 139];

                        //V_1124 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 140] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 140] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 140];

                        //V_1115 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 141] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 141] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 141];

                        //V_1106 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 142] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 142] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 142];

                        //V_197 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 143] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 143] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 143];

                        //V_188 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 144] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 144] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 144];

                        //V_179 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 145] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 145] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 145];

                        //V_1610 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 146] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 146] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 146];

                        //V_1511 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 147] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 147] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 147];

                        //V_1412 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 148] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 148] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 148];

                        //V_1313 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 149] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 149] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 149];

                        //V_1214 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 150] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 150] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 150];

                        //V_1115 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 151] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 135] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 135];

                        //V_1016 : STEP: x
                        AUX_INT__v_s_s_s[m * 171 + 152] = P_PA_x * AUX_INT__u_s_s_s[m * 153 + 152] - aop_PQ_x * AUX_INT__u_s_s_s[(m+1) * 153 + 152];

                        //V_0170 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 153] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 136] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 136]
                                      + 16 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  120] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 120] );

                        //V_0161 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 154] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 136] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 136];

                        //V_0152 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 155] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 137] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 137]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  120] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 120] );

                        //V_0143 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 156] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 138] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 138]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  121] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 121] );

                        //V_0134 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 157] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 139] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 139]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  122] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 122] );

                        //V_0125 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 158] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 140] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 140]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  123] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 123] );

                        //V_0116 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 159] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 141] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 141]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  124] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 124] );

                        //V_0107 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 160] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 142] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 142]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  125] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 125] );

                        //V_098 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 161] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 143] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 143]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  126] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 126] );

                        //V_089 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 162] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 145] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 145]
                                      + 7 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  129] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 129] );

                        //V_0710 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 163] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 146] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 146]
                                      + 6 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  130] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 130] );

                        //V_0611 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 164] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 147] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 147]
                                      + 5 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  131] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 131] );

                        //V_0512 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 165] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 148] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 148]
                                      + 4 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  132] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 132] );

                        //V_0413 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 166] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 149] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 149]
                                      + 3 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  133] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 133] );

                        //V_0314 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 167] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 150] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 150]
                                      + 2 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  134] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 134] );

                        //V_0215 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 168] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 151] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 151]
                                      + 1 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  135] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 135] );

                        //V_0116 : STEP: y
                        AUX_INT__v_s_s_s[m * 171 + 169] = P_PA_y * AUX_INT__u_s_s_s[m * 153 + 152] - aop_PQ_y * AUX_INT__u_s_s_s[(m+1) * 153 + 152];

                        //V_0017 : STEP: z
                        AUX_INT__v_s_s_s[m * 171 + 170] = P_PA_z * AUX_INT__u_s_s_s[m * 153 + 152] - aop_PQ_z * AUX_INT__u_s_s_s[(m+1) * 153 + 152]
                                      + 16 * one_over_2p * ( AUX_INT__t_s_s_s[m * 136 +  135] - a_over_p * AUX_INT__t_s_s_s[(m+1) * 136 + 135] );

                    }
}



// VRR to obtain AUX_INT__w_s_s_s
void VRR_AUX_INT__w_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__w_s_s_s,
                        double const * const restrict AUX_INT__v_s_s_s,
                        double const * const restrict AUX_INT__u_s_s_s)
{
                    // Forming AUX_INT__w_s_s_s[num_m * 190];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //W_1800 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 0] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 0] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 0]
                                      + 17 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  0] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 0] );

                        //W_1710 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 1] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 0] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 0];

                        //W_1701 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 2] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 0] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 0];

                        //W_1620 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 3] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 1] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  0] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 0] );

                        //W_1611 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 4] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 1] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 1];

                        //W_1602 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 5] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 2] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  0] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 0] );

                        //W_1530 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 6] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 3] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  1] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 1] );

                        //W_1521 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 7] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 3] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 3];

                        //W_1512 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 8] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 5] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 5];

                        //W_1503 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 9] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 5] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  2] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 2] );

                        //W_1440 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 10] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 6] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  3] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 3] );

                        //W_1431 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 11] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 6] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 6];

                        //W_1422 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 12] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 7] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  3] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 3] );

                        //W_1413 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 13] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 9] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 9];

                        //W_1404 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 14] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 9] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  5] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 5] );

                        //W_1350 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 15] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 10] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  6] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 6] );

                        //W_1341 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 16] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 10] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 10];

                        //W_1332 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 17] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 11] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  6] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 6] );

                        //W_1323 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 18] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 13] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  9] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 9] );

                        //W_1314 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 19] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 14] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 14];

                        //W_1305 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 20] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 14] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  9] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 9] );

                        //W_1260 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 21] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 15] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  10] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 10] );

                        //W_1251 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 22] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 15] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 15];

                        //W_1242 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 23] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 16] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  10] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 10] );

                        //W_1233 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 24] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 17] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  11] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 11] );

                        //W_1224 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 25] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 19] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  14] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 14] );

                        //W_1215 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 26] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 20] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 20];

                        //W_1206 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 27] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 20] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  14] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 14] );

                        //W_1170 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 28] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 21] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  15] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 15] );

                        //W_1161 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 29] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 21] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 21];

                        //W_1152 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 30] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 22] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  15] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 15] );

                        //W_1143 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 31] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 23] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  16] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 16] );

                        //W_1134 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 32] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 25] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  19] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 19] );

                        //W_1125 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 33] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 26] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  20] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 20] );

                        //W_1116 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 34] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 27] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 27];

                        //W_1107 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 35] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 27] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  20] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 20] );

                        //W_1080 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 36] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 28] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  21] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 21] );

                        //W_1071 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 37] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 28] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 28];

                        //W_1062 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 38] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 29] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  21] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 21] );

                        //W_1053 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 39] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 30] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  22] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 22] );

                        //W_1044 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 40] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 31] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  23] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 23] );

                        //W_1035 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 41] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 33] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  26] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 26] );

                        //W_1026 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 42] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 34] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  27] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 27] );

                        //W_1017 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 43] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 35] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 35];

                        //W_1008 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 44] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 35] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  27] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 27] );

                        //W_990 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 45] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 36] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 36]
                                      + 8 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  28] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 28] );

                        //W_981 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 46] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 36] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 36];

                        //W_972 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 47] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 37] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  28] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 28] );

                        //W_963 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 48] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 38] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  29] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 29] );

                        //W_954 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 49] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 39] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  30] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 30] );

                        //W_945 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 50] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 41] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  33] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 33] );

                        //W_936 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 51] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 42] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  34] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 34] );

                        //W_927 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 52] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 43] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  35] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 35] );

                        //W_918 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 53] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 44] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 44];

                        //W_909 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 54] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 44] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 44]
                                      + 8 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  35] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 35] );

                        //W_8100 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 55] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 55] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 55]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  55] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 55] );

                        //W_891 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 56] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 45] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 45];

                        //W_882 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 57] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 46] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  36] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 36] );

                        //W_873 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 58] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 47] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  37] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 37] );

                        //W_864 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 59] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 48] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  38] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 38] );

                        //W_855 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 60] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 49] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  39] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 39] );

                        //W_846 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 61] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 51] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  42] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 42] );

                        //W_837 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 62] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 52] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  43] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 43] );

                        //W_828 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 63] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 53] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  44] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 44] );

                        //W_819 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 64] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 54] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 54];

                        //W_8010 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 65] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 65] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 65]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  65] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 65] );

                        //W_7110 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 66] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 66] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 66]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  66] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 66] );

                        //W_7101 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 67] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 55] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 55];

                        //W_792 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 68] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 56] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  45] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 45] );

                        //W_783 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 69] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 57] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  46] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 46] );

                        //W_774 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 70] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 58] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  47] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 47] );

                        //W_765 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 71] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 59] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 59]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  48] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 48] );

                        //W_756 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 72] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 61] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 61]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  51] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 51] );

                        //W_747 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 73] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 62] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  52] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 52] );

                        //W_738 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 74] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 63] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  53] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 53] );

                        //W_729 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 75] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 64] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  54] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 54] );

                        //W_7110 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 76] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 65] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 65];

                        //W_7011 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 77] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 77] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 77]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  77] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 77] );

                        //W_6120 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 78] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 78] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 78]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  78] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 78] );

                        //W_6111 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 79] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 66] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 66];

                        //W_6102 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 80] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 67] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  55] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 55] );

                        //W_693 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 81] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 68] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  56] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 56] );

                        //W_684 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 82] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 69] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 69]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  57] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 57] );

                        //W_675 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 83] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 70] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 70]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  58] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 58] );

                        //W_666 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 84] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 71] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 71]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  59] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 59] );

                        //W_657 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 85] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 73] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 73]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  62] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 62] );

                        //W_648 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 86] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 74] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 74]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  63] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 63] );

                        //W_639 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 87] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 75] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  64] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 64] );

                        //W_6210 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 88] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 76] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  65] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 65] );

                        //W_6111 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 89] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 77] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 77];

                        //W_6012 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 90] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 90] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 90]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  90] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 90] );

                        //W_5130 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 91] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 91] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 91]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  91] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 91] );

                        //W_5121 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 92] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 78] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 78];

                        //W_5112 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 93] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 79] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  66] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 66] );

                        //W_5103 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 94] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 80] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 80]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  67] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 67] );

                        //W_594 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 95] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 81] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 81]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  68] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 68] );

                        //W_585 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 96] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 82] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 82]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  69] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 69] );

                        //W_576 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 97] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 97] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 97]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  97] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 97] );

                        //W_567 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 98] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 98] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 98]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  98] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 98] );

                        //W_558 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 99] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 86] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 86]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  74] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 74] );

                        //W_549 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 100] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 87] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 87]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  75] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 75] );

                        //W_5310 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 101] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 88] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 88]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  76] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 76] );

                        //W_5211 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 102] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 89] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  77] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 77] );

                        //W_5112 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 103] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 90] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 90];

                        //W_5013 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 104] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 104] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 104]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  104] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 104] );

                        //W_4140 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 105] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 105] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 105]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  105] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 105] );

                        //W_4131 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 106] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 91] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 91];

                        //W_4122 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 107] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 92] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 92]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  78] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 78] );

                        //W_4113 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 108] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 93] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 93]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  79] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 79] );

                        //W_4104 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 109] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 94] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 94]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  80] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 80] );

                        //W_495 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 110] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 110] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 110]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  110] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 110] );

                        //W_486 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 111] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 111] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 111]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  111] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 111] );

                        //W_477 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 112] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 112] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 112]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  112] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 112] );

                        //W_468 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 113] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 113] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 113]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  113] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 113] );

                        //W_459 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 114] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 114] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 114]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  114] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 114] );

                        //W_4410 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 115] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 101] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 101]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  88] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 88] );

                        //W_4311 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 116] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 102] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 102]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  89] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 89] );

                        //W_4212 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 117] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 103] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 103]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  90] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 90] );

                        //W_4113 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 118] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 104] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 104];

                        //W_4014 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 119] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 119] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 119]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  119] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 119] );

                        //W_3150 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 120] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 120] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 120]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  120] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 120] );

                        //W_3141 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 121] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 105] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 105];

                        //W_3132 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 122] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 106] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 106]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  91] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 91] );

                        //W_3123 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 123] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 107] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 107]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  92] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 92] );

                        //W_3114 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 124] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 124] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 124]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  124] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 124] );

                        //W_3105 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 125] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 125] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 125]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  125] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 125] );

                        //W_396 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 126] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 126] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 126]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  126] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 126] );

                        //W_387 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 127] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 127] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 127]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  127] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 127] );

                        //W_378 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 128] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 128] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 128]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  128] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 128] );

                        //W_369 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 129] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 129] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 129]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  129] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 129] );

                        //W_3510 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 130] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 130] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 130]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  130] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 130] );

                        //W_3411 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 131] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 131] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 131]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  131] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 131] );

                        //W_3312 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 132] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 117] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 117]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  103] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 103] );

                        //W_3213 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 133] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 118] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 118]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  104] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 104] );

                        //W_3114 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 134] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 119] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 119];

                        //W_3015 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 135] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 135] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 135]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  135] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 135] );

                        //W_2160 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 136] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 136] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 136]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  136] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 136] );

                        //W_2151 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 137] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 120] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 120];

                        //W_2142 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 138] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 121] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 121]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  105] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 105] );

                        //W_2133 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 139] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 139] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 139]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  139] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 139] );

                        //W_2124 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 140] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 140] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 140]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  140] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 140] );

                        //W_2115 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 141] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 141] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 141]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  141] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 141] );

                        //W_2106 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 142] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 142] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 142]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  142] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 142] );

                        //W_297 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 143] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 143] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 143]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  143] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 143] );

                        //W_288 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 144] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 144] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 144]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  144] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 144] );

                        //W_279 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 145] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 145] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 145]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  145] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 145] );

                        //W_2610 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 146] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 146] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 146]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  146] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 146] );

                        //W_2511 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 147] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 147] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 147]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  147] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 147] );

                        //W_2412 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 148] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 148] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 148]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  148] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 148] );

                        //W_2313 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 149] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 149] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 149]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  149] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 149] );

                        //W_2214 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 150] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 134] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 134]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  119] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 119] );

                        //W_2115 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 151] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 135] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 135];

                        //W_2016 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 152] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 152] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 152]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  152] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 152] );

                        //W_1170 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 153] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 153] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 153];

                        //W_1161 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 154] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 136] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 136];

                        //W_1152 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 155] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 155] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 155];

                        //W_1143 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 156] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 156] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 156];

                        //W_1134 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 157] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 157] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 157];

                        //W_1125 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 158] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 158] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 158];

                        //W_1116 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 159] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 159] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 159];

                        //W_1107 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 160] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 160] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 160];

                        //W_198 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 161] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 161] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 161];

                        //W_189 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 162] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 162] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 162];

                        //W_1710 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 163] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 163] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 163];

                        //W_1611 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 164] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 164] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 164];

                        //W_1512 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 165] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 165] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 165];

                        //W_1413 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 166] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 166] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 166];

                        //W_1314 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 167] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 167] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 167];

                        //W_1215 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 168] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 168] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 168];

                        //W_1116 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 169] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 152] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 152];

                        //W_1017 : STEP: x
                        AUX_INT__w_s_s_s[m * 190 + 170] = P_PA_x * AUX_INT__v_s_s_s[m * 171 + 170] - aop_PQ_x * AUX_INT__v_s_s_s[(m+1) * 171 + 170];

                        //W_0180 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 171] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 153] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 153]
                                      + 17 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  136] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 136] );

                        //W_0171 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 172] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 153] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 153];

                        //W_0162 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 173] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 154] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 154]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  136] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 136] );

                        //W_0153 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 174] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 155] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 155]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  137] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 137] );

                        //W_0144 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 175] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 156] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 156]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  138] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 138] );

                        //W_0135 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 176] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 157] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 157]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  139] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 139] );

                        //W_0126 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 177] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 158] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 158]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  140] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 140] );

                        //W_0117 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 178] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 159] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 159]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  141] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 141] );

                        //W_0108 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 179] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 160] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 160]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  142] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 142] );

                        //W_099 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 180] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 161] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 161]
                                      + 8 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  143] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 143] );

                        //W_0810 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 181] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 163] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 163]
                                      + 7 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  146] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 146] );

                        //W_0711 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 182] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 164] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 164]
                                      + 6 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  147] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 147] );

                        //W_0612 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 183] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 165] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 165]
                                      + 5 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  148] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 148] );

                        //W_0513 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 184] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 166] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 166]
                                      + 4 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  149] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 149] );

                        //W_0414 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 185] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 167] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 167]
                                      + 3 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  150] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 150] );

                        //W_0315 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 186] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 168] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 168]
                                      + 2 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  151] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 151] );

                        //W_0216 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 187] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 169] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 169]
                                      + 1 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  152] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 152] );

                        //W_0117 : STEP: y
                        AUX_INT__w_s_s_s[m * 190 + 188] = P_PA_y * AUX_INT__v_s_s_s[m * 171 + 170] - aop_PQ_y * AUX_INT__v_s_s_s[(m+1) * 171 + 170];

                        //W_0018 : STEP: z
                        AUX_INT__w_s_s_s[m * 190 + 189] = P_PA_z * AUX_INT__v_s_s_s[m * 171 + 170] - aop_PQ_z * AUX_INT__v_s_s_s[(m+1) * 171 + 170]
                                      + 17 * one_over_2p * ( AUX_INT__u_s_s_s[m * 153 +  152] - a_over_p * AUX_INT__u_s_s_s[(m+1) * 153 + 152] );

                    }
}



// VRR to obtain AUX_INT__x_s_s_s
void VRR_AUX_INT__x_s_s_s(const int num_m,
                        const double P_PA_x, const double P_PA_y, const double P_PA_z,
                        const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,
                        const double a_over_p, const double one_over_2p,
                        double * const restrict AUX_INT__x_s_s_s,
                        double const * const restrict AUX_INT__w_s_s_s,
                        double const * const restrict AUX_INT__v_s_s_s)
{
                    // Forming AUX_INT__x_s_s_s[num_m * 210];
                    for(m = 0; m < num_m; m++)  // loop over orders of auxiliary function
                    {
                        //X_1900 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 0] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 0] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 0]
                                      + 18 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  0] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 0] );

                        //X_1810 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 1] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 0] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 0];

                        //X_1801 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 2] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 0] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 0];

                        //X_1720 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 3] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 1] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  0] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 0] );

                        //X_1711 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 4] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 1] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 1];

                        //X_1702 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 5] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 2] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  0] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 0] );

                        //X_1630 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 6] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 3] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  1] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 1] );

                        //X_1621 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 7] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 3] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 3];

                        //X_1612 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 8] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 5] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 5];

                        //X_1603 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 9] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 5] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  2] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 2] );

                        //X_1540 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 10] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 6] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  3] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 3] );

                        //X_1531 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 11] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 6] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 6];

                        //X_1522 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 12] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 7] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  3] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 3] );

                        //X_1513 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 13] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 9] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 9];

                        //X_1504 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 14] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 9] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  5] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 5] );

                        //X_1450 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 15] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 10] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  6] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 6] );

                        //X_1441 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 16] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 10] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 10];

                        //X_1432 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 17] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 11] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  6] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 6] );

                        //X_1423 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 18] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 13] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  9] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 9] );

                        //X_1414 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 19] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 14] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 14];

                        //X_1405 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 20] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 14] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  9] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 9] );

                        //X_1360 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 21] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 15] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  10] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 10] );

                        //X_1351 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 22] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 15] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 15];

                        //X_1342 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 23] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 16] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  10] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 10] );

                        //X_1333 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 24] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 17] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  11] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 11] );

                        //X_1324 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 25] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 19] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  14] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 14] );

                        //X_1315 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 26] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 20] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 20];

                        //X_1306 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 27] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 20] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  14] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 14] );

                        //X_1270 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 28] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 21] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  15] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 15] );

                        //X_1261 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 29] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 21] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 21];

                        //X_1252 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 30] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 22] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  15] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 15] );

                        //X_1243 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 31] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 23] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  16] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 16] );

                        //X_1234 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 32] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 25] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  19] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 19] );

                        //X_1225 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 33] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 26] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  20] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 20] );

                        //X_1216 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 34] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 27] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 27];

                        //X_1207 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 35] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 27] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  20] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 20] );

                        //X_1180 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 36] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 28] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  21] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 21] );

                        //X_1171 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 37] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 28] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 28];

                        //X_1162 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 38] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 29] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  21] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 21] );

                        //X_1153 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 39] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 30] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  22] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 22] );

                        //X_1144 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 40] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 31] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  23] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 23] );

                        //X_1135 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 41] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 33] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  26] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 26] );

                        //X_1126 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 42] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 34] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  27] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 27] );

                        //X_1117 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 43] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 35] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 35];

                        //X_1108 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 44] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 35] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  27] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 27] );

                        //X_1090 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 45] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 36] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 36]
                                      + 8 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  28] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 28] );

                        //X_1081 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 46] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 36] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 36];

                        //X_1072 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 47] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 37] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 37]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  28] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 28] );

                        //X_1063 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 48] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 38] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 38]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  29] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 29] );

                        //X_1054 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 49] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 39] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 39]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  30] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 30] );

                        //X_1045 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 50] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 41] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 41]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  33] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 33] );

                        //X_1036 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 51] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 42] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 42]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  34] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 34] );

                        //X_1027 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 52] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 43] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 43]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  35] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 35] );

                        //X_1018 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 53] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 44] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 44];

                        //X_1009 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 54] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 44] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 44]
                                      + 8 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  35] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 35] );

                        //X_9100 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 55] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 55] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 55]
                                      + 8 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  55] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 55] );

                        //X_991 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 56] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 45] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 45];

                        //X_982 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 57] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 46] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 46]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  36] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 36] );

                        //X_973 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 58] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 47] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 47]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  37] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 37] );

                        //X_964 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 59] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 48] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 48]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  38] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 38] );

                        //X_955 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 60] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 49] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 49]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  39] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 39] );

                        //X_946 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 61] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 51] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 51]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  42] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 42] );

                        //X_937 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 62] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 52] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 52]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  43] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 43] );

                        //X_928 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 63] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 53] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 53]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  44] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 44] );

                        //X_919 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 64] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 54] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 54];

                        //X_9010 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 65] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 65] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 65]
                                      + 8 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  65] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 65] );

                        //X_8110 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 66] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 66] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 66]
                                      + 7 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  66] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 66] );

                        //X_8101 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 67] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 55] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 55];

                        //X_892 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 68] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 56] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 56]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  45] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 45] );

                        //X_883 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 69] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 57] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 57]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  46] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 46] );

                        //X_874 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 70] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 58] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 58]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  47] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 47] );

                        //X_865 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 71] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 59] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 59]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  48] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 48] );

                        //X_856 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 72] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 61] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 61]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  51] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 51] );

                        //X_847 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 73] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 62] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 62]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  52] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 52] );

                        //X_838 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 74] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 63] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 63]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  53] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 53] );

                        //X_829 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 75] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 64] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 64]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  54] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 54] );

                        //X_8110 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 76] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 65] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 65];

                        //X_8011 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 77] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 77] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 77]
                                      + 7 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  77] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 77] );

                        //X_7120 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 78] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 78] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 78]
                                      + 6 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  78] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 78] );

                        //X_7111 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 79] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 66] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 66];

                        //X_7102 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 80] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 67] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 67]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  55] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 55] );

                        //X_793 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 81] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 68] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 68]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  56] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 56] );

                        //X_784 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 82] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 69] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 69]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  57] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 57] );

                        //X_775 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 83] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 70] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 70]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  58] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 58] );

                        //X_766 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 84] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 71] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 71]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  59] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 59] );

                        //X_757 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 85] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 73] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 73]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  62] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 62] );

                        //X_748 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 86] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 74] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 74]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  63] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 63] );

                        //X_739 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 87] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 75] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 75]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  64] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 64] );

                        //X_7210 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 88] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 76] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 76]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  65] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 65] );

                        //X_7111 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 89] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 77] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 77];

                        //X_7012 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 90] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 90] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 90]
                                      + 6 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  90] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 90] );

                        //X_6130 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 91] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 91] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 91]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  91] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 91] );

                        //X_6121 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 92] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 78] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 78];

                        //X_6112 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 93] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 79] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 79]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  66] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 66] );

                        //X_6103 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 94] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 80] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 80]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  67] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 67] );

                        //X_694 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 95] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 81] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 81]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  68] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 68] );

                        //X_685 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 96] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 82] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 82]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  69] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 69] );

                        //X_676 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 97] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 83] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 83]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  70] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 70] );

                        //X_667 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 98] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 85] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 85]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  73] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 73] );

                        //X_658 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 99] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 86] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 86]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  74] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 74] );

                        //X_649 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 100] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 87] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 87]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  75] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 75] );

                        //X_6310 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 101] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 88] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 88]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  76] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 76] );

                        //X_6211 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 102] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 89] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 89]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  77] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 77] );

                        //X_6112 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 103] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 90] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 90];

                        //X_6013 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 104] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 104] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 104]
                                      + 5 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  104] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 104] );

                        //X_5140 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 105] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 105] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 105]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  105] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 105] );

                        //X_5131 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 106] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 91] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 91];

                        //X_5122 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 107] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 92] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 92]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  78] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 78] );

                        //X_5113 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 108] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 93] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 93]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  79] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 79] );

                        //X_5104 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 109] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 94] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 94]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  80] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 80] );

                        //X_595 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 110] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 95] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 95]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  81] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 81] );

                        //X_586 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 111] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 111] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 111]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  111] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 111] );

                        //X_577 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 112] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 112] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 112]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  112] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 112] );

                        //X_568 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 113] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 113] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 113]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  113] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 113] );

                        //X_559 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 114] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 100] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 100]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  87] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 87] );

                        //X_5410 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 115] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 101] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 101]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  88] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 88] );

                        //X_5311 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 116] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 102] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 102]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  89] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 89] );

                        //X_5212 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 117] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 103] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 103]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  90] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 90] );

                        //X_5113 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 118] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 104] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 104];

                        //X_5014 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 119] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 119] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 119]
                                      + 4 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  119] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 119] );

                        //X_4150 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 120] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 120] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 120]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  120] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 120] );

                        //X_4141 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 121] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 105] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 105];

                        //X_4132 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 122] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 106] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 106]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  91] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 91] );

                        //X_4123 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 123] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 107] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 107]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  92] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 92] );

                        //X_4114 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 124] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 108] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 108]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  93] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 93] );

                        //X_4105 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 125] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 125] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 125]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  125] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 125] );

                        //X_496 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 126] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 126] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 126]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  126] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 126] );

                        //X_487 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 127] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 127] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 127]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  127] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 127] );

                        //X_478 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 128] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 128] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 128]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  128] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 128] );

                        //X_469 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 129] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 129] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 129]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  129] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 129] );

                        //X_4510 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 130] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 130] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 130]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  130] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 130] );

                        //X_4411 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 131] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 116] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 116]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  102] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 102] );

                        //X_4312 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 132] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 117] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 117]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  103] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 103] );

                        //X_4213 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 133] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 118] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 118]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  104] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 104] );

                        //X_4114 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 134] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 119] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 119];

                        //X_4015 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 135] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 135] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 135]
                                      + 3 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  135] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 135] );

                        //X_3160 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 136] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 136] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 136]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  136] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 136] );

                        //X_3151 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 137] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 120] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 120];

                        //X_3142 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 138] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 121] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 121]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  105] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 105] );

                        //X_3133 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 139] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 122] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 122]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  106] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 106] );

                        //X_3124 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 140] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 140] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 140]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  140] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 140] );

                        //X_3115 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 141] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 141] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 141]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  141] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 141] );

                        //X_3106 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 142] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 142] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 142]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  142] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 142] );

                        //X_397 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 143] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 143] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 143]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  143] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 143] );

                        //X_388 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 144] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 144] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 144]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  144] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 144] );

                        //X_379 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 145] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 145] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 145]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  145] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 145] );

                        //X_3610 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 146] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 146] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 146]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  146] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 146] );

                        //X_3511 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 147] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 147] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 147]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  147] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 147] );

                        //X_3412 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 148] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 148] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 148]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  148] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 148] );

                        //X_3313 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 149] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 133] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 133]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  118] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 118] );

                        //X_3214 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 150] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 134] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 134]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  119] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 119] );

                        //X_3115 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 151] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 135] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 135];

                        //X_3016 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 152] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 152] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 152]
                                      + 2 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  152] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 152] );

                        //X_2170 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 153] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 153] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 153]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  153] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 153] );

                        //X_2161 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 154] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 136] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 136];

                        //X_2152 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 155] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 137] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 137]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  120] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 120] );

                        //X_2143 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 156] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 156] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 156]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  156] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 156] );

                        //X_2134 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 157] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 157] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 157]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  157] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 157] );

                        //X_2125 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 158] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 158] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 158]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  158] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 158] );

                        //X_2116 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 159] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 159] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 159]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  159] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 159] );

                        //X_2107 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 160] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 160] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 160]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  160] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 160] );

                        //X_298 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 161] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 161] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 161]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  161] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 161] );

                        //X_289 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 162] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 162] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 162]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  162] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 162] );

                        //X_2710 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 163] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 163] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 163]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  163] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 163] );

                        //X_2611 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 164] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 164] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 164]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  164] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 164] );

                        //X_2512 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 165] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 165] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 165]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  165] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 165] );

                        //X_2413 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 166] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 166] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 166]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  166] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 166] );

                        //X_2314 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 167] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 167] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 167]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  167] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 167] );

                        //X_2215 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 168] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 151] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 151]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  135] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 135] );

                        //X_2116 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 169] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 152] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 152];

                        //X_2017 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 170] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 170] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 170]
                                      + 1 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  170] - a_over_p * AUX_INT__v_s_s_s[(m+1) * 171 + 170] );

                        //X_1180 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 171] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 171] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 171];

                        //X_1171 : STEP: z
                        AUX_INT__x_s_s_s[m * 210 + 172] = P_PA_z * AUX_INT__w_s_s_s[m * 190 + 153] - aop_PQ_z * AUX_INT__w_s_s_s[(m+1) * 190 + 153];

                        //X_1162 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 173] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 173] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 173];

                        //X_1153 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 174] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 174] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 174];

                        //X_1144 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 175] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 175] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 175];

                        //X_1135 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 176] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 176] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 176];

                        //X_1126 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 177] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 177] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 177];

                        //X_1117 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 178] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 178] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 178];

                        //X_1108 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 179] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 179] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 179];

                        //X_199 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 180] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 180] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 180];

                        //X_1810 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 181] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 181] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 181];

                        //X_1711 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 182] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 182] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 182];

                        //X_1612 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 183] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 183] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 183];

                        //X_1513 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 184] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 184] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 184];

                        //X_1414 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 185] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 185] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 185];

                        //X_1315 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 186] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 186] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 186];

                        //X_1216 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 187] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 187] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 187];

                        //X_1117 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 188] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 170] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 170];

                        //X_1018 : STEP: x
                        AUX_INT__x_s_s_s[m * 210 + 189] = P_PA_x * AUX_INT__w_s_s_s[m * 190 + 189] - aop_PQ_x * AUX_INT__w_s_s_s[(m+1) * 190 + 189];

                        //X_0190 : STEP: y
                        AUX_INT__x_s_s_s[m * 210 + 190] = P_PA_y * AUX_INT__w_s_s_s[m * 190 + 171] - aop_PQ_y * AUX_INT__w_s_s_s[(m+1) * 190 + 171]
                                      + 18 * one_over_2p * ( AUX_INT__v_s_s_s[m * 171 +  153]