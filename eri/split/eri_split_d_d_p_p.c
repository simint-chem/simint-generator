#include <string.h>
#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"
#include "boys/boys_split.h"


int eri_split_d_d_p_p(struct multishell_pair const P,
                      struct multishell_pair const Q,
                      double * const restrict INT__d_d_p_p)
{

    ASSUME_ALIGN(P.x);
    ASSUME_ALIGN(P.y);
    ASSUME_ALIGN(P.z);
    ASSUME_ALIGN(P.PA_x);
    ASSUME_ALIGN(P.PA_y);
    ASSUME_ALIGN(P.PA_z);
    ASSUME_ALIGN(P.bAB_x);
    ASSUME_ALIGN(P.bAB_y);
    ASSUME_ALIGN(P.bAB_z);
    ASSUME_ALIGN(P.alpha);
    ASSUME_ALIGN(P.prefac);

    ASSUME_ALIGN(Q.x);
    ASSUME_ALIGN(Q.y);
    ASSUME_ALIGN(Q.z);
    ASSUME_ALIGN(Q.PA_x);
    ASSUME_ALIGN(Q.PA_y);
    ASSUME_ALIGN(Q.PA_z);
    ASSUME_ALIGN(Q.bAB_x);
    ASSUME_ALIGN(Q.bAB_y);
    ASSUME_ALIGN(Q.bAB_z);
    ASSUME_ALIGN(Q.alpha);
    ASSUME_ALIGN(Q.prefac);

    ASSUME_ALIGN(INT__d_d_p_p);

    const int nshell1234 = P.nshell12 * Q.nshell12;


    // Holds AB_{xyz} and CD_{xyz} in a flattened fashion for later
    double AB_x[nshell1234];  double CD_x[nshell1234];
    double AB_y[nshell1234];  double CD_y[nshell1234];
    double AB_z[nshell1234];  double CD_z[nshell1234];

    int ab, cd, abcd;
    int i, j;
    int m;
    int n;
    int iket;
    int ibra;

    // Workspace for contracted integrals
    double * const contwork = malloc(nshell1234 * 4824);
    memset(contwork, 0, nshell1234 * 4824);

    // partition workspace into shells
    double * const INT__d_s_p_s = contwork + (nshell1234 * 0);
    double * const INT__d_s_d_s = contwork + (nshell1234 * 18);
    double * const INT__d_d_p_s = contwork + (nshell1234 * 54);
    double * const INT__d_d_d_s = contwork + (nshell1234 * 162);
    double * const INT__f_s_p_s = contwork + (nshell1234 * 378);
    double * const INT__f_s_d_s = contwork + (nshell1234 * 408);
    double * const INT__g_s_p_s = contwork + (nshell1234 * 468);
    double * const INT__g_s_d_s = contwork + (nshell1234 * 513);


    ////////////////////////////////////////
    // Loop over shells and primitives
    ////////////////////////////////////////
    for(ab = 0, abcd = 0; ab < P.nshell12; ++ab)
    {
        const int abstart = P.primstart[ab];
        const int abend = P.primend[ab];

        // this should have been set/aligned in fill_multishell_pair or something else
        ASSUME(abstart%SIMD_ALIGN_DBL == 0);

        for(cd = 0; cd < Q.nshell12; ++cd, ++abcd)
        {
            // set up pointers to the contracted integrals - VRR
            // set up pointers to the contracted integrals - Electron Transfer
            double * const restrict PRIM_INT__d_s_p_s = INT__d_s_p_s + (abcd * 18);
            double * const restrict PRIM_INT__d_s_d_s = INT__d_s_d_s + (abcd * 36);
            double * const restrict PRIM_INT__f_s_p_s = INT__f_s_p_s + (abcd * 30);
            double * const restrict PRIM_INT__f_s_d_s = INT__f_s_d_s + (abcd * 60);
            double * const restrict PRIM_INT__g_s_p_s = INT__g_s_p_s + (abcd * 45);
            double * const restrict PRIM_INT__g_s_d_s = INT__g_s_d_s + (abcd * 90);

            const int cdstart = Q.primstart[cd];
            const int cdend = Q.primend[cd];

            // this should have been set/aligned in fill_multishell_pair or something else
            ASSUME(cdstart%SIMD_ALIGN_DBL == 0);

            // Store for later
            AB_x[abcd] = P.AB_x[ab];  CD_x[abcd] = Q.AB_x[cd];
            AB_y[abcd] = P.AB_y[ab];  CD_y[abcd] = Q.AB_y[cd];
            AB_z[abcd] = P.AB_z[ab];  CD_z[abcd] = Q.AB_z[cd];

            for(i = abstart; i < abend; ++i)
            {

                // Load these one per loop over i
                const double P_alpha = P.alpha[i];
                const double P_prefac = P.prefac[i];
                const double P_x = P.x[i];
                const double P_y = P.y[i];
                const double P_z = P.z[i];
                const double P_PA_x = P.PA_x[i];
                const double P_PA_y = P.PA_y[i];
                const double P_PA_z = P.PA_z[i];
                const double P_bAB_x = P.bAB_x[i];
                const double P_bAB_y = P.bAB_y[i];
                const double P_bAB_z = P.bAB_z[i];

                for(j = cdstart; j < cdend; ++j)
                {

                    // Holds the auxiliary integrals ( i 0 | 0 0 )^m in the primitive basis
                    // with m as the slowest index
                    // AM = 0: Needed from this AM: 1
                    double AUX_INT__s_s_s_s[7 * 1];

                    // AM = 1: Needed from this AM: 3
                    double AUX_INT__p_s_s_s[6 * 3];

                    // AM = 2: Needed from this AM: 6
                    double AUX_INT__d_s_s_s[5 * 6];

                    // AM = 3: Needed from this AM: 10
                    double AUX_INT__f_s_s_s[4 * 10];

                    // AM = 4: Needed from this AM: 15
                    double AUX_INT__g_s_s_s[3 * 15];

                    // AM = 5: Needed from this AM: 21
                    double AUX_INT__h_s_s_s[2 * 21];

                    // AM = 6: Needed from this AM: 28
                    double AUX_INT__i_s_s_s[1 * 28];



                    // Holds temporary integrals for electron transfer
                    double AUX_INT__p_s_p_s[9];
                    double AUX_INT__d_s_p_s[18];
                    double AUX_INT__d_s_d_s[36];
                    double AUX_INT__f_s_p_s[30];
                    double AUX_INT__f_s_d_s[60];
                    double AUX_INT__g_s_p_s[45];
                    double AUX_INT__g_s_d_s[90];
                    double AUX_INT__h_s_p_s[63];


                    const double PQalpha_mul = P_alpha * Q.alpha[j];
                    const double PQalpha_sum = P_alpha + Q.alpha[j];

                    const double pfac = TWO_PI_52 / (PQalpha_mul * sqrt(PQalpha_sum));

                    /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
                    const double PQ_x = P_x - Q.x[j];
                    const double PQ_y = P_y - Q.y[j];
                    const double PQ_z = P_z - Q.z[j];
                    const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;

                    // collected prefactors
                    const double allprefac =  pfac * P_prefac * Q.prefac[j];

                    // various factors
                    const double alpha = PQalpha_mul/PQalpha_sum;   // alpha from MEST
                    // for VRR
                    const double one_over_p = 1.0 / P_alpha;
                    const double a_over_p =  alpha * one_over_p;     // a/p from MEST
                    const double one_over_2p = 0.5 * one_over_p;  // gets multiplied by i in VRR
                    // for electron transfer
                    const double one_over_q = 1.0 / Q.alpha[j];
                    const double one_over_2q = 0.5 * one_over_q;
                    const double p_over_q = P_alpha * one_over_q;

                    const double etfac[3] = {
                                             -(P_bAB_x + Q.bAB_x[j]) * one_over_q,
                                             -(P_bAB_y + Q.bAB_y[j]) * one_over_q,
                                             -(P_bAB_z + Q.bAB_z[j]) * one_over_q,
                                            };


                    //////////////////////////////////////////////
                    // Boys function section
                    // Maximum v value: 6
                    //////////////////////////////////////////////
                    // The paremeter to the boys function
                    const double F_x = R2 * alpha;


                    Boys_F_split(AUX_INT__s_s_s_s, 6, F_x);
                    AUX_INT__s_s_s_s[0] *= allprefac;
                    AUX_INT__s_s_s_s[1] *= allprefac;
                    AUX_INT__s_s_s_s[2] *= allprefac;
                    AUX_INT__s_s_s_s[3] *= allprefac;
                    AUX_INT__s_s_s_s[4] *= allprefac;
                    AUX_INT__s_s_s_s[5] *= allprefac;
                    AUX_INT__s_s_s_s[6] *= allprefac;

                    //////////////////////////////////////////////
                    // Primitive integrals: Vertical recurrance
                    //////////////////////////////////////////////

                    // Forming AUX_INT__p_s_s_s[6 * 3];
                    // Needed from this AM:
                    //    P_100
                    //    P_010
                    //    P_001
                    for(m = 0; m < 6; m++)  // loop over orders of boys function
                    {
                        //P_100 : STEP: x
                        AUX_INT__p_s_s_s[m * 3 + 0] = P_PA_x * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_x * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_010 : STEP: y
                        AUX_INT__p_s_s_s[m * 3 + 1] = P_PA_y * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_y * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_001 : STEP: z
                        AUX_INT__p_s_s_s[m * 3 + 2] = P_PA_z * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_z * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                    }


                    // Forming AUX_INT__d_s_s_s[5 * 6];
                    // Needed from this AM:
                    //    D_200
                    //    D_110
                    //    D_101
                    //    D_020
                    //    D_011
                    //    D_002
                    for(m = 0; m < 5; m++)  // loop over orders of boys function
                    {
                        //D_200 : STEP: x
                        AUX_INT__d_s_s_s[m * 6 + 0] = P_PA_x * AUX_INT__p_s_s_s[m * 3 + 0] - a_over_p * PQ_x * AUX_INT__p_s_s_s[(m+1) * 3 + 0]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                        //D_110 : STEP: y
                        AUX_INT__d_s_s_s[m * 6 + 1] = P_PA_y * AUX_INT__p_s_s_s[m * 3 + 0] - a_over_p * PQ_y * AUX_INT__p_s_s_s[(m+1) * 3 + 0];

                        //D_101 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 2] = P_PA_z * AUX_INT__p_s_s_s[m * 3 + 0] - a_over_p * PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 0];

                        //D_020 : STEP: y
                        AUX_INT__d_s_s_s[m * 6 + 3] = P_PA_y * AUX_INT__p_s_s_s[m * 3 + 1] - a_over_p * PQ_y * AUX_INT__p_s_s_s[(m+1) * 3 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                        //D_011 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 4] = P_PA_z * AUX_INT__p_s_s_s[m * 3 + 1] - a_over_p * PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 1];

                        //D_002 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 5] = P_PA_z * AUX_INT__p_s_s_s[m * 3 + 2] - a_over_p * PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                    }


                    // Forming AUX_INT__f_s_s_s[4 * 10];
                    // Needed from this AM:
                    //    F_300
                    //    F_210
                    //    F_201
                    //    F_120
                    //    F_111
                    //    F_102
                    //    F_030
                    //    F_021
                    //    F_012
                    //    F_003
                    for(m = 0; m < 4; m++)  // loop over orders of boys function
                    {
                        //F_300 : STEP: x
                        AUX_INT__f_s_s_s[m * 10 + 0] = P_PA_x * AUX_INT__d_s_s_s[m * 6 + 0] - a_over_p * PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 0]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  0] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 0] );

                        //F_210 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 1] = P_PA_y * AUX_INT__d_s_s_s[m * 6 + 0] - a_over_p * PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 0];

                        //F_201 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 2] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 0] - a_over_p * PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 0];

                        //F_120 : STEP: x
                        AUX_INT__f_s_s_s[m * 10 + 3] = P_PA_x * AUX_INT__d_s_s_s[m * 6 + 3] - a_over_p * PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 3];

                        //F_111 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 4] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 1] - a_over_p * PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 1];

                        //F_102 : STEP: x
                        AUX_INT__f_s_s_s[m * 10 + 5] = P_PA_x * AUX_INT__d_s_s_s[m * 6 + 5] - a_over_p * PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 5];

                        //F_030 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 6] = P_PA_y * AUX_INT__d_s_s_s[m * 6 + 3] - a_over_p * PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  1] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 1] );

                        //F_021 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 7] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 3] - a_over_p * PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 3];

                        //F_012 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 8] = P_PA_y * AUX_INT__d_s_s_s[m * 6 + 5] - a_over_p * PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 5];

                        //F_003 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 9] = P_PA_z * AUX_INT__d_s_s_s[m * 6 + 5] - a_over_p * PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  2] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 2] );

                    }


                    // Forming AUX_INT__g_s_s_s[3 * 15];
                    // Needed from this AM:
                    //    G_400
                    //    G_310
                    //    G_301
                    //    G_220
                    //    G_211
                    //    G_202
                    //    G_130
                    //    G_121
                    //    G_112
                    //    G_103
                    //    G_040
                    //    G_031
                    //    G_022
                    //    G_013
                    //    G_004
                    for(m = 0; m < 3; m++)  // loop over orders of boys function
                    {
                        //G_400 : STEP: x
                        AUX_INT__g_s_s_s[m * 15 + 0] = P_PA_x * AUX_INT__f_s_s_s[m * 10 + 0] - a_over_p * PQ_x * AUX_INT__f_s_s_s[(m+1) * 10 + 0]
                                      + 3 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  0] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 0] );

                        //G_310 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 1] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 0] - a_over_p * PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 0];

                        //G_301 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 2] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 0] - a_over_p * PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 0];

                        //G_220 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 3] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 1] - a_over_p * PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  0] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 0] );

                        //G_211 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 4] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 1] - a_over_p * PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 1];

                        //G_202 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 5] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 2] - a_over_p * PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  0] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 0] );

                        //G_130 : STEP: x
                        AUX_INT__g_s_s_s[m * 15 + 6] = P_PA_x * AUX_INT__f_s_s_s[m * 10 + 6] - a_over_p * PQ_x * AUX_INT__f_s_s_s[(m+1) * 10 + 6];

                        //G_121 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 7] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 3] - a_over_p * PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 3];

                        //G_112 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 8] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 5] - a_over_p * PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 5];

                        //G_103 : STEP: x
                        AUX_INT__g_s_s_s[m * 15 + 9] = P_PA_x * AUX_INT__f_s_s_s[m * 10 + 9] - a_over_p * PQ_x * AUX_INT__f_s_s_s[(m+1) * 10 + 9];

                        //G_040 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 10] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 6] - a_over_p * PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  3] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 3] );

                        //G_031 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 11] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 6] - a_over_p * PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 6];

                        //G_022 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 12] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 7] - a_over_p * PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  3] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 3] );

                        //G_013 : STEP: y
                        AUX_INT__g_s_s_s[m * 15 + 13] = P_PA_y * AUX_INT__f_s_s_s[m * 10 + 9] - a_over_p * PQ_y * AUX_INT__f_s_s_s[(m+1) * 10 + 9];

                        //G_004 : STEP: z
                        AUX_INT__g_s_s_s[m * 15 + 14] = P_PA_z * AUX_INT__f_s_s_s[m * 10 + 9] - a_over_p * PQ_z * AUX_INT__f_s_s_s[(m+1) * 10 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__d_s_s_s[m * 6 +  5] - a_over_p * AUX_INT__d_s_s_s[(m+1) * 6 + 5] );

                    }


                    // Forming AUX_INT__h_s_s_s[2 * 21];
                    // Needed from this AM:
                    //    H_500
                    //    H_410
                    //    H_401
                    //    H_320
                    //    H_311
                    //    H_302
                    //    H_230
                    //    H_221
                    //    H_212
                    //    H_203
                    //    H_140
                    //    H_131
                    //    H_122
                    //    H_113
                    //    H_104
                    //    H_050
                    //    H_041
                    //    H_032
                    //    H_023
                    //    H_014
                    //    H_005
                    for(m = 0; m < 2; m++)  // loop over orders of boys function
                    {
                        //H_500 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 0] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 0] - a_over_p * PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 0]
                                      + 4 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  0] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 0] );

                        //H_410 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 1] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 0] - a_over_p * PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 0];

                        //H_401 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 2] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 0] - a_over_p * PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 0];

                        //H_320 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 3] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 1] - a_over_p * PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  0] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 0] );

                        //H_311 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 4] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 1] - a_over_p * PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 1];

                        //H_302 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 5] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 2] - a_over_p * PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  0] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 0] );

                        //H_230 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 6] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 6] - a_over_p * PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 6]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  6] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 6] );

                        //H_221 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 7] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 3] - a_over_p * PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 3];

                        //H_212 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 8] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 5] - a_over_p * PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 5];

                        //H_203 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 9] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 9] - a_over_p * PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 9]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  9] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 9] );

                        //H_140 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 10] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 10] - a_over_p * PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 10];

                        //H_131 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 11] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 6] - a_over_p * PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 6];

                        //H_122 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 12] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 12] - a_over_p * PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 12];

                        //H_113 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 13] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 9] - a_over_p * PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 9];

                        //H_104 : STEP: x
                        AUX_INT__h_s_s_s[m * 21 + 14] = P_PA_x * AUX_INT__g_s_s_s[m * 15 + 14] - a_over_p * PQ_x * AUX_INT__g_s_s_s[(m+1) * 15 + 14];

                        //H_050 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 15] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 10] - a_over_p * PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 10]
                                      + 4 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  6] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 6] );

                        //H_041 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 16] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 10] - a_over_p * PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 10];

                        //H_032 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 17] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 11] - a_over_p * PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  6] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 6] );

                        //H_023 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 18] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 13] - a_over_p * PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  9] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 9] );

                        //H_014 : STEP: y
                        AUX_INT__h_s_s_s[m * 21 + 19] = P_PA_y * AUX_INT__g_s_s_s[m * 15 + 14] - a_over_p * PQ_y * AUX_INT__g_s_s_s[(m+1) * 15 + 14];

                        //H_005 : STEP: z
                        AUX_INT__h_s_s_s[m * 21 + 20] = P_PA_z * AUX_INT__g_s_s_s[m * 15 + 14] - a_over_p * PQ_z * AUX_INT__g_s_s_s[(m+1) * 15 + 14]
                                      + 4 * one_over_2p * ( AUX_INT__f_s_s_s[m * 10 +  9] - a_over_p * AUX_INT__f_s_s_s[(m+1) * 10 + 9] );

                    }


                    // Forming AUX_INT__i_s_s_s[1 * 28];
                    // Needed from this AM:
                    //    I_600
                    //    I_510
                    //    I_501
                    //    I_420
                    //    I_411
                    //    I_402
                    //    I_330
                    //    I_321
                    //    I_312
                    //    I_303
                    //    I_240
                    //    I_231
                    //    I_222
                    //    I_213
                    //    I_204
                    //    I_150
                    //    I_141
                    //    I_132
                    //    I_123
                    //    I_114
                    //    I_105
                    //    I_060
                    //    I_051
                    //    I_042
                    //    I_033
                    //    I_024
                    //    I_015
                    //    I_006
                    for(m = 0; m < 1; m++)  // loop over orders of boys function
                    {
                        //I_600 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 0] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 0] - a_over_p * PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 0]
                                      + 5 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  0] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 0] );

                        //I_510 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 1] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 0] - a_over_p * PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 0];

                        //I_501 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 2] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 0] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 0];

                        //I_420 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 3] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 1] - a_over_p * PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  0] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 0] );

                        //I_411 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 4] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 1] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 1];

                        //I_402 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 5] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 2] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  0] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 0] );

                        //I_330 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 6] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 3] - a_over_p * PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  1] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 1] );

                        //I_321 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 7] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 3] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 3];

                        //I_312 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 8] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 5] - a_over_p * PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 5];

                        //I_303 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 9] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 5] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  2] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 2] );

                        //I_240 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 10] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 10] - a_over_p * PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 10]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  10] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 10] );

                        //I_231 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 11] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 6] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 6];

                        //I_222 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 12] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 7] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  3] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 3] );

                        //I_213 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 13] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 9] - a_over_p * PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 9];

                        //I_204 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 14] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 14] - a_over_p * PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 14]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  14] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 14] );

                        //I_150 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 15] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 15] - a_over_p * PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 15];

                        //I_141 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 16] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 10] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 10];

                        //I_132 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 17] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 17] - a_over_p * PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 17];

                        //I_123 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 18] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 18] - a_over_p * PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 18];

                        //I_114 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 19] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 14] - a_over_p * PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 14];

                        //I_105 : STEP: x
                        AUX_INT__i_s_s_s[m * 28 + 20] = P_PA_x * AUX_INT__h_s_s_s[m * 21 + 20] - a_over_p * PQ_x * AUX_INT__h_s_s_s[(m+1) * 21 + 20];

                        //I_060 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 21] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 15] - a_over_p * PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 15]
                                      + 5 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  10] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 10] );

                        //I_051 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 22] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 15] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 15];

                        //I_042 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 23] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 16] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  10] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 10] );

                        //I_033 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 24] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 17] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 17]
                                      + 2 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  11] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 11] );

                        //I_024 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 25] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 19] - a_over_p * PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  14] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 14] );

                        //I_015 : STEP: y
                        AUX_INT__i_s_s_s[m * 28 + 26] = P_PA_y * AUX_INT__h_s_s_s[m * 21 + 20] - a_over_p * PQ_y * AUX_INT__h_s_s_s[(m+1) * 21 + 20];

                        //I_006 : STEP: z
                        AUX_INT__i_s_s_s[m * 28 + 27] = P_PA_z * AUX_INT__h_s_s_s[m * 21 + 20] - a_over_p * PQ_z * AUX_INT__h_s_s_s[(m+1) * 21 + 20]
                                      + 5 * one_over_2p * ( AUX_INT__g_s_s_s[m * 15 +  14] - a_over_p * AUX_INT__g_s_s_s[(m+1) * 15 + 14] );

                    }




                    //////////////////////////////////////////////
                    // Primitive integrals: Electron transfer
                    //////////////////////////////////////////////

                    // ( D_200 S_000 | P_100 S_000 )^0_{t} = x * ( D_200 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_300 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[0] = etfac[0] * AUX_INT__d_s_s_s[0] + 2 * one_over_2q * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__f_s_s_s[0];

                    // ( D_200 S_000 | P_010 S_000 )^0_{t} = y * ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[1] = etfac[1] * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__f_s_s_s[1];

                    // ( D_200 S_000 | P_001 S_000 )^0_{t} = z * ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[2] = etfac[2] * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__f_s_s_s[2];

                    // ( D_110 S_000 | P_100 S_000 )^0_{t} = x * ( D_110 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[3] = etfac[0] * AUX_INT__d_s_s_s[1] + 1 * one_over_2q * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__f_s_s_s[1];

                    // ( D_110 S_000 | P_010 S_000 )^0_{t} = y * ( D_110 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[4] = etfac[1] * AUX_INT__d_s_s_s[1] + 1 * one_over_2q * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__f_s_s_s[3];

                    // ( D_110 S_000 | P_001 S_000 )^0_{t} = z * ( D_110 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[5] = etfac[2] * AUX_INT__d_s_s_s[1] - p_over_q * AUX_INT__f_s_s_s[4];

                    // ( D_101 S_000 | P_100 S_000 )^0_{t} = x * ( D_101 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[6] = etfac[0] * AUX_INT__d_s_s_s[2] + 1 * one_over_2q * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__f_s_s_s[2];

                    // ( D_101 S_000 | P_010 S_000 )^0_{t} = y * ( D_101 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[7] = etfac[1] * AUX_INT__d_s_s_s[2] - p_over_q * AUX_INT__f_s_s_s[4];

                    // ( D_101 S_000 | P_001 S_000 )^0_{t} = z * ( D_101 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[8] = etfac[2] * AUX_INT__d_s_s_s[2] + 1 * one_over_2q * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__f_s_s_s[5];

                    // ( D_020 S_000 | P_100 S_000 )^0_{t} = x * ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[9] = etfac[0] * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__f_s_s_s[3];

                    // ( D_020 S_000 | P_010 S_000 )^0_{t} = y * ( D_020 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_030 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[10] = etfac[1] * AUX_INT__d_s_s_s[3] + 2 * one_over_2q * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__f_s_s_s[6];

                    // ( D_020 S_000 | P_001 S_000 )^0_{t} = z * ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[11] = etfac[2] * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__f_s_s_s[7];

                    // ( D_011 S_000 | P_100 S_000 )^0_{t} = x * ( D_011 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[12] = etfac[0] * AUX_INT__d_s_s_s[4] - p_over_q * AUX_INT__f_s_s_s[4];

                    // ( D_011 S_000 | P_010 S_000 )^0_{t} = y * ( D_011 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[13] = etfac[1] * AUX_INT__d_s_s_s[4] + 1 * one_over_2q * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__f_s_s_s[7];

                    // ( D_011 S_000 | P_001 S_000 )^0_{t} = z * ( D_011 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[14] = etfac[2] * AUX_INT__d_s_s_s[4] + 1 * one_over_2q * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__f_s_s_s[8];

                    // ( D_002 S_000 | P_100 S_000 )^0_{t} = x * ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[15] = etfac[0] * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__f_s_s_s[5];

                    // ( D_002 S_000 | P_010 S_000 )^0_{t} = y * ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[16] = etfac[1] * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__f_s_s_s[8];

                    // ( D_002 S_000 | P_001 S_000 )^0_{t} = z * ( D_002 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_003 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[17] = etfac[2] * AUX_INT__d_s_s_s[5] + 2 * one_over_2q * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__f_s_s_s[9];

                    // ( F_300 S_000 | P_100 S_000 )^0_{t} = x * ( F_300 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_400 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[0] = etfac[0] * AUX_INT__f_s_s_s[0] + 3 * one_over_2q * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__g_s_s_s[0];

                    // ( F_300 S_000 | P_010 S_000 )^0_{t} = y * ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[1] = etfac[1] * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__g_s_s_s[1];

                    // ( F_300 S_000 | P_001 S_000 )^0_{t} = z * ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[2] = etfac[2] * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__g_s_s_s[2];

                    // ( F_210 S_000 | P_100 S_000 )^0_{t} = x * ( F_210 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[3] = etfac[0] * AUX_INT__f_s_s_s[1] + 2 * one_over_2q * AUX_INT__d_s_s_s[1] - p_over_q * AUX_INT__g_s_s_s[1];

                    // ( F_210 S_000 | P_010 S_000 )^0_{t} = y * ( F_210 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[4] = etfac[1] * AUX_INT__f_s_s_s[1] + 1 * one_over_2q * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__g_s_s_s[3];

                    // ( F_210 S_000 | P_001 S_000 )^0_{t} = z * ( F_210 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[5] = etfac[2] * AUX_INT__f_s_s_s[1] - p_over_q * AUX_INT__g_s_s_s[4];

                    // ( F_201 S_000 | P_100 S_000 )^0_{t} = x * ( F_201 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[6] = etfac[0] * AUX_INT__f_s_s_s[2] + 2 * one_over_2q * AUX_INT__d_s_s_s[2] - p_over_q * AUX_INT__g_s_s_s[2];

                    // ( F_201 S_000 | P_010 S_000 )^0_{t} = y * ( F_201 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[7] = etfac[1] * AUX_INT__f_s_s_s[2] - p_over_q * AUX_INT__g_s_s_s[4];

                    // ( F_201 S_000 | P_001 S_000 )^0_{t} = z * ( F_201 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[8] = etfac[2] * AUX_INT__f_s_s_s[2] + 1 * one_over_2q * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__g_s_s_s[5];

                    // ( F_120 S_000 | P_100 S_000 )^0_{t} = x * ( F_120 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[9] = etfac[0] * AUX_INT__f_s_s_s[3] + 1 * one_over_2q * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__g_s_s_s[3];

                    // ( F_120 S_000 | P_010 S_000 )^0_{t} = y * ( F_120 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[10] = etfac[1] * AUX_INT__f_s_s_s[3] + 2 * one_over_2q * AUX_INT__d_s_s_s[1] - p_over_q * AUX_INT__g_s_s_s[6];

                    // ( F_120 S_000 | P_001 S_000 )^0_{t} = z * ( F_120 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[11] = etfac[2] * AUX_INT__f_s_s_s[3] - p_over_q * AUX_INT__g_s_s_s[7];

                    // ( F_111 S_000 | P_100 S_000 )^0_{t} = x * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[12] = etfac[0] * AUX_INT__f_s_s_s[4] + 1 * one_over_2q * AUX_INT__d_s_s_s[4] - p_over_q * AUX_INT__g_s_s_s[4];

                    // ( F_111 S_000 | P_010 S_000 )^0_{t} = y * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[13] = etfac[1] * AUX_INT__f_s_s_s[4] + 1 * one_over_2q * AUX_INT__d_s_s_s[2] - p_over_q * AUX_INT__g_s_s_s[7];

                    // ( F_111 S_000 | P_001 S_000 )^0_{t} = z * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[14] = etfac[2] * AUX_INT__f_s_s_s[4] + 1 * one_over_2q * AUX_INT__d_s_s_s[1] - p_over_q * AUX_INT__g_s_s_s[8];

                    // ( F_102 S_000 | P_100 S_000 )^0_{t} = x * ( F_102 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[15] = etfac[0] * AUX_INT__f_s_s_s[5] + 1 * one_over_2q * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__g_s_s_s[5];

                    // ( F_102 S_000 | P_010 S_000 )^0_{t} = y * ( F_102 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[16] = etfac[1] * AUX_INT__f_s_s_s[5] - p_over_q * AUX_INT__g_s_s_s[8];

                    // ( F_102 S_000 | P_001 S_000 )^0_{t} = z * ( F_102 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[17] = etfac[2] * AUX_INT__f_s_s_s[5] + 2 * one_over_2q * AUX_INT__d_s_s_s[2] - p_over_q * AUX_INT__g_s_s_s[9];

                    // ( F_030 S_000 | P_100 S_000 )^0_{t} = x * ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[18] = etfac[0] * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__g_s_s_s[6];

                    // ( F_030 S_000 | P_010 S_000 )^0_{t} = y * ( F_030 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_040 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[19] = etfac[1] * AUX_INT__f_s_s_s[6] + 3 * one_over_2q * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__g_s_s_s[10];

                    // ( F_030 S_000 | P_001 S_000 )^0_{t} = z * ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[20] = etfac[2] * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__g_s_s_s[11];

                    // ( F_021 S_000 | P_100 S_000 )^0_{t} = x * ( F_021 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[21] = etfac[0] * AUX_INT__f_s_s_s[7] - p_over_q * AUX_INT__g_s_s_s[7];

                    // ( F_021 S_000 | P_010 S_000 )^0_{t} = y * ( F_021 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[22] = etfac[1] * AUX_INT__f_s_s_s[7] + 2 * one_over_2q * AUX_INT__d_s_s_s[4] - p_over_q * AUX_INT__g_s_s_s[11];

                    // ( F_021 S_000 | P_001 S_000 )^0_{t} = z * ( F_021 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[23] = etfac[2] * AUX_INT__f_s_s_s[7] + 1 * one_over_2q * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__g_s_s_s[12];

                    // ( F_012 S_000 | P_100 S_000 )^0_{t} = x * ( F_012 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[24] = etfac[0] * AUX_INT__f_s_s_s[8] - p_over_q * AUX_INT__g_s_s_s[8];

                    // ( F_012 S_000 | P_010 S_000 )^0_{t} = y * ( F_012 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[25] = etfac[1] * AUX_INT__f_s_s_s[8] + 1 * one_over_2q * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__g_s_s_s[12];

                    // ( F_012 S_000 | P_001 S_000 )^0_{t} = z * ( F_012 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[26] = etfac[2] * AUX_INT__f_s_s_s[8] + 2 * one_over_2q * AUX_INT__d_s_s_s[4] - p_over_q * AUX_INT__g_s_s_s[13];

                    // ( F_003 S_000 | P_100 S_000 )^0_{t} = x * ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[27] = etfac[0] * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__g_s_s_s[9];

                    // ( F_003 S_000 | P_010 S_000 )^0_{t} = y * ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[28] = etfac[1] * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__g_s_s_s[13];

                    // ( F_003 S_000 | P_001 S_000 )^0_{t} = z * ( F_003 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_004 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[29] = etfac[2] * AUX_INT__f_s_s_s[9] + 3 * one_over_2q * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__g_s_s_s[14];

                    // ( P_100 S_000 | P_100 S_000 )^0 = x * ( P_100 S_000 | S_000 S_000 )^0_{e} + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( D_200 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[0] = etfac[0] * AUX_INT__p_s_s_s[0] + 1 * one_over_2q * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__d_s_s_s[0];

                    // ( P_100 S_000 | P_010 S_000 )^0 = y * ( P_100 S_000 | S_000 S_000 )^0_{e} - ( D_110 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[1] = etfac[1] * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__d_s_s_s[1];

                    // ( P_100 S_000 | P_001 S_000 )^0 = z * ( P_100 S_000 | S_000 S_000 )^0_{e} - ( D_101 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[2] = etfac[2] * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__d_s_s_s[2];

                    // ( P_010 S_000 | P_100 S_000 )^0 = x * ( P_010 S_000 | S_000 S_000 )^0_{e} - ( D_110 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[3] = etfac[0] * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__d_s_s_s[1];

                    // ( P_010 S_000 | P_010 S_000 )^0 = y * ( P_010 S_000 | S_000 S_000 )^0_{e} + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( D_020 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[4] = etfac[1] * AUX_INT__p_s_s_s[1] + 1 * one_over_2q * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__d_s_s_s[3];

                    // ( P_010 S_000 | P_001 S_000 )^0 = z * ( P_010 S_000 | S_000 S_000 )^0_{e} - ( D_011 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[5] = etfac[2] * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__d_s_s_s[4];

                    // ( P_001 S_000 | P_100 S_000 )^0 = x * ( P_001 S_000 | S_000 S_000 )^0_{e} - ( D_101 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[6] = etfac[0] * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__d_s_s_s[2];

                    // ( P_001 S_000 | P_010 S_000 )^0 = y * ( P_001 S_000 | S_000 S_000 )^0_{e} - ( D_011 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[7] = etfac[1] * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__d_s_s_s[4];

                    // ( P_001 S_000 | P_001 S_000 )^0 = z * ( P_001 S_000 | S_000 S_000 )^0_{e} + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( D_002 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[8] = etfac[2] * AUX_INT__p_s_s_s[2] + 1 * one_over_2q * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__d_s_s_s[5];

                    // ( D_200 S_000 | D_200 S_000 )^0_{t} = x * ( D_200 S_000 | P_100 S_000 )^0 + ( P_100 S_000 | P_100 S_000 )^0 + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_300 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[0] = etfac[0] * AUX_INT__d_s_p_s[0] + 2 * one_over_2q * AUX_INT__p_s_p_s[0] + 1 * one_over_2q * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__f_s_p_s[0];

                    // ( D_200 S_000 | D_110 S_000 )^0_{t} = y * ( D_200 S_000 | P_100 S_000 )^0 - ( F_210 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[1] = etfac[1] * AUX_INT__d_s_p_s[0] - p_over_q * AUX_INT__f_s_p_s[3];

                    // ( D_200 S_000 | D_101 S_000 )^0_{t} = z * ( D_200 S_000 | P_100 S_000 )^0 - ( F_201 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[2] = etfac[2] * AUX_INT__d_s_p_s[0] - p_over_q * AUX_INT__f_s_p_s[6];

                    // ( D_200 S_000 | D_020 S_000 )^0_{t} = y * ( D_200 S_000 | P_010 S_000 )^0 + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[3] = etfac[1] * AUX_INT__d_s_p_s[1] + 1 * one_over_2q * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__f_s_p_s[4];

                    // ( D_200 S_000 | D_011 S_000 )^0_{t} = z * ( D_200 S_000 | P_010 S_000 )^0 - ( F_201 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[4] = etfac[2] * AUX_INT__d_s_p_s[1] - p_over_q * AUX_INT__f_s_p_s[7];

                    // ( D_200 S_000 | D_002 S_000 )^0_{t} = z * ( D_200 S_000 | P_001 S_000 )^0 + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | P_001 S_000 )^0
                    AUX_INT__d_s_d_s[5] = etfac[2] * AUX_INT__d_s_p_s[2] + 1 * one_over_2q * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__f_s_p_s[8];

                    // ( D_110 S_000 | D_200 S_000 )^0_{t} = x * ( D_110 S_000 | P_100 S_000 )^0 + ( P_010 S_000 | P_100 S_000 )^0 + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[6] = etfac[0] * AUX_INT__d_s_p_s[3] + 1 * one_over_2q * AUX_INT__p_s_p_s[3] + 1 * one_over_2q * AUX_INT__d_s_s_s[1] - p_over_q * AUX_INT__f_s_p_s[3];

                    // ( D_110 S_000 | D_110 S_000 )^0_{t} = y * ( D_110 S_000 | P_100 S_000 )^0 + ( P_100 S_000 | P_100 S_000 )^0 - ( F_120 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[7] = etfac[1] * AUX_INT__d_s_p_s[3] + 1 * one_over_2q * AUX_INT__p_s_p_s[0] - p_over_q * AUX_INT__f_s_p_s[9];

                    // ( D_110 S_000 | D_101 S_000 )^0_{t} = z * ( D_110 S_000 | P_100 S_000 )^0 - ( F_111 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[8] = etfac[2] * AUX_INT__d_s_p_s[3] - p_over_q * AUX_INT__f_s_p_s[12];

                    // ( D_110 S_000 | D_020 S_000 )^0_{t} = y * ( D_110 S_000 | P_010 S_000 )^0 + ( P_100 S_000 | P_010 S_000 )^0 + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[9] = etfac[1] * AUX_INT__d_s_p_s[4] + 1 * one_over_2q * AUX_INT__p_s_p_s[1] + 1 * one_over_2q * AUX_INT__d_s_s_s[1] - p_over_q * AUX_INT__f_s_p_s[10];

                    // ( D_110 S_000 | D_011 S_000 )^0_{t} = z * ( D_110 S_000 | P_010 S_000 )^0 - ( F_111 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[10] = etfac[2] * AUX_INT__d_s_p_s[4] - p_over_q * AUX_INT__f_s_p_s[13];

                    // ( D_110 S_000 | D_002 S_000 )^0_{t} = z * ( D_110 S_000 | P_001 S_000 )^0 + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | P_001 S_000 )^0
                    AUX_INT__d_s_d_s[11] = etfac[2] * AUX_INT__d_s_p_s[5] + 1 * one_over_2q * AUX_INT__d_s_s_s[1] - p_over_q * AUX_INT__f_s_p_s[14];

                    // ( D_101 S_000 | D_200 S_000 )^0_{t} = x * ( D_101 S_000 | P_100 S_000 )^0 + ( P_001 S_000 | P_100 S_000 )^0 + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[12] = etfac[0] * AUX_INT__d_s_p_s[6] + 1 * one_over_2q * AUX_INT__p_s_p_s[6] + 1 * one_over_2q * AUX_INT__d_s_s_s[2] - p_over_q * AUX_INT__f_s_p_s[6];

                    // ( D_101 S_000 | D_110 S_000 )^0_{t} = y * ( D_101 S_000 | P_100 S_000 )^0 - ( F_111 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[13] = etfac[1] * AUX_INT__d_s_p_s[6] - p_over_q * AUX_INT__f_s_p_s[12];

                    // ( D_101 S_000 | D_101 S_000 )^0_{t} = z * ( D_101 S_000 | P_100 S_000 )^0 + ( P_100 S_000 | P_100 S_000 )^0 - ( F_102 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[14] = etfac[2] * AUX_INT__d_s_p_s[6] + 1 * one_over_2q * AUX_INT__p_s_p_s[0] - p_over_q * AUX_INT__f_s_p_s[15];

                    // ( D_101 S_000 | D_020 S_000 )^0_{t} = y * ( D_101 S_000 | P_010 S_000 )^0 + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[15] = etfac[1] * AUX_INT__d_s_p_s[7] + 1 * one_over_2q * AUX_INT__d_s_s_s[2] - p_over_q * AUX_INT__f_s_p_s[13];

                    // ( D_101 S_000 | D_011 S_000 )^0_{t} = z * ( D_101 S_000 | P_010 S_000 )^0 + ( P_100 S_000 | P_010 S_000 )^0 - ( F_102 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[16] = etfac[2] * AUX_INT__d_s_p_s[7] + 1 * one_over_2q * AUX_INT__p_s_p_s[1] - p_over_q * AUX_INT__f_s_p_s[16];

                    // ( D_101 S_000 | D_002 S_000 )^0_{t} = z * ( D_101 S_000 | P_001 S_000 )^0 + ( P_100 S_000 | P_001 S_000 )^0 + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | P_001 S_000 )^0
                    AUX_INT__d_s_d_s[17] = etfac[2] * AUX_INT__d_s_p_s[8] + 1 * one_over_2q * AUX_INT__p_s_p_s[2] + 1 * one_over_2q * AUX_INT__d_s_s_s[2] - p_over_q * AUX_INT__f_s_p_s[17];

                    // ( D_020 S_000 | D_200 S_000 )^0_{t} = x * ( D_020 S_000 | P_100 S_000 )^0 + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[18] = etfac[0] * AUX_INT__d_s_p_s[9] + 1 * one_over_2q * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__f_s_p_s[9];

                    // ( D_020 S_000 | D_110 S_000 )^0_{t} = y * ( D_020 S_000 | P_100 S_000 )^0 + ( P_010 S_000 | P_100 S_000 )^0 - ( F_030 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[19] = etfac[1] * AUX_INT__d_s_p_s[9] + 2 * one_over_2q * AUX_INT__p_s_p_s[3] - p_over_q * AUX_INT__f_s_p_s[18];

                    // ( D_020 S_000 | D_101 S_000 )^0_{t} = z * ( D_020 S_000 | P_100 S_000 )^0 - ( F_021 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[20] = etfac[2] * AUX_INT__d_s_p_s[9] - p_over_q * AUX_INT__f_s_p_s[21];

                    // ( D_020 S_000 | D_020 S_000 )^0_{t} = y * ( D_020 S_000 | P_010 S_000 )^0 + ( P_010 S_000 | P_010 S_000 )^0 + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_030 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[21] = etfac[1] * AUX_INT__d_s_p_s[10] + 2 * one_over_2q * AUX_INT__p_s_p_s[4] + 1 * one_over_2q * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__f_s_p_s[19];

                    // ( D_020 S_000 | D_011 S_000 )^0_{t} = z * ( D_020 S_000 | P_010 S_000 )^0 - ( F_021 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[22] = etfac[2] * AUX_INT__d_s_p_s[10] - p_over_q * AUX_INT__f_s_p_s[22];

                    // ( D_020 S_000 | D_002 S_000 )^0_{t} = z * ( D_020 S_000 | P_001 S_000 )^0 + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | P_001 S_000 )^0
                    AUX_INT__d_s_d_s[23] = etfac[2] * AUX_INT__d_s_p_s[11] + 1 * one_over_2q * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__f_s_p_s[23];

                    // ( D_011 S_000 | D_200 S_000 )^0_{t} = x * ( D_011 S_000 | P_100 S_000 )^0 + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[24] = etfac[0] * AUX_INT__d_s_p_s[12] + 1 * one_over_2q * AUX_INT__d_s_s_s[4] - p_over_q * AUX_INT__f_s_p_s[12];

                    // ( D_011 S_000 | D_110 S_000 )^0_{t} = y * ( D_011 S_000 | P_100 S_000 )^0 + ( P_001 S_000 | P_100 S_000 )^0 - ( F_021 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[25] = etfac[1] * AUX_INT__d_s_p_s[12] + 1 * one_over_2q * AUX_INT__p_s_p_s[6] - p_over_q * AUX_INT__f_s_p_s[21];

                    // ( D_011 S_000 | D_101 S_000 )^0_{t} = z * ( D_011 S_000 | P_100 S_000 )^0 + ( P_010 S_000 | P_100 S_000 )^0 - ( F_012 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[26] = etfac[2] * AUX_INT__d_s_p_s[12] + 1 * one_over_2q * AUX_INT__p_s_p_s[3] - p_over_q * AUX_INT__f_s_p_s[24];

                    // ( D_011 S_000 | D_020 S_000 )^0_{t} = y * ( D_011 S_000 | P_010 S_000 )^0 + ( P_001 S_000 | P_010 S_000 )^0 + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[27] = etfac[1] * AUX_INT__d_s_p_s[13] + 1 * one_over_2q * AUX_INT__p_s_p_s[7] + 1 * one_over_2q * AUX_INT__d_s_s_s[4] - p_over_q * AUX_INT__f_s_p_s[22];

                    // ( D_011 S_000 | D_011 S_000 )^0_{t} = z * ( D_011 S_000 | P_010 S_000 )^0 + ( P_010 S_000 | P_010 S_000 )^0 - ( F_012 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[28] = etfac[2] * AUX_INT__d_s_p_s[13] + 1 * one_over_2q * AUX_INT__p_s_p_s[4] - p_over_q * AUX_INT__f_s_p_s[25];

                    // ( D_011 S_000 | D_002 S_000 )^0_{t} = z * ( D_011 S_000 | P_001 S_000 )^0 + ( P_010 S_000 | P_001 S_000 )^0 + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | P_001 S_000 )^0
                    AUX_INT__d_s_d_s[29] = etfac[2] * AUX_INT__d_s_p_s[14] + 1 * one_over_2q * AUX_INT__p_s_p_s[5] + 1 * one_over_2q * AUX_INT__d_s_s_s[4] - p_over_q * AUX_INT__f_s_p_s[26];

                    // ( D_002 S_000 | D_200 S_000 )^0_{t} = x * ( D_002 S_000 | P_100 S_000 )^0 + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[30] = etfac[0] * AUX_INT__d_s_p_s[15] + 1 * one_over_2q * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__f_s_p_s[15];

                    // ( D_002 S_000 | D_110 S_000 )^0_{t} = y * ( D_002 S_000 | P_100 S_000 )^0 - ( F_012 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[31] = etfac[1] * AUX_INT__d_s_p_s[15] - p_over_q * AUX_INT__f_s_p_s[24];

                    // ( D_002 S_000 | D_101 S_000 )^0_{t} = z * ( D_002 S_000 | P_100 S_000 )^0 + ( P_001 S_000 | P_100 S_000 )^0 - ( F_003 S_000 | P_100 S_000 )^0
                    AUX_INT__d_s_d_s[32] = etfac[2] * AUX_INT__d_s_p_s[15] + 2 * one_over_2q * AUX_INT__p_s_p_s[6] - p_over_q * AUX_INT__f_s_p_s[27];

                    // ( D_002 S_000 | D_020 S_000 )^0_{t} = y * ( D_002 S_000 | P_010 S_000 )^0 + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[33] = etfac[1] * AUX_INT__d_s_p_s[16] + 1 * one_over_2q * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__f_s_p_s[25];

                    // ( D_002 S_000 | D_011 S_000 )^0_{t} = z * ( D_002 S_000 | P_010 S_000 )^0 + ( P_001 S_000 | P_010 S_000 )^0 - ( F_003 S_000 | P_010 S_000 )^0
                    AUX_INT__d_s_d_s[34] = etfac[2] * AUX_INT__d_s_p_s[16] + 2 * one_over_2q * AUX_INT__p_s_p_s[7] - p_over_q * AUX_INT__f_s_p_s[28];

                    // ( D_002 S_000 | D_002 S_000 )^0_{t} = z * ( D_002 S_000 | P_001 S_000 )^0 + ( P_001 S_000 | P_001 S_000 )^0 + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_003 S_000 | P_001 S_000 )^0
                    AUX_INT__d_s_d_s[35] = etfac[2] * AUX_INT__d_s_p_s[17] + 2 * one_over_2q * AUX_INT__p_s_p_s[8] + 1 * one_over_2q * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__f_s_p_s[29];

                    // ( G_400 S_000 | P_100 S_000 )^0_{t} = x * ( G_400 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_500 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[0] = etfac[0] * AUX_INT__g_s_s_s[0] + 4 * one_over_2q * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__h_s_s_s[0];

                    // ( G_400 S_000 | P_010 S_000 )^0_{t} = y * ( G_400 S_000 | S_000 S_000 )^0_{e} - ( H_410 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[1] = etfac[1] * AUX_INT__g_s_s_s[0] - p_over_q * AUX_INT__h_s_s_s[1];

                    // ( G_400 S_000 | P_001 S_000 )^0_{t} = z * ( G_400 S_000 | S_000 S_000 )^0_{e} - ( H_401 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[2] = etfac[2] * AUX_INT__g_s_s_s[0] - p_over_q * AUX_INT__h_s_s_s[2];

                    // ( G_310 S_000 | P_100 S_000 )^0_{t} = x * ( G_310 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_410 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[3] = etfac[0] * AUX_INT__g_s_s_s[1] + 3 * one_over_2q * AUX_INT__f_s_s_s[1] - p_over_q * AUX_INT__h_s_s_s[1];

                    // ( G_310 S_000 | P_010 S_000 )^0_{t} = y * ( G_310 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_320 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[4] = etfac[1] * AUX_INT__g_s_s_s[1] + 1 * one_over_2q * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__h_s_s_s[3];

                    // ( G_310 S_000 | P_001 S_000 )^0_{t} = z * ( G_310 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[5] = etfac[2] * AUX_INT__g_s_s_s[1] - p_over_q * AUX_INT__h_s_s_s[4];

                    // ( G_301 S_000 | P_100 S_000 )^0_{t} = x * ( G_301 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_401 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[6] = etfac[0] * AUX_INT__g_s_s_s[2] + 3 * one_over_2q * AUX_INT__f_s_s_s[2] - p_over_q * AUX_INT__h_s_s_s[2];

                    // ( G_301 S_000 | P_010 S_000 )^0_{t} = y * ( G_301 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[7] = etfac[1] * AUX_INT__g_s_s_s[2] - p_over_q * AUX_INT__h_s_s_s[4];

                    // ( G_301 S_000 | P_001 S_000 )^0_{t} = z * ( G_301 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_302 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[8] = etfac[2] * AUX_INT__g_s_s_s[2] + 1 * one_over_2q * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__h_s_s_s[5];

                    // ( G_220 S_000 | P_100 S_000 )^0_{t} = x * ( G_220 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_320 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[9] = etfac[0] * AUX_INT__g_s_s_s[3] + 2 * one_over_2q * AUX_INT__f_s_s_s[3] - p_over_q * AUX_INT__h_s_s_s[3];

                    // ( G_220 S_000 | P_010 S_000 )^0_{t} = y * ( G_220 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_230 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[10] = etfac[1] * AUX_INT__g_s_s_s[3] + 2 * one_over_2q * AUX_INT__f_s_s_s[1] - p_over_q * AUX_INT__h_s_s_s[6];

                    // ( G_220 S_000 | P_001 S_000 )^0_{t} = z * ( G_220 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[11] = etfac[2] * AUX_INT__g_s_s_s[3] - p_over_q * AUX_INT__h_s_s_s[7];

                    // ( G_211 S_000 | P_100 S_000 )^0_{t} = x * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[12] = etfac[0] * AUX_INT__g_s_s_s[4] + 2 * one_over_2q * AUX_INT__f_s_s_s[4] - p_over_q * AUX_INT__h_s_s_s[4];

                    // ( G_211 S_000 | P_010 S_000 )^0_{t} = y * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[13] = etfac[1] * AUX_INT__g_s_s_s[4] + 1 * one_over_2q * AUX_INT__f_s_s_s[2] - p_over_q * AUX_INT__h_s_s_s[7];

                    // ( G_211 S_000 | P_001 S_000 )^0_{t} = z * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[14] = etfac[2] * AUX_INT__g_s_s_s[4] + 1 * one_over_2q * AUX_INT__f_s_s_s[1] - p_over_q * AUX_INT__h_s_s_s[8];

                    // ( G_202 S_000 | P_100 S_000 )^0_{t} = x * ( G_202 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_302 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[15] = etfac[0] * AUX_INT__g_s_s_s[5] + 2 * one_over_2q * AUX_INT__f_s_s_s[5] - p_over_q * AUX_INT__h_s_s_s[5];

                    // ( G_202 S_000 | P_010 S_000 )^0_{t} = y * ( G_202 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[16] = etfac[1] * AUX_INT__g_s_s_s[5] - p_over_q * AUX_INT__h_s_s_s[8];

                    // ( G_202 S_000 | P_001 S_000 )^0_{t} = z * ( G_202 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_203 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[17] = etfac[2] * AUX_INT__g_s_s_s[5] + 2 * one_over_2q * AUX_INT__f_s_s_s[2] - p_over_q * AUX_INT__h_s_s_s[9];

                    // ( G_130 S_000 | P_100 S_000 )^0_{t} = x * ( G_130 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_230 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[18] = etfac[0] * AUX_INT__g_s_s_s[6] + 1 * one_over_2q * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__h_s_s_s[6];

                    // ( G_130 S_000 | P_010 S_000 )^0_{t} = y * ( G_130 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_140 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[19] = etfac[1] * AUX_INT__g_s_s_s[6] + 3 * one_over_2q * AUX_INT__f_s_s_s[3] - p_over_q * AUX_INT__h_s_s_s[10];

                    // ( G_130 S_000 | P_001 S_000 )^0_{t} = z * ( G_130 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[20] = etfac[2] * AUX_INT__g_s_s_s[6] - p_over_q * AUX_INT__h_s_s_s[11];

                    // ( G_121 S_000 | P_100 S_000 )^0_{t} = x * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[21] = etfac[0] * AUX_INT__g_s_s_s[7] + 1 * one_over_2q * AUX_INT__f_s_s_s[7] - p_over_q * AUX_INT__h_s_s_s[7];

                    // ( G_121 S_000 | P_010 S_000 )^0_{t} = y * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[22] = etfac[1] * AUX_INT__g_s_s_s[7] + 2 * one_over_2q * AUX_INT__f_s_s_s[4] - p_over_q * AUX_INT__h_s_s_s[11];

                    // ( G_121 S_000 | P_001 S_000 )^0_{t} = z * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[23] = etfac[2] * AUX_INT__g_s_s_s[7] + 1 * one_over_2q * AUX_INT__f_s_s_s[3] - p_over_q * AUX_INT__h_s_s_s[12];

                    // ( G_112 S_000 | P_100 S_000 )^0_{t} = x * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[24] = etfac[0] * AUX_INT__g_s_s_s[8] + 1 * one_over_2q * AUX_INT__f_s_s_s[8] - p_over_q * AUX_INT__h_s_s_s[8];

                    // ( G_112 S_000 | P_010 S_000 )^0_{t} = y * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[25] = etfac[1] * AUX_INT__g_s_s_s[8] + 1 * one_over_2q * AUX_INT__f_s_s_s[5] - p_over_q * AUX_INT__h_s_s_s[12];

                    // ( G_112 S_000 | P_001 S_000 )^0_{t} = z * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[26] = etfac[2] * AUX_INT__g_s_s_s[8] + 2 * one_over_2q * AUX_INT__f_s_s_s[4] - p_over_q * AUX_INT__h_s_s_s[13];

                    // ( G_103 S_000 | P_100 S_000 )^0_{t} = x * ( G_103 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_203 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[27] = etfac[0] * AUX_INT__g_s_s_s[9] + 1 * one_over_2q * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__h_s_s_s[9];

                    // ( G_103 S_000 | P_010 S_000 )^0_{t} = y * ( G_103 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[28] = etfac[1] * AUX_INT__g_s_s_s[9] - p_over_q * AUX_INT__h_s_s_s[13];

                    // ( G_103 S_000 | P_001 S_000 )^0_{t} = z * ( G_103 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_104 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[29] = etfac[2] * AUX_INT__g_s_s_s[9] + 3 * one_over_2q * AUX_INT__f_s_s_s[5] - p_over_q * AUX_INT__h_s_s_s[14];

                    // ( G_040 S_000 | P_100 S_000 )^0_{t} = x * ( G_040 S_000 | S_000 S_000 )^0_{e} - ( H_140 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[30] = etfac[0] * AUX_INT__g_s_s_s[10] - p_over_q * AUX_INT__h_s_s_s[10];

                    // ( G_040 S_000 | P_010 S_000 )^0_{t} = y * ( G_040 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_050 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[31] = etfac[1] * AUX_INT__g_s_s_s[10] + 4 * one_over_2q * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__h_s_s_s[15];

                    // ( G_040 S_000 | P_001 S_000 )^0_{t} = z * ( G_040 S_000 | S_000 S_000 )^0_{e} - ( H_041 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[32] = etfac[2] * AUX_INT__g_s_s_s[10] - p_over_q * AUX_INT__h_s_s_s[16];

                    // ( G_031 S_000 | P_100 S_000 )^0_{t} = x * ( G_031 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[33] = etfac[0] * AUX_INT__g_s_s_s[11] - p_over_q * AUX_INT__h_s_s_s[11];

                    // ( G_031 S_000 | P_010 S_000 )^0_{t} = y * ( G_031 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_041 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[34] = etfac[1] * AUX_INT__g_s_s_s[11] + 3 * one_over_2q * AUX_INT__f_s_s_s[7] - p_over_q * AUX_INT__h_s_s_s[16];

                    // ( G_031 S_000 | P_001 S_000 )^0_{t} = z * ( G_031 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_032 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[35] = etfac[2] * AUX_INT__g_s_s_s[11] + 1 * one_over_2q * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__h_s_s_s[17];

                    // ( G_022 S_000 | P_100 S_000 )^0_{t} = x * ( G_022 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[36] = etfac[0] * AUX_INT__g_s_s_s[12] - p_over_q * AUX_INT__h_s_s_s[12];

                    // ( G_022 S_000 | P_010 S_000 )^0_{t} = y * ( G_022 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_032 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[37] = etfac[1] * AUX_INT__g_s_s_s[12] + 2 * one_over_2q * AUX_INT__f_s_s_s[8] - p_over_q * AUX_INT__h_s_s_s[17];

                    // ( G_022 S_000 | P_001 S_000 )^0_{t} = z * ( G_022 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_023 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[38] = etfac[2] * AUX_INT__g_s_s_s[12] + 2 * one_over_2q * AUX_INT__f_s_s_s[7] - p_over_q * AUX_INT__h_s_s_s[18];

                    // ( G_013 S_000 | P_100 S_000 )^0_{t} = x * ( G_013 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[39] = etfac[0] * AUX_INT__g_s_s_s[13] - p_over_q * AUX_INT__h_s_s_s[13];

                    // ( G_013 S_000 | P_010 S_000 )^0_{t} = y * ( G_013 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_023 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[40] = etfac[1] * AUX_INT__g_s_s_s[13] + 1 * one_over_2q * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__h_s_s_s[18];

                    // ( G_013 S_000 | P_001 S_000 )^0_{t} = z * ( G_013 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_014 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[41] = etfac[2] * AUX_INT__g_s_s_s[13] + 3 * one_over_2q * AUX_INT__f_s_s_s[8] - p_over_q * AUX_INT__h_s_s_s[19];

                    // ( G_004 S_000 | P_100 S_000 )^0_{t} = x * ( G_004 S_000 | S_000 S_000 )^0_{e} - ( H_104 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[42] = etfac[0] * AUX_INT__g_s_s_s[14] - p_over_q * AUX_INT__h_s_s_s[14];

                    // ( G_004 S_000 | P_010 S_000 )^0_{t} = y * ( G_004 S_000 | S_000 S_000 )^0_{e} - ( H_014 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[43] = etfac[1] * AUX_INT__g_s_s_s[14] - p_over_q * AUX_INT__h_s_s_s[19];

                    // ( G_004 S_000 | P_001 S_000 )^0_{t} = z * ( G_004 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_005 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[44] = etfac[2] * AUX_INT__g_s_s_s[14] + 4 * one_over_2q * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__h_s_s_s[20];

                    // ( F_300 S_000 | D_200 S_000 )^0_{t} = x * ( F_300 S_000 | P_100 S_000 )^0 + ( D_200 S_000 | P_100 S_000 )^0 + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_400 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[0] = etfac[0] * AUX_INT__f_s_p_s[0] + 3 * one_over_2q * AUX_INT__d_s_p_s[0] + 1 * one_over_2q * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__g_s_p_s[0];

                    // ( F_300 S_000 | D_110 S_000 )^0_{t} = y * ( F_300 S_000 | P_100 S_000 )^0 - ( G_310 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[1] = etfac[1] * AUX_INT__f_s_p_s[0] - p_over_q * AUX_INT__g_s_p_s[3];

                    // ( F_300 S_000 | D_101 S_000 )^0_{t} = z * ( F_300 S_000 | P_100 S_000 )^0 - ( G_301 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[2] = etfac[2] * AUX_INT__f_s_p_s[0] - p_over_q * AUX_INT__g_s_p_s[6];

                    // ( F_300 S_000 | D_020 S_000 )^0_{t} = y * ( F_300 S_000 | P_010 S_000 )^0 + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[3] = etfac[1] * AUX_INT__f_s_p_s[1] + 1 * one_over_2q * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__g_s_p_s[4];

                    // ( F_300 S_000 | D_011 S_000 )^0_{t} = z * ( F_300 S_000 | P_010 S_000 )^0 - ( G_301 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[4] = etfac[2] * AUX_INT__f_s_p_s[1] - p_over_q * AUX_INT__g_s_p_s[7];

                    // ( F_300 S_000 | D_002 S_000 )^0_{t} = z * ( F_300 S_000 | P_001 S_000 )^0 + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | P_001 S_000 )^0
                    AUX_INT__f_s_d_s[5] = etfac[2] * AUX_INT__f_s_p_s[2] + 1 * one_over_2q * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__g_s_p_s[8];

                    // ( F_210 S_000 | D_200 S_000 )^0_{t} = x * ( F_210 S_000 | P_100 S_000 )^0 + ( D_110 S_000 | P_100 S_000 )^0 + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[6] = etfac[0] * AUX_INT__f_s_p_s[3] + 2 * one_over_2q * AUX_INT__d_s_p_s[3] + 1 * one_over_2q * AUX_INT__f_s_s_s[1] - p_over_q * AUX_INT__g_s_p_s[3];

                    // ( F_210 S_000 | D_110 S_000 )^0_{t} = y * ( F_210 S_000 | P_100 S_000 )^0 + ( D_200 S_000 | P_100 S_000 )^0 - ( G_220 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[7] = etfac[1] * AUX_INT__f_s_p_s[3] + 1 * one_over_2q * AUX_INT__d_s_p_s[0] - p_over_q * AUX_INT__g_s_p_s[9];

                    // ( F_210 S_000 | D_101 S_000 )^0_{t} = z * ( F_210 S_000 | P_100 S_000 )^0 - ( G_211 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[8] = etfac[2] * AUX_INT__f_s_p_s[3] - p_over_q * AUX_INT__g_s_p_s[12];

                    // ( F_210 S_000 | D_020 S_000 )^0_{t} = y * ( F_210 S_000 | P_010 S_000 )^0 + ( D_200 S_000 | P_010 S_000 )^0 + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[9] = etfac[1] * AUX_INT__f_s_p_s[4] + 1 * one_over_2q * AUX_INT__d_s_p_s[1] + 1 * one_over_2q * AUX_INT__f_s_s_s[1] - p_over_q * AUX_INT__g_s_p_s[10];

                    // ( F_210 S_000 | D_011 S_000 )^0_{t} = z * ( F_210 S_000 | P_010 S_000 )^0 - ( G_211 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[10] = etfac[2] * AUX_INT__f_s_p_s[4] - p_over_q * AUX_INT__g_s_p_s[13];

                    // ( F_210 S_000 | D_002 S_000 )^0_{t} = z * ( F_210 S_000 | P_001 S_000 )^0 + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | P_001 S_000 )^0
                    AUX_INT__f_s_d_s[11] = etfac[2] * AUX_INT__f_s_p_s[5] + 1 * one_over_2q * AUX_INT__f_s_s_s[1] - p_over_q * AUX_INT__g_s_p_s[14];

                    // ( F_201 S_000 | D_200 S_000 )^0_{t} = x * ( F_201 S_000 | P_100 S_000 )^0 + ( D_101 S_000 | P_100 S_000 )^0 + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[12] = etfac[0] * AUX_INT__f_s_p_s[6] + 2 * one_over_2q * AUX_INT__d_s_p_s[6] + 1 * one_over_2q * AUX_INT__f_s_s_s[2] - p_over_q * AUX_INT__g_s_p_s[6];

                    // ( F_201 S_000 | D_110 S_000 )^0_{t} = y * ( F_201 S_000 | P_100 S_000 )^0 - ( G_211 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[13] = etfac[1] * AUX_INT__f_s_p_s[6] - p_over_q * AUX_INT__g_s_p_s[12];

                    // ( F_201 S_000 | D_101 S_000 )^0_{t} = z * ( F_201 S_000 | P_100 S_000 )^0 + ( D_200 S_000 | P_100 S_000 )^0 - ( G_202 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[14] = etfac[2] * AUX_INT__f_s_p_s[6] + 1 * one_over_2q * AUX_INT__d_s_p_s[0] - p_over_q * AUX_INT__g_s_p_s[15];

                    // ( F_201 S_000 | D_020 S_000 )^0_{t} = y * ( F_201 S_000 | P_010 S_000 )^0 + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[15] = etfac[1] * AUX_INT__f_s_p_s[7] + 1 * one_over_2q * AUX_INT__f_s_s_s[2] - p_over_q * AUX_INT__g_s_p_s[13];

                    // ( F_201 S_000 | D_011 S_000 )^0_{t} = z * ( F_201 S_000 | P_010 S_000 )^0 + ( D_200 S_000 | P_010 S_000 )^0 - ( G_202 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[16] = etfac[2] * AUX_INT__f_s_p_s[7] + 1 * one_over_2q * AUX_INT__d_s_p_s[1] - p_over_q * AUX_INT__g_s_p_s[16];

                    // ( F_201 S_000 | D_002 S_000 )^0_{t} = z * ( F_201 S_000 | P_001 S_000 )^0 + ( D_200 S_000 | P_001 S_000 )^0 + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | P_001 S_000 )^0
                    AUX_INT__f_s_d_s[17] = etfac[2] * AUX_INT__f_s_p_s[8] + 1 * one_over_2q * AUX_INT__d_s_p_s[2] + 1 * one_over_2q * AUX_INT__f_s_s_s[2] - p_over_q * AUX_INT__g_s_p_s[17];

                    // ( F_120 S_000 | D_200 S_000 )^0_{t} = x * ( F_120 S_000 | P_100 S_000 )^0 + ( D_020 S_000 | P_100 S_000 )^0 + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[18] = etfac[0] * AUX_INT__f_s_p_s[9] + 1 * one_over_2q * AUX_INT__d_s_p_s[9] + 1 * one_over_2q * AUX_INT__f_s_s_s[3] - p_over_q * AUX_INT__g_s_p_s[9];

                    // ( F_120 S_000 | D_110 S_000 )^0_{t} = y * ( F_120 S_000 | P_100 S_000 )^0 + ( D_110 S_000 | P_100 S_000 )^0 - ( G_130 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[19] = etfac[1] * AUX_INT__f_s_p_s[9] + 2 * one_over_2q * AUX_INT__d_s_p_s[3] - p_over_q * AUX_INT__g_s_p_s[18];

                    // ( F_120 S_000 | D_101 S_000 )^0_{t} = z * ( F_120 S_000 | P_100 S_000 )^0 - ( G_121 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[20] = etfac[2] * AUX_INT__f_s_p_s[9] - p_over_q * AUX_INT__g_s_p_s[21];

                    // ( F_120 S_000 | D_020 S_000 )^0_{t} = y * ( F_120 S_000 | P_010 S_000 )^0 + ( D_110 S_000 | P_010 S_000 )^0 + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[21] = etfac[1] * AUX_INT__f_s_p_s[10] + 2 * one_over_2q * AUX_INT__d_s_p_s[4] + 1 * one_over_2q * AUX_INT__f_s_s_s[3] - p_over_q * AUX_INT__g_s_p_s[19];

                    // ( F_120 S_000 | D_011 S_000 )^0_{t} = z * ( F_120 S_000 | P_010 S_000 )^0 - ( G_121 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[22] = etfac[2] * AUX_INT__f_s_p_s[10] - p_over_q * AUX_INT__g_s_p_s[22];

                    // ( F_120 S_000 | D_002 S_000 )^0_{t} = z * ( F_120 S_000 | P_001 S_000 )^0 + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | P_001 S_000 )^0
                    AUX_INT__f_s_d_s[23] = etfac[2] * AUX_INT__f_s_p_s[11] + 1 * one_over_2q * AUX_INT__f_s_s_s[3] - p_over_q * AUX_INT__g_s_p_s[23];

                    // ( F_111 S_000 | D_200 S_000 )^0_{t} = x * ( F_111 S_000 | P_100 S_000 )^0 + ( D_011 S_000 | P_100 S_000 )^0 + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[24] = etfac[0] * AUX_INT__f_s_p_s[12] + 1 * one_over_2q * AUX_INT__d_s_p_s[12] + 1 * one_over_2q * AUX_INT__f_s_s_s[4] - p_over_q * AUX_INT__g_s_p_s[12];

                    // ( F_111 S_000 | D_110 S_000 )^0_{t} = y * ( F_111 S_000 | P_100 S_000 )^0 + ( D_101 S_000 | P_100 S_000 )^0 - ( G_121 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[25] = etfac[1] * AUX_INT__f_s_p_s[12] + 1 * one_over_2q * AUX_INT__d_s_p_s[6] - p_over_q * AUX_INT__g_s_p_s[21];

                    // ( F_111 S_000 | D_101 S_000 )^0_{t} = z * ( F_111 S_000 | P_100 S_000 )^0 + ( D_110 S_000 | P_100 S_000 )^0 - ( G_112 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[26] = etfac[2] * AUX_INT__f_s_p_s[12] + 1 * one_over_2q * AUX_INT__d_s_p_s[3] - p_over_q * AUX_INT__g_s_p_s[24];

                    // ( F_111 S_000 | D_020 S_000 )^0_{t} = y * ( F_111 S_000 | P_010 S_000 )^0 + ( D_101 S_000 | P_010 S_000 )^0 + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[27] = etfac[1] * AUX_INT__f_s_p_s[13] + 1 * one_over_2q * AUX_INT__d_s_p_s[7] + 1 * one_over_2q * AUX_INT__f_s_s_s[4] - p_over_q * AUX_INT__g_s_p_s[22];

                    // ( F_111 S_000 | D_011 S_000 )^0_{t} = z * ( F_111 S_000 | P_010 S_000 )^0 + ( D_110 S_000 | P_010 S_000 )^0 - ( G_112 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[28] = etfac[2] * AUX_INT__f_s_p_s[13] + 1 * one_over_2q * AUX_INT__d_s_p_s[4] - p_over_q * AUX_INT__g_s_p_s[25];

                    // ( F_111 S_000 | D_002 S_000 )^0_{t} = z * ( F_111 S_000 | P_001 S_000 )^0 + ( D_110 S_000 | P_001 S_000 )^0 + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | P_001 S_000 )^0
                    AUX_INT__f_s_d_s[29] = etfac[2] * AUX_INT__f_s_p_s[14] + 1 * one_over_2q * AUX_INT__d_s_p_s[5] + 1 * one_over_2q * AUX_INT__f_s_s_s[4] - p_over_q * AUX_INT__g_s_p_s[26];

                    // ( F_102 S_000 | D_200 S_000 )^0_{t} = x * ( F_102 S_000 | P_100 S_000 )^0 + ( D_002 S_000 | P_100 S_000 )^0 + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[30] = etfac[0] * AUX_INT__f_s_p_s[15] + 1 * one_over_2q * AUX_INT__d_s_p_s[15] + 1 * one_over_2q * AUX_INT__f_s_s_s[5] - p_over_q * AUX_INT__g_s_p_s[15];

                    // ( F_102 S_000 | D_110 S_000 )^0_{t} = y * ( F_102 S_000 | P_100 S_000 )^0 - ( G_112 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[31] = etfac[1] * AUX_INT__f_s_p_s[15] - p_over_q * AUX_INT__g_s_p_s[24];

                    // ( F_102 S_000 | D_101 S_000 )^0_{t} = z * ( F_102 S_000 | P_100 S_000 )^0 + ( D_101 S_000 | P_100 S_000 )^0 - ( G_103 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[32] = etfac[2] * AUX_INT__f_s_p_s[15] + 2 * one_over_2q * AUX_INT__d_s_p_s[6] - p_over_q * AUX_INT__g_s_p_s[27];

                    // ( F_102 S_000 | D_020 S_000 )^0_{t} = y * ( F_102 S_000 | P_010 S_000 )^0 + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[33] = etfac[1] * AUX_INT__f_s_p_s[16] + 1 * one_over_2q * AUX_INT__f_s_s_s[5] - p_over_q * AUX_INT__g_s_p_s[25];

                    // ( F_102 S_000 | D_011 S_000 )^0_{t} = z * ( F_102 S_000 | P_010 S_000 )^0 + ( D_101 S_000 | P_010 S_000 )^0 - ( G_103 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[34] = etfac[2] * AUX_INT__f_s_p_s[16] + 2 * one_over_2q * AUX_INT__d_s_p_s[7] - p_over_q * AUX_INT__g_s_p_s[28];

                    // ( F_102 S_000 | D_002 S_000 )^0_{t} = z * ( F_102 S_000 | P_001 S_000 )^0 + ( D_101 S_000 | P_001 S_000 )^0 + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | P_001 S_000 )^0
                    AUX_INT__f_s_d_s[35] = etfac[2] * AUX_INT__f_s_p_s[17] + 2 * one_over_2q * AUX_INT__d_s_p_s[8] + 1 * one_over_2q * AUX_INT__f_s_s_s[5] - p_over_q * AUX_INT__g_s_p_s[29];

                    // ( F_030 S_000 | D_200 S_000 )^0_{t} = x * ( F_030 S_000 | P_100 S_000 )^0 + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[36] = etfac[0] * AUX_INT__f_s_p_s[18] + 1 * one_over_2q * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__g_s_p_s[18];

                    // ( F_030 S_000 | D_110 S_000 )^0_{t} = y * ( F_030 S_000 | P_100 S_000 )^0 + ( D_020 S_000 | P_100 S_000 )^0 - ( G_040 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[37] = etfac[1] * AUX_INT__f_s_p_s[18] + 3 * one_over_2q * AUX_INT__d_s_p_s[9] - p_over_q * AUX_INT__g_s_p_s[30];

                    // ( F_030 S_000 | D_101 S_000 )^0_{t} = z * ( F_030 S_000 | P_100 S_000 )^0 - ( G_031 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[38] = etfac[2] * AUX_INT__f_s_p_s[18] - p_over_q * AUX_INT__g_s_p_s[33];

                    // ( F_030 S_000 | D_020 S_000 )^0_{t} = y * ( F_030 S_000 | P_010 S_000 )^0 + ( D_020 S_000 | P_010 S_000 )^0 + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_040 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[39] = etfac[1] * AUX_INT__f_s_p_s[19] + 3 * one_over_2q * AUX_INT__d_s_p_s[10] + 1 * one_over_2q * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__g_s_p_s[31];

                    // ( F_030 S_000 | D_011 S_000 )^0_{t} = z * ( F_030 S_000 | P_010 S_000 )^0 - ( G_031 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[40] = etfac[2] * AUX_INT__f_s_p_s[19] - p_over_q * AUX_INT__g_s_p_s[34];

                    // ( F_030 S_000 | D_002 S_000 )^0_{t} = z * ( F_030 S_000 | P_001 S_000 )^0 + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | P_001 S_000 )^0
                    AUX_INT__f_s_d_s[41] = etfac[2] * AUX_INT__f_s_p_s[20] + 1 * one_over_2q * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__g_s_p_s[35];

                    // ( F_021 S_000 | D_200 S_000 )^0_{t} = x * ( F_021 S_000 | P_100 S_000 )^0 + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[42] = etfac[0] * AUX_INT__f_s_p_s[21] + 1 * one_over_2q * AUX_INT__f_s_s_s[7] - p_over_q * AUX_INT__g_s_p_s[21];

                    // ( F_021 S_000 | D_110 S_000 )^0_{t} = y * ( F_021 S_000 | P_100 S_000 )^0 + ( D_011 S_000 | P_100 S_000 )^0 - ( G_031 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[43] = etfac[1] * AUX_INT__f_s_p_s[21] + 2 * one_over_2q * AUX_INT__d_s_p_s[12] - p_over_q * AUX_INT__g_s_p_s[33];

                    // ( F_021 S_000 | D_101 S_000 )^0_{t} = z * ( F_021 S_000 | P_100 S_000 )^0 + ( D_020 S_000 | P_100 S_000 )^0 - ( G_022 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[44] = etfac[2] * AUX_INT__f_s_p_s[21] + 1 * one_over_2q * AUX_INT__d_s_p_s[9] - p_over_q * AUX_INT__g_s_p_s[36];

                    // ( F_021 S_000 | D_020 S_000 )^0_{t} = y * ( F_021 S_000 | P_010 S_000 )^0 + ( D_011 S_000 | P_010 S_000 )^0 + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[45] = etfac[1] * AUX_INT__f_s_p_s[22] + 2 * one_over_2q * AUX_INT__d_s_p_s[13] + 1 * one_over_2q * AUX_INT__f_s_s_s[7] - p_over_q * AUX_INT__g_s_p_s[34];

                    // ( F_021 S_000 | D_011 S_000 )^0_{t} = z * ( F_021 S_000 | P_010 S_000 )^0 + ( D_020 S_000 | P_010 S_000 )^0 - ( G_022 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[46] = etfac[2] * AUX_INT__f_s_p_s[22] + 1 * one_over_2q * AUX_INT__d_s_p_s[10] - p_over_q * AUX_INT__g_s_p_s[37];

                    // ( F_021 S_000 | D_002 S_000 )^0_{t} = z * ( F_021 S_000 | P_001 S_000 )^0 + ( D_020 S_000 | P_001 S_000 )^0 + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | P_001 S_000 )^0
                    AUX_INT__f_s_d_s[47] = etfac[2] * AUX_INT__f_s_p_s[23] + 1 * one_over_2q * AUX_INT__d_s_p_s[11] + 1 * one_over_2q * AUX_INT__f_s_s_s[7] - p_over_q * AUX_INT__g_s_p_s[38];

                    // ( F_012 S_000 | D_200 S_000 )^0_{t} = x * ( F_012 S_000 | P_100 S_000 )^0 + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[48] = etfac[0] * AUX_INT__f_s_p_s[24] + 1 * one_over_2q * AUX_INT__f_s_s_s[8] - p_over_q * AUX_INT__g_s_p_s[24];

                    // ( F_012 S_000 | D_110 S_000 )^0_{t} = y * ( F_012 S_000 | P_100 S_000 )^0 + ( D_002 S_000 | P_100 S_000 )^0 - ( G_022 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[49] = etfac[1] * AUX_INT__f_s_p_s[24] + 1 * one_over_2q * AUX_INT__d_s_p_s[15] - p_over_q * AUX_INT__g_s_p_s[36];

                    // ( F_012 S_000 | D_101 S_000 )^0_{t} = z * ( F_012 S_000 | P_100 S_000 )^0 + ( D_011 S_000 | P_100 S_000 )^0 - ( G_013 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[50] = etfac[2] * AUX_INT__f_s_p_s[24] + 2 * one_over_2q * AUX_INT__d_s_p_s[12] - p_over_q * AUX_INT__g_s_p_s[39];

                    // ( F_012 S_000 | D_020 S_000 )^0_{t} = y * ( F_012 S_000 | P_010 S_000 )^0 + ( D_002 S_000 | P_010 S_000 )^0 + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[51] = etfac[1] * AUX_INT__f_s_p_s[25] + 1 * one_over_2q * AUX_INT__d_s_p_s[16] + 1 * one_over_2q * AUX_INT__f_s_s_s[8] - p_over_q * AUX_INT__g_s_p_s[37];

                    // ( F_012 S_000 | D_011 S_000 )^0_{t} = z * ( F_012 S_000 | P_010 S_000 )^0 + ( D_011 S_000 | P_010 S_000 )^0 - ( G_013 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[52] = etfac[2] * AUX_INT__f_s_p_s[25] + 2 * one_over_2q * AUX_INT__d_s_p_s[13] - p_over_q * AUX_INT__g_s_p_s[40];

                    // ( F_012 S_000 | D_002 S_000 )^0_{t} = z * ( F_012 S_000 | P_001 S_000 )^0 + ( D_011 S_000 | P_001 S_000 )^0 + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | P_001 S_000 )^0
                    AUX_INT__f_s_d_s[53] = etfac[2] * AUX_INT__f_s_p_s[26] + 2 * one_over_2q * AUX_INT__d_s_p_s[14] + 1 * one_over_2q * AUX_INT__f_s_s_s[8] - p_over_q * AUX_INT__g_s_p_s[41];

                    // ( F_003 S_000 | D_200 S_000 )^0_{t} = x * ( F_003 S_000 | P_100 S_000 )^0 + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[54] = etfac[0] * AUX_INT__f_s_p_s[27] + 1 * one_over_2q * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__g_s_p_s[27];

                    // ( F_003 S_000 | D_110 S_000 )^0_{t} = y * ( F_003 S_000 | P_100 S_000 )^0 - ( G_013 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[55] = etfac[1] * AUX_INT__f_s_p_s[27] - p_over_q * AUX_INT__g_s_p_s[39];

                    // ( F_003 S_000 | D_101 S_000 )^0_{t} = z * ( F_003 S_000 | P_100 S_000 )^0 + ( D_002 S_000 | P_100 S_000 )^0 - ( G_004 S_000 | P_100 S_000 )^0
                    AUX_INT__f_s_d_s[56] = etfac[2] * AUX_INT__f_s_p_s[27] + 3 * one_over_2q * AUX_INT__d_s_p_s[15] - p_over_q * AUX_INT__g_s_p_s[42];

                    // ( F_003 S_000 | D_020 S_000 )^0_{t} = y * ( F_003 S_000 | P_010 S_000 )^0 + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[57] = etfac[1] * AUX_INT__f_s_p_s[28] + 1 * one_over_2q * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__g_s_p_s[40];

                    // ( F_003 S_000 | D_011 S_000 )^0_{t} = z * ( F_003 S_000 | P_010 S_000 )^0 + ( D_002 S_000 | P_010 S_000 )^0 - ( G_004 S_000 | P_010 S_000 )^0
                    AUX_INT__f_s_d_s[58] = etfac[2] * AUX_INT__f_s_p_s[28] + 3 * one_over_2q * AUX_INT__d_s_p_s[16] - p_over_q * AUX_INT__g_s_p_s[43];

                    // ( F_003 S_000 | D_002 S_000 )^0_{t} = z * ( F_003 S_000 | P_001 S_000 )^0 + ( D_002 S_000 | P_001 S_000 )^0 + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_004 S_000 | P_001 S_000 )^0
                    AUX_INT__f_s_d_s[59] = etfac[2] * AUX_INT__f_s_p_s[29] + 3 * one_over_2q * AUX_INT__d_s_p_s[17] + 1 * one_over_2q * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__g_s_p_s[44];

                    // ( H_500 S_000 | P_100 S_000 )^0 = x * ( H_500 S_000 | S_000 S_000 )^0_{e} + ( G_400 S_000 | S_000 S_000 )^0_{e} - ( I_600 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[0] = etfac[0] * AUX_INT__h_s_s_s[0] + 5 * one_over_2q * AUX_INT__g_s_s_s[0] - p_over_q * AUX_INT__i_s_s_s[0];

                    // ( H_410 S_000 | P_100 S_000 )^0 = x * ( H_410 S_000 | S_000 S_000 )^0_{e} + ( G_310 S_000 | S_000 S_000 )^0_{e} - ( I_510 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[3] = etfac[0] * AUX_INT__h_s_s_s[1] + 4 * one_over_2q * AUX_INT__g_s_s_s[1] - p_over_q * AUX_INT__i_s_s_s[1];

                    // ( H_410 S_000 | P_010 S_000 )^0 = y * ( H_410 S_000 | S_000 S_000 )^0_{e} + ( G_400 S_000 | S_000 S_000 )^0_{e} - ( I_420 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[4] = etfac[1] * AUX_INT__h_s_s_s[1] + 1 * one_over_2q * AUX_INT__g_s_s_s[0] - p_over_q * AUX_INT__i_s_s_s[3];

                    // ( H_401 S_000 | P_100 S_000 )^0 = x * ( H_401 S_000 | S_000 S_000 )^0_{e} + ( G_301 S_000 | S_000 S_000 )^0_{e} - ( I_501 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[6] = etfac[0] * AUX_INT__h_s_s_s[2] + 4 * one_over_2q * AUX_INT__g_s_s_s[2] - p_over_q * AUX_INT__i_s_s_s[2];

                    // ( H_401 S_000 | P_010 S_000 )^0 = y * ( H_401 S_000 | S_000 S_000 )^0_{e} - ( I_411 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[7] = etfac[1] * AUX_INT__h_s_s_s[2] - p_over_q * AUX_INT__i_s_s_s[4];

                    // ( H_401 S_000 | P_001 S_000 )^0 = z * ( H_401 S_000 | S_000 S_000 )^0_{e} + ( G_400 S_000 | S_000 S_000 )^0_{e} - ( I_402 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[8] = etfac[2] * AUX_INT__h_s_s_s[2] + 1 * one_over_2q * AUX_INT__g_s_s_s[0] - p_over_q * AUX_INT__i_s_s_s[5];

                    // ( H_320 S_000 | P_100 S_000 )^0 = x * ( H_320 S_000 | S_000 S_000 )^0_{e} + ( G_220 S_000 | S_000 S_000 )^0_{e} - ( I_420 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[9] = etfac[0] * AUX_INT__h_s_s_s[3] + 3 * one_over_2q * AUX_INT__g_s_s_s[3] - p_over_q * AUX_INT__i_s_s_s[3];

                    // ( H_320 S_000 | P_010 S_000 )^0 = y * ( H_320 S_000 | S_000 S_000 )^0_{e} + ( G_310 S_000 | S_000 S_000 )^0_{e} - ( I_330 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[10] = etfac[1] * AUX_INT__h_s_s_s[3] + 2 * one_over_2q * AUX_INT__g_s_s_s[1] - p_over_q * AUX_INT__i_s_s_s[6];

                    // ( H_311 S_000 | P_100 S_000 )^0 = x * ( H_311 S_000 | S_000 S_000 )^0_{e} + ( G_211 S_000 | S_000 S_000 )^0_{e} - ( I_411 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[12] = etfac[0] * AUX_INT__h_s_s_s[4] + 3 * one_over_2q * AUX_INT__g_s_s_s[4] - p_over_q * AUX_INT__i_s_s_s[4];

                    // ( H_311 S_000 | P_010 S_000 )^0 = y * ( H_311 S_000 | S_000 S_000 )^0_{e} + ( G_301 S_000 | S_000 S_000 )^0_{e} - ( I_321 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[13] = etfac[1] * AUX_INT__h_s_s_s[4] + 1 * one_over_2q * AUX_INT__g_s_s_s[2] - p_over_q * AUX_INT__i_s_s_s[7];

                    // ( H_311 S_000 | P_001 S_000 )^0 = z * ( H_311 S_000 | S_000 S_000 )^0_{e} + ( G_310 S_000 | S_000 S_000 )^0_{e} - ( I_312 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[14] = etfac[2] * AUX_INT__h_s_s_s[4] + 1 * one_over_2q * AUX_INT__g_s_s_s[1] - p_over_q * AUX_INT__i_s_s_s[8];

                    // ( H_302 S_000 | P_100 S_000 )^0 = x * ( H_302 S_000 | S_000 S_000 )^0_{e} + ( G_202 S_000 | S_000 S_000 )^0_{e} - ( I_402 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[15] = etfac[0] * AUX_INT__h_s_s_s[5] + 3 * one_over_2q * AUX_INT__g_s_s_s[5] - p_over_q * AUX_INT__i_s_s_s[5];

                    // ( H_302 S_000 | P_010 S_000 )^0 = y * ( H_302 S_000 | S_000 S_000 )^0_{e} - ( I_312 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[16] = etfac[1] * AUX_INT__h_s_s_s[5] - p_over_q * AUX_INT__i_s_s_s[8];

                    // ( H_302 S_000 | P_001 S_000 )^0 = z * ( H_302 S_000 | S_000 S_000 )^0_{e} + ( G_301 S_000 | S_000 S_000 )^0_{e} - ( I_303 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[17] = etfac[2] * AUX_INT__h_s_s_s[5] + 2 * one_over_2q * AUX_INT__g_s_s_s[2] - p_over_q * AUX_INT__i_s_s_s[9];

                    // ( H_230 S_000 | P_100 S_000 )^0 = x * ( H_230 S_000 | S_000 S_000 )^0_{e} + ( G_130 S_000 | S_000 S_000 )^0_{e} - ( I_330 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[18] = etfac[0] * AUX_INT__h_s_s_s[6] + 2 * one_over_2q * AUX_INT__g_s_s_s[6] - p_over_q * AUX_INT__i_s_s_s[6];

                    // ( H_230 S_000 | P_010 S_000 )^0 = y * ( H_230 S_000 | S_000 S_000 )^0_{e} + ( G_220 S_000 | S_000 S_000 )^0_{e} - ( I_240 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[19] = etfac[1] * AUX_INT__h_s_s_s[6] + 3 * one_over_2q * AUX_INT__g_s_s_s[3] - p_over_q * AUX_INT__i_s_s_s[10];

                    // ( H_221 S_000 | P_100 S_000 )^0 = x * ( H_221 S_000 | S_000 S_000 )^0_{e} + ( G_121 S_000 | S_000 S_000 )^0_{e} - ( I_321 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[21] = etfac[0] * AUX_INT__h_s_s_s[7] + 2 * one_over_2q * AUX_INT__g_s_s_s[7] - p_over_q * AUX_INT__i_s_s_s[7];

                    // ( H_221 S_000 | P_010 S_000 )^0 = y * ( H_221 S_000 | S_000 S_000 )^0_{e} + ( G_211 S_000 | S_000 S_000 )^0_{e} - ( I_231 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[22] = etfac[1] * AUX_INT__h_s_s_s[7] + 2 * one_over_2q * AUX_INT__g_s_s_s[4] - p_over_q * AUX_INT__i_s_s_s[11];

                    // ( H_221 S_000 | P_001 S_000 )^0 = z * ( H_221 S_000 | S_000 S_000 )^0_{e} + ( G_220 S_000 | S_000 S_000 )^0_{e} - ( I_222 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[23] = etfac[2] * AUX_INT__h_s_s_s[7] + 1 * one_over_2q * AUX_INT__g_s_s_s[3] - p_over_q * AUX_INT__i_s_s_s[12];

                    // ( H_212 S_000 | P_100 S_000 )^0 = x * ( H_212 S_000 | S_000 S_000 )^0_{e} + ( G_112 S_000 | S_000 S_000 )^0_{e} - ( I_312 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[24] = etfac[0] * AUX_INT__h_s_s_s[8] + 2 * one_over_2q * AUX_INT__g_s_s_s[8] - p_over_q * AUX_INT__i_s_s_s[8];

                    // ( H_212 S_000 | P_010 S_000 )^0 = y * ( H_212 S_000 | S_000 S_000 )^0_{e} + ( G_202 S_000 | S_000 S_000 )^0_{e} - ( I_222 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[25] = etfac[1] * AUX_INT__h_s_s_s[8] + 1 * one_over_2q * AUX_INT__g_s_s_s[5] - p_over_q * AUX_INT__i_s_s_s[12];

                    // ( H_212 S_000 | P_001 S_000 )^0 = z * ( H_212 S_000 | S_000 S_000 )^0_{e} + ( G_211 S_000 | S_000 S_000 )^0_{e} - ( I_213 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[26] = etfac[2] * AUX_INT__h_s_s_s[8] + 2 * one_over_2q * AUX_INT__g_s_s_s[4] - p_over_q * AUX_INT__i_s_s_s[13];

                    // ( H_203 S_000 | P_100 S_000 )^0 = x * ( H_203 S_000 | S_000 S_000 )^0_{e} + ( G_103 S_000 | S_000 S_000 )^0_{e} - ( I_303 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[27] = etfac[0] * AUX_INT__h_s_s_s[9] + 2 * one_over_2q * AUX_INT__g_s_s_s[9] - p_over_q * AUX_INT__i_s_s_s[9];

                    // ( H_203 S_000 | P_010 S_000 )^0 = y * ( H_203 S_000 | S_000 S_000 )^0_{e} - ( I_213 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[28] = etfac[1] * AUX_INT__h_s_s_s[9] - p_over_q * AUX_INT__i_s_s_s[13];

                    // ( H_203 S_000 | P_001 S_000 )^0 = z * ( H_203 S_000 | S_000 S_000 )^0_{e} + ( G_202 S_000 | S_000 S_000 )^0_{e} - ( I_204 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[29] = etfac[2] * AUX_INT__h_s_s_s[9] + 3 * one_over_2q * AUX_INT__g_s_s_s[5] - p_over_q * AUX_INT__i_s_s_s[14];

                    // ( H_140 S_000 | P_100 S_000 )^0 = x * ( H_140 S_000 | S_000 S_000 )^0_{e} + ( G_040 S_000 | S_000 S_000 )^0_{e} - ( I_240 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[30] = etfac[0] * AUX_INT__h_s_s_s[10] + 1 * one_over_2q * AUX_INT__g_s_s_s[10] - p_over_q * AUX_INT__i_s_s_s[10];

                    // ( H_140 S_000 | P_010 S_000 )^0 = y * ( H_140 S_000 | S_000 S_000 )^0_{e} + ( G_130 S_000 | S_000 S_000 )^0_{e} - ( I_150 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[31] = etfac[1] * AUX_INT__h_s_s_s[10] + 4 * one_over_2q * AUX_INT__g_s_s_s[6] - p_over_q * AUX_INT__i_s_s_s[15];

                    // ( H_131 S_000 | P_100 S_000 )^0 = x * ( H_131 S_000 | S_000 S_000 )^0_{e} + ( G_031 S_000 | S_000 S_000 )^0_{e} - ( I_231 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[33] = etfac[0] * AUX_INT__h_s_s_s[11] + 1 * one_over_2q * AUX_INT__g_s_s_s[11] - p_over_q * AUX_INT__i_s_s_s[11];

                    // ( H_131 S_000 | P_010 S_000 )^0 = y * ( H_131 S_000 | S_000 S_000 )^0_{e} + ( G_121 S_000 | S_000 S_000 )^0_{e} - ( I_141 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[34] = etfac[1] * AUX_INT__h_s_s_s[11] + 3 * one_over_2q * AUX_INT__g_s_s_s[7] - p_over_q * AUX_INT__i_s_s_s[16];

                    // ( H_131 S_000 | P_001 S_000 )^0 = z * ( H_131 S_000 | S_000 S_000 )^0_{e} + ( G_130 S_000 | S_000 S_000 )^0_{e} - ( I_132 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[35] = etfac[2] * AUX_INT__h_s_s_s[11] + 1 * one_over_2q * AUX_INT__g_s_s_s[6] - p_over_q * AUX_INT__i_s_s_s[17];

                    // ( H_122 S_000 | P_100 S_000 )^0 = x * ( H_122 S_000 | S_000 S_000 )^0_{e} + ( G_022 S_000 | S_000 S_000 )^0_{e} - ( I_222 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[36] = etfac[0] * AUX_INT__h_s_s_s[12] + 1 * one_over_2q * AUX_INT__g_s_s_s[12] - p_over_q * AUX_INT__i_s_s_s[12];

                    // ( H_122 S_000 | P_010 S_000 )^0 = y * ( H_122 S_000 | S_000 S_000 )^0_{e} + ( G_112 S_000 | S_000 S_000 )^0_{e} - ( I_132 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[37] = etfac[1] * AUX_INT__h_s_s_s[12] + 2 * one_over_2q * AUX_INT__g_s_s_s[8] - p_over_q * AUX_INT__i_s_s_s[17];

                    // ( H_122 S_000 | P_001 S_000 )^0 = z * ( H_122 S_000 | S_000 S_000 )^0_{e} + ( G_121 S_000 | S_000 S_000 )^0_{e} - ( I_123 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[38] = etfac[2] * AUX_INT__h_s_s_s[12] + 2 * one_over_2q * AUX_INT__g_s_s_s[7] - p_over_q * AUX_INT__i_s_s_s[18];

                    // ( H_113 S_000 | P_100 S_000 )^0 = x * ( H_113 S_000 | S_000 S_000 )^0_{e} + ( G_013 S_000 | S_000 S_000 )^0_{e} - ( I_213 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[39] = etfac[0] * AUX_INT__h_s_s_s[13] + 1 * one_over_2q * AUX_INT__g_s_s_s[13] - p_over_q * AUX_INT__i_s_s_s[13];

                    // ( H_113 S_000 | P_010 S_000 )^0 = y * ( H_113 S_000 | S_000 S_000 )^0_{e} + ( G_103 S_000 | S_000 S_000 )^0_{e} - ( I_123 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[40] = etfac[1] * AUX_INT__h_s_s_s[13] + 1 * one_over_2q * AUX_INT__g_s_s_s[9] - p_over_q * AUX_INT__i_s_s_s[18];

                    // ( H_113 S_000 | P_001 S_000 )^0 = z * ( H_113 S_000 | S_000 S_000 )^0_{e} + ( G_112 S_000 | S_000 S_000 )^0_{e} - ( I_114 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[41] = etfac[2] * AUX_INT__h_s_s_s[13] + 3 * one_over_2q * AUX_INT__g_s_s_s[8] - p_over_q * AUX_INT__i_s_s_s[19];

                    // ( H_104 S_000 | P_100 S_000 )^0 = x * ( H_104 S_000 | S_000 S_000 )^0_{e} + ( G_004 S_000 | S_000 S_000 )^0_{e} - ( I_204 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[42] = etfac[0] * AUX_INT__h_s_s_s[14] + 1 * one_over_2q * AUX_INT__g_s_s_s[14] - p_over_q * AUX_INT__i_s_s_s[14];

                    // ( H_104 S_000 | P_010 S_000 )^0 = y * ( H_104 S_000 | S_000 S_000 )^0_{e} - ( I_114 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[43] = etfac[1] * AUX_INT__h_s_s_s[14] - p_over_q * AUX_INT__i_s_s_s[19];

                    // ( H_104 S_000 | P_001 S_000 )^0 = z * ( H_104 S_000 | S_000 S_000 )^0_{e} + ( G_103 S_000 | S_000 S_000 )^0_{e} - ( I_105 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[44] = etfac[2] * AUX_INT__h_s_s_s[14] + 4 * one_over_2q * AUX_INT__g_s_s_s[9] - p_over_q * AUX_INT__i_s_s_s[20];

                    // ( H_050 S_000 | P_100 S_000 )^0 = x * ( H_050 S_000 | S_000 S_000 )^0_{e} - ( I_150 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[45] = etfac[0] * AUX_INT__h_s_s_s[15] - p_over_q * AUX_INT__i_s_s_s[15];

                    // ( H_050 S_000 | P_010 S_000 )^0 = y * ( H_050 S_000 | S_000 S_000 )^0_{e} + ( G_040 S_000 | S_000 S_000 )^0_{e} - ( I_060 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[46] = etfac[1] * AUX_INT__h_s_s_s[15] + 5 * one_over_2q * AUX_INT__g_s_s_s[10] - p_over_q * AUX_INT__i_s_s_s[21];

                    // ( H_041 S_000 | P_100 S_000 )^0 = x * ( H_041 S_000 | S_000 S_000 )^0_{e} - ( I_141 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[48] = etfac[0] * AUX_INT__h_s_s_s[16] - p_over_q * AUX_INT__i_s_s_s[16];

                    // ( H_041 S_000 | P_010 S_000 )^0 = y * ( H_041 S_000 | S_000 S_000 )^0_{e} + ( G_031 S_000 | S_000 S_000 )^0_{e} - ( I_051 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[49] = etfac[1] * AUX_INT__h_s_s_s[16] + 4 * one_over_2q * AUX_INT__g_s_s_s[11] - p_over_q * AUX_INT__i_s_s_s[22];

                    // ( H_041 S_000 | P_001 S_000 )^0 = z * ( H_041 S_000 | S_000 S_000 )^0_{e} + ( G_040 S_000 | S_000 S_000 )^0_{e} - ( I_042 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[50] = etfac[2] * AUX_INT__h_s_s_s[16] + 1 * one_over_2q * AUX_INT__g_s_s_s[10] - p_over_q * AUX_INT__i_s_s_s[23];

                    // ( H_032 S_000 | P_100 S_000 )^0 = x * ( H_032 S_000 | S_000 S_000 )^0_{e} - ( I_132 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[51] = etfac[0] * AUX_INT__h_s_s_s[17] - p_over_q * AUX_INT__i_s_s_s[17];

                    // ( H_032 S_000 | P_010 S_000 )^0 = y * ( H_032 S_000 | S_000 S_000 )^0_{e} + ( G_022 S_000 | S_000 S_000 )^0_{e} - ( I_042 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[52] = etfac[1] * AUX_INT__h_s_s_s[17] + 3 * one_over_2q * AUX_INT__g_s_s_s[12] - p_over_q * AUX_INT__i_s_s_s[23];

                    // ( H_032 S_000 | P_001 S_000 )^0 = z * ( H_032 S_000 | S_000 S_000 )^0_{e} + ( G_031 S_000 | S_000 S_000 )^0_{e} - ( I_033 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[53] = etfac[2] * AUX_INT__h_s_s_s[17] + 2 * one_over_2q * AUX_INT__g_s_s_s[11] - p_over_q * AUX_INT__i_s_s_s[24];

                    // ( H_023 S_000 | P_100 S_000 )^0 = x * ( H_023 S_000 | S_000 S_000 )^0_{e} - ( I_123 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[54] = etfac[0] * AUX_INT__h_s_s_s[18] - p_over_q * AUX_INT__i_s_s_s[18];

                    // ( H_023 S_000 | P_010 S_000 )^0 = y * ( H_023 S_000 | S_000 S_000 )^0_{e} + ( G_013 S_000 | S_000 S_000 )^0_{e} - ( I_033 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[55] = etfac[1] * AUX_INT__h_s_s_s[18] + 2 * one_over_2q * AUX_INT__g_s_s_s[13] - p_over_q * AUX_INT__i_s_s_s[24];

                    // ( H_023 S_000 | P_001 S_000 )^0 = z * ( H_023 S_000 | S_000 S_000 )^0_{e} + ( G_022 S_000 | S_000 S_000 )^0_{e} - ( I_024 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[56] = etfac[2] * AUX_INT__h_s_s_s[18] + 3 * one_over_2q * AUX_INT__g_s_s_s[12] - p_over_q * AUX_INT__i_s_s_s[25];

                    // ( H_014 S_000 | P_100 S_000 )^0 = x * ( H_014 S_000 | S_000 S_000 )^0_{e} - ( I_114 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[57] = etfac[0] * AUX_INT__h_s_s_s[19] - p_over_q * AUX_INT__i_s_s_s[19];

                    // ( H_014 S_000 | P_010 S_000 )^0 = y * ( H_014 S_000 | S_000 S_000 )^0_{e} + ( G_004 S_000 | S_000 S_000 )^0_{e} - ( I_024 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[58] = etfac[1] * AUX_INT__h_s_s_s[19] + 1 * one_over_2q * AUX_INT__g_s_s_s[14] - p_over_q * AUX_INT__i_s_s_s[25];

                    // ( H_014 S_000 | P_001 S_000 )^0 = z * ( H_014 S_000 | S_000 S_000 )^0_{e} + ( G_013 S_000 | S_000 S_000 )^0_{e} - ( I_015 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[59] = etfac[2] * AUX_INT__h_s_s_s[19] + 4 * one_over_2q * AUX_INT__g_s_s_s[13] - p_over_q * AUX_INT__i_s_s_s[26];

                    // ( H_005 S_000 | P_100 S_000 )^0 = x * ( H_005 S_000 | S_000 S_000 )^0_{e} - ( I_105 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[60] = etfac[0] * AUX_INT__h_s_s_s[20] - p_over_q * AUX_INT__i_s_s_s[20];

                    // ( H_005 S_000 | P_010 S_000 )^0 = y * ( H_005 S_000 | S_000 S_000 )^0_{e} - ( I_015 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[61] = etfac[1] * AUX_INT__h_s_s_s[20] - p_over_q * AUX_INT__i_s_s_s[26];

                    // ( H_005 S_000 | P_001 S_000 )^0 = z * ( H_005 S_000 | S_000 S_000 )^0_{e} + ( G_004 S_000 | S_000 S_000 )^0_{e} - ( I_006 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[62] = etfac[2] * AUX_INT__h_s_s_s[20] + 5 * one_over_2q * AUX_INT__g_s_s_s[14] - p_over_q * AUX_INT__i_s_s_s[27];

                    // ( G_400 S_000 | D_200 S_000 )^0_{t} = x * ( G_400 S_000 | P_100 S_000 )^0 + ( F_300 S_000 | P_100 S_000 )^0 + ( G_400 S_000 | S_000 S_000 )^0_{e} - ( H_500 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[0] = etfac[0] * AUX_INT__g_s_p_s[0] + 4 * one_over_2q * AUX_INT__f_s_p_s[0] + 1 * one_over_2q * AUX_INT__g_s_s_s[0] - p_over_q * AUX_INT__h_s_p_s[0];

                    // ( G_400 S_000 | D_110 S_000 )^0_{t} = y * ( G_400 S_000 | P_100 S_000 )^0 - ( H_410 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[1] = etfac[1] * AUX_INT__g_s_p_s[0] - p_over_q * AUX_INT__h_s_p_s[3];

                    // ( G_400 S_000 | D_101 S_000 )^0_{t} = z * ( G_400 S_000 | P_100 S_000 )^0 - ( H_401 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[2] = etfac[2] * AUX_INT__g_s_p_s[0] - p_over_q * AUX_INT__h_s_p_s[6];

                    // ( G_400 S_000 | D_020 S_000 )^0_{t} = y * ( G_400 S_000 | P_010 S_000 )^0 + ( G_400 S_000 | S_000 S_000 )^0_{e} - ( H_410 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[3] = etfac[1] * AUX_INT__g_s_p_s[1] + 1 * one_over_2q * AUX_INT__g_s_s_s[0] - p_over_q * AUX_INT__h_s_p_s[4];

                    // ( G_400 S_000 | D_011 S_000 )^0_{t} = z * ( G_400 S_000 | P_010 S_000 )^0 - ( H_401 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[4] = etfac[2] * AUX_INT__g_s_p_s[1] - p_over_q * AUX_INT__h_s_p_s[7];

                    // ( G_400 S_000 | D_002 S_000 )^0_{t} = z * ( G_400 S_000 | P_001 S_000 )^0 + ( G_400 S_000 | S_000 S_000 )^0_{e} - ( H_401 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[5] = etfac[2] * AUX_INT__g_s_p_s[2] + 1 * one_over_2q * AUX_INT__g_s_s_s[0] - p_over_q * AUX_INT__h_s_p_s[8];

                    // ( G_310 S_000 | D_200 S_000 )^0_{t} = x * ( G_310 S_000 | P_100 S_000 )^0 + ( F_210 S_000 | P_100 S_000 )^0 + ( G_310 S_000 | S_000 S_000 )^0_{e} - ( H_410 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[6] = etfac[0] * AUX_INT__g_s_p_s[3] + 3 * one_over_2q * AUX_INT__f_s_p_s[3] + 1 * one_over_2q * AUX_INT__g_s_s_s[1] - p_over_q * AUX_INT__h_s_p_s[3];

                    // ( G_310 S_000 | D_110 S_000 )^0_{t} = y * ( G_310 S_000 | P_100 S_000 )^0 + ( F_300 S_000 | P_100 S_000 )^0 - ( H_320 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[7] = etfac[1] * AUX_INT__g_s_p_s[3] + 1 * one_over_2q * AUX_INT__f_s_p_s[0] - p_over_q * AUX_INT__h_s_p_s[9];

                    // ( G_310 S_000 | D_101 S_000 )^0_{t} = z * ( G_310 S_000 | P_100 S_000 )^0 - ( H_311 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[8] = etfac[2] * AUX_INT__g_s_p_s[3] - p_over_q * AUX_INT__h_s_p_s[12];

                    // ( G_310 S_000 | D_020 S_000 )^0_{t} = y * ( G_310 S_000 | P_010 S_000 )^0 + ( F_300 S_000 | P_010 S_000 )^0 + ( G_310 S_000 | S_000 S_000 )^0_{e} - ( H_320 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[9] = etfac[1] * AUX_INT__g_s_p_s[4] + 1 * one_over_2q * AUX_INT__f_s_p_s[1] + 1 * one_over_2q * AUX_INT__g_s_s_s[1] - p_over_q * AUX_INT__h_s_p_s[10];

                    // ( G_310 S_000 | D_011 S_000 )^0_{t} = z * ( G_310 S_000 | P_010 S_000 )^0 - ( H_311 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[10] = etfac[2] * AUX_INT__g_s_p_s[4] - p_over_q * AUX_INT__h_s_p_s[13];

                    // ( G_310 S_000 | D_002 S_000 )^0_{t} = z * ( G_310 S_000 | P_001 S_000 )^0 + ( G_310 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[11] = etfac[2] * AUX_INT__g_s_p_s[5] + 1 * one_over_2q * AUX_INT__g_s_s_s[1] - p_over_q * AUX_INT__h_s_p_s[14];

                    // ( G_301 S_000 | D_200 S_000 )^0_{t} = x * ( G_301 S_000 | P_100 S_000 )^0 + ( F_201 S_000 | P_100 S_000 )^0 + ( G_301 S_000 | S_000 S_000 )^0_{e} - ( H_401 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[12] = etfac[0] * AUX_INT__g_s_p_s[6] + 3 * one_over_2q * AUX_INT__f_s_p_s[6] + 1 * one_over_2q * AUX_INT__g_s_s_s[2] - p_over_q * AUX_INT__h_s_p_s[6];

                    // ( G_301 S_000 | D_110 S_000 )^0_{t} = y * ( G_301 S_000 | P_100 S_000 )^0 - ( H_311 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[13] = etfac[1] * AUX_INT__g_s_p_s[6] - p_over_q * AUX_INT__h_s_p_s[12];

                    // ( G_301 S_000 | D_101 S_000 )^0_{t} = z * ( G_301 S_000 | P_100 S_000 )^0 + ( F_300 S_000 | P_100 S_000 )^0 - ( H_302 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[14] = etfac[2] * AUX_INT__g_s_p_s[6] + 1 * one_over_2q * AUX_INT__f_s_p_s[0] - p_over_q * AUX_INT__h_s_p_s[15];

                    // ( G_301 S_000 | D_020 S_000 )^0_{t} = y * ( G_301 S_000 | P_010 S_000 )^0 + ( G_301 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[15] = etfac[1] * AUX_INT__g_s_p_s[7] + 1 * one_over_2q * AUX_INT__g_s_s_s[2] - p_over_q * AUX_INT__h_s_p_s[13];

                    // ( G_301 S_000 | D_011 S_000 )^0_{t} = z * ( G_301 S_000 | P_010 S_000 )^0 + ( F_300 S_000 | P_010 S_000 )^0 - ( H_302 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[16] = etfac[2] * AUX_INT__g_s_p_s[7] + 1 * one_over_2q * AUX_INT__f_s_p_s[1] - p_over_q * AUX_INT__h_s_p_s[16];

                    // ( G_301 S_000 | D_002 S_000 )^0_{t} = z * ( G_301 S_000 | P_001 S_000 )^0 + ( F_300 S_000 | P_001 S_000 )^0 + ( G_301 S_000 | S_000 S_000 )^0_{e} - ( H_302 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[17] = etfac[2] * AUX_INT__g_s_p_s[8] + 1 * one_over_2q * AUX_INT__f_s_p_s[2] + 1 * one_over_2q * AUX_INT__g_s_s_s[2] - p_over_q * AUX_INT__h_s_p_s[17];

                    // ( G_220 S_000 | D_200 S_000 )^0_{t} = x * ( G_220 S_000 | P_100 S_000 )^0 + ( F_120 S_000 | P_100 S_000 )^0 + ( G_220 S_000 | S_000 S_000 )^0_{e} - ( H_320 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[18] = etfac[0] * AUX_INT__g_s_p_s[9] + 2 * one_over_2q * AUX_INT__f_s_p_s[9] + 1 * one_over_2q * AUX_INT__g_s_s_s[3] - p_over_q * AUX_INT__h_s_p_s[9];

                    // ( G_220 S_000 | D_110 S_000 )^0_{t} = y * ( G_220 S_000 | P_100 S_000 )^0 + ( F_210 S_000 | P_100 S_000 )^0 - ( H_230 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[19] = etfac[1] * AUX_INT__g_s_p_s[9] + 2 * one_over_2q * AUX_INT__f_s_p_s[3] - p_over_q * AUX_INT__h_s_p_s[18];

                    // ( G_220 S_000 | D_101 S_000 )^0_{t} = z * ( G_220 S_000 | P_100 S_000 )^0 - ( H_221 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[20] = etfac[2] * AUX_INT__g_s_p_s[9] - p_over_q * AUX_INT__h_s_p_s[21];

                    // ( G_220 S_000 | D_020 S_000 )^0_{t} = y * ( G_220 S_000 | P_010 S_000 )^0 + ( F_210 S_000 | P_010 S_000 )^0 + ( G_220 S_000 | S_000 S_000 )^0_{e} - ( H_230 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[21] = etfac[1] * AUX_INT__g_s_p_s[10] + 2 * one_over_2q * AUX_INT__f_s_p_s[4] + 1 * one_over_2q * AUX_INT__g_s_s_s[3] - p_over_q * AUX_INT__h_s_p_s[19];

                    // ( G_220 S_000 | D_011 S_000 )^0_{t} = z * ( G_220 S_000 | P_010 S_000 )^0 - ( H_221 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[22] = etfac[2] * AUX_INT__g_s_p_s[10] - p_over_q * AUX_INT__h_s_p_s[22];

                    // ( G_220 S_000 | D_002 S_000 )^0_{t} = z * ( G_220 S_000 | P_001 S_000 )^0 + ( G_220 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[23] = etfac[2] * AUX_INT__g_s_p_s[11] + 1 * one_over_2q * AUX_INT__g_s_s_s[3] - p_over_q * AUX_INT__h_s_p_s[23];

                    // ( G_211 S_000 | D_200 S_000 )^0_{t} = x * ( G_211 S_000 | P_100 S_000 )^0 + ( F_111 S_000 | P_100 S_000 )^0 + ( G_211 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[24] = etfac[0] * AUX_INT__g_s_p_s[12] + 2 * one_over_2q * AUX_INT__f_s_p_s[12] + 1 * one_over_2q * AUX_INT__g_s_s_s[4] - p_over_q * AUX_INT__h_s_p_s[12];

                    // ( G_211 S_000 | D_110 S_000 )^0_{t} = y * ( G_211 S_000 | P_100 S_000 )^0 + ( F_201 S_000 | P_100 S_000 )^0 - ( H_221 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[25] = etfac[1] * AUX_INT__g_s_p_s[12] + 1 * one_over_2q * AUX_INT__f_s_p_s[6] - p_over_q * AUX_INT__h_s_p_s[21];

                    // ( G_211 S_000 | D_101 S_000 )^0_{t} = z * ( G_211 S_000 | P_100 S_000 )^0 + ( F_210 S_000 | P_100 S_000 )^0 - ( H_212 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[26] = etfac[2] * AUX_INT__g_s_p_s[12] + 1 * one_over_2q * AUX_INT__f_s_p_s[3] - p_over_q * AUX_INT__h_s_p_s[24];

                    // ( G_211 S_000 | D_020 S_000 )^0_{t} = y * ( G_211 S_000 | P_010 S_000 )^0 + ( F_201 S_000 | P_010 S_000 )^0 + ( G_211 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[27] = etfac[1] * AUX_INT__g_s_p_s[13] + 1 * one_over_2q * AUX_INT__f_s_p_s[7] + 1 * one_over_2q * AUX_INT__g_s_s_s[4] - p_over_q * AUX_INT__h_s_p_s[22];

                    // ( G_211 S_000 | D_011 S_000 )^0_{t} = z * ( G_211 S_000 | P_010 S_000 )^0 + ( F_210 S_000 | P_010 S_000 )^0 - ( H_212 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[28] = etfac[2] * AUX_INT__g_s_p_s[13] + 1 * one_over_2q * AUX_INT__f_s_p_s[4] - p_over_q * AUX_INT__h_s_p_s[25];

                    // ( G_211 S_000 | D_002 S_000 )^0_{t} = z * ( G_211 S_000 | P_001 S_000 )^0 + ( F_210 S_000 | P_001 S_000 )^0 + ( G_211 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[29] = etfac[2] * AUX_INT__g_s_p_s[14] + 1 * one_over_2q * AUX_INT__f_s_p_s[5] + 1 * one_over_2q * AUX_INT__g_s_s_s[4] - p_over_q * AUX_INT__h_s_p_s[26];

                    // ( G_202 S_000 | D_200 S_000 )^0_{t} = x * ( G_202 S_000 | P_100 S_000 )^0 + ( F_102 S_000 | P_100 S_000 )^0 + ( G_202 S_000 | S_000 S_000 )^0_{e} - ( H_302 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[30] = etfac[0] * AUX_INT__g_s_p_s[15] + 2 * one_over_2q * AUX_INT__f_s_p_s[15] + 1 * one_over_2q * AUX_INT__g_s_s_s[5] - p_over_q * AUX_INT__h_s_p_s[15];

                    // ( G_202 S_000 | D_110 S_000 )^0_{t} = y * ( G_202 S_000 | P_100 S_000 )^0 - ( H_212 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[31] = etfac[1] * AUX_INT__g_s_p_s[15] - p_over_q * AUX_INT__h_s_p_s[24];

                    // ( G_202 S_000 | D_101 S_000 )^0_{t} = z * ( G_202 S_000 | P_100 S_000 )^0 + ( F_201 S_000 | P_100 S_000 )^0 - ( H_203 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[32] = etfac[2] * AUX_INT__g_s_p_s[15] + 2 * one_over_2q * AUX_INT__f_s_p_s[6] - p_over_q * AUX_INT__h_s_p_s[27];

                    // ( G_202 S_000 | D_020 S_000 )^0_{t} = y * ( G_202 S_000 | P_010 S_000 )^0 + ( G_202 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[33] = etfac[1] * AUX_INT__g_s_p_s[16] + 1 * one_over_2q * AUX_INT__g_s_s_s[5] - p_over_q * AUX_INT__h_s_p_s[25];

                    // ( G_202 S_000 | D_011 S_000 )^0_{t} = z * ( G_202 S_000 | P_010 S_000 )^0 + ( F_201 S_000 | P_010 S_000 )^0 - ( H_203 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[34] = etfac[2] * AUX_INT__g_s_p_s[16] + 2 * one_over_2q * AUX_INT__f_s_p_s[7] - p_over_q * AUX_INT__h_s_p_s[28];

                    // ( G_202 S_000 | D_002 S_000 )^0_{t} = z * ( G_202 S_000 | P_001 S_000 )^0 + ( F_201 S_000 | P_001 S_000 )^0 + ( G_202 S_000 | S_000 S_000 )^0_{e} - ( H_203 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[35] = etfac[2] * AUX_INT__g_s_p_s[17] + 2 * one_over_2q * AUX_INT__f_s_p_s[8] + 1 * one_over_2q * AUX_INT__g_s_s_s[5] - p_over_q * AUX_INT__h_s_p_s[29];

                    // ( G_130 S_000 | D_200 S_000 )^0_{t} = x * ( G_130 S_000 | P_100 S_000 )^0 + ( F_030 S_000 | P_100 S_000 )^0 + ( G_130 S_000 | S_000 S_000 )^0_{e} - ( H_230 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[36] = etfac[0] * AUX_INT__g_s_p_s[18] + 1 * one_over_2q * AUX_INT__f_s_p_s[18] + 1 * one_over_2q * AUX_INT__g_s_s_s[6] - p_over_q * AUX_INT__h_s_p_s[18];

                    // ( G_130 S_000 | D_110 S_000 )^0_{t} = y * ( G_130 S_000 | P_100 S_000 )^0 + ( F_120 S_000 | P_100 S_000 )^0 - ( H_140 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[37] = etfac[1] * AUX_INT__g_s_p_s[18] + 3 * one_over_2q * AUX_INT__f_s_p_s[9] - p_over_q * AUX_INT__h_s_p_s[30];

                    // ( G_130 S_000 | D_101 S_000 )^0_{t} = z * ( G_130 S_000 | P_100 S_000 )^0 - ( H_131 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[38] = etfac[2] * AUX_INT__g_s_p_s[18] - p_over_q * AUX_INT__h_s_p_s[33];

                    // ( G_130 S_000 | D_020 S_000 )^0_{t} = y * ( G_130 S_000 | P_010 S_000 )^0 + ( F_120 S_000 | P_010 S_000 )^0 + ( G_130 S_000 | S_000 S_000 )^0_{e} - ( H_140 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[39] = etfac[1] * AUX_INT__g_s_p_s[19] + 3 * one_over_2q * AUX_INT__f_s_p_s[10] + 1 * one_over_2q * AUX_INT__g_s_s_s[6] - p_over_q * AUX_INT__h_s_p_s[31];

                    // ( G_130 S_000 | D_011 S_000 )^0_{t} = z * ( G_130 S_000 | P_010 S_000 )^0 - ( H_131 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[40] = etfac[2] * AUX_INT__g_s_p_s[19] - p_over_q * AUX_INT__h_s_p_s[34];

                    // ( G_130 S_000 | D_002 S_000 )^0_{t} = z * ( G_130 S_000 | P_001 S_000 )^0 + ( G_130 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[41] = etfac[2] * AUX_INT__g_s_p_s[20] + 1 * one_over_2q * AUX_INT__g_s_s_s[6] - p_over_q * AUX_INT__h_s_p_s[35];

                    // ( G_121 S_000 | D_200 S_000 )^0_{t} = x * ( G_121 S_000 | P_100 S_000 )^0 + ( F_021 S_000 | P_100 S_000 )^0 + ( G_121 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[42] = etfac[0] * AUX_INT__g_s_p_s[21] + 1 * one_over_2q * AUX_INT__f_s_p_s[21] + 1 * one_over_2q * AUX_INT__g_s_s_s[7] - p_over_q * AUX_INT__h_s_p_s[21];

                    // ( G_121 S_000 | D_110 S_000 )^0_{t} = y * ( G_121 S_000 | P_100 S_000 )^0 + ( F_111 S_000 | P_100 S_000 )^0 - ( H_131 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[43] = etfac[1] * AUX_INT__g_s_p_s[21] + 2 * one_over_2q * AUX_INT__f_s_p_s[12] - p_over_q * AUX_INT__h_s_p_s[33];

                    // ( G_121 S_000 | D_101 S_000 )^0_{t} = z * ( G_121 S_000 | P_100 S_000 )^0 + ( F_120 S_000 | P_100 S_000 )^0 - ( H_122 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[44] = etfac[2] * AUX_INT__g_s_p_s[21] + 1 * one_over_2q * AUX_INT__f_s_p_s[9] - p_over_q * AUX_INT__h_s_p_s[36];

                    // ( G_121 S_000 | D_020 S_000 )^0_{t} = y * ( G_121 S_000 | P_010 S_000 )^0 + ( F_111 S_000 | P_010 S_000 )^0 + ( G_121 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[45] = etfac[1] * AUX_INT__g_s_p_s[22] + 2 * one_over_2q * AUX_INT__f_s_p_s[13] + 1 * one_over_2q * AUX_INT__g_s_s_s[7] - p_over_q * AUX_INT__h_s_p_s[34];

                    // ( G_121 S_000 | D_011 S_000 )^0_{t} = z * ( G_121 S_000 | P_010 S_000 )^0 + ( F_120 S_000 | P_010 S_000 )^0 - ( H_122 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[46] = etfac[2] * AUX_INT__g_s_p_s[22] + 1 * one_over_2q * AUX_INT__f_s_p_s[10] - p_over_q * AUX_INT__h_s_p_s[37];

                    // ( G_121 S_000 | D_002 S_000 )^0_{t} = z * ( G_121 S_000 | P_001 S_000 )^0 + ( F_120 S_000 | P_001 S_000 )^0 + ( G_121 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[47] = etfac[2] * AUX_INT__g_s_p_s[23] + 1 * one_over_2q * AUX_INT__f_s_p_s[11] + 1 * one_over_2q * AUX_INT__g_s_s_s[7] - p_over_q * AUX_INT__h_s_p_s[38];

                    // ( G_112 S_000 | D_200 S_000 )^0_{t} = x * ( G_112 S_000 | P_100 S_000 )^0 + ( F_012 S_000 | P_100 S_000 )^0 + ( G_112 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[48] = etfac[0] * AUX_INT__g_s_p_s[24] + 1 * one_over_2q * AUX_INT__f_s_p_s[24] + 1 * one_over_2q * AUX_INT__g_s_s_s[8] - p_over_q * AUX_INT__h_s_p_s[24];

                    // ( G_112 S_000 | D_110 S_000 )^0_{t} = y * ( G_112 S_000 | P_100 S_000 )^0 + ( F_102 S_000 | P_100 S_000 )^0 - ( H_122 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[49] = etfac[1] * AUX_INT__g_s_p_s[24] + 1 * one_over_2q * AUX_INT__f_s_p_s[15] - p_over_q * AUX_INT__h_s_p_s[36];

                    // ( G_112 S_000 | D_101 S_000 )^0_{t} = z * ( G_112 S_000 | P_100 S_000 )^0 + ( F_111 S_000 | P_100 S_000 )^0 - ( H_113 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[50] = etfac[2] * AUX_INT__g_s_p_s[24] + 2 * one_over_2q * AUX_INT__f_s_p_s[12] - p_over_q * AUX_INT__h_s_p_s[39];

                    // ( G_112 S_000 | D_020 S_000 )^0_{t} = y * ( G_112 S_000 | P_010 S_000 )^0 + ( F_102 S_000 | P_010 S_000 )^0 + ( G_112 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[51] = etfac[1] * AUX_INT__g_s_p_s[25] + 1 * one_over_2q * AUX_INT__f_s_p_s[16] + 1 * one_over_2q * AUX_INT__g_s_s_s[8] - p_over_q * AUX_INT__h_s_p_s[37];

                    // ( G_112 S_000 | D_011 S_000 )^0_{t} = z * ( G_112 S_000 | P_010 S_000 )^0 + ( F_111 S_000 | P_010 S_000 )^0 - ( H_113 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[52] = etfac[2] * AUX_INT__g_s_p_s[25] + 2 * one_over_2q * AUX_INT__f_s_p_s[13] - p_over_q * AUX_INT__h_s_p_s[40];

                    // ( G_112 S_000 | D_002 S_000 )^0_{t} = z * ( G_112 S_000 | P_001 S_000 )^0 + ( F_111 S_000 | P_001 S_000 )^0 + ( G_112 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[53] = etfac[2] * AUX_INT__g_s_p_s[26] + 2 * one_over_2q * AUX_INT__f_s_p_s[14] + 1 * one_over_2q * AUX_INT__g_s_s_s[8] - p_over_q * AUX_INT__h_s_p_s[41];

                    // ( G_103 S_000 | D_200 S_000 )^0_{t} = x * ( G_103 S_000 | P_100 S_000 )^0 + ( F_003 S_000 | P_100 S_000 )^0 + ( G_103 S_000 | S_000 S_000 )^0_{e} - ( H_203 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[54] = etfac[0] * AUX_INT__g_s_p_s[27] + 1 * one_over_2q * AUX_INT__f_s_p_s[27] + 1 * one_over_2q * AUX_INT__g_s_s_s[9] - p_over_q * AUX_INT__h_s_p_s[27];

                    // ( G_103 S_000 | D_110 S_000 )^0_{t} = y * ( G_103 S_000 | P_100 S_000 )^0 - ( H_113 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[55] = etfac[1] * AUX_INT__g_s_p_s[27] - p_over_q * AUX_INT__h_s_p_s[39];

                    // ( G_103 S_000 | D_101 S_000 )^0_{t} = z * ( G_103 S_000 | P_100 S_000 )^0 + ( F_102 S_000 | P_100 S_000 )^0 - ( H_104 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[56] = etfac[2] * AUX_INT__g_s_p_s[27] + 3 * one_over_2q * AUX_INT__f_s_p_s[15] - p_over_q * AUX_INT__h_s_p_s[42];

                    // ( G_103 S_000 | D_020 S_000 )^0_{t} = y * ( G_103 S_000 | P_010 S_000 )^0 + ( G_103 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[57] = etfac[1] * AUX_INT__g_s_p_s[28] + 1 * one_over_2q * AUX_INT__g_s_s_s[9] - p_over_q * AUX_INT__h_s_p_s[40];

                    // ( G_103 S_000 | D_011 S_000 )^0_{t} = z * ( G_103 S_000 | P_010 S_000 )^0 + ( F_102 S_000 | P_010 S_000 )^0 - ( H_104 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[58] = etfac[2] * AUX_INT__g_s_p_s[28] + 3 * one_over_2q * AUX_INT__f_s_p_s[16] - p_over_q * AUX_INT__h_s_p_s[43];

                    // ( G_103 S_000 | D_002 S_000 )^0_{t} = z * ( G_103 S_000 | P_001 S_000 )^0 + ( F_102 S_000 | P_001 S_000 )^0 + ( G_103 S_000 | S_000 S_000 )^0_{e} - ( H_104 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[59] = etfac[2] * AUX_INT__g_s_p_s[29] + 3 * one_over_2q * AUX_INT__f_s_p_s[17] + 1 * one_over_2q * AUX_INT__g_s_s_s[9] - p_over_q * AUX_INT__h_s_p_s[44];

                    // ( G_040 S_000 | D_200 S_000 )^0_{t} = x * ( G_040 S_000 | P_100 S_000 )^0 + ( G_040 S_000 | S_000 S_000 )^0_{e} - ( H_140 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[60] = etfac[0] * AUX_INT__g_s_p_s[30] + 1 * one_over_2q * AUX_INT__g_s_s_s[10] - p_over_q * AUX_INT__h_s_p_s[30];

                    // ( G_040 S_000 | D_110 S_000 )^0_{t} = y * ( G_040 S_000 | P_100 S_000 )^0 + ( F_030 S_000 | P_100 S_000 )^0 - ( H_050 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[61] = etfac[1] * AUX_INT__g_s_p_s[30] + 4 * one_over_2q * AUX_INT__f_s_p_s[18] - p_over_q * AUX_INT__h_s_p_s[45];

                    // ( G_040 S_000 | D_101 S_000 )^0_{t} = z * ( G_040 S_000 | P_100 S_000 )^0 - ( H_041 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[62] = etfac[2] * AUX_INT__g_s_p_s[30] - p_over_q * AUX_INT__h_s_p_s[48];

                    // ( G_040 S_000 | D_020 S_000 )^0_{t} = y * ( G_040 S_000 | P_010 S_000 )^0 + ( F_030 S_000 | P_010 S_000 )^0 + ( G_040 S_000 | S_000 S_000 )^0_{e} - ( H_050 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[63] = etfac[1] * AUX_INT__g_s_p_s[31] + 4 * one_over_2q * AUX_INT__f_s_p_s[19] + 1 * one_over_2q * AUX_INT__g_s_s_s[10] - p_over_q * AUX_INT__h_s_p_s[46];

                    // ( G_040 S_000 | D_011 S_000 )^0_{t} = z * ( G_040 S_000 | P_010 S_000 )^0 - ( H_041 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[64] = etfac[2] * AUX_INT__g_s_p_s[31] - p_over_q * AUX_INT__h_s_p_s[49];

                    // ( G_040 S_000 | D_002 S_000 )^0_{t} = z * ( G_040 S_000 | P_001 S_000 )^0 + ( G_040 S_000 | S_000 S_000 )^0_{e} - ( H_041 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[65] = etfac[2] * AUX_INT__g_s_p_s[32] + 1 * one_over_2q * AUX_INT__g_s_s_s[10] - p_over_q * AUX_INT__h_s_p_s[50];

                    // ( G_031 S_000 | D_200 S_000 )^0_{t} = x * ( G_031 S_000 | P_100 S_000 )^0 + ( G_031 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[66] = etfac[0] * AUX_INT__g_s_p_s[33] + 1 * one_over_2q * AUX_INT__g_s_s_s[11] - p_over_q * AUX_INT__h_s_p_s[33];

                    // ( G_031 S_000 | D_110 S_000 )^0_{t} = y * ( G_031 S_000 | P_100 S_000 )^0 + ( F_021 S_000 | P_100 S_000 )^0 - ( H_041 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[67] = etfac[1] * AUX_INT__g_s_p_s[33] + 3 * one_over_2q * AUX_INT__f_s_p_s[21] - p_over_q * AUX_INT__h_s_p_s[48];

                    // ( G_031 S_000 | D_101 S_000 )^0_{t} = z * ( G_031 S_000 | P_100 S_000 )^0 + ( F_030 S_000 | P_100 S_000 )^0 - ( H_032 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[68] = etfac[2] * AUX_INT__g_s_p_s[33] + 1 * one_over_2q * AUX_INT__f_s_p_s[18] - p_over_q * AUX_INT__h_s_p_s[51];

                    // ( G_031 S_000 | D_020 S_000 )^0_{t} = y * ( G_031 S_000 | P_010 S_000 )^0 + ( F_021 S_000 | P_010 S_000 )^0 + ( G_031 S_000 | S_000 S_000 )^0_{e} - ( H_041 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[69] = etfac[1] * AUX_INT__g_s_p_s[34] + 3 * one_over_2q * AUX_INT__f_s_p_s[22] + 1 * one_over_2q * AUX_INT__g_s_s_s[11] - p_over_q * AUX_INT__h_s_p_s[49];

                    // ( G_031 S_000 | D_011 S_000 )^0_{t} = z * ( G_031 S_000 | P_010 S_000 )^0 + ( F_030 S_000 | P_010 S_000 )^0 - ( H_032 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[70] = etfac[2] * AUX_INT__g_s_p_s[34] + 1 * one_over_2q * AUX_INT__f_s_p_s[19] - p_over_q * AUX_INT__h_s_p_s[52];

                    // ( G_031 S_000 | D_002 S_000 )^0_{t} = z * ( G_031 S_000 | P_001 S_000 )^0 + ( F_030 S_000 | P_001 S_000 )^0 + ( G_031 S_000 | S_000 S_000 )^0_{e} - ( H_032 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[71] = etfac[2] * AUX_INT__g_s_p_s[35] + 1 * one_over_2q * AUX_INT__f_s_p_s[20] + 1 * one_over_2q * AUX_INT__g_s_s_s[11] - p_over_q * AUX_INT__h_s_p_s[53];

                    // ( G_022 S_000 | D_200 S_000 )^0_{t} = x * ( G_022 S_000 | P_100 S_000 )^0 + ( G_022 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[72] = etfac[0] * AUX_INT__g_s_p_s[36] + 1 * one_over_2q * AUX_INT__g_s_s_s[12] - p_over_q * AUX_INT__h_s_p_s[36];

                    // ( G_022 S_000 | D_110 S_000 )^0_{t} = y * ( G_022 S_000 | P_100 S_000 )^0 + ( F_012 S_000 | P_100 S_000 )^0 - ( H_032 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[73] = etfac[1] * AUX_INT__g_s_p_s[36] + 2 * one_over_2q * AUX_INT__f_s_p_s[24] - p_over_q * AUX_INT__h_s_p_s[51];

                    // ( G_022 S_000 | D_101 S_000 )^0_{t} = z * ( G_022 S_000 | P_100 S_000 )^0 + ( F_021 S_000 | P_100 S_000 )^0 - ( H_023 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[74] = etfac[2] * AUX_INT__g_s_p_s[36] + 2 * one_over_2q * AUX_INT__f_s_p_s[21] - p_over_q * AUX_INT__h_s_p_s[54];

                    // ( G_022 S_000 | D_020 S_000 )^0_{t} = y * ( G_022 S_000 | P_010 S_000 )^0 + ( F_012 S_000 | P_010 S_000 )^0 + ( G_022 S_000 | S_000 S_000 )^0_{e} - ( H_032 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[75] = etfac[1] * AUX_INT__g_s_p_s[37] + 2 * one_over_2q * AUX_INT__f_s_p_s[25] + 1 * one_over_2q * AUX_INT__g_s_s_s[12] - p_over_q * AUX_INT__h_s_p_s[52];

                    // ( G_022 S_000 | D_011 S_000 )^0_{t} = z * ( G_022 S_000 | P_010 S_000 )^0 + ( F_021 S_000 | P_010 S_000 )^0 - ( H_023 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[76] = etfac[2] * AUX_INT__g_s_p_s[37] + 2 * one_over_2q * AUX_INT__f_s_p_s[22] - p_over_q * AUX_INT__h_s_p_s[55];

                    // ( G_022 S_000 | D_002 S_000 )^0_{t} = z * ( G_022 S_000 | P_001 S_000 )^0 + ( F_021 S_000 | P_001 S_000 )^0 + ( G_022 S_000 | S_000 S_000 )^0_{e} - ( H_023 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[77] = etfac[2] * AUX_INT__g_s_p_s[38] + 2 * one_over_2q * AUX_INT__f_s_p_s[23] + 1 * one_over_2q * AUX_INT__g_s_s_s[12] - p_over_q * AUX_INT__h_s_p_s[56];

                    // ( G_013 S_000 | D_200 S_000 )^0_{t} = x * ( G_013 S_000 | P_100 S_000 )^0 + ( G_013 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[78] = etfac[0] * AUX_INT__g_s_p_s[39] + 1 * one_over_2q * AUX_INT__g_s_s_s[13] - p_over_q * AUX_INT__h_s_p_s[39];

                    // ( G_013 S_000 | D_110 S_000 )^0_{t} = y * ( G_013 S_000 | P_100 S_000 )^0 + ( F_003 S_000 | P_100 S_000 )^0 - ( H_023 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[79] = etfac[1] * AUX_INT__g_s_p_s[39] + 1 * one_over_2q * AUX_INT__f_s_p_s[27] - p_over_q * AUX_INT__h_s_p_s[54];

                    // ( G_013 S_000 | D_101 S_000 )^0_{t} = z * ( G_013 S_000 | P_100 S_000 )^0 + ( F_012 S_000 | P_100 S_000 )^0 - ( H_014 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[80] = etfac[2] * AUX_INT__g_s_p_s[39] + 3 * one_over_2q * AUX_INT__f_s_p_s[24] - p_over_q * AUX_INT__h_s_p_s[57];

                    // ( G_013 S_000 | D_020 S_000 )^0_{t} = y * ( G_013 S_000 | P_010 S_000 )^0 + ( F_003 S_000 | P_010 S_000 )^0 + ( G_013 S_000 | S_000 S_000 )^0_{e} - ( H_023 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[81] = etfac[1] * AUX_INT__g_s_p_s[40] + 1 * one_over_2q * AUX_INT__f_s_p_s[28] + 1 * one_over_2q * AUX_INT__g_s_s_s[13] - p_over_q * AUX_INT__h_s_p_s[55];

                    // ( G_013 S_000 | D_011 S_000 )^0_{t} = z * ( G_013 S_000 | P_010 S_000 )^0 + ( F_012 S_000 | P_010 S_000 )^0 - ( H_014 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[82] = etfac[2] * AUX_INT__g_s_p_s[40] + 3 * one_over_2q * AUX_INT__f_s_p_s[25] - p_over_q * AUX_INT__h_s_p_s[58];

                    // ( G_013 S_000 | D_002 S_000 )^0_{t} = z * ( G_013 S_000 | P_001 S_000 )^0 + ( F_012 S_000 | P_001 S_000 )^0 + ( G_013 S_000 | S_000 S_000 )^0_{e} - ( H_014 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[83] = etfac[2] * AUX_INT__g_s_p_s[41] + 3 * one_over_2q * AUX_INT__f_s_p_s[26] + 1 * one_over_2q * AUX_INT__g_s_s_s[13] - p_over_q * AUX_INT__h_s_p_s[59];

                    // ( G_004 S_000 | D_200 S_000 )^0_{t} = x * ( G_004 S_000 | P_100 S_000 )^0 + ( G_004 S_000 | S_000 S_000 )^0_{e} - ( H_104 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[84] = etfac[0] * AUX_INT__g_s_p_s[42] + 1 * one_over_2q * AUX_INT__g_s_s_s[14] - p_over_q * AUX_INT__h_s_p_s[42];

                    // ( G_004 S_000 | D_110 S_000 )^0_{t} = y * ( G_004 S_000 | P_100 S_000 )^0 - ( H_014 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[85] = etfac[1] * AUX_INT__g_s_p_s[42] - p_over_q * AUX_INT__h_s_p_s[57];

                    // ( G_004 S_000 | D_101 S_000 )^0_{t} = z * ( G_004 S_000 | P_100 S_000 )^0 + ( F_003 S_000 | P_100 S_000 )^0 - ( H_005 S_000 | P_100 S_000 )^0
                    AUX_INT__g_s_d_s[86] = etfac[2] * AUX_INT__g_s_p_s[42] + 4 * one_over_2q * AUX_INT__f_s_p_s[27] - p_over_q * AUX_INT__h_s_p_s[60];

                    // ( G_004 S_000 | D_020 S_000 )^0_{t} = y * ( G_004 S_000 | P_010 S_000 )^0 + ( G_004 S_000 | S_000 S_000 )^0_{e} - ( H_014 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[87] = etfac[1] * AUX_INT__g_s_p_s[43] + 1 * one_over_2q * AUX_INT__g_s_s_s[14] - p_over_q * AUX_INT__h_s_p_s[58];

                    // ( G_004 S_000 | D_011 S_000 )^0_{t} = z * ( G_004 S_000 | P_010 S_000 )^0 + ( F_003 S_000 | P_010 S_000 )^0 - ( H_005 S_000 | P_010 S_000 )^0
                    AUX_INT__g_s_d_s[88] = etfac[2] * AUX_INT__g_s_p_s[43] + 4 * one_over_2q * AUX_INT__f_s_p_s[28] - p_over_q * AUX_INT__h_s_p_s[61];

                    // ( G_004 S_000 | D_002 S_000 )^0_{t} = z * ( G_004 S_000 | P_001 S_000 )^0 + ( F_003 S_000 | P_001 S_000 )^0 + ( G_004 S_000 | S_000 S_000 )^0_{e} - ( H_005 S_000 | P_001 S_000 )^0
                    AUX_INT__g_s_d_s[89] = etfac[2] * AUX_INT__g_s_p_s[44] + 4 * one_over_2q * AUX_INT__f_s_p_s[29] + 1 * one_over_2q * AUX_INT__g_s_s_s[14] - p_over_q * AUX_INT__h_s_p_s[62];


                    // Accumulating in contracted workspace
                    for(n = 0; n < 18; n++)
                        PRIM_INT__d_s_p_s[n] += AUX_INT__d_s_p_s[n];

                    // Accumulating in contracted workspace
                    for(n = 0; n < 36; n++)
                        PRIM_INT__d_s_d_s[n] += AUX_INT__d_s_d_s[n];

                    // Accumulating in contracted workspace
                    for(n = 0; n < 30; n++)
                        PRIM_INT__f_s_p_s[n] += AUX_INT__f_s_p_s[n];

                    // Accumulating in contracted workspace
                    for(n = 0; n < 60; n++)
                        PRIM_INT__f_s_d_s[n] += AUX_INT__f_s_d_s[n];

                    // Accumulating in contracted workspace
                    for(n = 0; n < 45; n++)
                        PRIM_INT__g_s_p_s[n] += AUX_INT__g_s_p_s[n];

                    // Accumulating in contracted workspace
                    for(n = 0; n < 90; n++)
                        PRIM_INT__g_s_d_s[n] += AUX_INT__g_s_d_s[n];

                 }
            }
        }
    }



    //////////////////////////////////////////////
    // Contracted integrals: Horizontal recurrance
    // Bra part
    // Steps: 79
    //////////////////////////////////////////////

    for(abcd = 0; abcd < nshell1234; ++abcd)
    {
        // form INT__d_d_p_s
        for(iket = 0; iket < 3; ++iket)
        {
            // (D_200 P_100| = (F_300 S_000|_{t} + x_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_p_s_s_p = INT__f_s_p_s[abcd * 30 + 0 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 0 * 3 + iket] );

            // (D_200 P_010| = (F_210 S_000|_{t} + y_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_s_p_s_p = INT__f_s_p_s[abcd * 30 + 1 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 0 * 3 + iket] );

            // (D_200 P_001| = (F_201 S_000|_{t} + z_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_s_s_p_p = INT__f_s_p_s[abcd * 30 + 2 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 0 * 3 + iket] );

            // (D_110 P_100| = (F_210 S_000|_{t} + x_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_p_s_s_p = INT__f_s_p_s[abcd * 30 + 1 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 1 * 3 + iket] );

            // (D_110 P_010| = (F_120 S_000|_{t} + y_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_s_p_s_p = INT__f_s_p_s[abcd * 30 + 3 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 1 * 3 + iket] );

            // (D_110 P_001| = (F_111 S_000|_{t} + z_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_s_s_p_p = INT__f_s_p_s[abcd * 30 + 4 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 1 * 3 + iket] );

            // (D_101 P_100| = (F_201 S_000|_{t} + x_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_p_s_s_p = INT__f_s_p_s[abcd * 30 + 2 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 2 * 3 + iket] );

            // (D_101 P_010| = (F_111 S_000|_{t} + y_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_s_p_s_p = INT__f_s_p_s[abcd * 30 + 4 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 2 * 3 + iket] );

            // (D_101 P_001| = (F_102 S_000|_{t} + z_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_s_s_p_p = INT__f_s_p_s[abcd * 30 + 5 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 2 * 3 + iket] );

            // (D_020 P_100| = (F_120 S_000|_{t} + x_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_p_s_s_p = INT__f_s_p_s[abcd * 30 + 3 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 3 * 3 + iket] );

            // (D_020 P_010| = (F_030 S_000|_{t} + y_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_s_p_s_p = INT__f_s_p_s[abcd * 30 + 6 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 3 * 3 + iket] );

            // (D_020 P_001| = (F_021 S_000|_{t} + z_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_s_s_p_p = INT__f_s_p_s[abcd * 30 + 7 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 3 * 3 + iket] );

            // (D_011 P_100| = (F_111 S_000|_{t} + x_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_p_s_s_p = INT__f_s_p_s[abcd * 30 + 4 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 4 * 3 + iket] );

            // (D_011 P_010| = (F_021 S_000|_{t} + y_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_s_p_s_p = INT__f_s_p_s[abcd * 30 + 7 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 4 * 3 + iket] );

            // (D_011 P_001| = (F_012 S_000|_{t} + z_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_s_s_p_p = INT__f_s_p_s[abcd * 30 + 8 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 4 * 3 + iket] );

            // (D_002 P_100| = (F_102 S_000|_{t} + x_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_p_s_s_p = INT__f_s_p_s[abcd * 30 + 5 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 5 * 3 + iket] );

            // (D_002 P_010| = (F_012 S_000|_{t} + y_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_s_p_s_p = INT__f_s_p_s[abcd * 30 + 8 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 5 * 3 + iket] );

            // (D_002 P_001| = (F_003 S_000|_{t} + z_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_s_s_p_p = INT__f_s_p_s[abcd * 30 + 9 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 5 * 3 + iket] );

            // (F_300 P_100| = (G_400 S_000|_{t} + x_ab * (F_300 S_000|_{t}
            const double Q_f_s_s_p_s_s_p = INT__g_s_p_s[abcd * 45 + 0 * 3 + iket] + ( AB_x[abcd] * INT__f_s_p_s[abcd * 30 + 0 * 3 + iket] );

            // (F_210 P_100| = (G_310 S_000|_{t} + x_ab * (F_210 S_000|_{t}
            const double Q_d_p_s_p_s_s_p = INT__g_s_p_s[abcd * 45 + 1 * 3 + iket] + ( AB_x[abcd] * INT__f_s_p_s[abcd * 30 + 1 * 3 + iket] );

            // (F_210 P_010| = (G_220 S_000|_{t} + y_ab * (F_210 S_000|_{t}
            const double Q_d_p_s_s_p_s_p = INT__g_s_p_s[abcd * 45 + 3 * 3 + iket] + ( AB_y[abcd] * INT__f_s_p_s[abcd * 30 + 1 * 3 + iket] );

            // (F_201 P_100| = (G_301 S_000|_{t} + x_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_p_s_s_p = INT__g_s_p_s[abcd * 45 + 2 * 3 + iket] + ( AB_x[abcd] * INT__f_s_p_s[abcd * 30 + 2 * 3 + iket] );

            // (F_201 P_010| = (G_211 S_000|_{t} + y_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_s_p_s_p = INT__g_s_p_s[abcd * 45 + 4 * 3 + iket] + ( AB_y[abcd] * INT__f_s_p_s[abcd * 30 + 2 * 3 + iket] );

            // (F_201 P_001| = (G_202 S_000|_{t} + z_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_s_s_p_p = INT__g_s_p_s[abcd * 45 + 5 * 3 + iket] + ( AB_z[abcd] * INT__f_s_p_s[abcd * 30 + 2 * 3 + iket] );

            // (F_120 P_100| = (G_220 S_000|_{t} + x_ab * (F_120 S_000|_{t}
            const double Q_p_d_s_p_s_s_p = INT__g_s_p_s[abcd * 45 + 3 * 3 + iket] + ( AB_x[abcd] * INT__f_s_p_s[abcd * 30 + 3 * 3 + iket] );

            // (F_120 P_010| = (G_130 S_000|_{t} + y_ab * (F_120 S_000|_{t}
            const double Q_p_d_s_s_p_s_p = INT__g_s_p_s[abcd * 45 + 6 * 3 + iket] + ( AB_y[abcd] * INT__f_s_p_s[abcd * 30 + 3 * 3 + iket] );

            // (F_111 P_100| = (G_211 S_000|_{t} + x_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_p_s_s_p = INT__g_s_p_s[abcd * 45 + 4 * 3 + iket] + ( AB_x[abcd] * INT__f_s_p_s[abcd * 30 + 4 * 3 + iket] );

            // (F_111 P_010| = (G_121 S_000|_{t} + y_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_s_p_s_p = INT__g_s_p_s[abcd * 45 + 7 * 3 + iket] + ( AB_y[abcd] * INT__f_s_p_s[abcd * 30 + 4 * 3 + iket] );

            // (F_111 P_001| = (G_112 S_000|_{t} + z_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_s_s_p_p = INT__g_s_p_s[abcd * 45 + 8 * 3 + iket] + ( AB_z[abcd] * INT__f_s_p_s[abcd * 30 + 4 * 3 + iket] );

            // (F_102 P_100| = (G_202 S_000|_{t} + x_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_p_s_s_p = INT__g_s_p_s[abcd * 45 + 5 * 3 + iket] + ( AB_x[abcd] * INT__f_s_p_s[abcd * 30 + 5 * 3 + iket] );

            // (F_102 P_010| = (G_112 S_000|_{t} + y_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_s_p_s_p = INT__g_s_p_s[abcd * 45 + 8 * 3 + iket] + ( AB_y[abcd] * INT__f_s_p_s[abcd * 30 + 5 * 3 + iket] );

            // (F_102 P_001| = (G_103 S_000|_{t} + z_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_s_s_p_p = INT__g_s_p_s[abcd * 45 + 9 * 3 + iket] + ( AB_z[abcd] * INT__f_s_p_s[abcd * 30 + 5 * 3 + iket] );

            // (F_030 P_100| = (G_130 S_000|_{t} + x_ab * (F_030 S_000|_{t}
            const double Q_s_f_s_p_s_s_p = INT__g_s_p_s[abcd * 45 + 6 * 3 + iket] + ( AB_x[abcd] * INT__f_s_p_s[abcd * 30 + 6 * 3 + iket] );

            // (F_030 P_010| = (G_040 S_000|_{t} + y_ab * (F_030 S_000|_{t}
            const double Q_s_f_s_s_p_s_p = INT__g_s_p_s[abcd * 45 + 10 * 3 + iket] + ( AB_y[abcd] * INT__f_s_p_s[abcd * 30 + 6 * 3 + iket] );

            // (F_021 P_100| = (G_121 S_000|_{t} + x_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_p_s_s_p = INT__g_s_p_s[abcd * 45 + 7 * 3 + iket] + ( AB_x[abcd] * INT__f_s_p_s[abcd * 30 + 7 * 3 + iket] );

            // (F_021 P_010| = (G_031 S_000|_{t} + y_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_s_p_s_p = INT__g_s_p_s[abcd * 45 + 11 * 3 + iket] + ( AB_y[abcd] * INT__f_s_p_s[abcd * 30 + 7 * 3 + iket] );

            // (F_021 P_001| = (G_022 S_000|_{t} + z_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_s_s_p_p = INT__g_s_p_s[abcd * 45 + 12 * 3 + iket] + ( AB_z[abcd] * INT__f_s_p_s[abcd * 30 + 7 * 3 + iket] );

            // (F_012 P_100| = (G_112 S_000|_{t} + x_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_p_s_s_p = INT__g_s_p_s[abcd * 45 + 8 * 3 + iket] + ( AB_x[abcd] * INT__f_s_p_s[abcd * 30 + 8 * 3 + iket] );

            // (F_012 P_010| = (G_022 S_000|_{t} + y_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_s_p_s_p = INT__g_s_p_s[abcd * 45 + 12 * 3 + iket] + ( AB_y[abcd] * INT__f_s_p_s[abcd * 30 + 8 * 3 + iket] );

            // (F_012 P_001| = (G_013 S_000|_{t} + z_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_s_s_p_p = INT__g_s_p_s[abcd * 45 + 13 * 3 + iket] + ( AB_z[abcd] * INT__f_s_p_s[abcd * 30 + 8 * 3 + iket] );

            // (F_003 P_100| = (G_103 S_000|_{t} + x_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_p_s_s_p = INT__g_s_p_s[abcd * 45 + 9 * 3 + iket] + ( AB_x[abcd] * INT__f_s_p_s[abcd * 30 + 9 * 3 + iket] );

            // (F_003 P_010| = (G_013 S_000|_{t} + y_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_s_p_s_p = INT__g_s_p_s[abcd * 45 + 13 * 3 + iket] + ( AB_y[abcd] * INT__f_s_p_s[abcd * 30 + 9 * 3 + iket] );

            // (F_003 P_001| = (G_004 S_000|_{t} + z_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_s_s_p_p = INT__g_s_p_s[abcd * 45 + 14 * 3 + iket] + ( AB_z[abcd] * INT__f_s_p_s[abcd * 30 + 9 * 3 + iket] );

            // (D_200 D_200|_{i} = (F_300 P_100| + x_ab * (D_200 P_100|
            INT__d_d_p_s[abcd * 108 + 0 * 3 + iket] = Q_f_s_s_p_s_s_p + ( AB_x[abcd] * Q_d_s_s_p_s_s_p );

            // (D_200 D_110|_{i} = (F_210 P_100| + y_ab * (D_200 P_100|
            INT__d_d_p_s[abcd * 108 + 1 * 3 + iket] = Q_d_p_s_p_s_s_p + ( AB_y[abcd] * Q_d_s_s_p_s_s_p );

            // (D_200 D_101|_{i} = (F_201 P_100| + z_ab * (D_200 P_100|
            INT__d_d_p_s[abcd * 108 + 2 * 3 + iket] = Q_d_s_p_p_s_s_p + ( AB_z[abcd] * Q_d_s_s_p_s_s_p );

            // (D_200 D_020|_{i} = (F_210 P_010| + y_ab * (D_200 P_010|
            INT__d_d_p_s[abcd * 108 + 3 * 3 + iket] = Q_d_p_s_s_p_s_p + ( AB_y[abcd] * Q_d_s_s_s_p_s_p );

            // (D_200 D_011|_{i} = (F_201 P_010| + z_ab * (D_200 P_010|
            INT__d_d_p_s[abcd * 108 + 4 * 3 + iket] = Q_d_s_p_s_p_s_p + ( AB_z[abcd] * Q_d_s_s_s_p_s_p );

            // (D_200 D_002|_{i} = (F_201 P_001| + z_ab * (D_200 P_001|
            INT__d_d_p_s[abcd * 108 + 5 * 3 + iket] = Q_d_s_p_s_s_p_p + ( AB_z[abcd] * Q_d_s_s_s_s_p_p );

            // (D_110 D_200|_{i} = (F_210 P_100| + x_ab * (D_110 P_100|
            INT__d_d_p_s[abcd * 108 + 6 * 3 + iket] = Q_d_p_s_p_s_s_p + ( AB_x[abcd] * Q_p_p_s_p_s_s_p );

            // (D_110 D_110|_{i} = (F_120 P_100| + y_ab * (D_110 P_100|
            INT__d_d_p_s[abcd * 108 + 7 * 3 + iket] = Q_p_d_s_p_s_s_p + ( AB_y[abcd] * Q_p_p_s_p_s_s_p );

            // (D_110 D_101|_{i} = (F_111 P_100| + z_ab * (D_110 P_100|
            INT__d_d_p_s[abcd * 108 + 8 * 3 + iket] = Q_p_p_p_p_s_s_p + ( AB_z[abcd] * Q_p_p_s_p_s_s_p );

            // (D_110 D_020|_{i} = (F_120 P_010| + y_ab * (D_110 P_010|
            INT__d_d_p_s[abcd * 108 + 9 * 3 + iket] = Q_p_d_s_s_p_s_p + ( AB_y[abcd] * Q_p_p_s_s_p_s_p );

            // (D_110 D_011|_{i} = (F_111 P_010| + z_ab * (D_110 P_010|
            INT__d_d_p_s[abcd * 108 + 10 * 3 + iket] = Q_p_p_p_s_p_s_p + ( AB_z[abcd] * Q_p_p_s_s_p_s_p );

            // (D_110 D_002|_{i} = (F_111 P_001| + z_ab * (D_110 P_001|
            INT__d_d_p_s[abcd * 108 + 11 * 3 + iket] = Q_p_p_p_s_s_p_p + ( AB_z[abcd] * Q_p_p_s_s_s_p_p );

            // (D_101 D_200|_{i} = (F_201 P_100| + x_ab * (D_101 P_100|
            INT__d_d_p_s[abcd * 108 + 12 * 3 + iket] = Q_d_s_p_p_s_s_p + ( AB_x[abcd] * Q_p_s_p_p_s_s_p );

            // (D_101 D_110|_{i} = (F_111 P_100| + y_ab * (D_101 P_100|
            INT__d_d_p_s[abcd * 108 + 13 * 3 + iket] = Q_p_p_p_p_s_s_p + ( AB_y[abcd] * Q_p_s_p_p_s_s_p );

            // (D_101 D_101|_{i} = (F_102 P_100| + z_ab * (D_101 P_100|
            INT__d_d_p_s[abcd * 108 + 14 * 3 + iket] = Q_p_s_d_p_s_s_p + ( AB_z[abcd] * Q_p_s_p_p_s_s_p );

            // (D_101 D_020|_{i} = (F_111 P_010| + y_ab * (D_101 P_010|
            INT__d_d_p_s[abcd * 108 + 15 * 3 + iket] = Q_p_p_p_s_p_s_p + ( AB_y[abcd] * Q_p_s_p_s_p_s_p );

            // (D_101 D_011|_{i} = (F_102 P_010| + z_ab * (D_101 P_010|
            INT__d_d_p_s[abcd * 108 + 16 * 3 + iket] = Q_p_s_d_s_p_s_p + ( AB_z[abcd] * Q_p_s_p_s_p_s_p );

            // (D_101 D_002|_{i} = (F_102 P_001| + z_ab * (D_101 P_001|
            INT__d_d_p_s[abcd * 108 + 17 * 3 + iket] = Q_p_s_d_s_s_p_p + ( AB_z[abcd] * Q_p_s_p_s_s_p_p );

            // (D_020 D_200|_{i} = (F_120 P_100| + x_ab * (D_020 P_100|
            INT__d_d_p_s[abcd * 108 + 18 * 3 + iket] = Q_p_d_s_p_s_s_p + ( AB_x[abcd] * Q_s_d_s_p_s_s_p );

            // (D_020 D_110|_{i} = (F_030 P_100| + y_ab * (D_020 P_100|
            INT__d_d_p_s[abcd * 108 + 19 * 3 + iket] = Q_s_f_s_p_s_s_p + ( AB_y[abcd] * Q_s_d_s_p_s_s_p );

            // (D_020 D_101|_{i} = (F_021 P_100| + z_ab * (D_020 P_100|
            INT__d_d_p_s[abcd * 108 + 20 * 3 + iket] = Q_s_d_p_p_s_s_p + ( AB_z[abcd] * Q_s_d_s_p_s_s_p );

            // (D_020 D_020|_{i} = (F_030 P_010| + y_ab * (D_020 P_010|
            INT__d_d_p_s[abcd * 108 + 21 * 3 + iket] = Q_s_f_s_s_p_s_p + ( AB_y[abcd] * Q_s_d_s_s_p_s_p );

            // (D_020 D_011|_{i} = (F_021 P_010| + z_ab * (D_020 P_010|
            INT__d_d_p_s[abcd * 108 + 22 * 3 + iket] = Q_s_d_p_s_p_s_p + ( AB_z[abcd] * Q_s_d_s_s_p_s_p );

            // (D_020 D_002|_{i} = (F_021 P_001| + z_ab * (D_020 P_001|
            INT__d_d_p_s[abcd * 108 + 23 * 3 + iket] = Q_s_d_p_s_s_p_p + ( AB_z[abcd] * Q_s_d_s_s_s_p_p );

            // (D_011 D_200|_{i} = (F_111 P_100| + x_ab * (D_011 P_100|
            INT__d_d_p_s[abcd * 108 + 24 * 3 + iket] = Q_p_p_p_p_s_s_p + ( AB_x[abcd] * Q_s_p_p_p_s_s_p );

            // (D_011 D_110|_{i} = (F_021 P_100| + y_ab * (D_011 P_100|
            INT__d_d_p_s[abcd * 108 + 25 * 3 + iket] = Q_s_d_p_p_s_s_p + ( AB_y[abcd] * Q_s_p_p_p_s_s_p );

            // (D_011 D_101|_{i} = (F_012 P_100| + z_ab * (D_011 P_100|
            INT__d_d_p_s[abcd * 108 + 26 * 3 + iket] = Q_s_p_d_p_s_s_p + ( AB_z[abcd] * Q_s_p_p_p_s_s_p );

            // (D_011 D_020|_{i} = (F_021 P_010| + y_ab * (D_011 P_010|
            INT__d_d_p_s[abcd * 108 + 27 * 3 + iket] = Q_s_d_p_s_p_s_p + ( AB_y[abcd] * Q_s_p_p_s_p_s_p );

            // (D_011 D_011|_{i} = (F_012 P_010| + z_ab * (D_011 P_010|
            INT__d_d_p_s[abcd * 108 + 28 * 3 + iket] = Q_s_p_d_s_p_s_p + ( AB_z[abcd] * Q_s_p_p_s_p_s_p );

            // (D_011 D_002|_{i} = (F_012 P_001| + z_ab * (D_011 P_001|
            INT__d_d_p_s[abcd * 108 + 29 * 3 + iket] = Q_s_p_d_s_s_p_p + ( AB_z[abcd] * Q_s_p_p_s_s_p_p );

            // (D_002 D_200|_{i} = (F_102 P_100| + x_ab * (D_002 P_100|
            INT__d_d_p_s[abcd * 108 + 30 * 3 + iket] = Q_p_s_d_p_s_s_p + ( AB_x[abcd] * Q_s_s_d_p_s_s_p );

            // (D_002 D_110|_{i} = (F_012 P_100| + y_ab * (D_002 P_100|
            INT__d_d_p_s[abcd * 108 + 31 * 3 + iket] = Q_s_p_d_p_s_s_p + ( AB_y[abcd] * Q_s_s_d_p_s_s_p );

            // (D_002 D_101|_{i} = (F_003 P_100| + z_ab * (D_002 P_100|
            INT__d_d_p_s[abcd * 108 + 32 * 3 + iket] = Q_s_s_f_p_s_s_p + ( AB_z[abcd] * Q_s_s_d_p_s_s_p );

            // (D_002 D_020|_{i} = (F_012 P_010| + y_ab * (D_002 P_010|
            INT__d_d_p_s[abcd * 108 + 33 * 3 + iket] = Q_s_p_d_s_p_s_p + ( AB_y[abcd] * Q_s_s_d_s_p_s_p );

            // (D_002 D_011|_{i} = (F_003 P_010| + z_ab * (D_002 P_010|
            INT__d_d_p_s[abcd * 108 + 34 * 3 + iket] = Q_s_s_f_s_p_s_p + ( AB_z[abcd] * Q_s_s_d_s_p_s_p );

            // (D_002 D_002|_{i} = (F_003 P_001| + z_ab * (D_002 P_001|
            INT__d_d_p_s[abcd * 108 + 35 * 3 + iket] = Q_s_s_f_s_s_p_p + ( AB_z[abcd] * Q_s_s_d_s_s_p_p );

        }

        // form INT__d_d_d_s
        for(iket = 0; iket < 6; ++iket)
        {
            // (D_200 P_100| = (F_300 S_000|_{t} + x_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_p_s_s_d = INT__f_s_d_s[abcd * 60 + 0 * 6 + iket] + ( AB_x[abcd] * INT__d_s_d_s[abcd * 36 + 0 * 6 + iket] );

            // (D_200 P_010| = (F_210 S_000|_{t} + y_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_s_p_s_d = INT__f_s_d_s[abcd * 60 + 1 * 6 + iket] + ( AB_y[abcd] * INT__d_s_d_s[abcd * 36 + 0 * 6 + iket] );

            // (D_200 P_001| = (F_201 S_000|_{t} + z_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_s_s_p_d = INT__f_s_d_s[abcd * 60 + 2 * 6 + iket] + ( AB_z[abcd] * INT__d_s_d_s[abcd * 36 + 0 * 6 + iket] );

            // (D_110 P_100| = (F_210 S_000|_{t} + x_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_p_s_s_d = INT__f_s_d_s[abcd * 60 + 1 * 6 + iket] + ( AB_x[abcd] * INT__d_s_d_s[abcd * 36 + 1 * 6 + iket] );

            // (D_110 P_010| = (F_120 S_000|_{t} + y_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_s_p_s_d = INT__f_s_d_s[abcd * 60 + 3 * 6 + iket] + ( AB_y[abcd] * INT__d_s_d_s[abcd * 36 + 1 * 6 + iket] );

            // (D_110 P_001| = (F_111 S_000|_{t} + z_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_s_s_p_d = INT__f_s_d_s[abcd * 60 + 4 * 6 + iket] + ( AB_z[abcd] * INT__d_s_d_s[abcd * 36 + 1 * 6 + iket] );

            // (D_101 P_100| = (F_201 S_000|_{t} + x_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_p_s_s_d = INT__f_s_d_s[abcd * 60 + 2 * 6 + iket] + ( AB_x[abcd] * INT__d_s_d_s[abcd * 36 + 2 * 6 + iket] );

            // (D_101 P_010| = (F_111 S_000|_{t} + y_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_s_p_s_d = INT__f_s_d_s[abcd * 60 + 4 * 6 + iket] + ( AB_y[abcd] * INT__d_s_d_s[abcd * 36 + 2 * 6 + iket] );

            // (D_101 P_001| = (F_102 S_000|_{t} + z_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_s_s_p_d = INT__f_s_d_s[abcd * 60 + 5 * 6 + iket] + ( AB_z[abcd] * INT__d_s_d_s[abcd * 36 + 2 * 6 + iket] );

            // (D_020 P_100| = (F_120 S_000|_{t} + x_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_p_s_s_d = INT__f_s_d_s[abcd * 60 + 3 * 6 + iket] + ( AB_x[abcd] * INT__d_s_d_s[abcd * 36 + 3 * 6 + iket] );

            // (D_020 P_010| = (F_030 S_000|_{t} + y_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_s_p_s_d = INT__f_s_d_s[abcd * 60 + 6 * 6 + iket] + ( AB_y[abcd] * INT__d_s_d_s[abcd * 36 + 3 * 6 + iket] );

            // (D_020 P_001| = (F_021 S_000|_{t} + z_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_s_s_p_d = INT__f_s_d_s[abcd * 60 + 7 * 6 + iket] + ( AB_z[abcd] * INT__d_s_d_s[abcd * 36 + 3 * 6 + iket] );

            // (D_011 P_100| = (F_111 S_000|_{t} + x_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_p_s_s_d = INT__f_s_d_s[abcd * 60 + 4 * 6 + iket] + ( AB_x[abcd] * INT__d_s_d_s[abcd * 36 + 4 * 6 + iket] );

            // (D_011 P_010| = (F_021 S_000|_{t} + y_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_s_p_s_d = INT__f_s_d_s[abcd * 60 + 7 * 6 + iket] + ( AB_y[abcd] * INT__d_s_d_s[abcd * 36 + 4 * 6 + iket] );

            // (D_011 P_001| = (F_012 S_000|_{t} + z_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_s_s_p_d = INT__f_s_d_s[abcd * 60 + 8 * 6 + iket] + ( AB_z[abcd] * INT__d_s_d_s[abcd * 36 + 4 * 6 + iket] );

            // (D_002 P_100| = (F_102 S_000|_{t} + x_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_p_s_s_d = INT__f_s_d_s[abcd * 60 + 5 * 6 + iket] + ( AB_x[abcd] * INT__d_s_d_s[abcd * 36 + 5 * 6 + iket] );

            // (D_002 P_010| = (F_012 S_000|_{t} + y_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_s_p_s_d = INT__f_s_d_s[abcd * 60 + 8 * 6 + iket] + ( AB_y[abcd] * INT__d_s_d_s[abcd * 36 + 5 * 6 + iket] );

            // (D_002 P_001| = (F_003 S_000|_{t} + z_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_s_s_p_d = INT__f_s_d_s[abcd * 60 + 9 * 6 + iket] + ( AB_z[abcd] * INT__d_s_d_s[abcd * 36 + 5 * 6 + iket] );

            // (F_300 P_100| = (G_400 S_000|_{t} + x_ab * (F_300 S_000|_{t}
            const double Q_f_s_s_p_s_s_d = INT__g_s_d_s[abcd * 90 + 0 * 6 + iket] + ( AB_x[abcd] * INT__f_s_d_s[abcd * 60 + 0 * 6 + iket] );

            // (F_210 P_100| = (G_310 S_000|_{t} + x_ab * (F_210 S_000|_{t}
            const double Q_d_p_s_p_s_s_d = INT__g_s_d_s[abcd * 90 + 1 * 6 + iket] + ( AB_x[abcd] * INT__f_s_d_s[abcd * 60 + 1 * 6 + iket] );

            // (F_210 P_010| = (G_220 S_000|_{t} + y_ab * (F_210 S_000|_{t}
            const double Q_d_p_s_s_p_s_d = INT__g_s_d_s[abcd * 90 + 3 * 6 + iket] + ( AB_y[abcd] * INT__f_s_d_s[abcd * 60 + 1 * 6 + iket] );

            // (F_201 P_100| = (G_301 S_000|_{t} + x_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_p_s_s_d = INT__g_s_d_s[abcd * 90 + 2 * 6 + iket] + ( AB_x[abcd] * INT__f_s_d_s[abcd * 60 + 2 * 6 + iket] );

            // (F_201 P_010| = (G_211 S_000|_{t} + y_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_s_p_s_d = INT__g_s_d_s[abcd * 90 + 4 * 6 + iket] + ( AB_y[abcd] * INT__f_s_d_s[abcd * 60 + 2 * 6 + iket] );

            // (F_201 P_001| = (G_202 S_000|_{t} + z_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_s_s_p_d = INT__g_s_d_s[abcd * 90 + 5 * 6 + iket] + ( AB_z[abcd] * INT__f_s_d_s[abcd * 60 + 2 * 6 + iket] );

            // (F_120 P_100| = (G_220 S_000|_{t} + x_ab * (F_120 S_000|_{t}
            const double Q_p_d_s_p_s_s_d = INT__g_s_d_s[abcd * 90 + 3 * 6 + iket] + ( AB_x[abcd] * INT__f_s_d_s[abcd * 60 + 3 * 6 + iket] );

            // (F_120 P_010| = (G_130 S_000|_{t} + y_ab * (F_120 S_000|_{t}
            const double Q_p_d_s_s_p_s_d = INT__g_s_d_s[abcd * 90 + 6 * 6 + iket] + ( AB_y[abcd] * INT__f_s_d_s[abcd * 60 + 3 * 6 + iket] );

            // (F_111 P_100| = (G_211 S_000|_{t} + x_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_p_s_s_d = INT__g_s_d_s[abcd * 90 + 4 * 6 + iket] + ( AB_x[abcd] * INT__f_s_d_s[abcd * 60 + 4 * 6 + iket] );

            // (F_111 P_010| = (G_121 S_000|_{t} + y_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_s_p_s_d = INT__g_s_d_s[abcd * 90 + 7 * 6 + iket] + ( AB_y[abcd] * INT__f_s_d_s[abcd * 60 + 4 * 6 + iket] );

            // (F_111 P_001| = (G_112 S_000|_{t} + z_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_s_s_p_d = INT__g_s_d_s[abcd * 90 + 8 * 6 + iket] + ( AB_z[abcd] * INT__f_s_d_s[abcd * 60 + 4 * 6 + iket] );

            // (F_102 P_100| = (G_202 S_000|_{t} + x_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_p_s_s_d = INT__g_s_d_s[abcd * 90 + 5 * 6 + iket] + ( AB_x[abcd] * INT__f_s_d_s[abcd * 60 + 5 * 6 + iket] );

            // (F_102 P_010| = (G_112 S_000|_{t} + y_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_s_p_s_d = INT__g_s_d_s[abcd * 90 + 8 * 6 + iket] + ( AB_y[abcd] * INT__f_s_d_s[abcd * 60 + 5 * 6 + iket] );

            // (F_102 P_001| = (G_103 S_000|_{t} + z_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_s_s_p_d = INT__g_s_d_s[abcd * 90 + 9 * 6 + iket] + ( AB_z[abcd] * INT__f_s_d_s[abcd * 60 + 5 * 6 + iket] );

            // (F_030 P_100| = (G_130 S_000|_{t} + x_ab * (F_030 S_000|_{t}
            const double Q_s_f_s_p_s_s_d = INT__g_s_d_s[abcd * 90 + 6 * 6 + iket] + ( AB_x[abcd] * INT__f_s_d_s[abcd * 60 + 6 * 6 + iket] );

            // (F_030 P_010| = (G_040 S_000|_{t} + y_ab * (F_030 S_000|_{t}
            const double Q_s_f_s_s_p_s_d = INT__g_s_d_s[abcd * 90 + 10 * 6 + iket] + ( AB_y[abcd] * INT__f_s_d_s[abcd * 60 + 6 * 6 + iket] );

            // (F_021 P_100| = (G_121 S_000|_{t} + x_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_p_s_s_d = INT__g_s_d_s[abcd * 90 + 7 * 6 + iket] + ( AB_x[abcd] * INT__f_s_d_s[abcd * 60 + 7 * 6 + iket] );

            // (F_021 P_010| = (G_031 S_000|_{t} + y_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_s_p_s_d = INT__g_s_d_s[abcd * 90 + 11 * 6 + iket] + ( AB_y[abcd] * INT__f_s_d_s[abcd * 60 + 7 * 6 + iket] );

            // (F_021 P_001| = (G_022 S_000|_{t} + z_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_s_s_p_d = INT__g_s_d_s[abcd * 90 + 12 * 6 + iket] + ( AB_z[abcd] * INT__f_s_d_s[abcd * 60 + 7 * 6 + iket] );

            // (F_012 P_100| = (G_112 S_000|_{t} + x_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_p_s_s_d = INT__g_s_d_s[abcd * 90 + 8 * 6 + iket] + ( AB_x[abcd] * INT__f_s_d_s[abcd * 60 + 8 * 6 + iket] );

            // (F_012 P_010| = (G_022 S_000|_{t} + y_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_s_p_s_d = INT__g_s_d_s[abcd * 90 + 12 * 6 + iket] + ( AB_y[abcd] * INT__f_s_d_s[abcd * 60 + 8 * 6 + iket] );

            // (F_012 P_001| = (G_013 S_000|_{t} + z_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_s_s_p_d = INT__g_s_d_s[abcd * 90 + 13 * 6 + iket] + ( AB_z[abcd] * INT__f_s_d_s[abcd * 60 + 8 * 6 + iket] );

            // (F_003 P_100| = (G_103 S_000|_{t} + x_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_p_s_s_d = INT__g_s_d_s[abcd * 90 + 9 * 6 + iket] + ( AB_x[abcd] * INT__f_s_d_s[abcd * 60 + 9 * 6 + iket] );

            // (F_003 P_010| = (G_013 S_000|_{t} + y_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_s_p_s_d = INT__g_s_d_s[abcd * 90 + 13 * 6 + iket] + ( AB_y[abcd] * INT__f_s_d_s[abcd * 60 + 9 * 6 + iket] );

            // (F_003 P_001| = (G_004 S_000|_{t} + z_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_s_s_p_d = INT__g_s_d_s[abcd * 90 + 14 * 6 + iket] + ( AB_z[abcd] * INT__f_s_d_s[abcd * 60 + 9 * 6 + iket] );

            // (D_200 D_200|_{i} = (F_300 P_100| + x_ab * (D_200 P_100|
            INT__d_d_d_s[abcd * 216 + 0 * 6 + iket] = Q_f_s_s_p_s_s_d + ( AB_x[abcd] * Q_d_s_s_p_s_s_d );

            // (D_200 D_110|_{i} = (F_210 P_100| + y_ab * (D_200 P_100|
            INT__d_d_d_s[abcd * 216 + 1 * 6 + iket] = Q_d_p_s_p_s_s_d + ( AB_y[abcd] * Q_d_s_s_p_s_s_d );

            // (D_200 D_101|_{i} = (F_201 P_100| + z_ab * (D_200 P_100|
            INT__d_d_d_s[abcd * 216 + 2 * 6 + iket] = Q_d_s_p_p_s_s_d + ( AB_z[abcd] * Q_d_s_s_p_s_s_d );

            // (D_200 D_020|_{i} = (F_210 P_010| + y_ab * (D_200 P_010|
            INT__d_d_d_s[abcd * 216 + 3 * 6 + iket] = Q_d_p_s_s_p_s_d + ( AB_y[abcd] * Q_d_s_s_s_p_s_d );

            // (D_200 D_011|_{i} = (F_201 P_010| + z_ab * (D_200 P_010|
            INT__d_d_d_s[abcd * 216 + 4 * 6 + iket] = Q_d_s_p_s_p_s_d + ( AB_z[abcd] * Q_d_s_s_s_p_s_d );

            // (D_200 D_002|_{i} = (F_201 P_001| + z_ab * (D_200 P_001|
            INT__d_d_d_s[abcd * 216 + 5 * 6 + iket] = Q_d_s_p_s_s_p_d + ( AB_z[abcd] * Q_d_s_s_s_s_p_d );

            // (D_110 D_200|_{i} = (F_210 P_100| + x_ab * (D_110 P_100|
            INT__d_d_d_s[abcd * 216 + 6 * 6 + iket] = Q_d_p_s_p_s_s_d + ( AB_x[abcd] * Q_p_p_s_p_s_s_d );

            // (D_110 D_110|_{i} = (F_120 P_100| + y_ab * (D_110 P_100|
            INT__d_d_d_s[abcd * 216 + 7 * 6 + iket] = Q_p_d_s_p_s_s_d + ( AB_y[abcd] * Q_p_p_s_p_s_s_d );

            // (D_110 D_101|_{i} = (F_111 P_100| + z_ab * (D_110 P_100|
            INT__d_d_d_s[abcd * 216 + 8 * 6 + iket] = Q_p_p_p_p_s_s_d + ( AB_z[abcd] * Q_p_p_s_p_s_s_d );

            // (D_110 D_020|_{i} = (F_120 P_010| + y_ab * (D_110 P_010|
            INT__d_d_d_s[abcd * 216 + 9 * 6 + iket] = Q_p_d_s_s_p_s_d + ( AB_y[abcd] * Q_p_p_s_s_p_s_d );

            // (D_110 D_011|_{i} = (F_111 P_010| + z_ab * (D_110 P_010|
            INT__d_d_d_s[abcd * 216 + 10 * 6 + iket] = Q_p_p_p_s_p_s_d + ( AB_z[abcd] * Q_p_p_s_s_p_s_d );

            // (D_110 D_002|_{i} = (F_111 P_001| + z_ab * (D_110 P_001|
            INT__d_d_d_s[abcd * 216 + 11 * 6 + iket] = Q_p_p_p_s_s_p_d + ( AB_z[abcd] * Q_p_p_s_s_s_p_d );

            // (D_101 D_200|_{i} = (F_201 P_100| + x_ab * (D_101 P_100|
            INT__d_d_d_s[abcd * 216 + 12 * 6 + iket] = Q_d_s_p_p_s_s_d + ( AB_x[abcd] * Q_p_s_p_p_s_s_d );

            // (D_101 D_110|_{i} = (F_111 P_100| + y_ab * (D_101 P_100|
            INT__d_d_d_s[abcd * 216 + 13 * 6 + iket] = Q_p_p_p_p_s_s_d + ( AB_y[abcd] * Q_p_s_p_p_s_s_d );

            // (D_101 D_101|_{i} = (F_102 P_100| + z_ab * (D_101 P_100|
            INT__d_d_d_s[abcd * 216 + 14 * 6 + iket] = Q_p_s_d_p_s_s_d + ( AB_z[abcd] * Q_p_s_p_p_s_s_d );

            // (D_101 D_020|_{i} = (F_111 P_010| + y_ab * (D_101 P_010|
            INT__d_d_d_s[abcd * 216 + 15 * 6 + iket] = Q_p_p_p_s_p_s_d + ( AB_y[abcd] * Q_p_s_p_s_p_s_d );

            // (D_101 D_011|_{i} = (F_102 P_010| + z_ab * (D_101 P_010|
            INT__d_d_d_s[abcd * 216 + 16 * 6 + iket] = Q_p_s_d_s_p_s_d + ( AB_z[abcd] * Q_p_s_p_s_p_s_d );

            // (D_101 D_002|_{i} = (F_102 P_001| + z_ab * (D_101 P_001|
            INT__d_d_d_s[abcd * 216 + 17 * 6 + iket] = Q_p_s_d_s_s_p_d + ( AB_z[abcd] * Q_p_s_p_s_s_p_d );

            // (D_020 D_200|_{i} = (F_120 P_100| + x_ab * (D_020 P_100|
            INT__d_d_d_s[abcd * 216 + 18 * 6 + iket] = Q_p_d_s_p_s_s_d + ( AB_x[abcd] * Q_s_d_s_p_s_s_d );

            // (D_020 D_110|_{i} = (F_030 P_100| + y_ab * (D_020 P_100|
            INT__d_d_d_s[abcd * 216 + 19 * 6 + iket] = Q_s_f_s_p_s_s_d + ( AB_y[abcd] * Q_s_d_s_p_s_s_d );

            // (D_020 D_101|_{i} = (F_021 P_100| + z_ab * (D_020 P_100|
            INT__d_d_d_s[abcd * 216 + 20 * 6 + iket] = Q_s_d_p_p_s_s_d + ( AB_z[abcd] * Q_s_d_s_p_s_s_d );

            // (D_020 D_020|_{i} = (F_030 P_010| + y_ab * (D_020 P_010|
            INT__d_d_d_s[abcd * 216 + 21 * 6 + iket] = Q_s_f_s_s_p_s_d + ( AB_y[abcd] * Q_s_d_s_s_p_s_d );

            // (D_020 D_011|_{i} = (F_021 P_010| + z_ab * (D_020 P_010|
            INT__d_d_d_s[abcd * 216 + 22 * 6 + iket] = Q_s_d_p_s_p_s_d + ( AB_z[abcd] * Q_s_d_s_s_p_s_d );

            // (D_020 D_002|_{i} = (F_021 P_001| + z_ab * (D_020 P_001|
            INT__d_d_d_s[abcd * 216 + 23 * 6 + iket] = Q_s_d_p_s_s_p_d + ( AB_z[abcd] * Q_s_d_s_s_s_p_d );

            // (D_011 D_200|_{i} = (F_111 P_100| + x_ab * (D_011 P_100|
            INT__d_d_d_s[abcd * 216 + 24 * 6 + iket] = Q_p_p_p_p_s_s_d + ( AB_x[abcd] * Q_s_p_p_p_s_s_d );

            // (D_011 D_110|_{i} = (F_021 P_100| + y_ab * (D_011 P_100|
            INT__d_d_d_s[abcd * 216 + 25 * 6 + iket] = Q_s_d_p_p_s_s_d + ( AB_y[abcd] * Q_s_p_p_p_s_s_d );

            // (D_011 D_101|_{i} = (F_012 P_100| + z_ab * (D_011 P_100|
            INT__d_d_d_s[abcd * 216 + 26 * 6 + iket] = Q_s_p_d_p_s_s_d + ( AB_z[abcd] * Q_s_p_p_p_s_s_d );

            // (D_011 D_020|_{i} = (F_021 P_010| + y_ab * (D_011 P_010|
            INT__d_d_d_s[abcd * 216 + 27 * 6 + iket] = Q_s_d_p_s_p_s_d + ( AB_y[abcd] * Q_s_p_p_s_p_s_d );

            // (D_011 D_011|_{i} = (F_012 P_010| + z_ab * (D_011 P_010|
            INT__d_d_d_s[abcd * 216 + 28 * 6 + iket] = Q_s_p_d_s_p_s_d + ( AB_z[abcd] * Q_s_p_p_s_p_s_d );

            // (D_011 D_002|_{i} = (F_012 P_001| + z_ab * (D_011 P_001|
            INT__d_d_d_s[abcd * 216 + 29 * 6 + iket] = Q_s_p_d_s_s_p_d + ( AB_z[abcd] * Q_s_p_p_s_s_p_d );

            // (D_002 D_200|_{i} = (F_102 P_100| + x_ab * (D_002 P_100|
            INT__d_d_d_s[abcd * 216 + 30 * 6 + iket] = Q_p_s_d_p_s_s_d + ( AB_x[abcd] * Q_s_s_d_p_s_s_d );

            // (D_002 D_110|_{i} = (F_012 P_100| + y_ab * (D_002 P_100|
            INT__d_d_d_s[abcd * 216 + 31 * 6 + iket] = Q_s_p_d_p_s_s_d + ( AB_y[abcd] * Q_s_s_d_p_s_s_d );

            // (D_002 D_101|_{i} = (F_003 P_100| + z_ab * (D_002 P_100|
            INT__d_d_d_s[abcd * 216 + 32 * 6 + iket] = Q_s_s_f_p_s_s_d + ( AB_z[abcd] * Q_s_s_d_p_s_s_d );

            // (D_002 D_020|_{i} = (F_012 P_010| + y_ab * (D_002 P_010|
            INT__d_d_d_s[abcd * 216 + 33 * 6 + iket] = Q_s_p_d_s_p_s_d + ( AB_y[abcd] * Q_s_s_d_s_p_s_d );

            // (D_002 D_011|_{i} = (F_003 P_010| + z_ab * (D_002 P_010|
            INT__d_d_d_s[abcd * 216 + 34 * 6 + iket] = Q_s_s_f_s_p_s_d + ( AB_z[abcd] * Q_s_s_d_s_p_s_d );

            // (D_002 D_002|_{i} = (F_003 P_001| + z_ab * (D_002 P_001|
            INT__d_d_d_s[abcd * 216 + 35 * 6 + iket] = Q_s_s_f_s_s_p_d + ( AB_z[abcd] * Q_s_s_d_s_s_p_d );

        }


    }


    //////////////////////////////////////////////
    // Contracted integrals: Horizontal recurrance
    // Ket part
    // Steps: 9
    //////////////////////////////////////////////

    for(abcd = 0; abcd < nshell1234; ++abcd)
    {
        for(ibra = 0; ibra < 36; ++ibra)
        {
            // |P_100 P_100)_{i} = |D_200 S_000)_{t} + x_cd * |P_100 S_000)_{t}
            INT__d_d_p_p[abcd * 324 + ibra * 9 + 0] = INT__d_d_d_s[abcd * 216 + ibra * 6 + 0] + ( CD_x[abcd] * INT__d_d_p_s[abcd * 108 + ibra * 3 + 0] );

            // |P_100 P_010)_{i} = |D_110 S_000)_{t} + y_cd * |P_100 S_000)_{t}
            INT__d_d_p_p[abcd * 324 + ibra * 9 + 1] = INT__d_d_d_s[abcd * 216 + ibra * 6 + 1] + ( CD_y[abcd] * INT__d_d_p_s[abcd * 108 + ibra * 3 + 0] );

            // |P_100 P_001)_{i} = |D_101 S_000)_{t} + z_cd * |P_100 S_000)_{t}
            INT__d_d_p_p[abcd * 324 + ibra * 9 + 2] = INT__d_d_d_s[abcd * 216 + ibra * 6 + 2] + ( CD_z[abcd] * INT__d_d_p_s[abcd * 108 + ibra * 3 + 0] );

            // |P_010 P_100)_{i} = |D_110 S_000)_{t} + x_cd * |P_010 S_000)_{t}
            INT__d_d_p_p[abcd * 324 + ibra * 9 + 3] = INT__d_d_d_s[abcd * 216 + ibra * 6 + 1] + ( CD_x[abcd] * INT__d_d_p_s[abcd * 108 + ibra * 3 + 1] );

            // |P_010 P_010)_{i} = |D_020 S_000)_{t} + y_cd * |P_010 S_000)_{t}
            INT__d_d_p_p[abcd * 324 + ibra * 9 + 4] = INT__d_d_d_s[abcd * 216 + ibra * 6 + 3] + ( CD_y[abcd] * INT__d_d_p_s[abcd * 108 + ibra * 3 + 1] );

            // |P_010 P_001)_{i} = |D_011 S_000)_{t} + z_cd * |P_010 S_000)_{t}
            INT__d_d_p_p[abcd * 324 + ibra * 9 + 5] = INT__d_d_d_s[abcd * 216 + ibra * 6 + 4] + ( CD_z[abcd] * INT__d_d_p_s[abcd * 108 + ibra * 3 + 1] );

            // |P_001 P_100)_{i} = |D_101 S_000)_{t} + x_cd * |P_001 S_000)_{t}
            INT__d_d_p_p[abcd * 324 + ibra * 9 + 6] = INT__d_d_d_s[abcd * 216 + ibra * 6 + 2] + ( CD_x[abcd] * INT__d_d_p_s[abcd * 108 + ibra * 3 + 2] );

            // |P_001 P_010)_{i} = |D_011 S_000)_{t} + y_cd * |P_001 S_000)_{t}
            INT__d_d_p_p[abcd * 324 + ibra * 9 + 7] = INT__d_d_d_s[abcd * 216 + ibra * 6 + 4] + ( CD_y[abcd] * INT__d_d_p_s[abcd * 108 + ibra * 3 + 2] );

            // |P_001 P_001)_{i} = |D_002 S_000)_{t} + z_cd * |P_001 S_000)_{t}
            INT__d_d_p_p[abcd * 324 + ibra * 9 + 8] = INT__d_d_d_s[abcd * 216 + ibra * 6 + 5] + ( CD_z[abcd] * INT__d_d_p_s[abcd * 108 + ibra * 3 + 2] );

        }
    }


    // Free contracted work space
    free(contwork);

    return nshell1234;
}

