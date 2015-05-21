#include <string.h>
#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"
#include "boys/boys_split.h"


int eri_split_d_p_d_s(struct multishell_pair const P,
                      struct multishell_pair const Q,
                      double * const restrict S_2_1_2_0)
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

    ASSUME_ALIGN(integrals)

    const int nshell1234 = P.nshell12 * Q.nshell12;


    // Holds AB_{xyz} and CD_{xyz} in a flattened fashion for later
    double AB_x[nshell1234];  double CD_x[nshell1234];
    double AB_y[nshell1234];  double CD_y[nshell1234];
    double AB_z[nshell1234];  double CD_z[nshell1234];

    int ab, cd, abcd;
    int i, j;

    // Workspace for contracted integrals
    double * const contwork = malloc(nshell1234 * 768);
    memset(contwork, 0, nshell1234 * 768);

    // partition workspace into shells
    double * const S_2_0_2_0 = contwork + (nshell1234 * 0);
    double * const S_3_0_2_0 = contwork + (nshell1234 * 36);


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
            double * const restrict PRIM_S_2_0_2_0 = S_2_0_2_0 + (abcd * 36);
            double * const restrict PRIM_S_3_0_2_0 = S_3_0_2_0 + (abcd * 60);

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
                for(j = cdstart; j < cdend; ++j)
                {

                    // Holds the auxiliary integrals ( i 0 | 0 0 )^m in the primitive basis
                    // with m as the slowest index
                    // AM = 0: Needed from this AM: 1
                    double AUX_S_0_0_0_0[6 * 1];

                    // AM = 1: Needed from this AM: 3
                    double AUX_S_1_0_0_0[5 * 3];

                    // AM = 2: Needed from this AM: 6
                    double AUX_S_2_0_0_0[4 * 6];

                    // AM = 3: Needed from this AM: 10
                    double AUX_S_3_0_0_0[3 * 10];

                    // AM = 4: Needed from this AM: 15
                    double AUX_S_4_0_0_0[2 * 15];

                    // AM = 5: Needed from this AM: 21
                    double AUX_S_5_0_0_0[1 * 21];



                    // Holds temporary integrals for electron transfer
                    double AUX_S_1_0_1_0[9];
                    double AUX_S_2_0_1_0[18];
                    double AUX_S_2_0_2_0[36];
                    double AUX_S_3_0_1_0[30];
                    double AUX_S_3_0_2_0[60];
                    double AUX_S_4_0_1_0[45];


                    const double PQalpha_mul = P.alpha[i] * Q.alpha[j];
                    const double PQalpha_sum = P.alpha[i] + Q.alpha[j];

                    const double pfac = TWO_PI_52 / (PQalpha_mul * sqrt(PQalpha_sum));

                    /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
                    const double PQ_x = P.x[i] - Q.x[j];
                    const double PQ_y = P.y[i] - Q.y[j];
                    const double PQ_z = P.z[i] - Q.z[j];
                    const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;

                    // collected prefactors
                    const double allprefac =  pfac * P.prefac[i] * Q.prefac[j];

                    // various factors
                    const double alpha = PQalpha_mul/PQalpha_sum;   // alpha from MEST
                    // for VRR
                    const double one_over_p = 1.0 / P.alpha[i];
                    const double a_over_p =  alpha * one_over_p;     // a/p from MEST
                    const double one_over_2p = 0.5 * one_over_p;  // gets multiplied by i in VRR
                    // for electron transfer
                    const double one_over_q = 1.0 / Q.alpha[j];
                    const double one_over_2q = 0.5 * one_over_q;
                    const double p_over_q = P.alpha[i] * one_over_q;

                    const double etfac[3] = {
                                             -(P.bAB_x[i] + Q.bAB_x[j]) * one_over_q,
                                             -(P.bAB_y[i] + Q.bAB_y[j]) * one_over_q,
                                             -(P.bAB_z[i] + Q.bAB_z[j]) * one_over_q,
                                            };


                    //////////////////////////////////////////////
                    // Boys function section
                    // Maximum v value: 5
                    //////////////////////////////////////////////
                    // The paremeter to the boys function
                    const double F_x = R2 * alpha;


                    Boys_F_split(AUX_S_0_0_0_0, 5, F_x);
                    AUX_S_0_0_0_0[0] *= allprefac;
                    AUX_S_0_0_0_0[1] *= allprefac;
                    AUX_S_0_0_0_0[2] *= allprefac;
                    AUX_S_0_0_0_0[3] *= allprefac;
                    AUX_S_0_0_0_0[4] *= allprefac;
                    AUX_S_0_0_0_0[5] *= allprefac;

                    //////////////////////////////////////////////
                    // Primitive integrals: Vertical recurrance
                    //////////////////////////////////////////////

                    // Forming AUX_S_1_0_0_0[5 * 3];
                    // Needed from this AM:
                    //    P_100
                    //    P_010
                    //    P_001
                    for(int m = 0; m < 5; m++)  // loop over orders of boys function
                    {
                        //P_100 : STEP: x
                        AUX_S_1_0_0_0[m * 3 + 0] = P.PA_x[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_x * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                        //P_010 : STEP: y
                        AUX_S_1_0_0_0[m * 3 + 1] = P.PA_y[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_y * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                        //P_001 : STEP: z
                        AUX_S_1_0_0_0[m * 3 + 2] = P.PA_z[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_z * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                    }


                    // Forming AUX_S_2_0_0_0[4 * 6];
                    // Needed from this AM:
                    //    D_200
                    //    D_110
                    //    D_101
                    //    D_020
                    //    D_011
                    //    D_002
                    for(int m = 0; m < 4; m++)  // loop over orders of boys function
                    {
                        //D_200 : STEP: x
                        AUX_S_2_0_0_0[m * 6 + 0] = P.PA_x[i] * AUX_S_1_0_0_0[m * 3 + 0] - a_over_p * PQ_x * AUX_S_1_0_0_0[(m+1) * 3 + 0]
                                      + 1 * one_over_2p * ( AUX_S_0_0_0_0[m * 1 +  0] - a_over_p * AUX_S_0_0_0_0[(m+1) * 1 + 0] );

                        //D_110 : STEP: y
                        AUX_S_2_0_0_0[m * 6 + 1] = P.PA_y[i] * AUX_S_1_0_0_0[m * 3 + 0] - a_over_p * PQ_y * AUX_S_1_0_0_0[(m+1) * 3 + 0];

                        //D_101 : STEP: z
                        AUX_S_2_0_0_0[m * 6 + 2] = P.PA_z[i] * AUX_S_1_0_0_0[m * 3 + 0] - a_over_p * PQ_z * AUX_S_1_0_0_0[(m+1) * 3 + 0];

                        //D_020 : STEP: y
                        AUX_S_2_0_0_0[m * 6 + 3] = P.PA_y[i] * AUX_S_1_0_0_0[m * 3 + 1] - a_over_p * PQ_y * AUX_S_1_0_0_0[(m+1) * 3 + 1]
                                      + 1 * one_over_2p * ( AUX_S_0_0_0_0[m * 1 +  0] - a_over_p * AUX_S_0_0_0_0[(m+1) * 1 + 0] );

                        //D_011 : STEP: z
                        AUX_S_2_0_0_0[m * 6 + 4] = P.PA_z[i] * AUX_S_1_0_0_0[m * 3 + 1] - a_over_p * PQ_z * AUX_S_1_0_0_0[(m+1) * 3 + 1];

                        //D_002 : STEP: z
                        AUX_S_2_0_0_0[m * 6 + 5] = P.PA_z[i] * AUX_S_1_0_0_0[m * 3 + 2] - a_over_p * PQ_z * AUX_S_1_0_0_0[(m+1) * 3 + 2]
                                      + 1 * one_over_2p * ( AUX_S_0_0_0_0[m * 1 +  0] - a_over_p * AUX_S_0_0_0_0[(m+1) * 1 + 0] );

                    }


                    // Forming AUX_S_3_0_0_0[3 * 10];
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
                    for(int m = 0; m < 3; m++)  // loop over orders of boys function
                    {
                        //F_300 : STEP: x
                        AUX_S_3_0_0_0[m * 10 + 0] = P.PA_x[i] * AUX_S_2_0_0_0[m * 6 + 0] - a_over_p * PQ_x * AUX_S_2_0_0_0[(m+1) * 6 + 0]
                                      + 2 * one_over_2p * ( AUX_S_1_0_0_0[m * 3 +  0] - a_over_p * AUX_S_1_0_0_0[(m+1) * 3 + 0] );

                        //F_210 : STEP: y
                        AUX_S_3_0_0_0[m * 10 + 1] = P.PA_y[i] * AUX_S_2_0_0_0[m * 6 + 0] - a_over_p * PQ_y * AUX_S_2_0_0_0[(m+1) * 6 + 0];

                        //F_201 : STEP: z
                        AUX_S_3_0_0_0[m * 10 + 2] = P.PA_z[i] * AUX_S_2_0_0_0[m * 6 + 0] - a_over_p * PQ_z * AUX_S_2_0_0_0[(m+1) * 6 + 0];

                        //F_120 : STEP: x
                        AUX_S_3_0_0_0[m * 10 + 3] = P.PA_x[i] * AUX_S_2_0_0_0[m * 6 + 3] - a_over_p * PQ_x * AUX_S_2_0_0_0[(m+1) * 6 + 3];

                        //F_111 : STEP: z
                        AUX_S_3_0_0_0[m * 10 + 4] = P.PA_z[i] * AUX_S_2_0_0_0[m * 6 + 1] - a_over_p * PQ_z * AUX_S_2_0_0_0[(m+1) * 6 + 1];

                        //F_102 : STEP: x
                        AUX_S_3_0_0_0[m * 10 + 5] = P.PA_x[i] * AUX_S_2_0_0_0[m * 6 + 5] - a_over_p * PQ_x * AUX_S_2_0_0_0[(m+1) * 6 + 5];

                        //F_030 : STEP: y
                        AUX_S_3_0_0_0[m * 10 + 6] = P.PA_y[i] * AUX_S_2_0_0_0[m * 6 + 3] - a_over_p * PQ_y * AUX_S_2_0_0_0[(m+1) * 6 + 3]
                                      + 2 * one_over_2p * ( AUX_S_1_0_0_0[m * 3 +  1] - a_over_p * AUX_S_1_0_0_0[(m+1) * 3 + 1] );

                        //F_021 : STEP: z
                        AUX_S_3_0_0_0[m * 10 + 7] = P.PA_z[i] * AUX_S_2_0_0_0[m * 6 + 3] - a_over_p * PQ_z * AUX_S_2_0_0_0[(m+1) * 6 + 3];

                        //F_012 : STEP: y
                        AUX_S_3_0_0_0[m * 10 + 8] = P.PA_y[i] * AUX_S_2_0_0_0[m * 6 + 5] - a_over_p * PQ_y * AUX_S_2_0_0_0[(m+1) * 6 + 5];

                        //F_003 : STEP: z
                        AUX_S_3_0_0_0[m * 10 + 9] = P.PA_z[i] * AUX_S_2_0_0_0[m * 6 + 5] - a_over_p * PQ_z * AUX_S_2_0_0_0[(m+1) * 6 + 5]
                                      + 2 * one_over_2p * ( AUX_S_1_0_0_0[m * 3 +  2] - a_over_p * AUX_S_1_0_0_0[(m+1) * 3 + 2] );

                    }


                    // Forming AUX_S_4_0_0_0[2 * 15];
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
                    for(int m = 0; m < 2; m++)  // loop over orders of boys function
                    {
                        //G_400 : STEP: x
                        AUX_S_4_0_0_0[m * 15 + 0] = P.PA_x[i] * AUX_S_3_0_0_0[m * 10 + 0] - a_over_p * PQ_x * AUX_S_3_0_0_0[(m+1) * 10 + 0]
                                      + 3 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  0] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 0] );

                        //G_310 : STEP: y
                        AUX_S_4_0_0_0[m * 15 + 1] = P.PA_y[i] * AUX_S_3_0_0_0[m * 10 + 0] - a_over_p * PQ_y * AUX_S_3_0_0_0[(m+1) * 10 + 0];

                        //G_301 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 2] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 0] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 0];

                        //G_220 : STEP: y
                        AUX_S_4_0_0_0[m * 15 + 3] = P.PA_y[i] * AUX_S_3_0_0_0[m * 10 + 1] - a_over_p * PQ_y * AUX_S_3_0_0_0[(m+1) * 10 + 1]
                                      + 1 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  0] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 0] );

                        //G_211 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 4] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 1] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 1];

                        //G_202 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 5] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 2] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 2]
                                      + 1 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  0] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 0] );

                        //G_130 : STEP: x
                        AUX_S_4_0_0_0[m * 15 + 6] = P.PA_x[i] * AUX_S_3_0_0_0[m * 10 + 6] - a_over_p * PQ_x * AUX_S_3_0_0_0[(m+1) * 10 + 6];

                        //G_121 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 7] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 3] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 3];

                        //G_112 : STEP: y
                        AUX_S_4_0_0_0[m * 15 + 8] = P.PA_y[i] * AUX_S_3_0_0_0[m * 10 + 5] - a_over_p * PQ_y * AUX_S_3_0_0_0[(m+1) * 10 + 5];

                        //G_103 : STEP: x
                        AUX_S_4_0_0_0[m * 15 + 9] = P.PA_x[i] * AUX_S_3_0_0_0[m * 10 + 9] - a_over_p * PQ_x * AUX_S_3_0_0_0[(m+1) * 10 + 9];

                        //G_040 : STEP: y
                        AUX_S_4_0_0_0[m * 15 + 10] = P.PA_y[i] * AUX_S_3_0_0_0[m * 10 + 6] - a_over_p * PQ_y * AUX_S_3_0_0_0[(m+1) * 10 + 6]
                                      + 3 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  3] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 3] );

                        //G_031 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 11] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 6] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 6];

                        //G_022 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 12] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 7] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 7]
                                      + 1 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  3] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 3] );

                        //G_013 : STEP: y
                        AUX_S_4_0_0_0[m * 15 + 13] = P.PA_y[i] * AUX_S_3_0_0_0[m * 10 + 9] - a_over_p * PQ_y * AUX_S_3_0_0_0[(m+1) * 10 + 9];

                        //G_004 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 14] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 9] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 9]
                                      + 3 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  5] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 5] );

                    }


                    // Forming AUX_S_5_0_0_0[1 * 21];
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
                    for(int m = 0; m < 1; m++)  // loop over orders of boys function
                    {
                        //H_500 : STEP: x
                        AUX_S_5_0_0_0[m * 21 + 0] = P.PA_x[i] * AUX_S_4_0_0_0[m * 15 + 0] - a_over_p * PQ_x * AUX_S_4_0_0_0[(m+1) * 15 + 0]
                                      + 4 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  0] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 0] );

                        //H_410 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 1] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 0] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 0];

                        //H_401 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 2] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 0] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 0];

                        //H_320 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 3] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 1] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 1]
                                      + 1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  0] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 0] );

                        //H_311 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 4] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 1] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 1];

                        //H_302 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 5] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 2] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 2]
                                      + 1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  0] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 0] );

                        //H_230 : STEP: x
                        AUX_S_5_0_0_0[m * 21 + 6] = P.PA_x[i] * AUX_S_4_0_0_0[m * 15 + 6] - a_over_p * PQ_x * AUX_S_4_0_0_0[(m+1) * 15 + 6]
                                      + 1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  6] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 6] );

                        //H_221 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 7] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 3] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 3];

                        //H_212 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 8] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 5] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 5];

                        //H_203 : STEP: x
                        AUX_S_5_0_0_0[m * 21 + 9] = P.PA_x[i] * AUX_S_4_0_0_0[m * 15 + 9] - a_over_p * PQ_x * AUX_S_4_0_0_0[(m+1) * 15 + 9]
                                      + 1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  9] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 9] );

                        //H_140 : STEP: x
                        AUX_S_5_0_0_0[m * 21 + 10] = P.PA_x[i] * AUX_S_4_0_0_0[m * 15 + 10] - a_over_p * PQ_x * AUX_S_4_0_0_0[(m+1) * 15 + 10];

                        //H_131 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 11] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 6] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 6];

                        //H_122 : STEP: x
                        AUX_S_5_0_0_0[m * 21 + 12] = P.PA_x[i] * AUX_S_4_0_0_0[m * 15 + 12] - a_over_p * PQ_x * AUX_S_4_0_0_0[(m+1) * 15 + 12];

                        //H_113 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 13] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 9] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 9];

                        //H_104 : STEP: x
                        AUX_S_5_0_0_0[m * 21 + 14] = P.PA_x[i] * AUX_S_4_0_0_0[m * 15 + 14] - a_over_p * PQ_x * AUX_S_4_0_0_0[(m+1) * 15 + 14];

                        //H_050 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 15] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 10] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 10]
                                      + 4 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  6] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 6] );

                        //H_041 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 16] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 10] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 10];

                        //H_032 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 17] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 11] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 11]
                                      + 1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  6] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 6] );

                        //H_023 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 18] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 13] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 13]
                                      + 1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  9] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 9] );

                        //H_014 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 19] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 14] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 14];

                        //H_005 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 20] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 14] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 14]
                                      + 4 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  9] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 9] );

                    }




                    //////////////////////////////////////////////
                    // Primitive integrals: Electron transfer
                    //////////////////////////////////////////////

                    // ( D_200 S_000 | P_100 S_000 )^0 = x * ( D_200 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_300 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[0] = etfac[0] * AUX_S_2_0_0_0[0] + 2 * one_over_2q * AUX_S_1_0_0_0[0] - p_over_q * AUX_S_3_0_0_0[0];

                    // ( D_200 S_000 | P_010 S_000 )^0 = y * ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[1] = etfac[1] * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_3_0_0_0[1];

                    // ( D_200 S_000 | P_001 S_000 )^0 = z * ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[2] = etfac[2] * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_3_0_0_0[2];

                    // ( D_110 S_000 | P_100 S_000 )^0 = x * ( D_110 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[3] = etfac[0] * AUX_S_2_0_0_0[1] + 1 * one_over_2q * AUX_S_1_0_0_0[1] - p_over_q * AUX_S_3_0_0_0[1];

                    // ( D_110 S_000 | P_010 S_000 )^0 = y * ( D_110 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[4] = etfac[1] * AUX_S_2_0_0_0[1] + 1 * one_over_2q * AUX_S_1_0_0_0[0] - p_over_q * AUX_S_3_0_0_0[3];

                    // ( D_110 S_000 | P_001 S_000 )^0 = z * ( D_110 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[5] = etfac[2] * AUX_S_2_0_0_0[1] - p_over_q * AUX_S_3_0_0_0[4];

                    // ( D_101 S_000 | P_100 S_000 )^0 = x * ( D_101 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[6] = etfac[0] * AUX_S_2_0_0_0[2] + 1 * one_over_2q * AUX_S_1_0_0_0[2] - p_over_q * AUX_S_3_0_0_0[2];

                    // ( D_101 S_000 | P_010 S_000 )^0 = y * ( D_101 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[7] = etfac[1] * AUX_S_2_0_0_0[2] - p_over_q * AUX_S_3_0_0_0[4];

                    // ( D_101 S_000 | P_001 S_000 )^0 = z * ( D_101 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[8] = etfac[2] * AUX_S_2_0_0_0[2] + 1 * one_over_2q * AUX_S_1_0_0_0[0] - p_over_q * AUX_S_3_0_0_0[5];

                    // ( D_020 S_000 | P_100 S_000 )^0 = x * ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[9] = etfac[0] * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_3_0_0_0[3];

                    // ( D_020 S_000 | P_010 S_000 )^0 = y * ( D_020 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_030 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[10] = etfac[1] * AUX_S_2_0_0_0[3] + 2 * one_over_2q * AUX_S_1_0_0_0[1] - p_over_q * AUX_S_3_0_0_0[6];

                    // ( D_020 S_000 | P_001 S_000 )^0 = z * ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[11] = etfac[2] * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_3_0_0_0[7];

                    // ( D_011 S_000 | P_100 S_000 )^0 = x * ( D_011 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[12] = etfac[0] * AUX_S_2_0_0_0[4] - p_over_q * AUX_S_3_0_0_0[4];

                    // ( D_011 S_000 | P_010 S_000 )^0 = y * ( D_011 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[13] = etfac[1] * AUX_S_2_0_0_0[4] + 1 * one_over_2q * AUX_S_1_0_0_0[2] - p_over_q * AUX_S_3_0_0_0[7];

                    // ( D_011 S_000 | P_001 S_000 )^0 = z * ( D_011 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[14] = etfac[2] * AUX_S_2_0_0_0[4] + 1 * one_over_2q * AUX_S_1_0_0_0[1] - p_over_q * AUX_S_3_0_0_0[8];

                    // ( D_002 S_000 | P_100 S_000 )^0 = x * ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[15] = etfac[0] * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_3_0_0_0[5];

                    // ( D_002 S_000 | P_010 S_000 )^0 = y * ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[16] = etfac[1] * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_3_0_0_0[8];

                    // ( D_002 S_000 | P_001 S_000 )^0 = z * ( D_002 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_003 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[17] = etfac[2] * AUX_S_2_0_0_0[5] + 2 * one_over_2q * AUX_S_1_0_0_0[2] - p_over_q * AUX_S_3_0_0_0[9];

                    // ( F_300 S_000 | P_100 S_000 )^0 = x * ( F_300 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_400 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[0] = etfac[0] * AUX_S_3_0_0_0[0] + 3 * one_over_2q * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_4_0_0_0[0];

                    // ( F_300 S_000 | P_010 S_000 )^0 = y * ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[1] = etfac[1] * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_4_0_0_0[1];

                    // ( F_300 S_000 | P_001 S_000 )^0 = z * ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[2] = etfac[2] * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_4_0_0_0[2];

                    // ( F_210 S_000 | P_100 S_000 )^0 = x * ( F_210 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[3] = etfac[0] * AUX_S_3_0_0_0[1] + 2 * one_over_2q * AUX_S_2_0_0_0[1] - p_over_q * AUX_S_4_0_0_0[1];

                    // ( F_210 S_000 | P_010 S_000 )^0 = y * ( F_210 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[4] = etfac[1] * AUX_S_3_0_0_0[1] + 1 * one_over_2q * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_4_0_0_0[3];

                    // ( F_210 S_000 | P_001 S_000 )^0 = z * ( F_210 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[5] = etfac[2] * AUX_S_3_0_0_0[1] - p_over_q * AUX_S_4_0_0_0[4];

                    // ( F_201 S_000 | P_100 S_000 )^0 = x * ( F_201 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[6] = etfac[0] * AUX_S_3_0_0_0[2] + 2 * one_over_2q * AUX_S_2_0_0_0[2] - p_over_q * AUX_S_4_0_0_0[2];

                    // ( F_201 S_000 | P_010 S_000 )^0 = y * ( F_201 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[7] = etfac[1] * AUX_S_3_0_0_0[2] - p_over_q * AUX_S_4_0_0_0[4];

                    // ( F_201 S_000 | P_001 S_000 )^0 = z * ( F_201 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[8] = etfac[2] * AUX_S_3_0_0_0[2] + 1 * one_over_2q * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_4_0_0_0[5];

                    // ( F_120 S_000 | P_100 S_000 )^0 = x * ( F_120 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[9] = etfac[0] * AUX_S_3_0_0_0[3] + 1 * one_over_2q * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_4_0_0_0[3];

                    // ( F_120 S_000 | P_010 S_000 )^0 = y * ( F_120 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[10] = etfac[1] * AUX_S_3_0_0_0[3] + 2 * one_over_2q * AUX_S_2_0_0_0[1] - p_over_q * AUX_S_4_0_0_0[6];

                    // ( F_120 S_000 | P_001 S_000 )^0 = z * ( F_120 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[11] = etfac[2] * AUX_S_3_0_0_0[3] - p_over_q * AUX_S_4_0_0_0[7];

                    // ( F_111 S_000 | P_100 S_000 )^0 = x * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[12] = etfac[0] * AUX_S_3_0_0_0[4] + 1 * one_over_2q * AUX_S_2_0_0_0[4] - p_over_q * AUX_S_4_0_0_0[4];

                    // ( F_111 S_000 | P_010 S_000 )^0 = y * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[13] = etfac[1] * AUX_S_3_0_0_0[4] + 1 * one_over_2q * AUX_S_2_0_0_0[2] - p_over_q * AUX_S_4_0_0_0[7];

                    // ( F_111 S_000 | P_001 S_000 )^0 = z * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[14] = etfac[2] * AUX_S_3_0_0_0[4] + 1 * one_over_2q * AUX_S_2_0_0_0[1] - p_over_q * AUX_S_4_0_0_0[8];

                    // ( F_102 S_000 | P_100 S_000 )^0 = x * ( F_102 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[15] = etfac[0] * AUX_S_3_0_0_0[5] + 1 * one_over_2q * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_4_0_0_0[5];

                    // ( F_102 S_000 | P_010 S_000 )^0 = y * ( F_102 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[16] = etfac[1] * AUX_S_3_0_0_0[5] - p_over_q * AUX_S_4_0_0_0[8];

                    // ( F_102 S_000 | P_001 S_000 )^0 = z * ( F_102 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[17] = etfac[2] * AUX_S_3_0_0_0[5] + 2 * one_over_2q * AUX_S_2_0_0_0[2] - p_over_q * AUX_S_4_0_0_0[9];

                    // ( F_030 S_000 | P_100 S_000 )^0 = x * ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[18] = etfac[0] * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_4_0_0_0[6];

                    // ( F_030 S_000 | P_010 S_000 )^0 = y * ( F_030 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_040 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[19] = etfac[1] * AUX_S_3_0_0_0[6] + 3 * one_over_2q * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_4_0_0_0[10];

                    // ( F_030 S_000 | P_001 S_000 )^0 = z * ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[20] = etfac[2] * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_4_0_0_0[11];

                    // ( F_021 S_000 | P_100 S_000 )^0 = x * ( F_021 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[21] = etfac[0] * AUX_S_3_0_0_0[7] - p_over_q * AUX_S_4_0_0_0[7];

                    // ( F_021 S_000 | P_010 S_000 )^0 = y * ( F_021 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[22] = etfac[1] * AUX_S_3_0_0_0[7] + 2 * one_over_2q * AUX_S_2_0_0_0[4] - p_over_q * AUX_S_4_0_0_0[11];

                    // ( F_021 S_000 | P_001 S_000 )^0 = z * ( F_021 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[23] = etfac[2] * AUX_S_3_0_0_0[7] + 1 * one_over_2q * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_4_0_0_0[12];

                    // ( F_012 S_000 | P_100 S_000 )^0 = x * ( F_012 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[24] = etfac[0] * AUX_S_3_0_0_0[8] - p_over_q * AUX_S_4_0_0_0[8];

                    // ( F_012 S_000 | P_010 S_000 )^0 = y * ( F_012 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[25] = etfac[1] * AUX_S_3_0_0_0[8] + 1 * one_over_2q * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_4_0_0_0[12];

                    // ( F_012 S_000 | P_001 S_000 )^0 = z * ( F_012 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[26] = etfac[2] * AUX_S_3_0_0_0[8] + 2 * one_over_2q * AUX_S_2_0_0_0[4] - p_over_q * AUX_S_4_0_0_0[13];

                    // ( F_003 S_000 | P_100 S_000 )^0 = x * ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[27] = etfac[0] * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_4_0_0_0[9];

                    // ( F_003 S_000 | P_010 S_000 )^0 = y * ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[28] = etfac[1] * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_4_0_0_0[13];

                    // ( F_003 S_000 | P_001 S_000 )^0 = z * ( F_003 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_004 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[29] = etfac[2] * AUX_S_3_0_0_0[9] + 3 * one_over_2q * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_4_0_0_0[14];

                    // ( P_100 S_000 | P_100 S_000 )^0 = x * ( P_100 S_000 | S_000 S_000 )^0_{e} + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( D_200 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_1_0_1_0[0] = etfac[0] * AUX_S_1_0_0_0[0] + 1 * one_over_2q * AUX_S_0_0_0_0[0] - p_over_q * AUX_S_2_0_0_0[0];

                    // ( P_100 S_000 | P_010 S_000 )^0 = y * ( P_100 S_000 | S_000 S_000 )^0_{e} - ( D_110 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_1_0_1_0[1] = etfac[1] * AUX_S_1_0_0_0[0] - p_over_q * AUX_S_2_0_0_0[1];

                    // ( P_100 S_000 | P_001 S_000 )^0 = z * ( P_100 S_000 | S_000 S_000 )^0_{e} - ( D_101 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_1_0_1_0[2] = etfac[2] * AUX_S_1_0_0_0[0] - p_over_q * AUX_S_2_0_0_0[2];

                    // ( P_010 S_000 | P_100 S_000 )^0 = x * ( P_010 S_000 | S_000 S_000 )^0_{e} - ( D_110 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_1_0_1_0[3] = etfac[0] * AUX_S_1_0_0_0[1] - p_over_q * AUX_S_2_0_0_0[1];

                    // ( P_010 S_000 | P_010 S_000 )^0 = y * ( P_010 S_000 | S_000 S_000 )^0_{e} + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( D_020 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_1_0_1_0[4] = etfac[1] * AUX_S_1_0_0_0[1] + 1 * one_over_2q * AUX_S_0_0_0_0[0] - p_over_q * AUX_S_2_0_0_0[3];

                    // ( P_010 S_000 | P_001 S_000 )^0 = z * ( P_010 S_000 | S_000 S_000 )^0_{e} - ( D_011 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_1_0_1_0[5] = etfac[2] * AUX_S_1_0_0_0[1] - p_over_q * AUX_S_2_0_0_0[4];

                    // ( P_001 S_000 | P_100 S_000 )^0 = x * ( P_001 S_000 | S_000 S_000 )^0_{e} - ( D_101 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_1_0_1_0[6] = etfac[0] * AUX_S_1_0_0_0[2] - p_over_q * AUX_S_2_0_0_0[2];

                    // ( P_001 S_000 | P_010 S_000 )^0 = y * ( P_001 S_000 | S_000 S_000 )^0_{e} - ( D_011 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_1_0_1_0[7] = etfac[1] * AUX_S_1_0_0_0[2] - p_over_q * AUX_S_2_0_0_0[4];

                    // ( P_001 S_000 | P_001 S_000 )^0 = z * ( P_001 S_000 | S_000 S_000 )^0_{e} + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( D_002 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_1_0_1_0[8] = etfac[2] * AUX_S_1_0_0_0[2] + 1 * one_over_2q * AUX_S_0_0_0_0[0] - p_over_q * AUX_S_2_0_0_0[5];

                    // ( D_200 S_000 | D_200 S_000 )^0_{t} = x * ( D_200 S_000 | P_100 S_000 )^0 + ( P_100 S_000 | P_100 S_000 )^0 + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_300 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[0] = etfac[0] * AUX_S_2_0_1_0[0] + 2 * one_over_2q * AUX_S_1_0_1_0[0] + 1 * one_over_2q * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_3_0_1_0[0];

                    // ( D_200 S_000 | D_110 S_000 )^0_{t} = y * ( D_200 S_000 | P_100 S_000 )^0 - ( F_210 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[1] = etfac[1] * AUX_S_2_0_1_0[0] - p_over_q * AUX_S_3_0_1_0[3];

                    // ( D_200 S_000 | D_101 S_000 )^0_{t} = z * ( D_200 S_000 | P_100 S_000 )^0 - ( F_201 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[2] = etfac[2] * AUX_S_2_0_1_0[0] - p_over_q * AUX_S_3_0_1_0[6];

                    // ( D_200 S_000 | D_020 S_000 )^0_{t} = y * ( D_200 S_000 | P_010 S_000 )^0 + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[3] = etfac[1] * AUX_S_2_0_1_0[1] + 1 * one_over_2q * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_3_0_1_0[4];

                    // ( D_200 S_000 | D_011 S_000 )^0_{t} = z * ( D_200 S_000 | P_010 S_000 )^0 - ( F_201 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[4] = etfac[2] * AUX_S_2_0_1_0[1] - p_over_q * AUX_S_3_0_1_0[7];

                    // ( D_200 S_000 | D_002 S_000 )^0_{t} = z * ( D_200 S_000 | P_001 S_000 )^0 + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | P_001 S_000 )^0
                    AUX_S_2_0_2_0[5] = etfac[2] * AUX_S_2_0_1_0[2] + 1 * one_over_2q * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_3_0_1_0[8];

                    // ( D_110 S_000 | D_200 S_000 )^0_{t} = x * ( D_110 S_000 | P_100 S_000 )^0 + ( P_010 S_000 | P_100 S_000 )^0 + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[6] = etfac[0] * AUX_S_2_0_1_0[3] + 1 * one_over_2q * AUX_S_1_0_1_0[3] + 1 * one_over_2q * AUX_S_2_0_0_0[1] - p_over_q * AUX_S_3_0_1_0[3];

                    // ( D_110 S_000 | D_110 S_000 )^0_{t} = y * ( D_110 S_000 | P_100 S_000 )^0 + ( P_100 S_000 | P_100 S_000 )^0 - ( F_120 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[7] = etfac[1] * AUX_S_2_0_1_0[3] + 1 * one_over_2q * AUX_S_1_0_1_0[0] - p_over_q * AUX_S_3_0_1_0[9];

                    // ( D_110 S_000 | D_101 S_000 )^0_{t} = z * ( D_110 S_000 | P_100 S_000 )^0 - ( F_111 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[8] = etfac[2] * AUX_S_2_0_1_0[3] - p_over_q * AUX_S_3_0_1_0[12];

                    // ( D_110 S_000 | D_020 S_000 )^0_{t} = y * ( D_110 S_000 | P_010 S_000 )^0 + ( P_100 S_000 | P_010 S_000 )^0 + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[9] = etfac[1] * AUX_S_2_0_1_0[4] + 1 * one_over_2q * AUX_S_1_0_1_0[1] + 1 * one_over_2q * AUX_S_2_0_0_0[1] - p_over_q * AUX_S_3_0_1_0[10];

                    // ( D_110 S_000 | D_011 S_000 )^0_{t} = z * ( D_110 S_000 | P_010 S_000 )^0 - ( F_111 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[10] = etfac[2] * AUX_S_2_0_1_0[4] - p_over_q * AUX_S_3_0_1_0[13];

                    // ( D_110 S_000 | D_002 S_000 )^0_{t} = z * ( D_110 S_000 | P_001 S_000 )^0 + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | P_001 S_000 )^0
                    AUX_S_2_0_2_0[11] = etfac[2] * AUX_S_2_0_1_0[5] + 1 * one_over_2q * AUX_S_2_0_0_0[1] - p_over_q * AUX_S_3_0_1_0[14];

                    // ( D_101 S_000 | D_200 S_000 )^0_{t} = x * ( D_101 S_000 | P_100 S_000 )^0 + ( P_001 S_000 | P_100 S_000 )^0 + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[12] = etfac[0] * AUX_S_2_0_1_0[6] + 1 * one_over_2q * AUX_S_1_0_1_0[6] + 1 * one_over_2q * AUX_S_2_0_0_0[2] - p_over_q * AUX_S_3_0_1_0[6];

                    // ( D_101 S_000 | D_110 S_000 )^0_{t} = y * ( D_101 S_000 | P_100 S_000 )^0 - ( F_111 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[13] = etfac[1] * AUX_S_2_0_1_0[6] - p_over_q * AUX_S_3_0_1_0[12];

                    // ( D_101 S_000 | D_101 S_000 )^0_{t} = z * ( D_101 S_000 | P_100 S_000 )^0 + ( P_100 S_000 | P_100 S_000 )^0 - ( F_102 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[14] = etfac[2] * AUX_S_2_0_1_0[6] + 1 * one_over_2q * AUX_S_1_0_1_0[0] - p_over_q * AUX_S_3_0_1_0[15];

                    // ( D_101 S_000 | D_020 S_000 )^0_{t} = y * ( D_101 S_000 | P_010 S_000 )^0 + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[15] = etfac[1] * AUX_S_2_0_1_0[7] + 1 * one_over_2q * AUX_S_2_0_0_0[2] - p_over_q * AUX_S_3_0_1_0[13];

                    // ( D_101 S_000 | D_011 S_000 )^0_{t} = z * ( D_101 S_000 | P_010 S_000 )^0 + ( P_100 S_000 | P_010 S_000 )^0 - ( F_102 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[16] = etfac[2] * AUX_S_2_0_1_0[7] + 1 * one_over_2q * AUX_S_1_0_1_0[1] - p_over_q * AUX_S_3_0_1_0[16];

                    // ( D_101 S_000 | D_002 S_000 )^0_{t} = z * ( D_101 S_000 | P_001 S_000 )^0 + ( P_100 S_000 | P_001 S_000 )^0 + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | P_001 S_000 )^0
                    AUX_S_2_0_2_0[17] = etfac[2] * AUX_S_2_0_1_0[8] + 1 * one_over_2q * AUX_S_1_0_1_0[2] + 1 * one_over_2q * AUX_S_2_0_0_0[2] - p_over_q * AUX_S_3_0_1_0[17];

                    // ( D_020 S_000 | D_200 S_000 )^0_{t} = x * ( D_020 S_000 | P_100 S_000 )^0 + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[18] = etfac[0] * AUX_S_2_0_1_0[9] + 1 * one_over_2q * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_3_0_1_0[9];

                    // ( D_020 S_000 | D_110 S_000 )^0_{t} = y * ( D_020 S_000 | P_100 S_000 )^0 + ( P_010 S_000 | P_100 S_000 )^0 - ( F_030 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[19] = etfac[1] * AUX_S_2_0_1_0[9] + 2 * one_over_2q * AUX_S_1_0_1_0[3] - p_over_q * AUX_S_3_0_1_0[18];

                    // ( D_020 S_000 | D_101 S_000 )^0_{t} = z * ( D_020 S_000 | P_100 S_000 )^0 - ( F_021 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[20] = etfac[2] * AUX_S_2_0_1_0[9] - p_over_q * AUX_S_3_0_1_0[21];

                    // ( D_020 S_000 | D_020 S_000 )^0_{t} = y * ( D_020 S_000 | P_010 S_000 )^0 + ( P_010 S_000 | P_010 S_000 )^0 + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_030 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[21] = etfac[1] * AUX_S_2_0_1_0[10] + 2 * one_over_2q * AUX_S_1_0_1_0[4] + 1 * one_over_2q * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_3_0_1_0[19];

                    // ( D_020 S_000 | D_011 S_000 )^0_{t} = z * ( D_020 S_000 | P_010 S_000 )^0 - ( F_021 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[22] = etfac[2] * AUX_S_2_0_1_0[10] - p_over_q * AUX_S_3_0_1_0[22];

                    // ( D_020 S_000 | D_002 S_000 )^0_{t} = z * ( D_020 S_000 | P_001 S_000 )^0 + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | P_001 S_000 )^0
                    AUX_S_2_0_2_0[23] = etfac[2] * AUX_S_2_0_1_0[11] + 1 * one_over_2q * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_3_0_1_0[23];

                    // ( D_011 S_000 | D_200 S_000 )^0_{t} = x * ( D_011 S_000 | P_100 S_000 )^0 + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[24] = etfac[0] * AUX_S_2_0_1_0[12] + 1 * one_over_2q * AUX_S_2_0_0_0[4] - p_over_q * AUX_S_3_0_1_0[12];

                    // ( D_011 S_000 | D_110 S_000 )^0_{t} = y * ( D_011 S_000 | P_100 S_000 )^0 + ( P_001 S_000 | P_100 S_000 )^0 - ( F_021 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[25] = etfac[1] * AUX_S_2_0_1_0[12] + 1 * one_over_2q * AUX_S_1_0_1_0[6] - p_over_q * AUX_S_3_0_1_0[21];

                    // ( D_011 S_000 | D_101 S_000 )^0_{t} = z * ( D_011 S_000 | P_100 S_000 )^0 + ( P_010 S_000 | P_100 S_000 )^0 - ( F_012 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[26] = etfac[2] * AUX_S_2_0_1_0[12] + 1 * one_over_2q * AUX_S_1_0_1_0[3] - p_over_q * AUX_S_3_0_1_0[24];

                    // ( D_011 S_000 | D_020 S_000 )^0_{t} = y * ( D_011 S_000 | P_010 S_000 )^0 + ( P_001 S_000 | P_010 S_000 )^0 + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[27] = etfac[1] * AUX_S_2_0_1_0[13] + 1 * one_over_2q * AUX_S_1_0_1_0[7] + 1 * one_over_2q * AUX_S_2_0_0_0[4] - p_over_q * AUX_S_3_0_1_0[22];

                    // ( D_011 S_000 | D_011 S_000 )^0_{t} = z * ( D_011 S_000 | P_010 S_000 )^0 + ( P_010 S_000 | P_010 S_000 )^0 - ( F_012 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[28] = etfac[2] * AUX_S_2_0_1_0[13] + 1 * one_over_2q * AUX_S_1_0_1_0[4] - p_over_q * AUX_S_3_0_1_0[25];

                    // ( D_011 S_000 | D_002 S_000 )^0_{t} = z * ( D_011 S_000 | P_001 S_000 )^0 + ( P_010 S_000 | P_001 S_000 )^0 + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | P_001 S_000 )^0
                    AUX_S_2_0_2_0[29] = etfac[2] * AUX_S_2_0_1_0[14] + 1 * one_over_2q * AUX_S_1_0_1_0[5] + 1 * one_over_2q * AUX_S_2_0_0_0[4] - p_over_q * AUX_S_3_0_1_0[26];

                    // ( D_002 S_000 | D_200 S_000 )^0_{t} = x * ( D_002 S_000 | P_100 S_000 )^0 + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[30] = etfac[0] * AUX_S_2_0_1_0[15] + 1 * one_over_2q * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_3_0_1_0[15];

                    // ( D_002 S_000 | D_110 S_000 )^0_{t} = y * ( D_002 S_000 | P_100 S_000 )^0 - ( F_012 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[31] = etfac[1] * AUX_S_2_0_1_0[15] - p_over_q * AUX_S_3_0_1_0[24];

                    // ( D_002 S_000 | D_101 S_000 )^0_{t} = z * ( D_002 S_000 | P_100 S_000 )^0 + ( P_001 S_000 | P_100 S_000 )^0 - ( F_003 S_000 | P_100 S_000 )^0
                    AUX_S_2_0_2_0[32] = etfac[2] * AUX_S_2_0_1_0[15] + 2 * one_over_2q * AUX_S_1_0_1_0[6] - p_over_q * AUX_S_3_0_1_0[27];

                    // ( D_002 S_000 | D_020 S_000 )^0_{t} = y * ( D_002 S_000 | P_010 S_000 )^0 + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[33] = etfac[1] * AUX_S_2_0_1_0[16] + 1 * one_over_2q * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_3_0_1_0[25];

                    // ( D_002 S_000 | D_011 S_000 )^0_{t} = z * ( D_002 S_000 | P_010 S_000 )^0 + ( P_001 S_000 | P_010 S_000 )^0 - ( F_003 S_000 | P_010 S_000 )^0
                    AUX_S_2_0_2_0[34] = etfac[2] * AUX_S_2_0_1_0[16] + 2 * one_over_2q * AUX_S_1_0_1_0[7] - p_over_q * AUX_S_3_0_1_0[28];

                    // ( D_002 S_000 | D_002 S_000 )^0_{t} = z * ( D_002 S_000 | P_001 S_000 )^0 + ( P_001 S_000 | P_001 S_000 )^0 + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_003 S_000 | P_001 S_000 )^0
                    AUX_S_2_0_2_0[35] = etfac[2] * AUX_S_2_0_1_0[17] + 2 * one_over_2q * AUX_S_1_0_1_0[8] + 1 * one_over_2q * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_3_0_1_0[29];

                    // ( G_400 S_000 | P_100 S_000 )^0 = x * ( G_400 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_500 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[0] = etfac[0] * AUX_S_4_0_0_0[0] + 4 * one_over_2q * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_5_0_0_0[0];

                    // ( G_310 S_000 | P_100 S_000 )^0 = x * ( G_310 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_410 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[3] = etfac[0] * AUX_S_4_0_0_0[1] + 3 * one_over_2q * AUX_S_3_0_0_0[1] - p_over_q * AUX_S_5_0_0_0[1];

                    // ( G_310 S_000 | P_010 S_000 )^0 = y * ( G_310 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_320 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[4] = etfac[1] * AUX_S_4_0_0_0[1] + 1 * one_over_2q * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_5_0_0_0[3];

                    // ( G_301 S_000 | P_100 S_000 )^0 = x * ( G_301 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_401 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[6] = etfac[0] * AUX_S_4_0_0_0[2] + 3 * one_over_2q * AUX_S_3_0_0_0[2] - p_over_q * AUX_S_5_0_0_0[2];

                    // ( G_301 S_000 | P_010 S_000 )^0 = y * ( G_301 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[7] = etfac[1] * AUX_S_4_0_0_0[2] - p_over_q * AUX_S_5_0_0_0[4];

                    // ( G_301 S_000 | P_001 S_000 )^0 = z * ( G_301 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_302 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[8] = etfac[2] * AUX_S_4_0_0_0[2] + 1 * one_over_2q * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_5_0_0_0[5];

                    // ( G_220 S_000 | P_100 S_000 )^0 = x * ( G_220 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_320 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[9] = etfac[0] * AUX_S_4_0_0_0[3] + 2 * one_over_2q * AUX_S_3_0_0_0[3] - p_over_q * AUX_S_5_0_0_0[3];

                    // ( G_220 S_000 | P_010 S_000 )^0 = y * ( G_220 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_230 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[10] = etfac[1] * AUX_S_4_0_0_0[3] + 2 * one_over_2q * AUX_S_3_0_0_0[1] - p_over_q * AUX_S_5_0_0_0[6];

                    // ( G_211 S_000 | P_100 S_000 )^0 = x * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[12] = etfac[0] * AUX_S_4_0_0_0[4] + 2 * one_over_2q * AUX_S_3_0_0_0[4] - p_over_q * AUX_S_5_0_0_0[4];

                    // ( G_211 S_000 | P_010 S_000 )^0 = y * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[13] = etfac[1] * AUX_S_4_0_0_0[4] + 1 * one_over_2q * AUX_S_3_0_0_0[2] - p_over_q * AUX_S_5_0_0_0[7];

                    // ( G_211 S_000 | P_001 S_000 )^0 = z * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[14] = etfac[2] * AUX_S_4_0_0_0[4] + 1 * one_over_2q * AUX_S_3_0_0_0[1] - p_over_q * AUX_S_5_0_0_0[8];

                    // ( G_202 S_000 | P_100 S_000 )^0 = x * ( G_202 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_302 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[15] = etfac[0] * AUX_S_4_0_0_0[5] + 2 * one_over_2q * AUX_S_3_0_0_0[5] - p_over_q * AUX_S_5_0_0_0[5];

                    // ( G_202 S_000 | P_010 S_000 )^0 = y * ( G_202 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[16] = etfac[1] * AUX_S_4_0_0_0[5] - p_over_q * AUX_S_5_0_0_0[8];

                    // ( G_202 S_000 | P_001 S_000 )^0 = z * ( G_202 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_203 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[17] = etfac[2] * AUX_S_4_0_0_0[5] + 2 * one_over_2q * AUX_S_3_0_0_0[2] - p_over_q * AUX_S_5_0_0_0[9];

                    // ( G_130 S_000 | P_100 S_000 )^0 = x * ( G_130 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_230 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[18] = etfac[0] * AUX_S_4_0_0_0[6] + 1 * one_over_2q * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_5_0_0_0[6];

                    // ( G_130 S_000 | P_010 S_000 )^0 = y * ( G_130 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_140 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[19] = etfac[1] * AUX_S_4_0_0_0[6] + 3 * one_over_2q * AUX_S_3_0_0_0[3] - p_over_q * AUX_S_5_0_0_0[10];

                    // ( G_121 S_000 | P_100 S_000 )^0 = x * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[21] = etfac[0] * AUX_S_4_0_0_0[7] + 1 * one_over_2q * AUX_S_3_0_0_0[7] - p_over_q * AUX_S_5_0_0_0[7];

                    // ( G_121 S_000 | P_010 S_000 )^0 = y * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[22] = etfac[1] * AUX_S_4_0_0_0[7] + 2 * one_over_2q * AUX_S_3_0_0_0[4] - p_over_q * AUX_S_5_0_0_0[11];

                    // ( G_121 S_000 | P_001 S_000 )^0 = z * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[23] = etfac[2] * AUX_S_4_0_0_0[7] + 1 * one_over_2q * AUX_S_3_0_0_0[3] - p_over_q * AUX_S_5_0_0_0[12];

                    // ( G_112 S_000 | P_100 S_000 )^0 = x * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[24] = etfac[0] * AUX_S_4_0_0_0[8] + 1 * one_over_2q * AUX_S_3_0_0_0[8] - p_over_q * AUX_S_5_0_0_0[8];

                    // ( G_112 S_000 | P_010 S_000 )^0 = y * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[25] = etfac[1] * AUX_S_4_0_0_0[8] + 1 * one_over_2q * AUX_S_3_0_0_0[5] - p_over_q * AUX_S_5_0_0_0[12];

                    // ( G_112 S_000 | P_001 S_000 )^0 = z * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[26] = etfac[2] * AUX_S_4_0_0_0[8] + 2 * one_over_2q * AUX_S_3_0_0_0[4] - p_over_q * AUX_S_5_0_0_0[13];

                    // ( G_103 S_000 | P_100 S_000 )^0 = x * ( G_103 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_203 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[27] = etfac[0] * AUX_S_4_0_0_0[9] + 1 * one_over_2q * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_5_0_0_0[9];

                    // ( G_103 S_000 | P_010 S_000 )^0 = y * ( G_103 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[28] = etfac[1] * AUX_S_4_0_0_0[9] - p_over_q * AUX_S_5_0_0_0[13];

                    // ( G_103 S_000 | P_001 S_000 )^0 = z * ( G_103 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_104 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[29] = etfac[2] * AUX_S_4_0_0_0[9] + 3 * one_over_2q * AUX_S_3_0_0_0[5] - p_over_q * AUX_S_5_0_0_0[14];

                    // ( G_040 S_000 | P_100 S_000 )^0 = x * ( G_040 S_000 | S_000 S_000 )^0_{e} - ( H_140 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[30] = etfac[0] * AUX_S_4_0_0_0[10] - p_over_q * AUX_S_5_0_0_0[10];

                    // ( G_040 S_000 | P_010 S_000 )^0 = y * ( G_040 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_050 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[31] = etfac[1] * AUX_S_4_0_0_0[10] + 4 * one_over_2q * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_5_0_0_0[15];

                    // ( G_031 S_000 | P_100 S_000 )^0 = x * ( G_031 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[33] = etfac[0] * AUX_S_4_0_0_0[11] - p_over_q * AUX_S_5_0_0_0[11];

                    // ( G_031 S_000 | P_010 S_000 )^0 = y * ( G_031 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_041 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[34] = etfac[1] * AUX_S_4_0_0_0[11] + 3 * one_over_2q * AUX_S_3_0_0_0[7] - p_over_q * AUX_S_5_0_0_0[16];

                    // ( G_031 S_000 | P_001 S_000 )^0 = z * ( G_031 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_032 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[35] = etfac[2] * AUX_S_4_0_0_0[11] + 1 * one_over_2q * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_5_0_0_0[17];

                    // ( G_022 S_000 | P_100 S_000 )^0 = x * ( G_022 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[36] = etfac[0] * AUX_S_4_0_0_0[12] - p_over_q * AUX_S_5_0_0_0[12];

                    // ( G_022 S_000 | P_010 S_000 )^0 = y * ( G_022 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_032 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[37] = etfac[1] * AUX_S_4_0_0_0[12] + 2 * one_over_2q * AUX_S_3_0_0_0[8] - p_over_q * AUX_S_5_0_0_0[17];

                    // ( G_022 S_000 | P_001 S_000 )^0 = z * ( G_022 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_023 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[38] = etfac[2] * AUX_S_4_0_0_0[12] + 2 * one_over_2q * AUX_S_3_0_0_0[7] - p_over_q * AUX_S_5_0_0_0[18];

                    // ( G_013 S_000 | P_100 S_000 )^0 = x * ( G_013 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[39] = etfac[0] * AUX_S_4_0_0_0[13] - p_over_q * AUX_S_5_0_0_0[13];

                    // ( G_013 S_000 | P_010 S_000 )^0 = y * ( G_013 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_023 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[40] = etfac[1] * AUX_S_4_0_0_0[13] + 1 * one_over_2q * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_5_0_0_0[18];

                    // ( G_013 S_000 | P_001 S_000 )^0 = z * ( G_013 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_014 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[41] = etfac[2] * AUX_S_4_0_0_0[13] + 3 * one_over_2q * AUX_S_3_0_0_0[8] - p_over_q * AUX_S_5_0_0_0[19];

                    // ( G_004 S_000 | P_100 S_000 )^0 = x * ( G_004 S_000 | S_000 S_000 )^0_{e} - ( H_104 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[42] = etfac[0] * AUX_S_4_0_0_0[14] - p_over_q * AUX_S_5_0_0_0[14];

                    // ( G_004 S_000 | P_010 S_000 )^0 = y * ( G_004 S_000 | S_000 S_000 )^0_{e} - ( H_014 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[43] = etfac[1] * AUX_S_4_0_0_0[14] - p_over_q * AUX_S_5_0_0_0[19];

                    // ( G_004 S_000 | P_001 S_000 )^0 = z * ( G_004 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_005 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[44] = etfac[2] * AUX_S_4_0_0_0[14] + 4 * one_over_2q * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_5_0_0_0[20];

                    // ( F_300 S_000 | D_200 S_000 )^0_{t} = x * ( F_300 S_000 | P_100 S_000 )^0 + ( D_200 S_000 | P_100 S_000 )^0 + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_400 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[0] = etfac[0] * AUX_S_3_0_1_0[0] + 3 * one_over_2q * AUX_S_2_0_1_0[0] + 1 * one_over_2q * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_4_0_1_0[0];

                    // ( F_300 S_000 | D_110 S_000 )^0_{t} = y * ( F_300 S_000 | P_100 S_000 )^0 - ( G_310 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[1] = etfac[1] * AUX_S_3_0_1_0[0] - p_over_q * AUX_S_4_0_1_0[3];

                    // ( F_300 S_000 | D_101 S_000 )^0_{t} = z * ( F_300 S_000 | P_100 S_000 )^0 - ( G_301 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[2] = etfac[2] * AUX_S_3_0_1_0[0] - p_over_q * AUX_S_4_0_1_0[6];

                    // ( F_300 S_000 | D_020 S_000 )^0_{t} = y * ( F_300 S_000 | P_010 S_000 )^0 + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[3] = etfac[1] * AUX_S_3_0_1_0[1] + 1 * one_over_2q * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_4_0_1_0[4];

                    // ( F_300 S_000 | D_011 S_000 )^0_{t} = z * ( F_300 S_000 | P_010 S_000 )^0 - ( G_301 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[4] = etfac[2] * AUX_S_3_0_1_0[1] - p_over_q * AUX_S_4_0_1_0[7];

                    // ( F_300 S_000 | D_002 S_000 )^0_{t} = z * ( F_300 S_000 | P_001 S_000 )^0 + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | P_001 S_000 )^0
                    AUX_S_3_0_2_0[5] = etfac[2] * AUX_S_3_0_1_0[2] + 1 * one_over_2q * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_4_0_1_0[8];

                    // ( F_210 S_000 | D_200 S_000 )^0_{t} = x * ( F_210 S_000 | P_100 S_000 )^0 + ( D_110 S_000 | P_100 S_000 )^0 + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[6] = etfac[0] * AUX_S_3_0_1_0[3] + 2 * one_over_2q * AUX_S_2_0_1_0[3] + 1 * one_over_2q * AUX_S_3_0_0_0[1] - p_over_q * AUX_S_4_0_1_0[3];

                    // ( F_210 S_000 | D_110 S_000 )^0_{t} = y * ( F_210 S_000 | P_100 S_000 )^0 + ( D_200 S_000 | P_100 S_000 )^0 - ( G_220 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[7] = etfac[1] * AUX_S_3_0_1_0[3] + 1 * one_over_2q * AUX_S_2_0_1_0[0] - p_over_q * AUX_S_4_0_1_0[9];

                    // ( F_210 S_000 | D_101 S_000 )^0_{t} = z * ( F_210 S_000 | P_100 S_000 )^0 - ( G_211 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[8] = etfac[2] * AUX_S_3_0_1_0[3] - p_over_q * AUX_S_4_0_1_0[12];

                    // ( F_210 S_000 | D_020 S_000 )^0_{t} = y * ( F_210 S_000 | P_010 S_000 )^0 + ( D_200 S_000 | P_010 S_000 )^0 + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[9] = etfac[1] * AUX_S_3_0_1_0[4] + 1 * one_over_2q * AUX_S_2_0_1_0[1] + 1 * one_over_2q * AUX_S_3_0_0_0[1] - p_over_q * AUX_S_4_0_1_0[10];

                    // ( F_210 S_000 | D_011 S_000 )^0_{t} = z * ( F_210 S_000 | P_010 S_000 )^0 - ( G_211 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[10] = etfac[2] * AUX_S_3_0_1_0[4] - p_over_q * AUX_S_4_0_1_0[13];

                    // ( F_210 S_000 | D_002 S_000 )^0_{t} = z * ( F_210 S_000 | P_001 S_000 )^0 + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | P_001 S_000 )^0
                    AUX_S_3_0_2_0[11] = etfac[2] * AUX_S_3_0_1_0[5] + 1 * one_over_2q * AUX_S_3_0_0_0[1] - p_over_q * AUX_S_4_0_1_0[14];

                    // ( F_201 S_000 | D_200 S_000 )^0_{t} = x * ( F_201 S_000 | P_100 S_000 )^0 + ( D_101 S_000 | P_100 S_000 )^0 + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[12] = etfac[0] * AUX_S_3_0_1_0[6] + 2 * one_over_2q * AUX_S_2_0_1_0[6] + 1 * one_over_2q * AUX_S_3_0_0_0[2] - p_over_q * AUX_S_4_0_1_0[6];

                    // ( F_201 S_000 | D_110 S_000 )^0_{t} = y * ( F_201 S_000 | P_100 S_000 )^0 - ( G_211 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[13] = etfac[1] * AUX_S_3_0_1_0[6] - p_over_q * AUX_S_4_0_1_0[12];

                    // ( F_201 S_000 | D_101 S_000 )^0_{t} = z * ( F_201 S_000 | P_100 S_000 )^0 + ( D_200 S_000 | P_100 S_000 )^0 - ( G_202 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[14] = etfac[2] * AUX_S_3_0_1_0[6] + 1 * one_over_2q * AUX_S_2_0_1_0[0] - p_over_q * AUX_S_4_0_1_0[15];

                    // ( F_201 S_000 | D_020 S_000 )^0_{t} = y * ( F_201 S_000 | P_010 S_000 )^0 + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[15] = etfac[1] * AUX_S_3_0_1_0[7] + 1 * one_over_2q * AUX_S_3_0_0_0[2] - p_over_q * AUX_S_4_0_1_0[13];

                    // ( F_201 S_000 | D_011 S_000 )^0_{t} = z * ( F_201 S_000 | P_010 S_000 )^0 + ( D_200 S_000 | P_010 S_000 )^0 - ( G_202 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[16] = etfac[2] * AUX_S_3_0_1_0[7] + 1 * one_over_2q * AUX_S_2_0_1_0[1] - p_over_q * AUX_S_4_0_1_0[16];

                    // ( F_201 S_000 | D_002 S_000 )^0_{t} = z * ( F_201 S_000 | P_001 S_000 )^0 + ( D_200 S_000 | P_001 S_000 )^0 + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | P_001 S_000 )^0
                    AUX_S_3_0_2_0[17] = etfac[2] * AUX_S_3_0_1_0[8] + 1 * one_over_2q * AUX_S_2_0_1_0[2] + 1 * one_over_2q * AUX_S_3_0_0_0[2] - p_over_q * AUX_S_4_0_1_0[17];

                    // ( F_120 S_000 | D_200 S_000 )^0_{t} = x * ( F_120 S_000 | P_100 S_000 )^0 + ( D_020 S_000 | P_100 S_000 )^0 + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[18] = etfac[0] * AUX_S_3_0_1_0[9] + 1 * one_over_2q * AUX_S_2_0_1_0[9] + 1 * one_over_2q * AUX_S_3_0_0_0[3] - p_over_q * AUX_S_4_0_1_0[9];

                    // ( F_120 S_000 | D_110 S_000 )^0_{t} = y * ( F_120 S_000 | P_100 S_000 )^0 + ( D_110 S_000 | P_100 S_000 )^0 - ( G_130 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[19] = etfac[1] * AUX_S_3_0_1_0[9] + 2 * one_over_2q * AUX_S_2_0_1_0[3] - p_over_q * AUX_S_4_0_1_0[18];

                    // ( F_120 S_000 | D_101 S_000 )^0_{t} = z * ( F_120 S_000 | P_100 S_000 )^0 - ( G_121 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[20] = etfac[2] * AUX_S_3_0_1_0[9] - p_over_q * AUX_S_4_0_1_0[21];

                    // ( F_120 S_000 | D_020 S_000 )^0_{t} = y * ( F_120 S_000 | P_010 S_000 )^0 + ( D_110 S_000 | P_010 S_000 )^0 + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[21] = etfac[1] * AUX_S_3_0_1_0[10] + 2 * one_over_2q * AUX_S_2_0_1_0[4] + 1 * one_over_2q * AUX_S_3_0_0_0[3] - p_over_q * AUX_S_4_0_1_0[19];

                    // ( F_120 S_000 | D_011 S_000 )^0_{t} = z * ( F_120 S_000 | P_010 S_000 )^0 - ( G_121 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[22] = etfac[2] * AUX_S_3_0_1_0[10] - p_over_q * AUX_S_4_0_1_0[22];

                    // ( F_120 S_000 | D_002 S_000 )^0_{t} = z * ( F_120 S_000 | P_001 S_000 )^0 + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | P_001 S_000 )^0
                    AUX_S_3_0_2_0[23] = etfac[2] * AUX_S_3_0_1_0[11] + 1 * one_over_2q * AUX_S_3_0_0_0[3] - p_over_q * AUX_S_4_0_1_0[23];

                    // ( F_111 S_000 | D_200 S_000 )^0_{t} = x * ( F_111 S_000 | P_100 S_000 )^0 + ( D_011 S_000 | P_100 S_000 )^0 + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[24] = etfac[0] * AUX_S_3_0_1_0[12] + 1 * one_over_2q * AUX_S_2_0_1_0[12] + 1 * one_over_2q * AUX_S_3_0_0_0[4] - p_over_q * AUX_S_4_0_1_0[12];

                    // ( F_111 S_000 | D_110 S_000 )^0_{t} = y * ( F_111 S_000 | P_100 S_000 )^0 + ( D_101 S_000 | P_100 S_000 )^0 - ( G_121 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[25] = etfac[1] * AUX_S_3_0_1_0[12] + 1 * one_over_2q * AUX_S_2_0_1_0[6] - p_over_q * AUX_S_4_0_1_0[21];

                    // ( F_111 S_000 | D_101 S_000 )^0_{t} = z * ( F_111 S_000 | P_100 S_000 )^0 + ( D_110 S_000 | P_100 S_000 )^0 - ( G_112 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[26] = etfac[2] * AUX_S_3_0_1_0[12] + 1 * one_over_2q * AUX_S_2_0_1_0[3] - p_over_q * AUX_S_4_0_1_0[24];

                    // ( F_111 S_000 | D_020 S_000 )^0_{t} = y * ( F_111 S_000 | P_010 S_000 )^0 + ( D_101 S_000 | P_010 S_000 )^0 + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[27] = etfac[1] * AUX_S_3_0_1_0[13] + 1 * one_over_2q * AUX_S_2_0_1_0[7] + 1 * one_over_2q * AUX_S_3_0_0_0[4] - p_over_q * AUX_S_4_0_1_0[22];

                    // ( F_111 S_000 | D_011 S_000 )^0_{t} = z * ( F_111 S_000 | P_010 S_000 )^0 + ( D_110 S_000 | P_010 S_000 )^0 - ( G_112 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[28] = etfac[2] * AUX_S_3_0_1_0[13] + 1 * one_over_2q * AUX_S_2_0_1_0[4] - p_over_q * AUX_S_4_0_1_0[25];

                    // ( F_111 S_000 | D_002 S_000 )^0_{t} = z * ( F_111 S_000 | P_001 S_000 )^0 + ( D_110 S_000 | P_001 S_000 )^0 + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | P_001 S_000 )^0
                    AUX_S_3_0_2_0[29] = etfac[2] * AUX_S_3_0_1_0[14] + 1 * one_over_2q * AUX_S_2_0_1_0[5] + 1 * one_over_2q * AUX_S_3_0_0_0[4] - p_over_q * AUX_S_4_0_1_0[26];

                    // ( F_102 S_000 | D_200 S_000 )^0_{t} = x * ( F_102 S_000 | P_100 S_000 )^0 + ( D_002 S_000 | P_100 S_000 )^0 + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[30] = etfac[0] * AUX_S_3_0_1_0[15] + 1 * one_over_2q * AUX_S_2_0_1_0[15] + 1 * one_over_2q * AUX_S_3_0_0_0[5] - p_over_q * AUX_S_4_0_1_0[15];

                    // ( F_102 S_000 | D_110 S_000 )^0_{t} = y * ( F_102 S_000 | P_100 S_000 )^0 - ( G_112 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[31] = etfac[1] * AUX_S_3_0_1_0[15] - p_over_q * AUX_S_4_0_1_0[24];

                    // ( F_102 S_000 | D_101 S_000 )^0_{t} = z * ( F_102 S_000 | P_100 S_000 )^0 + ( D_101 S_000 | P_100 S_000 )^0 - ( G_103 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[32] = etfac[2] * AUX_S_3_0_1_0[15] + 2 * one_over_2q * AUX_S_2_0_1_0[6] - p_over_q * AUX_S_4_0_1_0[27];

                    // ( F_102 S_000 | D_020 S_000 )^0_{t} = y * ( F_102 S_000 | P_010 S_000 )^0 + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[33] = etfac[1] * AUX_S_3_0_1_0[16] + 1 * one_over_2q * AUX_S_3_0_0_0[5] - p_over_q * AUX_S_4_0_1_0[25];

                    // ( F_102 S_000 | D_011 S_000 )^0_{t} = z * ( F_102 S_000 | P_010 S_000 )^0 + ( D_101 S_000 | P_010 S_000 )^0 - ( G_103 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[34] = etfac[2] * AUX_S_3_0_1_0[16] + 2 * one_over_2q * AUX_S_2_0_1_0[7] - p_over_q * AUX_S_4_0_1_0[28];

                    // ( F_102 S_000 | D_002 S_000 )^0_{t} = z * ( F_102 S_000 | P_001 S_000 )^0 + ( D_101 S_000 | P_001 S_000 )^0 + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | P_001 S_000 )^0
                    AUX_S_3_0_2_0[35] = etfac[2] * AUX_S_3_0_1_0[17] + 2 * one_over_2q * AUX_S_2_0_1_0[8] + 1 * one_over_2q * AUX_S_3_0_0_0[5] - p_over_q * AUX_S_4_0_1_0[29];

                    // ( F_030 S_000 | D_200 S_000 )^0_{t} = x * ( F_030 S_000 | P_100 S_000 )^0 + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[36] = etfac[0] * AUX_S_3_0_1_0[18] + 1 * one_over_2q * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_4_0_1_0[18];

                    // ( F_030 S_000 | D_110 S_000 )^0_{t} = y * ( F_030 S_000 | P_100 S_000 )^0 + ( D_020 S_000 | P_100 S_000 )^0 - ( G_040 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[37] = etfac[1] * AUX_S_3_0_1_0[18] + 3 * one_over_2q * AUX_S_2_0_1_0[9] - p_over_q * AUX_S_4_0_1_0[30];

                    // ( F_030 S_000 | D_101 S_000 )^0_{t} = z * ( F_030 S_000 | P_100 S_000 )^0 - ( G_031 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[38] = etfac[2] * AUX_S_3_0_1_0[18] - p_over_q * AUX_S_4_0_1_0[33];

                    // ( F_030 S_000 | D_020 S_000 )^0_{t} = y * ( F_030 S_000 | P_010 S_000 )^0 + ( D_020 S_000 | P_010 S_000 )^0 + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_040 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[39] = etfac[1] * AUX_S_3_0_1_0[19] + 3 * one_over_2q * AUX_S_2_0_1_0[10] + 1 * one_over_2q * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_4_0_1_0[31];

                    // ( F_030 S_000 | D_011 S_000 )^0_{t} = z * ( F_030 S_000 | P_010 S_000 )^0 - ( G_031 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[40] = etfac[2] * AUX_S_3_0_1_0[19] - p_over_q * AUX_S_4_0_1_0[34];

                    // ( F_030 S_000 | D_002 S_000 )^0_{t} = z * ( F_030 S_000 | P_001 S_000 )^0 + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | P_001 S_000 )^0
                    AUX_S_3_0_2_0[41] = etfac[2] * AUX_S_3_0_1_0[20] + 1 * one_over_2q * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_4_0_1_0[35];

                    // ( F_021 S_000 | D_200 S_000 )^0_{t} = x * ( F_021 S_000 | P_100 S_000 )^0 + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[42] = etfac[0] * AUX_S_3_0_1_0[21] + 1 * one_over_2q * AUX_S_3_0_0_0[7] - p_over_q * AUX_S_4_0_1_0[21];

                    // ( F_021 S_000 | D_110 S_000 )^0_{t} = y * ( F_021 S_000 | P_100 S_000 )^0 + ( D_011 S_000 | P_100 S_000 )^0 - ( G_031 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[43] = etfac[1] * AUX_S_3_0_1_0[21] + 2 * one_over_2q * AUX_S_2_0_1_0[12] - p_over_q * AUX_S_4_0_1_0[33];

                    // ( F_021 S_000 | D_101 S_000 )^0_{t} = z * ( F_021 S_000 | P_100 S_000 )^0 + ( D_020 S_000 | P_100 S_000 )^0 - ( G_022 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[44] = etfac[2] * AUX_S_3_0_1_0[21] + 1 * one_over_2q * AUX_S_2_0_1_0[9] - p_over_q * AUX_S_4_0_1_0[36];

                    // ( F_021 S_000 | D_020 S_000 )^0_{t} = y * ( F_021 S_000 | P_010 S_000 )^0 + ( D_011 S_000 | P_010 S_000 )^0 + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[45] = etfac[1] * AUX_S_3_0_1_0[22] + 2 * one_over_2q * AUX_S_2_0_1_0[13] + 1 * one_over_2q * AUX_S_3_0_0_0[7] - p_over_q * AUX_S_4_0_1_0[34];

                    // ( F_021 S_000 | D_011 S_000 )^0_{t} = z * ( F_021 S_000 | P_010 S_000 )^0 + ( D_020 S_000 | P_010 S_000 )^0 - ( G_022 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[46] = etfac[2] * AUX_S_3_0_1_0[22] + 1 * one_over_2q * AUX_S_2_0_1_0[10] - p_over_q * AUX_S_4_0_1_0[37];

                    // ( F_021 S_000 | D_002 S_000 )^0_{t} = z * ( F_021 S_000 | P_001 S_000 )^0 + ( D_020 S_000 | P_001 S_000 )^0 + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | P_001 S_000 )^0
                    AUX_S_3_0_2_0[47] = etfac[2] * AUX_S_3_0_1_0[23] + 1 * one_over_2q * AUX_S_2_0_1_0[11] + 1 * one_over_2q * AUX_S_3_0_0_0[7] - p_over_q * AUX_S_4_0_1_0[38];

                    // ( F_012 S_000 | D_200 S_000 )^0_{t} = x * ( F_012 S_000 | P_100 S_000 )^0 + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[48] = etfac[0] * AUX_S_3_0_1_0[24] + 1 * one_over_2q * AUX_S_3_0_0_0[8] - p_over_q * AUX_S_4_0_1_0[24];

                    // ( F_012 S_000 | D_110 S_000 )^0_{t} = y * ( F_012 S_000 | P_100 S_000 )^0 + ( D_002 S_000 | P_100 S_000 )^0 - ( G_022 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[49] = etfac[1] * AUX_S_3_0_1_0[24] + 1 * one_over_2q * AUX_S_2_0_1_0[15] - p_over_q * AUX_S_4_0_1_0[36];

                    // ( F_012 S_000 | D_101 S_000 )^0_{t} = z * ( F_012 S_000 | P_100 S_000 )^0 + ( D_011 S_000 | P_100 S_000 )^0 - ( G_013 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[50] = etfac[2] * AUX_S_3_0_1_0[24] + 2 * one_over_2q * AUX_S_2_0_1_0[12] - p_over_q * AUX_S_4_0_1_0[39];

                    // ( F_012 S_000 | D_020 S_000 )^0_{t} = y * ( F_012 S_000 | P_010 S_000 )^0 + ( D_002 S_000 | P_010 S_000 )^0 + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[51] = etfac[1] * AUX_S_3_0_1_0[25] + 1 * one_over_2q * AUX_S_2_0_1_0[16] + 1 * one_over_2q * AUX_S_3_0_0_0[8] - p_over_q * AUX_S_4_0_1_0[37];

                    // ( F_012 S_000 | D_011 S_000 )^0_{t} = z * ( F_012 S_000 | P_010 S_000 )^0 + ( D_011 S_000 | P_010 S_000 )^0 - ( G_013 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[52] = etfac[2] * AUX_S_3_0_1_0[25] + 2 * one_over_2q * AUX_S_2_0_1_0[13] - p_over_q * AUX_S_4_0_1_0[40];

                    // ( F_012 S_000 | D_002 S_000 )^0_{t} = z * ( F_012 S_000 | P_001 S_000 )^0 + ( D_011 S_000 | P_001 S_000 )^0 + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | P_001 S_000 )^0
                    AUX_S_3_0_2_0[53] = etfac[2] * AUX_S_3_0_1_0[26] + 2 * one_over_2q * AUX_S_2_0_1_0[14] + 1 * one_over_2q * AUX_S_3_0_0_0[8] - p_over_q * AUX_S_4_0_1_0[41];

                    // ( F_003 S_000 | D_200 S_000 )^0_{t} = x * ( F_003 S_000 | P_100 S_000 )^0 + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[54] = etfac[0] * AUX_S_3_0_1_0[27] + 1 * one_over_2q * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_4_0_1_0[27];

                    // ( F_003 S_000 | D_110 S_000 )^0_{t} = y * ( F_003 S_000 | P_100 S_000 )^0 - ( G_013 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[55] = etfac[1] * AUX_S_3_0_1_0[27] - p_over_q * AUX_S_4_0_1_0[39];

                    // ( F_003 S_000 | D_101 S_000 )^0_{t} = z * ( F_003 S_000 | P_100 S_000 )^0 + ( D_002 S_000 | P_100 S_000 )^0 - ( G_004 S_000 | P_100 S_000 )^0
                    AUX_S_3_0_2_0[56] = etfac[2] * AUX_S_3_0_1_0[27] + 3 * one_over_2q * AUX_S_2_0_1_0[15] - p_over_q * AUX_S_4_0_1_0[42];

                    // ( F_003 S_000 | D_020 S_000 )^0_{t} = y * ( F_003 S_000 | P_010 S_000 )^0 + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[57] = etfac[1] * AUX_S_3_0_1_0[28] + 1 * one_over_2q * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_4_0_1_0[40];

                    // ( F_003 S_000 | D_011 S_000 )^0_{t} = z * ( F_003 S_000 | P_010 S_000 )^0 + ( D_002 S_000 | P_010 S_000 )^0 - ( G_004 S_000 | P_010 S_000 )^0
                    AUX_S_3_0_2_0[58] = etfac[2] * AUX_S_3_0_1_0[28] + 3 * one_over_2q * AUX_S_2_0_1_0[16] - p_over_q * AUX_S_4_0_1_0[43];

                    // ( F_003 S_000 | D_002 S_000 )^0_{t} = z * ( F_003 S_000 | P_001 S_000 )^0 + ( D_002 S_000 | P_001 S_000 )^0 + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_004 S_000 | P_001 S_000 )^0
                    AUX_S_3_0_2_0[59] = etfac[2] * AUX_S_3_0_1_0[29] + 3 * one_over_2q * AUX_S_2_0_1_0[17] + 1 * one_over_2q * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_4_0_1_0[44];


                    // Accumulating in contracted workspace
                    for(int i = 0; i < 36; i++)
                        PRIM_S_2_0_2_0[i] += AUX_S_2_0_2_0[i];

                    // Accumulating in contracted workspace
                    for(int i = 0; i < 60; i++)
                        PRIM_S_3_0_2_0[i] += AUX_S_3_0_2_0[i];

                 }
            }
        }
    }


    //////////////////////////////////////////////
    // Contracted integrals: Horizontal recurrance
    // Bra part
    // Steps: 18
    //////////////////////////////////////////////

    for(abcd = 0; abcd < nshell1234; ++abcd)
    {
        // form S_2_1_2_0
        for(int iket = 0; iket < 6; ++iket)
        {
            // (D_200 P_100|_{i} = (F_300 S_000|_{t} + x_ab * (D_200 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 0 * 6 + iket] = S_3_0_2_0[abcd * 60 + 0 * 6 + iket] + ( AB_x[abcd] * S_2_0_2_0[abcd * 36 + 0 * 6 + iket] );

            // (D_200 P_010|_{i} = (F_210 S_000|_{t} + y_ab * (D_200 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 1 * 6 + iket] = S_3_0_2_0[abcd * 60 + 1 * 6 + iket] + ( AB_y[abcd] * S_2_0_2_0[abcd * 36 + 0 * 6 + iket] );

            // (D_200 P_001|_{i} = (F_201 S_000|_{t} + z_ab * (D_200 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 2 * 6 + iket] = S_3_0_2_0[abcd * 60 + 2 * 6 + iket] + ( AB_z[abcd] * S_2_0_2_0[abcd * 36 + 0 * 6 + iket] );

            // (D_110 P_100|_{i} = (F_210 S_000|_{t} + x_ab * (D_110 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 3 * 6 + iket] = S_3_0_2_0[abcd * 60 + 1 * 6 + iket] + ( AB_x[abcd] * S_2_0_2_0[abcd * 36 + 1 * 6 + iket] );

            // (D_110 P_010|_{i} = (F_120 S_000|_{t} + y_ab * (D_110 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 4 * 6 + iket] = S_3_0_2_0[abcd * 60 + 3 * 6 + iket] + ( AB_y[abcd] * S_2_0_2_0[abcd * 36 + 1 * 6 + iket] );

            // (D_110 P_001|_{i} = (F_111 S_000|_{t} + z_ab * (D_110 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 5 * 6 + iket] = S_3_0_2_0[abcd * 60 + 4 * 6 + iket] + ( AB_z[abcd] * S_2_0_2_0[abcd * 36 + 1 * 6 + iket] );

            // (D_101 P_100|_{i} = (F_201 S_000|_{t} + x_ab * (D_101 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 6 * 6 + iket] = S_3_0_2_0[abcd * 60 + 2 * 6 + iket] + ( AB_x[abcd] * S_2_0_2_0[abcd * 36 + 2 * 6 + iket] );

            // (D_101 P_010|_{i} = (F_111 S_000|_{t} + y_ab * (D_101 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 7 * 6 + iket] = S_3_0_2_0[abcd * 60 + 4 * 6 + iket] + ( AB_y[abcd] * S_2_0_2_0[abcd * 36 + 2 * 6 + iket] );

            // (D_101 P_001|_{i} = (F_102 S_000|_{t} + z_ab * (D_101 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 8 * 6 + iket] = S_3_0_2_0[abcd * 60 + 5 * 6 + iket] + ( AB_z[abcd] * S_2_0_2_0[abcd * 36 + 2 * 6 + iket] );

            // (D_020 P_100|_{i} = (F_120 S_000|_{t} + x_ab * (D_020 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 9 * 6 + iket] = S_3_0_2_0[abcd * 60 + 3 * 6 + iket] + ( AB_x[abcd] * S_2_0_2_0[abcd * 36 + 3 * 6 + iket] );

            // (D_020 P_010|_{i} = (F_030 S_000|_{t} + y_ab * (D_020 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 10 * 6 + iket] = S_3_0_2_0[abcd * 60 + 6 * 6 + iket] + ( AB_y[abcd] * S_2_0_2_0[abcd * 36 + 3 * 6 + iket] );

            // (D_020 P_001|_{i} = (F_021 S_000|_{t} + z_ab * (D_020 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 11 * 6 + iket] = S_3_0_2_0[abcd * 60 + 7 * 6 + iket] + ( AB_z[abcd] * S_2_0_2_0[abcd * 36 + 3 * 6 + iket] );

            // (D_011 P_100|_{i} = (F_111 S_000|_{t} + x_ab * (D_011 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 12 * 6 + iket] = S_3_0_2_0[abcd * 60 + 4 * 6 + iket] + ( AB_x[abcd] * S_2_0_2_0[abcd * 36 + 4 * 6 + iket] );

            // (D_011 P_010|_{i} = (F_021 S_000|_{t} + y_ab * (D_011 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 13 * 6 + iket] = S_3_0_2_0[abcd * 60 + 7 * 6 + iket] + ( AB_y[abcd] * S_2_0_2_0[abcd * 36 + 4 * 6 + iket] );

            // (D_011 P_001|_{i} = (F_012 S_000|_{t} + z_ab * (D_011 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 14 * 6 + iket] = S_3_0_2_0[abcd * 60 + 8 * 6 + iket] + ( AB_z[abcd] * S_2_0_2_0[abcd * 36 + 4 * 6 + iket] );

            // (D_002 P_100|_{i} = (F_102 S_000|_{t} + x_ab * (D_002 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 15 * 6 + iket] = S_3_0_2_0[abcd * 60 + 5 * 6 + iket] + ( AB_x[abcd] * S_2_0_2_0[abcd * 36 + 5 * 6 + iket] );

            // (D_002 P_010|_{i} = (F_012 S_000|_{t} + y_ab * (D_002 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 16 * 6 + iket] = S_3_0_2_0[abcd * 60 + 8 * 6 + iket] + ( AB_y[abcd] * S_2_0_2_0[abcd * 36 + 5 * 6 + iket] );

            // (D_002 P_001|_{i} = (F_003 S_000|_{t} + z_ab * (D_002 S_000|_{t}
            S_2_1_2_0[abcd * 108 + 17 * 6 + iket] = S_3_0_2_0[abcd * 60 + 9 * 6 + iket] + ( AB_z[abcd] * S_2_0_2_0[abcd * 36 + 5 * 6 + iket] );

        }


    }


    //////////////////////////////////////////////
    // Contracted integrals: Horizontal recurrance
    // Ket part
    // Steps: 0
    // Forming final integrals
    //////////////////////////////////////////////

    //Nothing to do.....


    // Free contracted work space
    free(contwork);

    return nshell1234;
}

