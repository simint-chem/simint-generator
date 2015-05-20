#include <string.h>
#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"
#include "boys/boys_split.h"


int eri_split_d_d_p_s(struct multishell_pair const P,
                      struct multishell_pair const Q,
                      double * const restrict S_2_2_1_0)
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

    memset(S_2_2_1_0, 0, nshell1234*108*sizeof(double));

    // Holds AB_{xyz} and CD_{xyz} in a flattened fashion for later
    double AB_x[nshell1234];  double CD_x[nshell1234];
    double AB_y[nshell1234];  double CD_y[nshell1234];
    double AB_z[nshell1234];  double CD_z[nshell1234];

    int ab, cd, abcd;
    int i, j;

    // Workspace for contracted integrals
    double * const contwork = malloc(nshell1234 * 744);
    memset(contwork, 0, nshell1234 * 744);

    // partition workspace into shells
    double * const restrict S_2_0_1_0 = contwork + (nshell1234 * 0);
    double * const restrict S_3_0_1_0 = contwork + (nshell1234 * 18);
    double * const restrict S_4_0_1_0 = contwork + (nshell1234 * 48);


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
                    double AUX_S_2_0_1_0[18];
                    double AUX_S_3_0_1_0[30];
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
                                      +1 * one_over_2p * ( AUX_S_0_0_0_0[m * 1 +  0] - a_over_p * AUX_S_0_0_0_0[(m+1) * 1 + 0] );

                        //D_110 : STEP: y
                        AUX_S_2_0_0_0[m * 6 + 1] = P.PA_y[i] * AUX_S_1_0_0_0[m * 3 + 0] - a_over_p * PQ_y * AUX_S_1_0_0_0[(m+1) * 3 + 0];

                        //D_101 : STEP: z
                        AUX_S_2_0_0_0[m * 6 + 2] = P.PA_z[i] * AUX_S_1_0_0_0[m * 3 + 0] - a_over_p * PQ_z * AUX_S_1_0_0_0[(m+1) * 3 + 0];

                        //D_020 : STEP: y
                        AUX_S_2_0_0_0[m * 6 + 3] = P.PA_y[i] * AUX_S_1_0_0_0[m * 3 + 1] - a_over_p * PQ_y * AUX_S_1_0_0_0[(m+1) * 3 + 1]
                                      +1 * one_over_2p * ( AUX_S_0_0_0_0[m * 1 +  0] - a_over_p * AUX_S_0_0_0_0[(m+1) * 1 + 0] );

                        //D_011 : STEP: z
                        AUX_S_2_0_0_0[m * 6 + 4] = P.PA_z[i] * AUX_S_1_0_0_0[m * 3 + 1] - a_over_p * PQ_z * AUX_S_1_0_0_0[(m+1) * 3 + 1];

                        //D_002 : STEP: z
                        AUX_S_2_0_0_0[m * 6 + 5] = P.PA_z[i] * AUX_S_1_0_0_0[m * 3 + 2] - a_over_p * PQ_z * AUX_S_1_0_0_0[(m+1) * 3 + 2]
                                      +1 * one_over_2p * ( AUX_S_0_0_0_0[m * 1 +  0] - a_over_p * AUX_S_0_0_0_0[(m+1) * 1 + 0] );

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
                                      +2 * one_over_2p * ( AUX_S_1_0_0_0[m * 3 +  0] - a_over_p * AUX_S_1_0_0_0[(m+1) * 3 + 0] );

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
                                      +2 * one_over_2p * ( AUX_S_1_0_0_0[m * 3 +  1] - a_over_p * AUX_S_1_0_0_0[(m+1) * 3 + 1] );

                        //F_021 : STEP: z
                        AUX_S_3_0_0_0[m * 10 + 7] = P.PA_z[i] * AUX_S_2_0_0_0[m * 6 + 3] - a_over_p * PQ_z * AUX_S_2_0_0_0[(m+1) * 6 + 3];

                        //F_012 : STEP: y
                        AUX_S_3_0_0_0[m * 10 + 8] = P.PA_y[i] * AUX_S_2_0_0_0[m * 6 + 5] - a_over_p * PQ_y * AUX_S_2_0_0_0[(m+1) * 6 + 5];

                        //F_003 : STEP: z
                        AUX_S_3_0_0_0[m * 10 + 9] = P.PA_z[i] * AUX_S_2_0_0_0[m * 6 + 5] - a_over_p * PQ_z * AUX_S_2_0_0_0[(m+1) * 6 + 5]
                                      +2 * one_over_2p * ( AUX_S_1_0_0_0[m * 3 +  2] - a_over_p * AUX_S_1_0_0_0[(m+1) * 3 + 2] );

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
                                      +3 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  0] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 0] );

                        //G_310 : STEP: y
                        AUX_S_4_0_0_0[m * 15 + 1] = P.PA_y[i] * AUX_S_3_0_0_0[m * 10 + 0] - a_over_p * PQ_y * AUX_S_3_0_0_0[(m+1) * 10 + 0];

                        //G_301 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 2] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 0] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 0];

                        //G_220 : STEP: y
                        AUX_S_4_0_0_0[m * 15 + 3] = P.PA_y[i] * AUX_S_3_0_0_0[m * 10 + 1] - a_over_p * PQ_y * AUX_S_3_0_0_0[(m+1) * 10 + 1]
                                      +1 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  0] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 0] );

                        //G_211 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 4] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 1] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 1];

                        //G_202 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 5] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 2] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 2]
                                      +1 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  0] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 0] );

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
                                      +3 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  3] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 3] );

                        //G_031 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 11] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 6] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 6];

                        //G_022 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 12] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 7] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 7]
                                      +1 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  3] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 3] );

                        //G_013 : STEP: y
                        AUX_S_4_0_0_0[m * 15 + 13] = P.PA_y[i] * AUX_S_3_0_0_0[m * 10 + 9] - a_over_p * PQ_y * AUX_S_3_0_0_0[(m+1) * 10 + 9];

                        //G_004 : STEP: z
                        AUX_S_4_0_0_0[m * 15 + 14] = P.PA_z[i] * AUX_S_3_0_0_0[m * 10 + 9] - a_over_p * PQ_z * AUX_S_3_0_0_0[(m+1) * 10 + 9]
                                      +3 * one_over_2p * ( AUX_S_2_0_0_0[m * 6 +  5] - a_over_p * AUX_S_2_0_0_0[(m+1) * 6 + 5] );

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
                                      +4 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  0] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 0] );

                        //H_410 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 1] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 0] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 0];

                        //H_401 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 2] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 0] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 0];

                        //H_320 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 3] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 1] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 1]
                                      +1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  0] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 0] );

                        //H_311 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 4] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 1] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 1];

                        //H_302 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 5] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 2] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 2]
                                      +1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  0] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 0] );

                        //H_230 : STEP: x
                        AUX_S_5_0_0_0[m * 21 + 6] = P.PA_x[i] * AUX_S_4_0_0_0[m * 15 + 6] - a_over_p * PQ_x * AUX_S_4_0_0_0[(m+1) * 15 + 6]
                                      +1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  6] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 6] );

                        //H_221 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 7] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 3] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 3];

                        //H_212 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 8] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 5] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 5];

                        //H_203 : STEP: x
                        AUX_S_5_0_0_0[m * 21 + 9] = P.PA_x[i] * AUX_S_4_0_0_0[m * 15 + 9] - a_over_p * PQ_x * AUX_S_4_0_0_0[(m+1) * 15 + 9]
                                      +1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  9] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 9] );

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
                                      +4 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  6] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 6] );

                        //H_041 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 16] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 10] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 10];

                        //H_032 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 17] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 11] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 11]
                                      +1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  6] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 6] );

                        //H_023 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 18] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 13] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 13]
                                      +1 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  9] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 9] );

                        //H_014 : STEP: y
                        AUX_S_5_0_0_0[m * 21 + 19] = P.PA_y[i] * AUX_S_4_0_0_0[m * 15 + 14] - a_over_p * PQ_y * AUX_S_4_0_0_0[(m+1) * 15 + 14];

                        //H_005 : STEP: z
                        AUX_S_5_0_0_0[m * 21 + 20] = P.PA_z[i] * AUX_S_4_0_0_0[m * 15 + 14] - a_over_p * PQ_z * AUX_S_4_0_0_0[(m+1) * 15 + 14]
                                      +4 * one_over_2p * ( AUX_S_3_0_0_0[m * 10 +  9] - a_over_p * AUX_S_3_0_0_0[(m+1) * 10 + 9] );

                    }




                    //////////////////////////////////////////////
                    // Primitive integrals: Electron transfer
                    //////////////////////////////////////////////

                    // ( D_200 S_000 | P_100 S_000 )^0_{t} = x * ( D_200 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_300 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[0] = etfac[0] * AUX_S_2_0_0_0[0] + 2 * one_over_2q * AUX_S_1_0_0_0[0] - p_over_q * AUX_S_3_0_0_0[0];

                    // ( D_200 S_000 | P_010 S_000 )^0_{t} = y * ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[1] = etfac[1] * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_3_0_0_0[1];

                    // ( D_200 S_000 | P_001 S_000 )^0_{t} = z * ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[2] = etfac[2] * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_3_0_0_0[2];

                    // ( D_110 S_000 | P_100 S_000 )^0_{t} = x * ( D_110 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[3] = etfac[0] * AUX_S_2_0_0_0[1] + 1 * one_over_2q * AUX_S_1_0_0_0[1] - p_over_q * AUX_S_3_0_0_0[1];

                    // ( D_110 S_000 | P_010 S_000 )^0_{t} = y * ( D_110 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[4] = etfac[1] * AUX_S_2_0_0_0[1] + 1 * one_over_2q * AUX_S_1_0_0_0[0] - p_over_q * AUX_S_3_0_0_0[3];

                    // ( D_110 S_000 | P_001 S_000 )^0_{t} = z * ( D_110 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[5] = etfac[2] * AUX_S_2_0_0_0[1] - p_over_q * AUX_S_3_0_0_0[4];

                    // ( D_101 S_000 | P_100 S_000 )^0_{t} = x * ( D_101 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[6] = etfac[0] * AUX_S_2_0_0_0[2] + 1 * one_over_2q * AUX_S_1_0_0_0[2] - p_over_q * AUX_S_3_0_0_0[2];

                    // ( D_101 S_000 | P_010 S_000 )^0_{t} = y * ( D_101 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[7] = etfac[1] * AUX_S_2_0_0_0[2] - p_over_q * AUX_S_3_0_0_0[4];

                    // ( D_101 S_000 | P_001 S_000 )^0_{t} = z * ( D_101 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[8] = etfac[2] * AUX_S_2_0_0_0[2] + 1 * one_over_2q * AUX_S_1_0_0_0[0] - p_over_q * AUX_S_3_0_0_0[5];

                    // ( D_020 S_000 | P_100 S_000 )^0_{t} = x * ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[9] = etfac[0] * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_3_0_0_0[3];

                    // ( D_020 S_000 | P_010 S_000 )^0_{t} = y * ( D_020 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_030 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[10] = etfac[1] * AUX_S_2_0_0_0[3] + 2 * one_over_2q * AUX_S_1_0_0_0[1] - p_over_q * AUX_S_3_0_0_0[6];

                    // ( D_020 S_000 | P_001 S_000 )^0_{t} = z * ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[11] = etfac[2] * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_3_0_0_0[7];

                    // ( D_011 S_000 | P_100 S_000 )^0_{t} = x * ( D_011 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[12] = etfac[0] * AUX_S_2_0_0_0[4] - p_over_q * AUX_S_3_0_0_0[4];

                    // ( D_011 S_000 | P_010 S_000 )^0_{t} = y * ( D_011 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[13] = etfac[1] * AUX_S_2_0_0_0[4] + 1 * one_over_2q * AUX_S_1_0_0_0[2] - p_over_q * AUX_S_3_0_0_0[7];

                    // ( D_011 S_000 | P_001 S_000 )^0_{t} = z * ( D_011 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[14] = etfac[2] * AUX_S_2_0_0_0[4] + 1 * one_over_2q * AUX_S_1_0_0_0[1] - p_over_q * AUX_S_3_0_0_0[8];

                    // ( D_002 S_000 | P_100 S_000 )^0_{t} = x * ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[15] = etfac[0] * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_3_0_0_0[5];

                    // ( D_002 S_000 | P_010 S_000 )^0_{t} = y * ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[16] = etfac[1] * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_3_0_0_0[8];

                    // ( D_002 S_000 | P_001 S_000 )^0_{t} = z * ( D_002 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_003 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_2_0_1_0[17] = etfac[2] * AUX_S_2_0_0_0[5] + 2 * one_over_2q * AUX_S_1_0_0_0[2] - p_over_q * AUX_S_3_0_0_0[9];

                    // ( F_300 S_000 | P_100 S_000 )^0_{t} = x * ( F_300 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_400 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[0] = etfac[0] * AUX_S_3_0_0_0[0] + 3 * one_over_2q * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_4_0_0_0[0];

                    // ( F_300 S_000 | P_010 S_000 )^0_{t} = y * ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[1] = etfac[1] * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_4_0_0_0[1];

                    // ( F_300 S_000 | P_001 S_000 )^0_{t} = z * ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[2] = etfac[2] * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_4_0_0_0[2];

                    // ( F_210 S_000 | P_100 S_000 )^0_{t} = x * ( F_210 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[3] = etfac[0] * AUX_S_3_0_0_0[1] + 2 * one_over_2q * AUX_S_2_0_0_0[1] - p_over_q * AUX_S_4_0_0_0[1];

                    // ( F_210 S_000 | P_010 S_000 )^0_{t} = y * ( F_210 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[4] = etfac[1] * AUX_S_3_0_0_0[1] + 1 * one_over_2q * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_4_0_0_0[3];

                    // ( F_210 S_000 | P_001 S_000 )^0_{t} = z * ( F_210 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[5] = etfac[2] * AUX_S_3_0_0_0[1] - p_over_q * AUX_S_4_0_0_0[4];

                    // ( F_201 S_000 | P_100 S_000 )^0_{t} = x * ( F_201 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[6] = etfac[0] * AUX_S_3_0_0_0[2] + 2 * one_over_2q * AUX_S_2_0_0_0[2] - p_over_q * AUX_S_4_0_0_0[2];

                    // ( F_201 S_000 | P_010 S_000 )^0_{t} = y * ( F_201 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[7] = etfac[1] * AUX_S_3_0_0_0[2] - p_over_q * AUX_S_4_0_0_0[4];

                    // ( F_201 S_000 | P_001 S_000 )^0_{t} = z * ( F_201 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[8] = etfac[2] * AUX_S_3_0_0_0[2] + 1 * one_over_2q * AUX_S_2_0_0_0[0] - p_over_q * AUX_S_4_0_0_0[5];

                    // ( F_120 S_000 | P_100 S_000 )^0_{t} = x * ( F_120 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[9] = etfac[0] * AUX_S_3_0_0_0[3] + 1 * one_over_2q * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_4_0_0_0[3];

                    // ( F_120 S_000 | P_010 S_000 )^0_{t} = y * ( F_120 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[10] = etfac[1] * AUX_S_3_0_0_0[3] + 2 * one_over_2q * AUX_S_2_0_0_0[1] - p_over_q * AUX_S_4_0_0_0[6];

                    // ( F_120 S_000 | P_001 S_000 )^0_{t} = z * ( F_120 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[11] = etfac[2] * AUX_S_3_0_0_0[3] - p_over_q * AUX_S_4_0_0_0[7];

                    // ( F_111 S_000 | P_100 S_000 )^0_{t} = x * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[12] = etfac[0] * AUX_S_3_0_0_0[4] + 1 * one_over_2q * AUX_S_2_0_0_0[4] - p_over_q * AUX_S_4_0_0_0[4];

                    // ( F_111 S_000 | P_010 S_000 )^0_{t} = y * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[13] = etfac[1] * AUX_S_3_0_0_0[4] + 1 * one_over_2q * AUX_S_2_0_0_0[2] - p_over_q * AUX_S_4_0_0_0[7];

                    // ( F_111 S_000 | P_001 S_000 )^0_{t} = z * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[14] = etfac[2] * AUX_S_3_0_0_0[4] + 1 * one_over_2q * AUX_S_2_0_0_0[1] - p_over_q * AUX_S_4_0_0_0[8];

                    // ( F_102 S_000 | P_100 S_000 )^0_{t} = x * ( F_102 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[15] = etfac[0] * AUX_S_3_0_0_0[5] + 1 * one_over_2q * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_4_0_0_0[5];

                    // ( F_102 S_000 | P_010 S_000 )^0_{t} = y * ( F_102 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[16] = etfac[1] * AUX_S_3_0_0_0[5] - p_over_q * AUX_S_4_0_0_0[8];

                    // ( F_102 S_000 | P_001 S_000 )^0_{t} = z * ( F_102 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[17] = etfac[2] * AUX_S_3_0_0_0[5] + 2 * one_over_2q * AUX_S_2_0_0_0[2] - p_over_q * AUX_S_4_0_0_0[9];

                    // ( F_030 S_000 | P_100 S_000 )^0_{t} = x * ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[18] = etfac[0] * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_4_0_0_0[6];

                    // ( F_030 S_000 | P_010 S_000 )^0_{t} = y * ( F_030 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_040 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[19] = etfac[1] * AUX_S_3_0_0_0[6] + 3 * one_over_2q * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_4_0_0_0[10];

                    // ( F_030 S_000 | P_001 S_000 )^0_{t} = z * ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[20] = etfac[2] * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_4_0_0_0[11];

                    // ( F_021 S_000 | P_100 S_000 )^0_{t} = x * ( F_021 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[21] = etfac[0] * AUX_S_3_0_0_0[7] - p_over_q * AUX_S_4_0_0_0[7];

                    // ( F_021 S_000 | P_010 S_000 )^0_{t} = y * ( F_021 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[22] = etfac[1] * AUX_S_3_0_0_0[7] + 2 * one_over_2q * AUX_S_2_0_0_0[4] - p_over_q * AUX_S_4_0_0_0[11];

                    // ( F_021 S_000 | P_001 S_000 )^0_{t} = z * ( F_021 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[23] = etfac[2] * AUX_S_3_0_0_0[7] + 1 * one_over_2q * AUX_S_2_0_0_0[3] - p_over_q * AUX_S_4_0_0_0[12];

                    // ( F_012 S_000 | P_100 S_000 )^0_{t} = x * ( F_012 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[24] = etfac[0] * AUX_S_3_0_0_0[8] - p_over_q * AUX_S_4_0_0_0[8];

                    // ( F_012 S_000 | P_010 S_000 )^0_{t} = y * ( F_012 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[25] = etfac[1] * AUX_S_3_0_0_0[8] + 1 * one_over_2q * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_4_0_0_0[12];

                    // ( F_012 S_000 | P_001 S_000 )^0_{t} = z * ( F_012 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[26] = etfac[2] * AUX_S_3_0_0_0[8] + 2 * one_over_2q * AUX_S_2_0_0_0[4] - p_over_q * AUX_S_4_0_0_0[13];

                    // ( F_003 S_000 | P_100 S_000 )^0_{t} = x * ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[27] = etfac[0] * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_4_0_0_0[9];

                    // ( F_003 S_000 | P_010 S_000 )^0_{t} = y * ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[28] = etfac[1] * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_4_0_0_0[13];

                    // ( F_003 S_000 | P_001 S_000 )^0_{t} = z * ( F_003 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_004 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_3_0_1_0[29] = etfac[2] * AUX_S_3_0_0_0[9] + 3 * one_over_2q * AUX_S_2_0_0_0[5] - p_over_q * AUX_S_4_0_0_0[14];

                    // ( G_400 S_000 | P_100 S_000 )^0_{t} = x * ( G_400 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_500 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[0] = etfac[0] * AUX_S_4_0_0_0[0] + 4 * one_over_2q * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_5_0_0_0[0];

                    // ( G_400 S_000 | P_010 S_000 )^0_{t} = y * ( G_400 S_000 | S_000 S_000 )^0_{e} - ( H_410 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[1] = etfac[1] * AUX_S_4_0_0_0[0] - p_over_q * AUX_S_5_0_0_0[1];

                    // ( G_400 S_000 | P_001 S_000 )^0_{t} = z * ( G_400 S_000 | S_000 S_000 )^0_{e} - ( H_401 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[2] = etfac[2] * AUX_S_4_0_0_0[0] - p_over_q * AUX_S_5_0_0_0[2];

                    // ( G_310 S_000 | P_100 S_000 )^0_{t} = x * ( G_310 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_410 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[3] = etfac[0] * AUX_S_4_0_0_0[1] + 3 * one_over_2q * AUX_S_3_0_0_0[1] - p_over_q * AUX_S_5_0_0_0[1];

                    // ( G_310 S_000 | P_010 S_000 )^0_{t} = y * ( G_310 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_320 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[4] = etfac[1] * AUX_S_4_0_0_0[1] + 1 * one_over_2q * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_5_0_0_0[3];

                    // ( G_310 S_000 | P_001 S_000 )^0_{t} = z * ( G_310 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[5] = etfac[2] * AUX_S_4_0_0_0[1] - p_over_q * AUX_S_5_0_0_0[4];

                    // ( G_301 S_000 | P_100 S_000 )^0_{t} = x * ( G_301 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_401 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[6] = etfac[0] * AUX_S_4_0_0_0[2] + 3 * one_over_2q * AUX_S_3_0_0_0[2] - p_over_q * AUX_S_5_0_0_0[2];

                    // ( G_301 S_000 | P_010 S_000 )^0_{t} = y * ( G_301 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[7] = etfac[1] * AUX_S_4_0_0_0[2] - p_over_q * AUX_S_5_0_0_0[4];

                    // ( G_301 S_000 | P_001 S_000 )^0_{t} = z * ( G_301 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_302 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[8] = etfac[2] * AUX_S_4_0_0_0[2] + 1 * one_over_2q * AUX_S_3_0_0_0[0] - p_over_q * AUX_S_5_0_0_0[5];

                    // ( G_220 S_000 | P_100 S_000 )^0_{t} = x * ( G_220 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_320 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[9] = etfac[0] * AUX_S_4_0_0_0[3] + 2 * one_over_2q * AUX_S_3_0_0_0[3] - p_over_q * AUX_S_5_0_0_0[3];

                    // ( G_220 S_000 | P_010 S_000 )^0_{t} = y * ( G_220 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_230 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[10] = etfac[1] * AUX_S_4_0_0_0[3] + 2 * one_over_2q * AUX_S_3_0_0_0[1] - p_over_q * AUX_S_5_0_0_0[6];

                    // ( G_220 S_000 | P_001 S_000 )^0_{t} = z * ( G_220 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[11] = etfac[2] * AUX_S_4_0_0_0[3] - p_over_q * AUX_S_5_0_0_0[7];

                    // ( G_211 S_000 | P_100 S_000 )^0_{t} = x * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[12] = etfac[0] * AUX_S_4_0_0_0[4] + 2 * one_over_2q * AUX_S_3_0_0_0[4] - p_over_q * AUX_S_5_0_0_0[4];

                    // ( G_211 S_000 | P_010 S_000 )^0_{t} = y * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[13] = etfac[1] * AUX_S_4_0_0_0[4] + 1 * one_over_2q * AUX_S_3_0_0_0[2] - p_over_q * AUX_S_5_0_0_0[7];

                    // ( G_211 S_000 | P_001 S_000 )^0_{t} = z * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[14] = etfac[2] * AUX_S_4_0_0_0[4] + 1 * one_over_2q * AUX_S_3_0_0_0[1] - p_over_q * AUX_S_5_0_0_0[8];

                    // ( G_202 S_000 | P_100 S_000 )^0_{t} = x * ( G_202 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_302 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[15] = etfac[0] * AUX_S_4_0_0_0[5] + 2 * one_over_2q * AUX_S_3_0_0_0[5] - p_over_q * AUX_S_5_0_0_0[5];

                    // ( G_202 S_000 | P_010 S_000 )^0_{t} = y * ( G_202 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[16] = etfac[1] * AUX_S_4_0_0_0[5] - p_over_q * AUX_S_5_0_0_0[8];

                    // ( G_202 S_000 | P_001 S_000 )^0_{t} = z * ( G_202 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_203 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[17] = etfac[2] * AUX_S_4_0_0_0[5] + 2 * one_over_2q * AUX_S_3_0_0_0[2] - p_over_q * AUX_S_5_0_0_0[9];

                    // ( G_130 S_000 | P_100 S_000 )^0_{t} = x * ( G_130 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_230 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[18] = etfac[0] * AUX_S_4_0_0_0[6] + 1 * one_over_2q * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_5_0_0_0[6];

                    // ( G_130 S_000 | P_010 S_000 )^0_{t} = y * ( G_130 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_140 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[19] = etfac[1] * AUX_S_4_0_0_0[6] + 3 * one_over_2q * AUX_S_3_0_0_0[3] - p_over_q * AUX_S_5_0_0_0[10];

                    // ( G_130 S_000 | P_001 S_000 )^0_{t} = z * ( G_130 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[20] = etfac[2] * AUX_S_4_0_0_0[6] - p_over_q * AUX_S_5_0_0_0[11];

                    // ( G_121 S_000 | P_100 S_000 )^0_{t} = x * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[21] = etfac[0] * AUX_S_4_0_0_0[7] + 1 * one_over_2q * AUX_S_3_0_0_0[7] - p_over_q * AUX_S_5_0_0_0[7];

                    // ( G_121 S_000 | P_010 S_000 )^0_{t} = y * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[22] = etfac[1] * AUX_S_4_0_0_0[7] + 2 * one_over_2q * AUX_S_3_0_0_0[4] - p_over_q * AUX_S_5_0_0_0[11];

                    // ( G_121 S_000 | P_001 S_000 )^0_{t} = z * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[23] = etfac[2] * AUX_S_4_0_0_0[7] + 1 * one_over_2q * AUX_S_3_0_0_0[3] - p_over_q * AUX_S_5_0_0_0[12];

                    // ( G_112 S_000 | P_100 S_000 )^0_{t} = x * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[24] = etfac[0] * AUX_S_4_0_0_0[8] + 1 * one_over_2q * AUX_S_3_0_0_0[8] - p_over_q * AUX_S_5_0_0_0[8];

                    // ( G_112 S_000 | P_010 S_000 )^0_{t} = y * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[25] = etfac[1] * AUX_S_4_0_0_0[8] + 1 * one_over_2q * AUX_S_3_0_0_0[5] - p_over_q * AUX_S_5_0_0_0[12];

                    // ( G_112 S_000 | P_001 S_000 )^0_{t} = z * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[26] = etfac[2] * AUX_S_4_0_0_0[8] + 2 * one_over_2q * AUX_S_3_0_0_0[4] - p_over_q * AUX_S_5_0_0_0[13];

                    // ( G_103 S_000 | P_100 S_000 )^0_{t} = x * ( G_103 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_203 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[27] = etfac[0] * AUX_S_4_0_0_0[9] + 1 * one_over_2q * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_5_0_0_0[9];

                    // ( G_103 S_000 | P_010 S_000 )^0_{t} = y * ( G_103 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[28] = etfac[1] * AUX_S_4_0_0_0[9] - p_over_q * AUX_S_5_0_0_0[13];

                    // ( G_103 S_000 | P_001 S_000 )^0_{t} = z * ( G_103 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_104 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[29] = etfac[2] * AUX_S_4_0_0_0[9] + 3 * one_over_2q * AUX_S_3_0_0_0[5] - p_over_q * AUX_S_5_0_0_0[14];

                    // ( G_040 S_000 | P_100 S_000 )^0_{t} = x * ( G_040 S_000 | S_000 S_000 )^0_{e} - ( H_140 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[30] = etfac[0] * AUX_S_4_0_0_0[10] - p_over_q * AUX_S_5_0_0_0[10];

                    // ( G_040 S_000 | P_010 S_000 )^0_{t} = y * ( G_040 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_050 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[31] = etfac[1] * AUX_S_4_0_0_0[10] + 4 * one_over_2q * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_5_0_0_0[15];

                    // ( G_040 S_000 | P_001 S_000 )^0_{t} = z * ( G_040 S_000 | S_000 S_000 )^0_{e} - ( H_041 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[32] = etfac[2] * AUX_S_4_0_0_0[10] - p_over_q * AUX_S_5_0_0_0[16];

                    // ( G_031 S_000 | P_100 S_000 )^0_{t} = x * ( G_031 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[33] = etfac[0] * AUX_S_4_0_0_0[11] - p_over_q * AUX_S_5_0_0_0[11];

                    // ( G_031 S_000 | P_010 S_000 )^0_{t} = y * ( G_031 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_041 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[34] = etfac[1] * AUX_S_4_0_0_0[11] + 3 * one_over_2q * AUX_S_3_0_0_0[7] - p_over_q * AUX_S_5_0_0_0[16];

                    // ( G_031 S_000 | P_001 S_000 )^0_{t} = z * ( G_031 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_032 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[35] = etfac[2] * AUX_S_4_0_0_0[11] + 1 * one_over_2q * AUX_S_3_0_0_0[6] - p_over_q * AUX_S_5_0_0_0[17];

                    // ( G_022 S_000 | P_100 S_000 )^0_{t} = x * ( G_022 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[36] = etfac[0] * AUX_S_4_0_0_0[12] - p_over_q * AUX_S_5_0_0_0[12];

                    // ( G_022 S_000 | P_010 S_000 )^0_{t} = y * ( G_022 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_032 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[37] = etfac[1] * AUX_S_4_0_0_0[12] + 2 * one_over_2q * AUX_S_3_0_0_0[8] - p_over_q * AUX_S_5_0_0_0[17];

                    // ( G_022 S_000 | P_001 S_000 )^0_{t} = z * ( G_022 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_023 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[38] = etfac[2] * AUX_S_4_0_0_0[12] + 2 * one_over_2q * AUX_S_3_0_0_0[7] - p_over_q * AUX_S_5_0_0_0[18];

                    // ( G_013 S_000 | P_100 S_000 )^0_{t} = x * ( G_013 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[39] = etfac[0] * AUX_S_4_0_0_0[13] - p_over_q * AUX_S_5_0_0_0[13];

                    // ( G_013 S_000 | P_010 S_000 )^0_{t} = y * ( G_013 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_023 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[40] = etfac[1] * AUX_S_4_0_0_0[13] + 1 * one_over_2q * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_5_0_0_0[18];

                    // ( G_013 S_000 | P_001 S_000 )^0_{t} = z * ( G_013 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_014 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[41] = etfac[2] * AUX_S_4_0_0_0[13] + 3 * one_over_2q * AUX_S_3_0_0_0[8] - p_over_q * AUX_S_5_0_0_0[19];

                    // ( G_004 S_000 | P_100 S_000 )^0_{t} = x * ( G_004 S_000 | S_000 S_000 )^0_{e} - ( H_104 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[42] = etfac[0] * AUX_S_4_0_0_0[14] - p_over_q * AUX_S_5_0_0_0[14];

                    // ( G_004 S_000 | P_010 S_000 )^0_{t} = y * ( G_004 S_000 | S_000 S_000 )^0_{e} - ( H_014 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[43] = etfac[1] * AUX_S_4_0_0_0[14] - p_over_q * AUX_S_5_0_0_0[19];

                    // ( G_004 S_000 | P_001 S_000 )^0_{t} = z * ( G_004 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_005 S_000 | S_000 S_000 )^0_{e}
                    AUX_S_4_0_1_0[44] = etfac[2] * AUX_S_4_0_0_0[14] + 4 * one_over_2q * AUX_S_3_0_0_0[9] - p_over_q * AUX_S_5_0_0_0[20];


                    // Accumulating in contracted workspace
                    S_2_0_1_0[abcd * 18 + 0] += AUX_S_2_0_1_0[0];
                    S_2_0_1_0[abcd * 18 + 1] += AUX_S_2_0_1_0[1];
                    S_2_0_1_0[abcd * 18 + 2] += AUX_S_2_0_1_0[2];
                    S_2_0_1_0[abcd * 18 + 3] += AUX_S_2_0_1_0[3];
                    S_2_0_1_0[abcd * 18 + 4] += AUX_S_2_0_1_0[4];
                    S_2_0_1_0[abcd * 18 + 5] += AUX_S_2_0_1_0[5];
                    S_2_0_1_0[abcd * 18 + 6] += AUX_S_2_0_1_0[6];
                    S_2_0_1_0[abcd * 18 + 7] += AUX_S_2_0_1_0[7];
                    S_2_0_1_0[abcd * 18 + 8] += AUX_S_2_0_1_0[8];
                    S_2_0_1_0[abcd * 18 + 9] += AUX_S_2_0_1_0[9];
                    S_2_0_1_0[abcd * 18 + 10] += AUX_S_2_0_1_0[10];
                    S_2_0_1_0[abcd * 18 + 11] += AUX_S_2_0_1_0[11];
                    S_2_0_1_0[abcd * 18 + 12] += AUX_S_2_0_1_0[12];
                    S_2_0_1_0[abcd * 18 + 13] += AUX_S_2_0_1_0[13];
                    S_2_0_1_0[abcd * 18 + 14] += AUX_S_2_0_1_0[14];
                    S_2_0_1_0[abcd * 18 + 15] += AUX_S_2_0_1_0[15];
                    S_2_0_1_0[abcd * 18 + 16] += AUX_S_2_0_1_0[16];
                    S_2_0_1_0[abcd * 18 + 17] += AUX_S_2_0_1_0[17];

                    // Accumulating in contracted workspace
                    S_3_0_1_0[abcd * 30 + 0] += AUX_S_3_0_1_0[0];
                    S_3_0_1_0[abcd * 30 + 1] += AUX_S_3_0_1_0[1];
                    S_3_0_1_0[abcd * 30 + 2] += AUX_S_3_0_1_0[2];
                    S_3_0_1_0[abcd * 30 + 3] += AUX_S_3_0_1_0[3];
                    S_3_0_1_0[abcd * 30 + 4] += AUX_S_3_0_1_0[4];
                    S_3_0_1_0[abcd * 30 + 5] += AUX_S_3_0_1_0[5];
                    S_3_0_1_0[abcd * 30 + 6] += AUX_S_3_0_1_0[6];
                    S_3_0_1_0[abcd * 30 + 7] += AUX_S_3_0_1_0[7];
                    S_3_0_1_0[abcd * 30 + 8] += AUX_S_3_0_1_0[8];
                    S_3_0_1_0[abcd * 30 + 9] += AUX_S_3_0_1_0[9];
                    S_3_0_1_0[abcd * 30 + 10] += AUX_S_3_0_1_0[10];
                    S_3_0_1_0[abcd * 30 + 11] += AUX_S_3_0_1_0[11];
                    S_3_0_1_0[abcd * 30 + 12] += AUX_S_3_0_1_0[12];
                    S_3_0_1_0[abcd * 30 + 13] += AUX_S_3_0_1_0[13];
                    S_3_0_1_0[abcd * 30 + 14] += AUX_S_3_0_1_0[14];
                    S_3_0_1_0[abcd * 30 + 15] += AUX_S_3_0_1_0[15];
                    S_3_0_1_0[abcd * 30 + 16] += AUX_S_3_0_1_0[16];
                    S_3_0_1_0[abcd * 30 + 17] += AUX_S_3_0_1_0[17];
                    S_3_0_1_0[abcd * 30 + 18] += AUX_S_3_0_1_0[18];
                    S_3_0_1_0[abcd * 30 + 19] += AUX_S_3_0_1_0[19];
                    S_3_0_1_0[abcd * 30 + 20] += AUX_S_3_0_1_0[20];
                    S_3_0_1_0[abcd * 30 + 21] += AUX_S_3_0_1_0[21];
                    S_3_0_1_0[abcd * 30 + 22] += AUX_S_3_0_1_0[22];
                    S_3_0_1_0[abcd * 30 + 23] += AUX_S_3_0_1_0[23];
                    S_3_0_1_0[abcd * 30 + 24] += AUX_S_3_0_1_0[24];
                    S_3_0_1_0[abcd * 30 + 25] += AUX_S_3_0_1_0[25];
                    S_3_0_1_0[abcd * 30 + 26] += AUX_S_3_0_1_0[26];
                    S_3_0_1_0[abcd * 30 + 27] += AUX_S_3_0_1_0[27];
                    S_3_0_1_0[abcd * 30 + 28] += AUX_S_3_0_1_0[28];
                    S_3_0_1_0[abcd * 30 + 29] += AUX_S_3_0_1_0[29];

                    // Accumulating in contracted workspace
                    S_4_0_1_0[abcd * 45 + 0] += AUX_S_4_0_1_0[0];
                    S_4_0_1_0[abcd * 45 + 1] += AUX_S_4_0_1_0[1];
                    S_4_0_1_0[abcd * 45 + 2] += AUX_S_4_0_1_0[2];
                    S_4_0_1_0[abcd * 45 + 3] += AUX_S_4_0_1_0[3];
                    S_4_0_1_0[abcd * 45 + 4] += AUX_S_4_0_1_0[4];
                    S_4_0_1_0[abcd * 45 + 5] += AUX_S_4_0_1_0[5];
                    S_4_0_1_0[abcd * 45 + 6] += AUX_S_4_0_1_0[6];
                    S_4_0_1_0[abcd * 45 + 7] += AUX_S_4_0_1_0[7];
                    S_4_0_1_0[abcd * 45 + 8] += AUX_S_4_0_1_0[8];
                    S_4_0_1_0[abcd * 45 + 9] += AUX_S_4_0_1_0[9];
                    S_4_0_1_0[abcd * 45 + 10] += AUX_S_4_0_1_0[10];
                    S_4_0_1_0[abcd * 45 + 11] += AUX_S_4_0_1_0[11];
                    S_4_0_1_0[abcd * 45 + 12] += AUX_S_4_0_1_0[12];
                    S_4_0_1_0[abcd * 45 + 13] += AUX_S_4_0_1_0[13];
                    S_4_0_1_0[abcd * 45 + 14] += AUX_S_4_0_1_0[14];
                    S_4_0_1_0[abcd * 45 + 15] += AUX_S_4_0_1_0[15];
                    S_4_0_1_0[abcd * 45 + 16] += AUX_S_4_0_1_0[16];
                    S_4_0_1_0[abcd * 45 + 17] += AUX_S_4_0_1_0[17];
                    S_4_0_1_0[abcd * 45 + 18] += AUX_S_4_0_1_0[18];
                    S_4_0_1_0[abcd * 45 + 19] += AUX_S_4_0_1_0[19];
                    S_4_0_1_0[abcd * 45 + 20] += AUX_S_4_0_1_0[20];
                    S_4_0_1_0[abcd * 45 + 21] += AUX_S_4_0_1_0[21];
                    S_4_0_1_0[abcd * 45 + 22] += AUX_S_4_0_1_0[22];
                    S_4_0_1_0[abcd * 45 + 23] += AUX_S_4_0_1_0[23];
                    S_4_0_1_0[abcd * 45 + 24] += AUX_S_4_0_1_0[24];
                    S_4_0_1_0[abcd * 45 + 25] += AUX_S_4_0_1_0[25];
                    S_4_0_1_0[abcd * 45 + 26] += AUX_S_4_0_1_0[26];
                    S_4_0_1_0[abcd * 45 + 27] += AUX_S_4_0_1_0[27];
                    S_4_0_1_0[abcd * 45 + 28] += AUX_S_4_0_1_0[28];
                    S_4_0_1_0[abcd * 45 + 29] += AUX_S_4_0_1_0[29];
                    S_4_0_1_0[abcd * 45 + 30] += AUX_S_4_0_1_0[30];
                    S_4_0_1_0[abcd * 45 + 31] += AUX_S_4_0_1_0[31];
                    S_4_0_1_0[abcd * 45 + 32] += AUX_S_4_0_1_0[32];
                    S_4_0_1_0[abcd * 45 + 33] += AUX_S_4_0_1_0[33];
                    S_4_0_1_0[abcd * 45 + 34] += AUX_S_4_0_1_0[34];
                    S_4_0_1_0[abcd * 45 + 35] += AUX_S_4_0_1_0[35];
                    S_4_0_1_0[abcd * 45 + 36] += AUX_S_4_0_1_0[36];
                    S_4_0_1_0[abcd * 45 + 37] += AUX_S_4_0_1_0[37];
                    S_4_0_1_0[abcd * 45 + 38] += AUX_S_4_0_1_0[38];
                    S_4_0_1_0[abcd * 45 + 39] += AUX_S_4_0_1_0[39];
                    S_4_0_1_0[abcd * 45 + 40] += AUX_S_4_0_1_0[40];
                    S_4_0_1_0[abcd * 45 + 41] += AUX_S_4_0_1_0[41];
                    S_4_0_1_0[abcd * 45 + 42] += AUX_S_4_0_1_0[42];
                    S_4_0_1_0[abcd * 45 + 43] += AUX_S_4_0_1_0[43];
                    S_4_0_1_0[abcd * 45 + 44] += AUX_S_4_0_1_0[44];

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
        // form S_2_2_1_0
        for(int iket = 0; iket < 3; ++iket)
        {
            // (D_200 P_100| = (F_300 S_000|_{t} + x_ab * (D_200 S_000|_{t}
            const double Q_2_0_0_1_0_0_1 = S_3_0_1_0[abcd * 30 + 0 * 3 + iket] + ( AB_x[abcd] * S_2_0_1_0[abcd * 18 + 0 * 3 + iket] );

            // (D_200 P_010| = (F_210 S_000|_{t} + y_ab * (D_200 S_000|_{t}
            const double Q_2_0_0_0_1_0_1 = S_3_0_1_0[abcd * 30 + 1 * 3 + iket] + ( AB_y[abcd] * S_2_0_1_0[abcd * 18 + 0 * 3 + iket] );

            // (D_200 P_001| = (F_201 S_000|_{t} + z_ab * (D_200 S_000|_{t}
            const double Q_2_0_0_0_0_1_1 = S_3_0_1_0[abcd * 30 + 2 * 3 + iket] + ( AB_z[abcd] * S_2_0_1_0[abcd * 18 + 0 * 3 + iket] );

            // (D_110 P_100| = (F_210 S_000|_{t} + x_ab * (D_110 S_000|_{t}
            const double Q_1_1_0_1_0_0_1 = S_3_0_1_0[abcd * 30 + 1 * 3 + iket] + ( AB_x[abcd] * S_2_0_1_0[abcd * 18 + 1 * 3 + iket] );

            // (D_110 P_010| = (F_120 S_000|_{t} + y_ab * (D_110 S_000|_{t}
            const double Q_1_1_0_0_1_0_1 = S_3_0_1_0[abcd * 30 + 3 * 3 + iket] + ( AB_y[abcd] * S_2_0_1_0[abcd * 18 + 1 * 3 + iket] );

            // (D_110 P_001| = (F_111 S_000|_{t} + z_ab * (D_110 S_000|_{t}
            const double Q_1_1_0_0_0_1_1 = S_3_0_1_0[abcd * 30 + 4 * 3 + iket] + ( AB_z[abcd] * S_2_0_1_0[abcd * 18 + 1 * 3 + iket] );

            // (D_101 P_100| = (F_201 S_000|_{t} + x_ab * (D_101 S_000|_{t}
            const double Q_1_0_1_1_0_0_1 = S_3_0_1_0[abcd * 30 + 2 * 3 + iket] + ( AB_x[abcd] * S_2_0_1_0[abcd * 18 + 2 * 3 + iket] );

            // (D_101 P_010| = (F_111 S_000|_{t} + y_ab * (D_101 S_000|_{t}
            const double Q_1_0_1_0_1_0_1 = S_3_0_1_0[abcd * 30 + 4 * 3 + iket] + ( AB_y[abcd] * S_2_0_1_0[abcd * 18 + 2 * 3 + iket] );

            // (D_101 P_001| = (F_102 S_000|_{t} + z_ab * (D_101 S_000|_{t}
            const double Q_1_0_1_0_0_1_1 = S_3_0_1_0[abcd * 30 + 5 * 3 + iket] + ( AB_z[abcd] * S_2_0_1_0[abcd * 18 + 2 * 3 + iket] );

            // (D_020 P_100| = (F_120 S_000|_{t} + x_ab * (D_020 S_000|_{t}
            const double Q_0_2_0_1_0_0_1 = S_3_0_1_0[abcd * 30 + 3 * 3 + iket] + ( AB_x[abcd] * S_2_0_1_0[abcd * 18 + 3 * 3 + iket] );

            // (D_020 P_010| = (F_030 S_000|_{t} + y_ab * (D_020 S_000|_{t}
            const double Q_0_2_0_0_1_0_1 = S_3_0_1_0[abcd * 30 + 6 * 3 + iket] + ( AB_y[abcd] * S_2_0_1_0[abcd * 18 + 3 * 3 + iket] );

            // (D_020 P_001| = (F_021 S_000|_{t} + z_ab * (D_020 S_000|_{t}
            const double Q_0_2_0_0_0_1_1 = S_3_0_1_0[abcd * 30 + 7 * 3 + iket] + ( AB_z[abcd] * S_2_0_1_0[abcd * 18 + 3 * 3 + iket] );

            // (D_011 P_100| = (F_111 S_000|_{t} + x_ab * (D_011 S_000|_{t}
            const double Q_0_1_1_1_0_0_1 = S_3_0_1_0[abcd * 30 + 4 * 3 + iket] + ( AB_x[abcd] * S_2_0_1_0[abcd * 18 + 4 * 3 + iket] );

            // (D_011 P_010| = (F_021 S_000|_{t} + y_ab * (D_011 S_000|_{t}
            const double Q_0_1_1_0_1_0_1 = S_3_0_1_0[abcd * 30 + 7 * 3 + iket] + ( AB_y[abcd] * S_2_0_1_0[abcd * 18 + 4 * 3 + iket] );

            // (D_011 P_001| = (F_012 S_000|_{t} + z_ab * (D_011 S_000|_{t}
            const double Q_0_1_1_0_0_1_1 = S_3_0_1_0[abcd * 30 + 8 * 3 + iket] + ( AB_z[abcd] * S_2_0_1_0[abcd * 18 + 4 * 3 + iket] );

            // (D_002 P_100| = (F_102 S_000|_{t} + x_ab * (D_002 S_000|_{t}
            const double Q_0_0_2_1_0_0_1 = S_3_0_1_0[abcd * 30 + 5 * 3 + iket] + ( AB_x[abcd] * S_2_0_1_0[abcd * 18 + 5 * 3 + iket] );

            // (D_002 P_010| = (F_012 S_000|_{t} + y_ab * (D_002 S_000|_{t}
            const double Q_0_0_2_0_1_0_1 = S_3_0_1_0[abcd * 30 + 8 * 3 + iket] + ( AB_y[abcd] * S_2_0_1_0[abcd * 18 + 5 * 3 + iket] );

            // (D_002 P_001| = (F_003 S_000|_{t} + z_ab * (D_002 S_000|_{t}
            const double Q_0_0_2_0_0_1_1 = S_3_0_1_0[abcd * 30 + 9 * 3 + iket] + ( AB_z[abcd] * S_2_0_1_0[abcd * 18 + 5 * 3 + iket] );

            // (F_300 P_100| = (G_400 S_000|_{t} + x_ab * (F_300 S_000|_{t}
            const double Q_3_0_0_1_0_0_1 = S_4_0_1_0[abcd * 45 + 0 * 3 + iket] + ( AB_x[abcd] * S_3_0_1_0[abcd * 30 + 0 * 3 + iket] );

            // (F_210 P_100| = (G_310 S_000|_{t} + x_ab * (F_210 S_000|_{t}
            const double Q_2_1_0_1_0_0_1 = S_4_0_1_0[abcd * 45 + 1 * 3 + iket] + ( AB_x[abcd] * S_3_0_1_0[abcd * 30 + 1 * 3 + iket] );

            // (F_210 P_010| = (G_220 S_000|_{t} + y_ab * (F_210 S_000|_{t}
            const double Q_2_1_0_0_1_0_1 = S_4_0_1_0[abcd * 45 + 3 * 3 + iket] + ( AB_y[abcd] * S_3_0_1_0[abcd * 30 + 1 * 3 + iket] );

            // (F_201 P_100| = (G_301 S_000|_{t} + x_ab * (F_201 S_000|_{t}
            const double Q_2_0_1_1_0_0_1 = S_4_0_1_0[abcd * 45 + 2 * 3 + iket] + ( AB_x[abcd] * S_3_0_1_0[abcd * 30 + 2 * 3 + iket] );

            // (F_201 P_010| = (G_211 S_000|_{t} + y_ab * (F_201 S_000|_{t}
            const double Q_2_0_1_0_1_0_1 = S_4_0_1_0[abcd * 45 + 4 * 3 + iket] + ( AB_y[abcd] * S_3_0_1_0[abcd * 30 + 2 * 3 + iket] );

            // (F_201 P_001| = (G_202 S_000|_{t} + z_ab * (F_201 S_000|_{t}
            const double Q_2_0_1_0_0_1_1 = S_4_0_1_0[abcd * 45 + 5 * 3 + iket] + ( AB_z[abcd] * S_3_0_1_0[abcd * 30 + 2 * 3 + iket] );

            // (F_120 P_100| = (G_220 S_000|_{t} + x_ab * (F_120 S_000|_{t}
            const double Q_1_2_0_1_0_0_1 = S_4_0_1_0[abcd * 45 + 3 * 3 + iket] + ( AB_x[abcd] * S_3_0_1_0[abcd * 30 + 3 * 3 + iket] );

            // (F_120 P_010| = (G_130 S_000|_{t} + y_ab * (F_120 S_000|_{t}
            const double Q_1_2_0_0_1_0_1 = S_4_0_1_0[abcd * 45 + 6 * 3 + iket] + ( AB_y[abcd] * S_3_0_1_0[abcd * 30 + 3 * 3 + iket] );

            // (F_111 P_100| = (G_211 S_000|_{t} + x_ab * (F_111 S_000|_{t}
            const double Q_1_1_1_1_0_0_1 = S_4_0_1_0[abcd * 45 + 4 * 3 + iket] + ( AB_x[abcd] * S_3_0_1_0[abcd * 30 + 4 * 3 + iket] );

            // (F_111 P_010| = (G_121 S_000|_{t} + y_ab * (F_111 S_000|_{t}
            const double Q_1_1_1_0_1_0_1 = S_4_0_1_0[abcd * 45 + 7 * 3 + iket] + ( AB_y[abcd] * S_3_0_1_0[abcd * 30 + 4 * 3 + iket] );

            // (F_111 P_001| = (G_112 S_000|_{t} + z_ab * (F_111 S_000|_{t}
            const double Q_1_1_1_0_0_1_1 = S_4_0_1_0[abcd * 45 + 8 * 3 + iket] + ( AB_z[abcd] * S_3_0_1_0[abcd * 30 + 4 * 3 + iket] );

            // (F_102 P_100| = (G_202 S_000|_{t} + x_ab * (F_102 S_000|_{t}
            const double Q_1_0_2_1_0_0_1 = S_4_0_1_0[abcd * 45 + 5 * 3 + iket] + ( AB_x[abcd] * S_3_0_1_0[abcd * 30 + 5 * 3 + iket] );

            // (F_102 P_010| = (G_112 S_000|_{t} + y_ab * (F_102 S_000|_{t}
            const double Q_1_0_2_0_1_0_1 = S_4_0_1_0[abcd * 45 + 8 * 3 + iket] + ( AB_y[abcd] * S_3_0_1_0[abcd * 30 + 5 * 3 + iket] );

            // (F_102 P_001| = (G_103 S_000|_{t} + z_ab * (F_102 S_000|_{t}
            const double Q_1_0_2_0_0_1_1 = S_4_0_1_0[abcd * 45 + 9 * 3 + iket] + ( AB_z[abcd] * S_3_0_1_0[abcd * 30 + 5 * 3 + iket] );

            // (F_030 P_100| = (G_130 S_000|_{t} + x_ab * (F_030 S_000|_{t}
            const double Q_0_3_0_1_0_0_1 = S_4_0_1_0[abcd * 45 + 6 * 3 + iket] + ( AB_x[abcd] * S_3_0_1_0[abcd * 30 + 6 * 3 + iket] );

            // (F_030 P_010| = (G_040 S_000|_{t} + y_ab * (F_030 S_000|_{t}
            const double Q_0_3_0_0_1_0_1 = S_4_0_1_0[abcd * 45 + 10 * 3 + iket] + ( AB_y[abcd] * S_3_0_1_0[abcd * 30 + 6 * 3 + iket] );

            // (F_021 P_100| = (G_121 S_000|_{t} + x_ab * (F_021 S_000|_{t}
            const double Q_0_2_1_1_0_0_1 = S_4_0_1_0[abcd * 45 + 7 * 3 + iket] + ( AB_x[abcd] * S_3_0_1_0[abcd * 30 + 7 * 3 + iket] );

            // (F_021 P_010| = (G_031 S_000|_{t} + y_ab * (F_021 S_000|_{t}
            const double Q_0_2_1_0_1_0_1 = S_4_0_1_0[abcd * 45 + 11 * 3 + iket] + ( AB_y[abcd] * S_3_0_1_0[abcd * 30 + 7 * 3 + iket] );

            // (F_021 P_001| = (G_022 S_000|_{t} + z_ab * (F_021 S_000|_{t}
            const double Q_0_2_1_0_0_1_1 = S_4_0_1_0[abcd * 45 + 12 * 3 + iket] + ( AB_z[abcd] * S_3_0_1_0[abcd * 30 + 7 * 3 + iket] );

            // (F_012 P_100| = (G_112 S_000|_{t} + x_ab * (F_012 S_000|_{t}
            const double Q_0_1_2_1_0_0_1 = S_4_0_1_0[abcd * 45 + 8 * 3 + iket] + ( AB_x[abcd] * S_3_0_1_0[abcd * 30 + 8 * 3 + iket] );

            // (F_012 P_010| = (G_022 S_000|_{t} + y_ab * (F_012 S_000|_{t}
            const double Q_0_1_2_0_1_0_1 = S_4_0_1_0[abcd * 45 + 12 * 3 + iket] + ( AB_y[abcd] * S_3_0_1_0[abcd * 30 + 8 * 3 + iket] );

            // (F_012 P_001| = (G_013 S_000|_{t} + z_ab * (F_012 S_000|_{t}
            const double Q_0_1_2_0_0_1_1 = S_4_0_1_0[abcd * 45 + 13 * 3 + iket] + ( AB_z[abcd] * S_3_0_1_0[abcd * 30 + 8 * 3 + iket] );

            // (F_003 P_100| = (G_103 S_000|_{t} + x_ab * (F_003 S_000|_{t}
            const double Q_0_0_3_1_0_0_1 = S_4_0_1_0[abcd * 45 + 9 * 3 + iket] + ( AB_x[abcd] * S_3_0_1_0[abcd * 30 + 9 * 3 + iket] );

            // (F_003 P_010| = (G_013 S_000|_{t} + y_ab * (F_003 S_000|_{t}
            const double Q_0_0_3_0_1_0_1 = S_4_0_1_0[abcd * 45 + 13 * 3 + iket] + ( AB_y[abcd] * S_3_0_1_0[abcd * 30 + 9 * 3 + iket] );

            // (F_003 P_001| = (G_004 S_000|_{t} + z_ab * (F_003 S_000|_{t}
            const double Q_0_0_3_0_0_1_1 = S_4_0_1_0[abcd * 45 + 14 * 3 + iket] + ( AB_z[abcd] * S_3_0_1_0[abcd * 30 + 9 * 3 + iket] );

            // (D_200 D_200|_{i} = (F_300 P_100| + x_ab * (D_200 P_100|
            S_2_2_1_0[abcd * 108 + 0 * 3 + iket] = Q_3_0_0_1_0_0_1 + ( AB_x[abcd] * Q_2_0_0_1_0_0_1 );

            // (D_200 D_110|_{i} = (F_210 P_100| + y_ab * (D_200 P_100|
            S_2_2_1_0[abcd * 108 + 1 * 3 + iket] = Q_2_1_0_1_0_0_1 + ( AB_y[abcd] * Q_2_0_0_1_0_0_1 );

            // (D_200 D_101|_{i} = (F_201 P_100| + z_ab * (D_200 P_100|
            S_2_2_1_0[abcd * 108 + 2 * 3 + iket] = Q_2_0_1_1_0_0_1 + ( AB_z[abcd] * Q_2_0_0_1_0_0_1 );

            // (D_200 D_020|_{i} = (F_210 P_010| + y_ab * (D_200 P_010|
            S_2_2_1_0[abcd * 108 + 3 * 3 + iket] = Q_2_1_0_0_1_0_1 + ( AB_y[abcd] * Q_2_0_0_0_1_0_1 );

            // (D_200 D_011|_{i} = (F_201 P_010| + z_ab * (D_200 P_010|
            S_2_2_1_0[abcd * 108 + 4 * 3 + iket] = Q_2_0_1_0_1_0_1 + ( AB_z[abcd] * Q_2_0_0_0_1_0_1 );

            // (D_200 D_002|_{i} = (F_201 P_001| + z_ab * (D_200 P_001|
            S_2_2_1_0[abcd * 108 + 5 * 3 + iket] = Q_2_0_1_0_0_1_1 + ( AB_z[abcd] * Q_2_0_0_0_0_1_1 );

            // (D_110 D_200|_{i} = (F_210 P_100| + x_ab * (D_110 P_100|
            S_2_2_1_0[abcd * 108 + 6 * 3 + iket] = Q_2_1_0_1_0_0_1 + ( AB_x[abcd] * Q_1_1_0_1_0_0_1 );

            // (D_110 D_110|_{i} = (F_120 P_100| + y_ab * (D_110 P_100|
            S_2_2_1_0[abcd * 108 + 7 * 3 + iket] = Q_1_2_0_1_0_0_1 + ( AB_y[abcd] * Q_1_1_0_1_0_0_1 );

            // (D_110 D_101|_{i} = (F_111 P_100| + z_ab * (D_110 P_100|
            S_2_2_1_0[abcd * 108 + 8 * 3 + iket] = Q_1_1_1_1_0_0_1 + ( AB_z[abcd] * Q_1_1_0_1_0_0_1 );

            // (D_110 D_020|_{i} = (F_120 P_010| + y_ab * (D_110 P_010|
            S_2_2_1_0[abcd * 108 + 9 * 3 + iket] = Q_1_2_0_0_1_0_1 + ( AB_y[abcd] * Q_1_1_0_0_1_0_1 );

            // (D_110 D_011|_{i} = (F_111 P_010| + z_ab * (D_110 P_010|
            S_2_2_1_0[abcd * 108 + 10 * 3 + iket] = Q_1_1_1_0_1_0_1 + ( AB_z[abcd] * Q_1_1_0_0_1_0_1 );

            // (D_110 D_002|_{i} = (F_111 P_001| + z_ab * (D_110 P_001|
            S_2_2_1_0[abcd * 108 + 11 * 3 + iket] = Q_1_1_1_0_0_1_1 + ( AB_z[abcd] * Q_1_1_0_0_0_1_1 );

            // (D_101 D_200|_{i} = (F_201 P_100| + x_ab * (D_101 P_100|
            S_2_2_1_0[abcd * 108 + 12 * 3 + iket] = Q_2_0_1_1_0_0_1 + ( AB_x[abcd] * Q_1_0_1_1_0_0_1 );

            // (D_101 D_110|_{i} = (F_111 P_100| + y_ab * (D_101 P_100|
            S_2_2_1_0[abcd * 108 + 13 * 3 + iket] = Q_1_1_1_1_0_0_1 + ( AB_y[abcd] * Q_1_0_1_1_0_0_1 );

            // (D_101 D_101|_{i} = (F_102 P_100| + z_ab * (D_101 P_100|
            S_2_2_1_0[abcd * 108 + 14 * 3 + iket] = Q_1_0_2_1_0_0_1 + ( AB_z[abcd] * Q_1_0_1_1_0_0_1 );

            // (D_101 D_020|_{i} = (F_111 P_010| + y_ab * (D_101 P_010|
            S_2_2_1_0[abcd * 108 + 15 * 3 + iket] = Q_1_1_1_0_1_0_1 + ( AB_y[abcd] * Q_1_0_1_0_1_0_1 );

            // (D_101 D_011|_{i} = (F_102 P_010| + z_ab * (D_101 P_010|
            S_2_2_1_0[abcd * 108 + 16 * 3 + iket] = Q_1_0_2_0_1_0_1 + ( AB_z[abcd] * Q_1_0_1_0_1_0_1 );

            // (D_101 D_002|_{i} = (F_102 P_001| + z_ab * (D_101 P_001|
            S_2_2_1_0[abcd * 108 + 17 * 3 + iket] = Q_1_0_2_0_0_1_1 + ( AB_z[abcd] * Q_1_0_1_0_0_1_1 );

            // (D_020 D_200|_{i} = (F_120 P_100| + x_ab * (D_020 P_100|
            S_2_2_1_0[abcd * 108 + 18 * 3 + iket] = Q_1_2_0_1_0_0_1 + ( AB_x[abcd] * Q_0_2_0_1_0_0_1 );

            // (D_020 D_110|_{i} = (F_030 P_100| + y_ab * (D_020 P_100|
            S_2_2_1_0[abcd * 108 + 19 * 3 + iket] = Q_0_3_0_1_0_0_1 + ( AB_y[abcd] * Q_0_2_0_1_0_0_1 );

            // (D_020 D_101|_{i} = (F_021 P_100| + z_ab * (D_020 P_100|
            S_2_2_1_0[abcd * 108 + 20 * 3 + iket] = Q_0_2_1_1_0_0_1 + ( AB_z[abcd] * Q_0_2_0_1_0_0_1 );

            // (D_020 D_020|_{i} = (F_030 P_010| + y_ab * (D_020 P_010|
            S_2_2_1_0[abcd * 108 + 21 * 3 + iket] = Q_0_3_0_0_1_0_1 + ( AB_y[abcd] * Q_0_2_0_0_1_0_1 );

            // (D_020 D_011|_{i} = (F_021 P_010| + z_ab * (D_020 P_010|
            S_2_2_1_0[abcd * 108 + 22 * 3 + iket] = Q_0_2_1_0_1_0_1 + ( AB_z[abcd] * Q_0_2_0_0_1_0_1 );

            // (D_020 D_002|_{i} = (F_021 P_001| + z_ab * (D_020 P_001|
            S_2_2_1_0[abcd * 108 + 23 * 3 + iket] = Q_0_2_1_0_0_1_1 + ( AB_z[abcd] * Q_0_2_0_0_0_1_1 );

            // (D_011 D_200|_{i} = (F_111 P_100| + x_ab * (D_011 P_100|
            S_2_2_1_0[abcd * 108 + 24 * 3 + iket] = Q_1_1_1_1_0_0_1 + ( AB_x[abcd] * Q_0_1_1_1_0_0_1 );

            // (D_011 D_110|_{i} = (F_021 P_100| + y_ab * (D_011 P_100|
            S_2_2_1_0[abcd * 108 + 25 * 3 + iket] = Q_0_2_1_1_0_0_1 + ( AB_y[abcd] * Q_0_1_1_1_0_0_1 );

            // (D_011 D_101|_{i} = (F_012 P_100| + z_ab * (D_011 P_100|
            S_2_2_1_0[abcd * 108 + 26 * 3 + iket] = Q_0_1_2_1_0_0_1 + ( AB_z[abcd] * Q_0_1_1_1_0_0_1 );

            // (D_011 D_020|_{i} = (F_021 P_010| + y_ab * (D_011 P_010|
            S_2_2_1_0[abcd * 108 + 27 * 3 + iket] = Q_0_2_1_0_1_0_1 + ( AB_y[abcd] * Q_0_1_1_0_1_0_1 );

            // (D_011 D_011|_{i} = (F_012 P_010| + z_ab * (D_011 P_010|
            S_2_2_1_0[abcd * 108 + 28 * 3 + iket] = Q_0_1_2_0_1_0_1 + ( AB_z[abcd] * Q_0_1_1_0_1_0_1 );

            // (D_011 D_002|_{i} = (F_012 P_001| + z_ab * (D_011 P_001|
            S_2_2_1_0[abcd * 108 + 29 * 3 + iket] = Q_0_1_2_0_0_1_1 + ( AB_z[abcd] * Q_0_1_1_0_0_1_1 );

            // (D_002 D_200|_{i} = (F_102 P_100| + x_ab * (D_002 P_100|
            S_2_2_1_0[abcd * 108 + 30 * 3 + iket] = Q_1_0_2_1_0_0_1 + ( AB_x[abcd] * Q_0_0_2_1_0_0_1 );

            // (D_002 D_110|_{i} = (F_012 P_100| + y_ab * (D_002 P_100|
            S_2_2_1_0[abcd * 108 + 31 * 3 + iket] = Q_0_1_2_1_0_0_1 + ( AB_y[abcd] * Q_0_0_2_1_0_0_1 );

            // (D_002 D_101|_{i} = (F_003 P_100| + z_ab * (D_002 P_100|
            S_2_2_1_0[abcd * 108 + 32 * 3 + iket] = Q_0_0_3_1_0_0_1 + ( AB_z[abcd] * Q_0_0_2_1_0_0_1 );

            // (D_002 D_020|_{i} = (F_012 P_010| + y_ab * (D_002 P_010|
            S_2_2_1_0[abcd * 108 + 33 * 3 + iket] = Q_0_1_2_0_1_0_1 + ( AB_y[abcd] * Q_0_0_2_0_1_0_1 );

            // (D_002 D_011|_{i} = (F_003 P_010| + z_ab * (D_002 P_010|
            S_2_2_1_0[abcd * 108 + 34 * 3 + iket] = Q_0_0_3_0_1_0_1 + ( AB_z[abcd] * Q_0_0_2_0_1_0_1 );

            // (D_002 D_002|_{i} = (F_003 P_001| + z_ab * (D_002 P_001|
            S_2_2_1_0[abcd * 108 + 35 * 3 + iket] = Q_0_0_3_0_0_1_1 + ( AB_z[abcd] * Q_0_0_2_0_0_1_1 );

        }


    }


    //////////////////////////////////////////////
    // Contracted integrals: Horizontal recurrance
    // Ket part
    // Steps: 0
    // Forming final integrals
    //////////////////////////////////////////////

    //Nothing to do.....


    return nshell1234;
}

