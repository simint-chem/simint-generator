#include <string.h>
#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"
#include "boys/boys_split.h"


int eri_split_d_d_s_s(struct multishell_pair const P,
                      struct multishell_pair const Q,
                      double * const restrict S_2_2_0_0)
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
    double * const contwork = malloc(nshell1234 * 248);
    memset(contwork, 0, nshell1234 * 248);

    // partition workspace into shells
    double * const S_2_0_0_0 = contwork + (nshell1234 * 0);
    double * const S_3_0_0_0 = contwork + (nshell1234 * 6);
    double * const S_4_0_0_0 = contwork + (nshell1234 * 16);


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
        double * const restrict PRIM_S_2_0_0_0 = S_2_0_0_0 + (abcd * 6);
        double * const restrict PRIM_S_3_0_0_0 = S_3_0_0_0 + (abcd * 10);
        double * const restrict PRIM_S_4_0_0_0 = S_4_0_0_0 + (abcd * 15);
            // set up pointers to the contracted integrals - Electron Transfer

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
                    double AUX_S_0_0_0_0[5 * 1];

                    // AM = 1: Needed from this AM: 3
                    double AUX_S_1_0_0_0[4 * 3];

                    // AM = 2: Needed from this AM: 6
                    double AUX_S_2_0_0_0[3 * 6];

                    // AM = 3: Needed from this AM: 10
                    double AUX_S_3_0_0_0[2 * 10];

                    // AM = 4: Needed from this AM: 15
                    double AUX_S_4_0_0_0[1 * 15];



                    // Holds temporary integrals for electron transfer


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


                    //////////////////////////////////////////////
                    // Boys function section
                    // Maximum v value: 4
                    //////////////////////////////////////////////
                    // The paremeter to the boys function
                    const double F_x = R2 * alpha;


                    Boys_F_split(AUX_S_0_0_0_0, 4, F_x);
                    AUX_S_0_0_0_0[0] *= allprefac;
                    AUX_S_0_0_0_0[1] *= allprefac;
                    AUX_S_0_0_0_0[2] *= allprefac;
                    AUX_S_0_0_0_0[3] *= allprefac;
                    AUX_S_0_0_0_0[4] *= allprefac;

                    //////////////////////////////////////////////
                    // Primitive integrals: Vertical recurrance
                    //////////////////////////////////////////////

                    // Forming AUX_S_1_0_0_0[4 * 3];
                    // Needed from this AM:
                    //    P_100
                    //    P_010
                    //    P_001
                    for(int m = 0; m < 4; m++)  // loop over orders of boys function
                    {
                        //P_100 : STEP: x
                        AUX_S_1_0_0_0[m * 3 + 0] = P.PA_x[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_x * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                        //P_010 : STEP: y
                        AUX_S_1_0_0_0[m * 3 + 1] = P.PA_y[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_y * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                        //P_001 : STEP: z
                        AUX_S_1_0_0_0[m * 3 + 2] = P.PA_z[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_z * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                    }


                    // Forming AUX_S_2_0_0_0[3 * 6];
                    // Needed from this AM:
                    //    D_200
                    //    D_110
                    //    D_101
                    //    D_020
                    //    D_011
                    //    D_002
                    for(int m = 0; m < 3; m++)  // loop over orders of boys function
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

                    // Accumulating in contracted workspace
                    for(int i = 0; i < 6; i++)
                        PRIM_S_2_0_0_0[i] += AUX_S_2_0_0_0[i];


                    // Forming AUX_S_3_0_0_0[2 * 10];
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
                    for(int m = 0; m < 2; m++)  // loop over orders of boys function
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

                    // Accumulating in contracted workspace
                    for(int i = 0; i < 10; i++)
                        PRIM_S_3_0_0_0[i] += AUX_S_3_0_0_0[i];


                    // Forming AUX_S_4_0_0_0[1 * 15];
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
                    for(int m = 0; m < 1; m++)  // loop over orders of boys function
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

                    // Accumulating in contracted workspace
                    for(int i = 0; i < 15; i++)
                        PRIM_S_4_0_0_0[i] += AUX_S_4_0_0_0[i];




                    //////////////////////////////////////////////
                    // Primitive integrals: Electron transfer
                    //////////////////////////////////////////////

                    //...nothing to do...

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
        // form S_2_2_0_0
        for(int iket = 0; iket < 1; ++iket)
        {
            // (D_200 P_100| = (F_300 S_000|_{t} + x_ab * (D_200 S_000|_{t}
            const double Q_2_0_0_1_0_0_0 = S_3_0_0_0[abcd * 10 + 0 * 1 + iket] + ( AB_x[abcd] * S_2_0_0_0[abcd * 6 + 0 * 1 + iket] );

            // (D_200 P_010| = (F_210 S_000|_{t} + y_ab * (D_200 S_000|_{t}
            const double Q_2_0_0_0_1_0_0 = S_3_0_0_0[abcd * 10 + 1 * 1 + iket] + ( AB_y[abcd] * S_2_0_0_0[abcd * 6 + 0 * 1 + iket] );

            // (D_200 P_001| = (F_201 S_000|_{t} + z_ab * (D_200 S_000|_{t}
            const double Q_2_0_0_0_0_1_0 = S_3_0_0_0[abcd * 10 + 2 * 1 + iket] + ( AB_z[abcd] * S_2_0_0_0[abcd * 6 + 0 * 1 + iket] );

            // (D_110 P_100| = (F_210 S_000|_{t} + x_ab * (D_110 S_000|_{t}
            const double Q_1_1_0_1_0_0_0 = S_3_0_0_0[abcd * 10 + 1 * 1 + iket] + ( AB_x[abcd] * S_2_0_0_0[abcd * 6 + 1 * 1 + iket] );

            // (D_110 P_010| = (F_120 S_000|_{t} + y_ab * (D_110 S_000|_{t}
            const double Q_1_1_0_0_1_0_0 = S_3_0_0_0[abcd * 10 + 3 * 1 + iket] + ( AB_y[abcd] * S_2_0_0_0[abcd * 6 + 1 * 1 + iket] );

            // (D_110 P_001| = (F_111 S_000|_{t} + z_ab * (D_110 S_000|_{t}
            const double Q_1_1_0_0_0_1_0 = S_3_0_0_0[abcd * 10 + 4 * 1 + iket] + ( AB_z[abcd] * S_2_0_0_0[abcd * 6 + 1 * 1 + iket] );

            // (D_101 P_100| = (F_201 S_000|_{t} + x_ab * (D_101 S_000|_{t}
            const double Q_1_0_1_1_0_0_0 = S_3_0_0_0[abcd * 10 + 2 * 1 + iket] + ( AB_x[abcd] * S_2_0_0_0[abcd * 6 + 2 * 1 + iket] );

            // (D_101 P_010| = (F_111 S_000|_{t} + y_ab * (D_101 S_000|_{t}
            const double Q_1_0_1_0_1_0_0 = S_3_0_0_0[abcd * 10 + 4 * 1 + iket] + ( AB_y[abcd] * S_2_0_0_0[abcd * 6 + 2 * 1 + iket] );

            // (D_101 P_001| = (F_102 S_000|_{t} + z_ab * (D_101 S_000|_{t}
            const double Q_1_0_1_0_0_1_0 = S_3_0_0_0[abcd * 10 + 5 * 1 + iket] + ( AB_z[abcd] * S_2_0_0_0[abcd * 6 + 2 * 1 + iket] );

            // (D_020 P_100| = (F_120 S_000|_{t} + x_ab * (D_020 S_000|_{t}
            const double Q_0_2_0_1_0_0_0 = S_3_0_0_0[abcd * 10 + 3 * 1 + iket] + ( AB_x[abcd] * S_2_0_0_0[abcd * 6 + 3 * 1 + iket] );

            // (D_020 P_010| = (F_030 S_000|_{t} + y_ab * (D_020 S_000|_{t}
            const double Q_0_2_0_0_1_0_0 = S_3_0_0_0[abcd * 10 + 6 * 1 + iket] + ( AB_y[abcd] * S_2_0_0_0[abcd * 6 + 3 * 1 + iket] );

            // (D_020 P_001| = (F_021 S_000|_{t} + z_ab * (D_020 S_000|_{t}
            const double Q_0_2_0_0_0_1_0 = S_3_0_0_0[abcd * 10 + 7 * 1 + iket] + ( AB_z[abcd] * S_2_0_0_0[abcd * 6 + 3 * 1 + iket] );

            // (D_011 P_100| = (F_111 S_000|_{t} + x_ab * (D_011 S_000|_{t}
            const double Q_0_1_1_1_0_0_0 = S_3_0_0_0[abcd * 10 + 4 * 1 + iket] + ( AB_x[abcd] * S_2_0_0_0[abcd * 6 + 4 * 1 + iket] );

            // (D_011 P_010| = (F_021 S_000|_{t} + y_ab * (D_011 S_000|_{t}
            const double Q_0_1_1_0_1_0_0 = S_3_0_0_0[abcd * 10 + 7 * 1 + iket] + ( AB_y[abcd] * S_2_0_0_0[abcd * 6 + 4 * 1 + iket] );

            // (D_011 P_001| = (F_012 S_000|_{t} + z_ab * (D_011 S_000|_{t}
            const double Q_0_1_1_0_0_1_0 = S_3_0_0_0[abcd * 10 + 8 * 1 + iket] + ( AB_z[abcd] * S_2_0_0_0[abcd * 6 + 4 * 1 + iket] );

            // (D_002 P_100| = (F_102 S_000|_{t} + x_ab * (D_002 S_000|_{t}
            const double Q_0_0_2_1_0_0_0 = S_3_0_0_0[abcd * 10 + 5 * 1 + iket] + ( AB_x[abcd] * S_2_0_0_0[abcd * 6 + 5 * 1 + iket] );

            // (D_002 P_010| = (F_012 S_000|_{t} + y_ab * (D_002 S_000|_{t}
            const double Q_0_0_2_0_1_0_0 = S_3_0_0_0[abcd * 10 + 8 * 1 + iket] + ( AB_y[abcd] * S_2_0_0_0[abcd * 6 + 5 * 1 + iket] );

            // (D_002 P_001| = (F_003 S_000|_{t} + z_ab * (D_002 S_000|_{t}
            const double Q_0_0_2_0_0_1_0 = S_3_0_0_0[abcd * 10 + 9 * 1 + iket] + ( AB_z[abcd] * S_2_0_0_0[abcd * 6 + 5 * 1 + iket] );

            // (F_300 P_100| = (G_400 S_000|_{t} + x_ab * (F_300 S_000|_{t}
            const double Q_3_0_0_1_0_0_0 = S_4_0_0_0[abcd * 15 + 0 * 1 + iket] + ( AB_x[abcd] * S_3_0_0_0[abcd * 10 + 0 * 1 + iket] );

            // (F_210 P_100| = (G_310 S_000|_{t} + x_ab * (F_210 S_000|_{t}
            const double Q_2_1_0_1_0_0_0 = S_4_0_0_0[abcd * 15 + 1 * 1 + iket] + ( AB_x[abcd] * S_3_0_0_0[abcd * 10 + 1 * 1 + iket] );

            // (F_210 P_010| = (G_220 S_000|_{t} + y_ab * (F_210 S_000|_{t}
            const double Q_2_1_0_0_1_0_0 = S_4_0_0_0[abcd * 15 + 3 * 1 + iket] + ( AB_y[abcd] * S_3_0_0_0[abcd * 10 + 1 * 1 + iket] );

            // (F_201 P_100| = (G_301 S_000|_{t} + x_ab * (F_201 S_000|_{t}
            const double Q_2_0_1_1_0_0_0 = S_4_0_0_0[abcd * 15 + 2 * 1 + iket] + ( AB_x[abcd] * S_3_0_0_0[abcd * 10 + 2 * 1 + iket] );

            // (F_201 P_010| = (G_211 S_000|_{t} + y_ab * (F_201 S_000|_{t}
            const double Q_2_0_1_0_1_0_0 = S_4_0_0_0[abcd * 15 + 4 * 1 + iket] + ( AB_y[abcd] * S_3_0_0_0[abcd * 10 + 2 * 1 + iket] );

            // (F_201 P_001| = (G_202 S_000|_{t} + z_ab * (F_201 S_000|_{t}
            const double Q_2_0_1_0_0_1_0 = S_4_0_0_0[abcd * 15 + 5 * 1 + iket] + ( AB_z[abcd] * S_3_0_0_0[abcd * 10 + 2 * 1 + iket] );

            // (F_120 P_100| = (G_220 S_000|_{t} + x_ab * (F_120 S_000|_{t}
            const double Q_1_2_0_1_0_0_0 = S_4_0_0_0[abcd * 15 + 3 * 1 + iket] + ( AB_x[abcd] * S_3_0_0_0[abcd * 10 + 3 * 1 + iket] );

            // (F_120 P_010| = (G_130 S_000|_{t} + y_ab * (F_120 S_000|_{t}
            const double Q_1_2_0_0_1_0_0 = S_4_0_0_0[abcd * 15 + 6 * 1 + iket] + ( AB_y[abcd] * S_3_0_0_0[abcd * 10 + 3 * 1 + iket] );

            // (F_111 P_100| = (G_211 S_000|_{t} + x_ab * (F_111 S_000|_{t}
            const double Q_1_1_1_1_0_0_0 = S_4_0_0_0[abcd * 15 + 4 * 1 + iket] + ( AB_x[abcd] * S_3_0_0_0[abcd * 10 + 4 * 1 + iket] );

            // (F_111 P_010| = (G_121 S_000|_{t} + y_ab * (F_111 S_000|_{t}
            const double Q_1_1_1_0_1_0_0 = S_4_0_0_0[abcd * 15 + 7 * 1 + iket] + ( AB_y[abcd] * S_3_0_0_0[abcd * 10 + 4 * 1 + iket] );

            // (F_111 P_001| = (G_112 S_000|_{t} + z_ab * (F_111 S_000|_{t}
            const double Q_1_1_1_0_0_1_0 = S_4_0_0_0[abcd * 15 + 8 * 1 + iket] + ( AB_z[abcd] * S_3_0_0_0[abcd * 10 + 4 * 1 + iket] );

            // (F_102 P_100| = (G_202 S_000|_{t} + x_ab * (F_102 S_000|_{t}
            const double Q_1_0_2_1_0_0_0 = S_4_0_0_0[abcd * 15 + 5 * 1 + iket] + ( AB_x[abcd] * S_3_0_0_0[abcd * 10 + 5 * 1 + iket] );

            // (F_102 P_010| = (G_112 S_000|_{t} + y_ab * (F_102 S_000|_{t}
            const double Q_1_0_2_0_1_0_0 = S_4_0_0_0[abcd * 15 + 8 * 1 + iket] + ( AB_y[abcd] * S_3_0_0_0[abcd * 10 + 5 * 1 + iket] );

            // (F_102 P_001| = (G_103 S_000|_{t} + z_ab * (F_102 S_000|_{t}
            const double Q_1_0_2_0_0_1_0 = S_4_0_0_0[abcd * 15 + 9 * 1 + iket] + ( AB_z[abcd] * S_3_0_0_0[abcd * 10 + 5 * 1 + iket] );

            // (F_030 P_100| = (G_130 S_000|_{t} + x_ab * (F_030 S_000|_{t}
            const double Q_0_3_0_1_0_0_0 = S_4_0_0_0[abcd * 15 + 6 * 1 + iket] + ( AB_x[abcd] * S_3_0_0_0[abcd * 10 + 6 * 1 + iket] );

            // (F_030 P_010| = (G_040 S_000|_{t} + y_ab * (F_030 S_000|_{t}
            const double Q_0_3_0_0_1_0_0 = S_4_0_0_0[abcd * 15 + 10 * 1 + iket] + ( AB_y[abcd] * S_3_0_0_0[abcd * 10 + 6 * 1 + iket] );

            // (F_021 P_100| = (G_121 S_000|_{t} + x_ab * (F_021 S_000|_{t}
            const double Q_0_2_1_1_0_0_0 = S_4_0_0_0[abcd * 15 + 7 * 1 + iket] + ( AB_x[abcd] * S_3_0_0_0[abcd * 10 + 7 * 1 + iket] );

            // (F_021 P_010| = (G_031 S_000|_{t} + y_ab * (F_021 S_000|_{t}
            const double Q_0_2_1_0_1_0_0 = S_4_0_0_0[abcd * 15 + 11 * 1 + iket] + ( AB_y[abcd] * S_3_0_0_0[abcd * 10 + 7 * 1 + iket] );

            // (F_021 P_001| = (G_022 S_000|_{t} + z_ab * (F_021 S_000|_{t}
            const double Q_0_2_1_0_0_1_0 = S_4_0_0_0[abcd * 15 + 12 * 1 + iket] + ( AB_z[abcd] * S_3_0_0_0[abcd * 10 + 7 * 1 + iket] );

            // (F_012 P_100| = (G_112 S_000|_{t} + x_ab * (F_012 S_000|_{t}
            const double Q_0_1_2_1_0_0_0 = S_4_0_0_0[abcd * 15 + 8 * 1 + iket] + ( AB_x[abcd] * S_3_0_0_0[abcd * 10 + 8 * 1 + iket] );

            // (F_012 P_010| = (G_022 S_000|_{t} + y_ab * (F_012 S_000|_{t}
            const double Q_0_1_2_0_1_0_0 = S_4_0_0_0[abcd * 15 + 12 * 1 + iket] + ( AB_y[abcd] * S_3_0_0_0[abcd * 10 + 8 * 1 + iket] );

            // (F_012 P_001| = (G_013 S_000|_{t} + z_ab * (F_012 S_000|_{t}
            const double Q_0_1_2_0_0_1_0 = S_4_0_0_0[abcd * 15 + 13 * 1 + iket] + ( AB_z[abcd] * S_3_0_0_0[abcd * 10 + 8 * 1 + iket] );

            // (F_003 P_100| = (G_103 S_000|_{t} + x_ab * (F_003 S_000|_{t}
            const double Q_0_0_3_1_0_0_0 = S_4_0_0_0[abcd * 15 + 9 * 1 + iket] + ( AB_x[abcd] * S_3_0_0_0[abcd * 10 + 9 * 1 + iket] );

            // (F_003 P_010| = (G_013 S_000|_{t} + y_ab * (F_003 S_000|_{t}
            const double Q_0_0_3_0_1_0_0 = S_4_0_0_0[abcd * 15 + 13 * 1 + iket] + ( AB_y[abcd] * S_3_0_0_0[abcd * 10 + 9 * 1 + iket] );

            // (F_003 P_001| = (G_004 S_000|_{t} + z_ab * (F_003 S_000|_{t}
            const double Q_0_0_3_0_0_1_0 = S_4_0_0_0[abcd * 15 + 14 * 1 + iket] + ( AB_z[abcd] * S_3_0_0_0[abcd * 10 + 9 * 1 + iket] );

            // (D_200 D_200|_{i} = (F_300 P_100| + x_ab * (D_200 P_100|
            S_2_2_0_0[abcd * 36 + 0 * 1 + iket] = Q_3_0_0_1_0_0_0 + ( AB_x[abcd] * Q_2_0_0_1_0_0_0 );

            // (D_200 D_110|_{i} = (F_210 P_100| + y_ab * (D_200 P_100|
            S_2_2_0_0[abcd * 36 + 1 * 1 + iket] = Q_2_1_0_1_0_0_0 + ( AB_y[abcd] * Q_2_0_0_1_0_0_0 );

            // (D_200 D_101|_{i} = (F_201 P_100| + z_ab * (D_200 P_100|
            S_2_2_0_0[abcd * 36 + 2 * 1 + iket] = Q_2_0_1_1_0_0_0 + ( AB_z[abcd] * Q_2_0_0_1_0_0_0 );

            // (D_200 D_020|_{i} = (F_210 P_010| + y_ab * (D_200 P_010|
            S_2_2_0_0[abcd * 36 + 3 * 1 + iket] = Q_2_1_0_0_1_0_0 + ( AB_y[abcd] * Q_2_0_0_0_1_0_0 );

            // (D_200 D_011|_{i} = (F_201 P_010| + z_ab * (D_200 P_010|
            S_2_2_0_0[abcd * 36 + 4 * 1 + iket] = Q_2_0_1_0_1_0_0 + ( AB_z[abcd] * Q_2_0_0_0_1_0_0 );

            // (D_200 D_002|_{i} = (F_201 P_001| + z_ab * (D_200 P_001|
            S_2_2_0_0[abcd * 36 + 5 * 1 + iket] = Q_2_0_1_0_0_1_0 + ( AB_z[abcd] * Q_2_0_0_0_0_1_0 );

            // (D_110 D_200|_{i} = (F_210 P_100| + x_ab * (D_110 P_100|
            S_2_2_0_0[abcd * 36 + 6 * 1 + iket] = Q_2_1_0_1_0_0_0 + ( AB_x[abcd] * Q_1_1_0_1_0_0_0 );

            // (D_110 D_110|_{i} = (F_120 P_100| + y_ab * (D_110 P_100|
            S_2_2_0_0[abcd * 36 + 7 * 1 + iket] = Q_1_2_0_1_0_0_0 + ( AB_y[abcd] * Q_1_1_0_1_0_0_0 );

            // (D_110 D_101|_{i} = (F_111 P_100| + z_ab * (D_110 P_100|
            S_2_2_0_0[abcd * 36 + 8 * 1 + iket] = Q_1_1_1_1_0_0_0 + ( AB_z[abcd] * Q_1_1_0_1_0_0_0 );

            // (D_110 D_020|_{i} = (F_120 P_010| + y_ab * (D_110 P_010|
            S_2_2_0_0[abcd * 36 + 9 * 1 + iket] = Q_1_2_0_0_1_0_0 + ( AB_y[abcd] * Q_1_1_0_0_1_0_0 );

            // (D_110 D_011|_{i} = (F_111 P_010| + z_ab * (D_110 P_010|
            S_2_2_0_0[abcd * 36 + 10 * 1 + iket] = Q_1_1_1_0_1_0_0 + ( AB_z[abcd] * Q_1_1_0_0_1_0_0 );

            // (D_110 D_002|_{i} = (F_111 P_001| + z_ab * (D_110 P_001|
            S_2_2_0_0[abcd * 36 + 11 * 1 + iket] = Q_1_1_1_0_0_1_0 + ( AB_z[abcd] * Q_1_1_0_0_0_1_0 );

            // (D_101 D_200|_{i} = (F_201 P_100| + x_ab * (D_101 P_100|
            S_2_2_0_0[abcd * 36 + 12 * 1 + iket] = Q_2_0_1_1_0_0_0 + ( AB_x[abcd] * Q_1_0_1_1_0_0_0 );

            // (D_101 D_110|_{i} = (F_111 P_100| + y_ab * (D_101 P_100|
            S_2_2_0_0[abcd * 36 + 13 * 1 + iket] = Q_1_1_1_1_0_0_0 + ( AB_y[abcd] * Q_1_0_1_1_0_0_0 );

            // (D_101 D_101|_{i} = (F_102 P_100| + z_ab * (D_101 P_100|
            S_2_2_0_0[abcd * 36 + 14 * 1 + iket] = Q_1_0_2_1_0_0_0 + ( AB_z[abcd] * Q_1_0_1_1_0_0_0 );

            // (D_101 D_020|_{i} = (F_111 P_010| + y_ab * (D_101 P_010|
            S_2_2_0_0[abcd * 36 + 15 * 1 + iket] = Q_1_1_1_0_1_0_0 + ( AB_y[abcd] * Q_1_0_1_0_1_0_0 );

            // (D_101 D_011|_{i} = (F_102 P_010| + z_ab * (D_101 P_010|
            S_2_2_0_0[abcd * 36 + 16 * 1 + iket] = Q_1_0_2_0_1_0_0 + ( AB_z[abcd] * Q_1_0_1_0_1_0_0 );

            // (D_101 D_002|_{i} = (F_102 P_001| + z_ab * (D_101 P_001|
            S_2_2_0_0[abcd * 36 + 17 * 1 + iket] = Q_1_0_2_0_0_1_0 + ( AB_z[abcd] * Q_1_0_1_0_0_1_0 );

            // (D_020 D_200|_{i} = (F_120 P_100| + x_ab * (D_020 P_100|
            S_2_2_0_0[abcd * 36 + 18 * 1 + iket] = Q_1_2_0_1_0_0_0 + ( AB_x[abcd] * Q_0_2_0_1_0_0_0 );

            // (D_020 D_110|_{i} = (F_030 P_100| + y_ab * (D_020 P_100|
            S_2_2_0_0[abcd * 36 + 19 * 1 + iket] = Q_0_3_0_1_0_0_0 + ( AB_y[abcd] * Q_0_2_0_1_0_0_0 );

            // (D_020 D_101|_{i} = (F_021 P_100| + z_ab * (D_020 P_100|
            S_2_2_0_0[abcd * 36 + 20 * 1 + iket] = Q_0_2_1_1_0_0_0 + ( AB_z[abcd] * Q_0_2_0_1_0_0_0 );

            // (D_020 D_020|_{i} = (F_030 P_010| + y_ab * (D_020 P_010|
            S_2_2_0_0[abcd * 36 + 21 * 1 + iket] = Q_0_3_0_0_1_0_0 + ( AB_y[abcd] * Q_0_2_0_0_1_0_0 );

            // (D_020 D_011|_{i} = (F_021 P_010| + z_ab * (D_020 P_010|
            S_2_2_0_0[abcd * 36 + 22 * 1 + iket] = Q_0_2_1_0_1_0_0 + ( AB_z[abcd] * Q_0_2_0_0_1_0_0 );

            // (D_020 D_002|_{i} = (F_021 P_001| + z_ab * (D_020 P_001|
            S_2_2_0_0[abcd * 36 + 23 * 1 + iket] = Q_0_2_1_0_0_1_0 + ( AB_z[abcd] * Q_0_2_0_0_0_1_0 );

            // (D_011 D_200|_{i} = (F_111 P_100| + x_ab * (D_011 P_100|
            S_2_2_0_0[abcd * 36 + 24 * 1 + iket] = Q_1_1_1_1_0_0_0 + ( AB_x[abcd] * Q_0_1_1_1_0_0_0 );

            // (D_011 D_110|_{i} = (F_021 P_100| + y_ab * (D_011 P_100|
            S_2_2_0_0[abcd * 36 + 25 * 1 + iket] = Q_0_2_1_1_0_0_0 + ( AB_y[abcd] * Q_0_1_1_1_0_0_0 );

            // (D_011 D_101|_{i} = (F_012 P_100| + z_ab * (D_011 P_100|
            S_2_2_0_0[abcd * 36 + 26 * 1 + iket] = Q_0_1_2_1_0_0_0 + ( AB_z[abcd] * Q_0_1_1_1_0_0_0 );

            // (D_011 D_020|_{i} = (F_021 P_010| + y_ab * (D_011 P_010|
            S_2_2_0_0[abcd * 36 + 27 * 1 + iket] = Q_0_2_1_0_1_0_0 + ( AB_y[abcd] * Q_0_1_1_0_1_0_0 );

            // (D_011 D_011|_{i} = (F_012 P_010| + z_ab * (D_011 P_010|
            S_2_2_0_0[abcd * 36 + 28 * 1 + iket] = Q_0_1_2_0_1_0_0 + ( AB_z[abcd] * Q_0_1_1_0_1_0_0 );

            // (D_011 D_002|_{i} = (F_012 P_001| + z_ab * (D_011 P_001|
            S_2_2_0_0[abcd * 36 + 29 * 1 + iket] = Q_0_1_2_0_0_1_0 + ( AB_z[abcd] * Q_0_1_1_0_0_1_0 );

            // (D_002 D_200|_{i} = (F_102 P_100| + x_ab * (D_002 P_100|
            S_2_2_0_0[abcd * 36 + 30 * 1 + iket] = Q_1_0_2_1_0_0_0 + ( AB_x[abcd] * Q_0_0_2_1_0_0_0 );

            // (D_002 D_110|_{i} = (F_012 P_100| + y_ab * (D_002 P_100|
            S_2_2_0_0[abcd * 36 + 31 * 1 + iket] = Q_0_1_2_1_0_0_0 + ( AB_y[abcd] * Q_0_0_2_1_0_0_0 );

            // (D_002 D_101|_{i} = (F_003 P_100| + z_ab * (D_002 P_100|
            S_2_2_0_0[abcd * 36 + 32 * 1 + iket] = Q_0_0_3_1_0_0_0 + ( AB_z[abcd] * Q_0_0_2_1_0_0_0 );

            // (D_002 D_020|_{i} = (F_012 P_010| + y_ab * (D_002 P_010|
            S_2_2_0_0[abcd * 36 + 33 * 1 + iket] = Q_0_1_2_0_1_0_0 + ( AB_y[abcd] * Q_0_0_2_0_1_0_0 );

            // (D_002 D_011|_{i} = (F_003 P_010| + z_ab * (D_002 P_010|
            S_2_2_0_0[abcd * 36 + 34 * 1 + iket] = Q_0_0_3_0_1_0_0 + ( AB_z[abcd] * Q_0_0_2_0_1_0_0 );

            // (D_002 D_002|_{i} = (F_003 P_001| + z_ab * (D_002 P_001|
            S_2_2_0_0[abcd * 36 + 35 * 1 + iket] = Q_0_0_3_0_0_1_0 + ( AB_z[abcd] * Q_0_0_2_0_0_1_0 );

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

