#include <string.h>
#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"
#include "boys/boys_split.h"


int eri_split_d_d_d_d(struct multishell_pair const P,
                      struct multishell_pair const Q,
                      double * const restrict INT__d_d_d_d)
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

    ASSUME_ALIGN(INT__d_d_d_d);

    const int nshell1234 = P.nshell12 * Q.nshell12;


    // Holds AB_{xyz} and CD_{xyz} in a flattened fashion for later
    double AB_x[nshell1234];  double CD_x[nshell1234];
    double AB_y[nshell1234];  double CD_y[nshell1234];
    double AB_z[nshell1234];  double CD_z[nshell1234];

    int ab, cd, abcd;
    int i, j;

    // Workspace for contracted integrals
    double * const contwork = malloc(nshell1234 * 16616);
    memset(contwork, 0, nshell1234 * 16616);

    // partition workspace into shells
    double * const INT__d_s_d_s = contwork + (nshell1234 * 0);
    double * const INT__d_s_f_s = contwork + (nshell1234 * 36);
    double * const INT__d_s_g_s = contwork + (nshell1234 * 96);
    double * const INT__d_d_d_s = contwork + (nshell1234 * 186);
    double * const INT__d_d_f_s = contwork + (nshell1234 * 402);
    double * const INT__d_d_g_s = contwork + (nshell1234 * 762);
    double * const INT__f_s_d_s = contwork + (nshell1234 * 1302);
    double * const INT__f_s_f_s = contwork + (nshell1234 * 1362);
    double * const INT__f_s_g_s = contwork + (nshell1234 * 1462);
    double * const INT__g_s_d_s = contwork + (nshell1234 * 1612);
    double * const INT__g_s_f_s = contwork + (nshell1234 * 1702);
    double * const INT__g_s_g_s = contwork + (nshell1234 * 1852);


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
            double * const restrict PRIM_INT__d_s_d_s = INT__d_s_d_s + (abcd * 36);
            double * const restrict PRIM_INT__d_s_f_s = INT__d_s_f_s + (abcd * 60);
            double * const restrict PRIM_INT__d_s_g_s = INT__d_s_g_s + (abcd * 90);
            double * const restrict PRIM_INT__f_s_d_s = INT__f_s_d_s + (abcd * 60);
            double * const restrict PRIM_INT__f_s_f_s = INT__f_s_f_s + (abcd * 100);
            double * const restrict PRIM_INT__f_s_g_s = INT__f_s_g_s + (abcd * 150);
            double * const restrict PRIM_INT__g_s_d_s = INT__g_s_d_s + (abcd * 90);
            double * const restrict PRIM_INT__g_s_f_s = INT__g_s_f_s + (abcd * 150);
            double * const restrict PRIM_INT__g_s_g_s = INT__g_s_g_s + (abcd * 225);

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
                    double AUX_INT__s_s_s_s[9 * 1];

                    // AM = 1: Needed from this AM: 3
                    double AUX_INT__p_s_s_s[8 * 3];

                    // AM = 2: Needed from this AM: 6
                    double AUX_INT__d_s_s_s[7 * 6];

                    // AM = 3: Needed from this AM: 10
                    double AUX_INT__f_s_s_s[6 * 10];

                    // AM = 4: Needed from this AM: 15
                    double AUX_INT__g_s_s_s[5 * 15];

                    // AM = 5: Needed from this AM: 21
                    double AUX_INT__h_s_s_s[4 * 21];

                    // AM = 6: Needed from this AM: 28
                    double AUX_INT__i_s_s_s[3 * 28];

                    // AM = 7: Needed from this AM: 36
                    double AUX_INT__j_s_s_s[2 * 36];

                    // AM = 8: Needed from this AM: 45
                    double AUX_INT__k_s_s_s[1 * 45];



                    // Holds temporary integrals for electron transfer
                    double AUX_INT__s_s_p_s[3];
                    double AUX_INT__s_s_d_s[6];
                    double AUX_INT__p_s_p_s[9];
                    double AUX_INT__p_s_d_s[18];
                    double AUX_INT__p_s_f_s[30];
                    double AUX_INT__d_s_p_s[18];
                    double AUX_INT__d_s_d_s[36];
                    double AUX_INT__d_s_f_s[60];
                    double AUX_INT__d_s_g_s[90];
                    double AUX_INT__f_s_p_s[30];
                    double AUX_INT__f_s_d_s[60];
                    double AUX_INT__f_s_f_s[100];
                    double AUX_INT__f_s_g_s[150];
                    double AUX_INT__g_s_p_s[45];
                    double AUX_INT__g_s_d_s[90];
                    double AUX_INT__g_s_f_s[150];
                    double AUX_INT__g_s_g_s[225];
                    double AUX_INT__h_s_p_s[63];
                    double AUX_INT__h_s_d_s[126];
                    double AUX_INT__h_s_f_s[210];
                    double AUX_INT__i_s_p_s[84];
                    double AUX_INT__i_s_d_s[168];
                    double AUX_INT__j_s_p_s[108];


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
                    // Maximum v value: 8
                    //////////////////////////////////////////////
                    // The paremeter to the boys function
                    const double F_x = R2 * alpha;


                    Boys_F_split(AUX_INT__s_s_s_s, 8, F_x);
                    AUX_INT__s_s_s_s[0] *= allprefac;
                    AUX_INT__s_s_s_s[1] *= allprefac;
                    AUX_INT__s_s_s_s[2] *= allprefac;
                    AUX_INT__s_s_s_s[3] *= allprefac;
                    AUX_INT__s_s_s_s[4] *= allprefac;
                    AUX_INT__s_s_s_s[5] *= allprefac;
                    AUX_INT__s_s_s_s[6] *= allprefac;
                    AUX_INT__s_s_s_s[7] *= allprefac;
                    AUX_INT__s_s_s_s[8] *= allprefac;

                    //////////////////////////////////////////////
                    // Primitive integrals: Vertical recurrance
                    //////////////////////////////////////////////

                    // Forming AUX_INT__p_s_s_s[8 * 3];
                    // Needed from this AM:
                    //    P_100
                    //    P_010
                    //    P_001
                    for(int m = 0; m < 8; m++)  // loop over orders of boys function
                    {
                        //P_100 : STEP: x
                        AUX_INT__p_s_s_s[m * 3 + 0] = P_PA_x * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_x * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_010 : STEP: y
                        AUX_INT__p_s_s_s[m * 3 + 1] = P_PA_y * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_y * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_001 : STEP: z
                        AUX_INT__p_s_s_s[m * 3 + 2] = P_PA_z * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_z * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                    }


                    // Forming AUX_INT__d_s_s_s[7 * 6];
                    // Needed from this AM:
                    //    D_200
                    //    D_110
                    //    D_101
                    //    D_020
                    //    D_011
                    //    D_002
                    for(int m = 0; m < 7; m++)  // loop over orders of boys function
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


                    // Forming AUX_INT__f_s_s_s[6 * 10];
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
                    for(int m = 0; m < 6; m++)  // loop over orders of boys function
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


                    // Forming AUX_INT__g_s_s_s[5 * 15];
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
                    for(int m = 0; m < 5; m++)  // loop over orders of boys function
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


                    // Forming AUX_INT__h_s_s_s[4 * 21];
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
                    for(int m = 0; m < 4; m++)  // loop over orders of boys function
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


                    // Forming AUX_INT__i_s_s_s[3 * 28];
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
                    for(int m = 0; m < 3; m++)  // loop over orders of boys function
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


                    // Forming AUX_INT__j_s_s_s[2 * 36];
                    // Needed from this AM:
                    //    J_700
                    //    J_610
                    //    J_601
                    //    J_520
                    //    J_511
                    //    J_502
                    //    J_430
                    //    J_421
                    //    J_412
                    //    J_403
                    //    J_340
                    //    J_331
                    //    J_322
                    //    J_313
                    //    J_304
                    //    J_250
                    //    J_241
                    //    J_232
                    //    J_223
                    //    J_214
                    //    J_205
                    //    J_160
                    //    J_151
                    //    J_142
                    //    J_133
                    //    J_124
                    //    J_115
                    //    J_106
                    //    J_070
                    //    J_061
                    //    J_052
                    //    J_043
                    //    J_034
                    //    J_025
                    //    J_016
                    //    J_007
                    for(int m = 0; m < 2; m++)  // loop over orders of boys function
                    {
                        //J_700 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 0] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 0] - a_over_p * PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 0]
                                      + 6 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  0] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 0] );

                        //J_610 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 1] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 0] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 0];

                        //J_601 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 2] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 0] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 0];

                        //J_520 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 3] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 1] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  0] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 0] );

                        //J_511 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 4] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 1] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 1];

                        //J_502 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 5] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 2] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  0] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 0] );

                        //J_430 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 6] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 3] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  1] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 1] );

                        //J_421 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 7] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 3] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 3];

                        //J_412 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 8] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 5] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 5];

                        //J_403 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 9] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 5] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  2] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 2] );

                        //J_340 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 10] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 10] - a_over_p * PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 10]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  10] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 10] );

                        //J_331 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 11] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 6] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 6];

                        //J_322 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 12] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 7] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  3] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 3] );

                        //J_313 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 13] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 9] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 9];

                        //J_304 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 14] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 14] - a_over_p * PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 14]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  14] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 14] );

                        //J_250 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 15] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 15] - a_over_p * PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 15]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  15] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 15] );

                        //J_241 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 16] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 10] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 10];

                        //J_232 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 17] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 11] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  6] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 6] );

                        //J_223 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 18] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 13] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  9] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 9] );

                        //J_214 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 19] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 14] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 14];

                        //J_205 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 20] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 20] - a_over_p * PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 20]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  20] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 20] );

                        //J_160 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 21] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 21] - a_over_p * PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 21];

                        //J_151 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 22] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 15] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 15];

                        //J_142 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 23] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 23] - a_over_p * PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 23];

                        //J_133 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 24] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 24] - a_over_p * PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 24];

                        //J_124 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 25] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 25] - a_over_p * PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 25];

                        //J_115 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 26] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 20] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 20];

                        //J_106 : STEP: x
                        AUX_INT__j_s_s_s[m * 36 + 27] = P_PA_x * AUX_INT__i_s_s_s[m * 28 + 27] - a_over_p * PQ_x * AUX_INT__i_s_s_s[(m+1) * 28 + 27];

                        //J_070 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 28] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 21] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 21]
                                      + 6 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  15] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 15] );

                        //J_061 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 29] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 21] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 21];

                        //J_052 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 30] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 22] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 22]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  15] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 15] );

                        //J_043 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 31] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 23] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 23]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  16] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 16] );

                        //J_034 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 32] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 25] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 25]
                                      + 2 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  19] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 19] );

                        //J_025 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 33] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 26] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 26]
                                      + 1 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  20] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 20] );

                        //J_016 : STEP: y
                        AUX_INT__j_s_s_s[m * 36 + 34] = P_PA_y * AUX_INT__i_s_s_s[m * 28 + 27] - a_over_p * PQ_y * AUX_INT__i_s_s_s[(m+1) * 28 + 27];

                        //J_007 : STEP: z
                        AUX_INT__j_s_s_s[m * 36 + 35] = P_PA_z * AUX_INT__i_s_s_s[m * 28 + 27] - a_over_p * PQ_z * AUX_INT__i_s_s_s[(m+1) * 28 + 27]
                                      + 6 * one_over_2p * ( AUX_INT__h_s_s_s[m * 21 +  20] - a_over_p * AUX_INT__h_s_s_s[(m+1) * 21 + 20] );

                    }


                    // Forming AUX_INT__k_s_s_s[1 * 45];
                    // Needed from this AM:
                    //    K_800
                    //    K_710
                    //    K_701
                    //    K_620
                    //    K_611
                    //    K_602
                    //    K_530
                    //    K_521
                    //    K_512
                    //    K_503
                    //    K_440
                    //    K_431
                    //    K_422
                    //    K_413
                    //    K_404
                    //    K_350
                    //    K_341
                    //    K_332
                    //    K_323
                    //    K_314
                    //    K_305
                    //    K_260
                    //    K_251
                    //    K_242
                    //    K_233
                    //    K_224
                    //    K_215
                    //    K_206
                    //    K_170
                    //    K_161
                    //    K_152
                    //    K_143
                    //    K_134
                    //    K_125
                    //    K_116
                    //    K_107
                    //    K_080
                    //    K_071
                    //    K_062
                    //    K_053
                    //    K_044
                    //    K_035
                    //    K_026
                    //    K_017
                    //    K_008
                    for(int m = 0; m < 1; m++)  // loop over orders of boys function
                    {
                        //K_800 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 0] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 0] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 0]
                                      + 7 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  0] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 0] );

                        //K_710 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 1] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 0] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 0];

                        //K_701 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 2] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 0] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 0];

                        //K_620 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 3] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 1] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  0] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 0] );

                        //K_611 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 4] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 1] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 1];

                        //K_602 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 5] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 2] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  0] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 0] );

                        //K_530 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 6] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 3] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  1] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 1] );

                        //K_521 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 7] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 3] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 3];

                        //K_512 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 8] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 5] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 5];

                        //K_503 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 9] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 5] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  2] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 2] );

                        //K_440 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 10] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 6] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 6]
                                      + 3 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  3] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 3] );

                        //K_431 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 11] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 6] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 6];

                        //K_422 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 12] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 7] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 7]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  3] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 3] );

                        //K_413 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 13] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 9] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 9];

                        //K_404 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 14] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 9] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 9]
                                      + 3 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  5] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 5] );

                        //K_350 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 15] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 15] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 15]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  15] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 15] );

                        //K_341 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 16] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 10] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 10];

                        //K_332 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 17] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 11] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 11]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  6] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 6] );

                        //K_323 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 18] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 13] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 13]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  9] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 9] );

                        //K_314 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 19] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 14] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 14];

                        //K_305 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 20] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 20] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 20]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  20] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 20] );

                        //K_260 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 21] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 21] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 21]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  21] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 21] );

                        //K_251 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 22] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 15] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 15];

                        //K_242 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 23] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 16] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 16]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  10] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 10] );

                        //K_233 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 24] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 24] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 24]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  24] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 24] );

                        //K_224 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 25] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 19] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 19]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  14] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 14] );

                        //K_215 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 26] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 20] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 20];

                        //K_206 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 27] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 27] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 27]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  27] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 27] );

                        //K_170 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 28] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 28] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 28];

                        //K_161 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 29] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 21] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 21];

                        //K_152 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 30] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 30] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 30];

                        //K_143 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 31] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 31] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 31];

                        //K_134 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 32] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 32] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 32];

                        //K_125 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 33] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 33] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 33];

                        //K_116 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 34] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 27] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 27];

                        //K_107 : STEP: x
                        AUX_INT__k_s_s_s[m * 45 + 35] = P_PA_x * AUX_INT__j_s_s_s[m * 36 + 35] - a_over_p * PQ_x * AUX_INT__j_s_s_s[(m+1) * 36 + 35];

                        //K_080 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 36] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 28] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 28]
                                      + 7 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  21] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 21] );

                        //K_071 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 37] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 28] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 28];

                        //K_062 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 38] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 29] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 29]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  21] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 21] );

                        //K_053 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 39] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 30] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 30]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  22] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 22] );

                        //K_044 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 40] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 31] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 31]
                                      + 3 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  23] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 23] );

                        //K_035 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 41] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 33] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 33]
                                      + 2 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  26] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 26] );

                        //K_026 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 42] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 34] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 34]
                                      + 1 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  27] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 27] );

                        //K_017 : STEP: y
                        AUX_INT__k_s_s_s[m * 45 + 43] = P_PA_y * AUX_INT__j_s_s_s[m * 36 + 35] - a_over_p * PQ_y * AUX_INT__j_s_s_s[(m+1) * 36 + 35];

                        //K_008 : STEP: z
                        AUX_INT__k_s_s_s[m * 45 + 44] = P_PA_z * AUX_INT__j_s_s_s[m * 36 + 35] - a_over_p * PQ_z * AUX_INT__j_s_s_s[(m+1) * 36 + 35]
                                      + 7 * one_over_2p * ( AUX_INT__i_s_s_s[m * 28 +  27] - a_over_p * AUX_INT__i_s_s_s[(m+1) * 28 + 27] );

                    }




                    //////////////////////////////////////////////
                    // Primitive integrals: Electron transfer
                    //////////////////////////////////////////////

                    // ( D_200 S_000 | P_100 S_000 )^0 = x * ( D_200 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_300 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[0] = etfac[0] * AUX_INT__d_s_s_s[0] + 2 * one_over_2q * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__f_s_s_s[0];

                    // ( D_200 S_000 | P_010 S_000 )^0 = y * ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[1] = etfac[1] * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__f_s_s_s[1];

                    // ( D_200 S_000 | P_001 S_000 )^0 = z * ( D_200 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[2] = etfac[2] * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__f_s_s_s[2];

                    // ( D_110 S_000 | P_100 S_000 )^0 = x * ( D_110 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_210 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[3] = etfac[0] * AUX_INT__d_s_s_s[1] + 1 * one_over_2q * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__f_s_s_s[1];

                    // ( D_110 S_000 | P_010 S_000 )^0 = y * ( D_110 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[4] = etfac[1] * AUX_INT__d_s_s_s[1] + 1 * one_over_2q * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__f_s_s_s[3];

                    // ( D_110 S_000 | P_001 S_000 )^0 = z * ( D_110 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[5] = etfac[2] * AUX_INT__d_s_s_s[1] - p_over_q * AUX_INT__f_s_s_s[4];

                    // ( D_101 S_000 | P_100 S_000 )^0 = x * ( D_101 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_201 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[6] = etfac[0] * AUX_INT__d_s_s_s[2] + 1 * one_over_2q * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__f_s_s_s[2];

                    // ( D_101 S_000 | P_010 S_000 )^0 = y * ( D_101 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[7] = etfac[1] * AUX_INT__d_s_s_s[2] - p_over_q * AUX_INT__f_s_s_s[4];

                    // ( D_101 S_000 | P_001 S_000 )^0 = z * ( D_101 S_000 | S_000 S_000 )^0_{e} + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[8] = etfac[2] * AUX_INT__d_s_s_s[2] + 1 * one_over_2q * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__f_s_s_s[5];

                    // ( D_020 S_000 | P_100 S_000 )^0 = x * ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_120 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[9] = etfac[0] * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__f_s_s_s[3];

                    // ( D_020 S_000 | P_010 S_000 )^0 = y * ( D_020 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_030 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[10] = etfac[1] * AUX_INT__d_s_s_s[3] + 2 * one_over_2q * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__f_s_s_s[6];

                    // ( D_020 S_000 | P_001 S_000 )^0 = z * ( D_020 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[11] = etfac[2] * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__f_s_s_s[7];

                    // ( D_011 S_000 | P_100 S_000 )^0 = x * ( D_011 S_000 | S_000 S_000 )^0_{e} - ( F_111 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[12] = etfac[0] * AUX_INT__d_s_s_s[4] - p_over_q * AUX_INT__f_s_s_s[4];

                    // ( D_011 S_000 | P_010 S_000 )^0 = y * ( D_011 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_021 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[13] = etfac[1] * AUX_INT__d_s_s_s[4] + 1 * one_over_2q * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__f_s_s_s[7];

                    // ( D_011 S_000 | P_001 S_000 )^0 = z * ( D_011 S_000 | S_000 S_000 )^0_{e} + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[14] = etfac[2] * AUX_INT__d_s_s_s[4] + 1 * one_over_2q * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__f_s_s_s[8];

                    // ( D_002 S_000 | P_100 S_000 )^0 = x * ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_102 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[15] = etfac[0] * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__f_s_s_s[5];

                    // ( D_002 S_000 | P_010 S_000 )^0 = y * ( D_002 S_000 | S_000 S_000 )^0_{e} - ( F_012 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[16] = etfac[1] * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__f_s_s_s[8];

                    // ( D_002 S_000 | P_001 S_000 )^0 = z * ( D_002 S_000 | S_000 S_000 )^0_{e} + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( F_003 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__d_s_p_s[17] = etfac[2] * AUX_INT__d_s_s_s[5] + 2 * one_over_2q * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__f_s_s_s[9];

                    // ( F_300 S_000 | P_100 S_000 )^0 = x * ( F_300 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_400 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[0] = etfac[0] * AUX_INT__f_s_s_s[0] + 3 * one_over_2q * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__g_s_s_s[0];

                    // ( F_300 S_000 | P_010 S_000 )^0 = y * ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[1] = etfac[1] * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__g_s_s_s[1];

                    // ( F_300 S_000 | P_001 S_000 )^0 = z * ( F_300 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[2] = etfac[2] * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__g_s_s_s[2];

                    // ( F_210 S_000 | P_100 S_000 )^0 = x * ( F_210 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_310 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[3] = etfac[0] * AUX_INT__f_s_s_s[1] + 2 * one_over_2q * AUX_INT__d_s_s_s[1] - p_over_q * AUX_INT__g_s_s_s[1];

                    // ( F_210 S_000 | P_010 S_000 )^0 = y * ( F_210 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[4] = etfac[1] * AUX_INT__f_s_s_s[1] + 1 * one_over_2q * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__g_s_s_s[3];

                    // ( F_210 S_000 | P_001 S_000 )^0 = z * ( F_210 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[5] = etfac[2] * AUX_INT__f_s_s_s[1] - p_over_q * AUX_INT__g_s_s_s[4];

                    // ( F_201 S_000 | P_100 S_000 )^0 = x * ( F_201 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_301 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[6] = etfac[0] * AUX_INT__f_s_s_s[2] + 2 * one_over_2q * AUX_INT__d_s_s_s[2] - p_over_q * AUX_INT__g_s_s_s[2];

                    // ( F_201 S_000 | P_010 S_000 )^0 = y * ( F_201 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[7] = etfac[1] * AUX_INT__f_s_s_s[2] - p_over_q * AUX_INT__g_s_s_s[4];

                    // ( F_201 S_000 | P_001 S_000 )^0 = z * ( F_201 S_000 | S_000 S_000 )^0_{e} + ( D_200 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[8] = etfac[2] * AUX_INT__f_s_s_s[2] + 1 * one_over_2q * AUX_INT__d_s_s_s[0] - p_over_q * AUX_INT__g_s_s_s[5];

                    // ( F_120 S_000 | P_100 S_000 )^0 = x * ( F_120 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_220 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[9] = etfac[0] * AUX_INT__f_s_s_s[3] + 1 * one_over_2q * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__g_s_s_s[3];

                    // ( F_120 S_000 | P_010 S_000 )^0 = y * ( F_120 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[10] = etfac[1] * AUX_INT__f_s_s_s[3] + 2 * one_over_2q * AUX_INT__d_s_s_s[1] - p_over_q * AUX_INT__g_s_s_s[6];

                    // ( F_120 S_000 | P_001 S_000 )^0 = z * ( F_120 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[11] = etfac[2] * AUX_INT__f_s_s_s[3] - p_over_q * AUX_INT__g_s_s_s[7];

                    // ( F_111 S_000 | P_100 S_000 )^0 = x * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_211 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[12] = etfac[0] * AUX_INT__f_s_s_s[4] + 1 * one_over_2q * AUX_INT__d_s_s_s[4] - p_over_q * AUX_INT__g_s_s_s[4];

                    // ( F_111 S_000 | P_010 S_000 )^0 = y * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[13] = etfac[1] * AUX_INT__f_s_s_s[4] + 1 * one_over_2q * AUX_INT__d_s_s_s[2] - p_over_q * AUX_INT__g_s_s_s[7];

                    // ( F_111 S_000 | P_001 S_000 )^0 = z * ( F_111 S_000 | S_000 S_000 )^0_{e} + ( D_110 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[14] = etfac[2] * AUX_INT__f_s_s_s[4] + 1 * one_over_2q * AUX_INT__d_s_s_s[1] - p_over_q * AUX_INT__g_s_s_s[8];

                    // ( F_102 S_000 | P_100 S_000 )^0 = x * ( F_102 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_202 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[15] = etfac[0] * AUX_INT__f_s_s_s[5] + 1 * one_over_2q * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__g_s_s_s[5];

                    // ( F_102 S_000 | P_010 S_000 )^0 = y * ( F_102 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[16] = etfac[1] * AUX_INT__f_s_s_s[5] - p_over_q * AUX_INT__g_s_s_s[8];

                    // ( F_102 S_000 | P_001 S_000 )^0 = z * ( F_102 S_000 | S_000 S_000 )^0_{e} + ( D_101 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[17] = etfac[2] * AUX_INT__f_s_s_s[5] + 2 * one_over_2q * AUX_INT__d_s_s_s[2] - p_over_q * AUX_INT__g_s_s_s[9];

                    // ( F_030 S_000 | P_100 S_000 )^0 = x * ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_130 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[18] = etfac[0] * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__g_s_s_s[6];

                    // ( F_030 S_000 | P_010 S_000 )^0 = y * ( F_030 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_040 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[19] = etfac[1] * AUX_INT__f_s_s_s[6] + 3 * one_over_2q * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__g_s_s_s[10];

                    // ( F_030 S_000 | P_001 S_000 )^0 = z * ( F_030 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[20] = etfac[2] * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__g_s_s_s[11];

                    // ( F_021 S_000 | P_100 S_000 )^0 = x * ( F_021 S_000 | S_000 S_000 )^0_{e} - ( G_121 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[21] = etfac[0] * AUX_INT__f_s_s_s[7] - p_over_q * AUX_INT__g_s_s_s[7];

                    // ( F_021 S_000 | P_010 S_000 )^0 = y * ( F_021 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_031 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[22] = etfac[1] * AUX_INT__f_s_s_s[7] + 2 * one_over_2q * AUX_INT__d_s_s_s[4] - p_over_q * AUX_INT__g_s_s_s[11];

                    // ( F_021 S_000 | P_001 S_000 )^0 = z * ( F_021 S_000 | S_000 S_000 )^0_{e} + ( D_020 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[23] = etfac[2] * AUX_INT__f_s_s_s[7] + 1 * one_over_2q * AUX_INT__d_s_s_s[3] - p_over_q * AUX_INT__g_s_s_s[12];

                    // ( F_012 S_000 | P_100 S_000 )^0 = x * ( F_012 S_000 | S_000 S_000 )^0_{e} - ( G_112 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[24] = etfac[0] * AUX_INT__f_s_s_s[8] - p_over_q * AUX_INT__g_s_s_s[8];

                    // ( F_012 S_000 | P_010 S_000 )^0 = y * ( F_012 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_022 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[25] = etfac[1] * AUX_INT__f_s_s_s[8] + 1 * one_over_2q * AUX_INT__d_s_s_s[5] - p_over_q * AUX_INT__g_s_s_s[12];

                    // ( F_012 S_000 | P_001 S_000 )^0 = z * ( F_012 S_000 | S_000 S_000 )^0_{e} + ( D_011 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[26] = etfac[2] * AUX_INT__f_s_s_s[8] + 2 * one_over_2q * AUX_INT__d_s_s_s[4] - p_over_q * AUX_INT__g_s_s_s[13];

                    // ( F_003 S_000 | P_100 S_000 )^0 = x * ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_103 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[27] = etfac[0] * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__g_s_s_s[9];

                    // ( F_003 S_000 | P_010 S_000 )^0 = y * ( F_003 S_000 | S_000 S_000 )^0_{e} - ( G_013 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__f_s_p_s[28] = etfac[1] * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__g_s_s_s[13];

                    // ( F_003 S_000 | P_001 S_000 )^0 = z * ( F_003 S_000 | S_000 S_000 )^0_{e} + ( D_002 S_000 | S_000 S_000 )^0_{e} - ( G_004 S_000 | S_000 S_000 )^0_{e}
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

                    // ( G_400 S_000 | P_100 S_000 )^0 = x * ( G_400 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_500 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[0] = etfac[0] * AUX_INT__g_s_s_s[0] + 4 * one_over_2q * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__h_s_s_s[0];

                    // ( G_400 S_000 | P_010 S_000 )^0 = y * ( G_400 S_000 | S_000 S_000 )^0_{e} - ( H_410 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[1] = etfac[1] * AUX_INT__g_s_s_s[0] - p_over_q * AUX_INT__h_s_s_s[1];

                    // ( G_400 S_000 | P_001 S_000 )^0 = z * ( G_400 S_000 | S_000 S_000 )^0_{e} - ( H_401 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[2] = etfac[2] * AUX_INT__g_s_s_s[0] - p_over_q * AUX_INT__h_s_s_s[2];

                    // ( G_310 S_000 | P_100 S_000 )^0 = x * ( G_310 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_410 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[3] = etfac[0] * AUX_INT__g_s_s_s[1] + 3 * one_over_2q * AUX_INT__f_s_s_s[1] - p_over_q * AUX_INT__h_s_s_s[1];

                    // ( G_310 S_000 | P_010 S_000 )^0 = y * ( G_310 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_320 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[4] = etfac[1] * AUX_INT__g_s_s_s[1] + 1 * one_over_2q * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__h_s_s_s[3];

                    // ( G_310 S_000 | P_001 S_000 )^0 = z * ( G_310 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[5] = etfac[2] * AUX_INT__g_s_s_s[1] - p_over_q * AUX_INT__h_s_s_s[4];

                    // ( G_301 S_000 | P_100 S_000 )^0 = x * ( G_301 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_401 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[6] = etfac[0] * AUX_INT__g_s_s_s[2] + 3 * one_over_2q * AUX_INT__f_s_s_s[2] - p_over_q * AUX_INT__h_s_s_s[2];

                    // ( G_301 S_000 | P_010 S_000 )^0 = y * ( G_301 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[7] = etfac[1] * AUX_INT__g_s_s_s[2] - p_over_q * AUX_INT__h_s_s_s[4];

                    // ( G_301 S_000 | P_001 S_000 )^0 = z * ( G_301 S_000 | S_000 S_000 )^0_{e} + ( F_300 S_000 | S_000 S_000 )^0_{e} - ( H_302 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[8] = etfac[2] * AUX_INT__g_s_s_s[2] + 1 * one_over_2q * AUX_INT__f_s_s_s[0] - p_over_q * AUX_INT__h_s_s_s[5];

                    // ( G_220 S_000 | P_100 S_000 )^0 = x * ( G_220 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_320 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[9] = etfac[0] * AUX_INT__g_s_s_s[3] + 2 * one_over_2q * AUX_INT__f_s_s_s[3] - p_over_q * AUX_INT__h_s_s_s[3];

                    // ( G_220 S_000 | P_010 S_000 )^0 = y * ( G_220 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_230 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[10] = etfac[1] * AUX_INT__g_s_s_s[3] + 2 * one_over_2q * AUX_INT__f_s_s_s[1] - p_over_q * AUX_INT__h_s_s_s[6];

                    // ( G_220 S_000 | P_001 S_000 )^0 = z * ( G_220 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[11] = etfac[2] * AUX_INT__g_s_s_s[3] - p_over_q * AUX_INT__h_s_s_s[7];

                    // ( G_211 S_000 | P_100 S_000 )^0 = x * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_311 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[12] = etfac[0] * AUX_INT__g_s_s_s[4] + 2 * one_over_2q * AUX_INT__f_s_s_s[4] - p_over_q * AUX_INT__h_s_s_s[4];

                    // ( G_211 S_000 | P_010 S_000 )^0 = y * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[13] = etfac[1] * AUX_INT__g_s_s_s[4] + 1 * one_over_2q * AUX_INT__f_s_s_s[2] - p_over_q * AUX_INT__h_s_s_s[7];

                    // ( G_211 S_000 | P_001 S_000 )^0 = z * ( G_211 S_000 | S_000 S_000 )^0_{e} + ( F_210 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[14] = etfac[2] * AUX_INT__g_s_s_s[4] + 1 * one_over_2q * AUX_INT__f_s_s_s[1] - p_over_q * AUX_INT__h_s_s_s[8];

                    // ( G_202 S_000 | P_100 S_000 )^0 = x * ( G_202 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_302 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[15] = etfac[0] * AUX_INT__g_s_s_s[5] + 2 * one_over_2q * AUX_INT__f_s_s_s[5] - p_over_q * AUX_INT__h_s_s_s[5];

                    // ( G_202 S_000 | P_010 S_000 )^0 = y * ( G_202 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[16] = etfac[1] * AUX_INT__g_s_s_s[5] - p_over_q * AUX_INT__h_s_s_s[8];

                    // ( G_202 S_000 | P_001 S_000 )^0 = z * ( G_202 S_000 | S_000 S_000 )^0_{e} + ( F_201 S_000 | S_000 S_000 )^0_{e} - ( H_203 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[17] = etfac[2] * AUX_INT__g_s_s_s[5] + 2 * one_over_2q * AUX_INT__f_s_s_s[2] - p_over_q * AUX_INT__h_s_s_s[9];

                    // ( G_130 S_000 | P_100 S_000 )^0 = x * ( G_130 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_230 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[18] = etfac[0] * AUX_INT__g_s_s_s[6] + 1 * one_over_2q * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__h_s_s_s[6];

                    // ( G_130 S_000 | P_010 S_000 )^0 = y * ( G_130 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_140 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[19] = etfac[1] * AUX_INT__g_s_s_s[6] + 3 * one_over_2q * AUX_INT__f_s_s_s[3] - p_over_q * AUX_INT__h_s_s_s[10];

                    // ( G_130 S_000 | P_001 S_000 )^0 = z * ( G_130 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[20] = etfac[2] * AUX_INT__g_s_s_s[6] - p_over_q * AUX_INT__h_s_s_s[11];

                    // ( G_121 S_000 | P_100 S_000 )^0 = x * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_221 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[21] = etfac[0] * AUX_INT__g_s_s_s[7] + 1 * one_over_2q * AUX_INT__f_s_s_s[7] - p_over_q * AUX_INT__h_s_s_s[7];

                    // ( G_121 S_000 | P_010 S_000 )^0 = y * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[22] = etfac[1] * AUX_INT__g_s_s_s[7] + 2 * one_over_2q * AUX_INT__f_s_s_s[4] - p_over_q * AUX_INT__h_s_s_s[11];

                    // ( G_121 S_000 | P_001 S_000 )^0 = z * ( G_121 S_000 | S_000 S_000 )^0_{e} + ( F_120 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[23] = etfac[2] * AUX_INT__g_s_s_s[7] + 1 * one_over_2q * AUX_INT__f_s_s_s[3] - p_over_q * AUX_INT__h_s_s_s[12];

                    // ( G_112 S_000 | P_100 S_000 )^0 = x * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_212 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[24] = etfac[0] * AUX_INT__g_s_s_s[8] + 1 * one_over_2q * AUX_INT__f_s_s_s[8] - p_over_q * AUX_INT__h_s_s_s[8];

                    // ( G_112 S_000 | P_010 S_000 )^0 = y * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[25] = etfac[1] * AUX_INT__g_s_s_s[8] + 1 * one_over_2q * AUX_INT__f_s_s_s[5] - p_over_q * AUX_INT__h_s_s_s[12];

                    // ( G_112 S_000 | P_001 S_000 )^0 = z * ( G_112 S_000 | S_000 S_000 )^0_{e} + ( F_111 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[26] = etfac[2] * AUX_INT__g_s_s_s[8] + 2 * one_over_2q * AUX_INT__f_s_s_s[4] - p_over_q * AUX_INT__h_s_s_s[13];

                    // ( G_103 S_000 | P_100 S_000 )^0 = x * ( G_103 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_203 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[27] = etfac[0] * AUX_INT__g_s_s_s[9] + 1 * one_over_2q * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__h_s_s_s[9];

                    // ( G_103 S_000 | P_010 S_000 )^0 = y * ( G_103 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[28] = etfac[1] * AUX_INT__g_s_s_s[9] - p_over_q * AUX_INT__h_s_s_s[13];

                    // ( G_103 S_000 | P_001 S_000 )^0 = z * ( G_103 S_000 | S_000 S_000 )^0_{e} + ( F_102 S_000 | S_000 S_000 )^0_{e} - ( H_104 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[29] = etfac[2] * AUX_INT__g_s_s_s[9] + 3 * one_over_2q * AUX_INT__f_s_s_s[5] - p_over_q * AUX_INT__h_s_s_s[14];

                    // ( G_040 S_000 | P_100 S_000 )^0 = x * ( G_040 S_000 | S_000 S_000 )^0_{e} - ( H_140 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[30] = etfac[0] * AUX_INT__g_s_s_s[10] - p_over_q * AUX_INT__h_s_s_s[10];

                    // ( G_040 S_000 | P_010 S_000 )^0 = y * ( G_040 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_050 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[31] = etfac[1] * AUX_INT__g_s_s_s[10] + 4 * one_over_2q * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__h_s_s_s[15];

                    // ( G_040 S_000 | P_001 S_000 )^0 = z * ( G_040 S_000 | S_000 S_000 )^0_{e} - ( H_041 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[32] = etfac[2] * AUX_INT__g_s_s_s[10] - p_over_q * AUX_INT__h_s_s_s[16];

                    // ( G_031 S_000 | P_100 S_000 )^0 = x * ( G_031 S_000 | S_000 S_000 )^0_{e} - ( H_131 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[33] = etfac[0] * AUX_INT__g_s_s_s[11] - p_over_q * AUX_INT__h_s_s_s[11];

                    // ( G_031 S_000 | P_010 S_000 )^0 = y * ( G_031 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_041 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[34] = etfac[1] * AUX_INT__g_s_s_s[11] + 3 * one_over_2q * AUX_INT__f_s_s_s[7] - p_over_q * AUX_INT__h_s_s_s[16];

                    // ( G_031 S_000 | P_001 S_000 )^0 = z * ( G_031 S_000 | S_000 S_000 )^0_{e} + ( F_030 S_000 | S_000 S_000 )^0_{e} - ( H_032 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[35] = etfac[2] * AUX_INT__g_s_s_s[11] + 1 * one_over_2q * AUX_INT__f_s_s_s[6] - p_over_q * AUX_INT__h_s_s_s[17];

                    // ( G_022 S_000 | P_100 S_000 )^0 = x * ( G_022 S_000 | S_000 S_000 )^0_{e} - ( H_122 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[36] = etfac[0] * AUX_INT__g_s_s_s[12] - p_over_q * AUX_INT__h_s_s_s[12];

                    // ( G_022 S_000 | P_010 S_000 )^0 = y * ( G_022 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_032 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[37] = etfac[1] * AUX_INT__g_s_s_s[12] + 2 * one_over_2q * AUX_INT__f_s_s_s[8] - p_over_q * AUX_INT__h_s_s_s[17];

                    // ( G_022 S_000 | P_001 S_000 )^0 = z * ( G_022 S_000 | S_000 S_000 )^0_{e} + ( F_021 S_000 | S_000 S_000 )^0_{e} - ( H_023 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[38] = etfac[2] * AUX_INT__g_s_s_s[12] + 2 * one_over_2q * AUX_INT__f_s_s_s[7] - p_over_q * AUX_INT__h_s_s_s[18];

                    // ( G_013 S_000 | P_100 S_000 )^0 = x * ( G_013 S_000 | S_000 S_000 )^0_{e} - ( H_113 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[39] = etfac[0] * AUX_INT__g_s_s_s[13] - p_over_q * AUX_INT__h_s_s_s[13];

                    // ( G_013 S_000 | P_010 S_000 )^0 = y * ( G_013 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_023 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[40] = etfac[1] * AUX_INT__g_s_s_s[13] + 1 * one_over_2q * AUX_INT__f_s_s_s[9] - p_over_q * AUX_INT__h_s_s_s[18];

                    // ( G_013 S_000 | P_001 S_000 )^0 = z * ( G_013 S_000 | S_000 S_000 )^0_{e} + ( F_012 S_000 | S_000 S_000 )^0_{e} - ( H_014 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[41] = etfac[2] * AUX_INT__g_s_s_s[13] + 3 * one_over_2q * AUX_INT__f_s_s_s[8] - p_over_q * AUX_INT__h_s_s_s[19];

                    // ( G_004 S_000 | P_100 S_000 )^0 = x * ( G_004 S_000 | S_000 S_000 )^0_{e} - ( H_104 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[42] = etfac[0] * AUX_INT__g_s_s_s[14] - p_over_q * AUX_INT__h_s_s_s[14];

                    // ( G_004 S_000 | P_010 S_000 )^0 = y * ( G_004 S_000 | S_000 S_000 )^0_{e} - ( H_014 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__g_s_p_s[43] = etfac[1] * AUX_INT__g_s_s_s[14] - p_over_q * AUX_INT__h_s_s_s[19];

                    // ( G_004 S_000 | P_001 S_000 )^0 = z * ( G_004 S_000 | S_000 S_000 )^0_{e} + ( F_003 S_000 | S_000 S_000 )^0_{e} - ( H_005 S_000 | S_000 S_000 )^0_{e}
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

                    // ( S_000 S_000 | P_100 S_000 )^0 = x * ( S_000 S_000 | S_000 S_000 )^0_{e} - ( P_100 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__s_s_p_s[0] = etfac[0] * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__p_s_s_s[0];

                    // ( S_000 S_000 | P_010 S_000 )^0 = y * ( S_000 S_000 | S_000 S_000 )^0_{e} - ( P_010 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__s_s_p_s[1] = etfac[1] * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__p_s_s_s[1];

                    // ( S_000 S_000 | P_001 S_000 )^0 = z * ( S_000 S_000 | S_000 S_000 )^0_{e} - ( P_001 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__s_s_p_s[2] = etfac[2] * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__p_s_s_s[2];

                    // ( P_100 S_000 | D_200 S_000 )^0 = x * ( P_100 S_000 | P_100 S_000 )^0 + ( S_000 S_000 | P_100 S_000 )^0 + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( D_200 S_000 | P_100 S_000 )^0
                    AUX_INT__p_s_d_s[0] = etfac[0] * AUX_INT__p_s_p_s[0] + 1 * one_over_2q * AUX_INT__s_s_p_s[0] + 1 * one_over_2q * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__d_s_p_s[0];

                    // ( P_100 S_000 | D_110 S_000 )^0 = y * ( P_100 S_000 | P_100 S_000 )^0 - ( D_110 S_000 | P_100 S_000 )^0
                    AUX_INT__p_s_d_s[1] = etfac[1] * AUX_INT__p_s_p_s[0] - p_over_q * AUX_INT__d_s_p_s[3];

                    // ( P_100 S_000 | D_020 S_000 )^0 = y * ( P_100 S_000 | P_010 S_000 )^0 + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( D_110 S_000 | P_010 S_000 )^0
                    AUX_INT__p_s_d_s[3] = etfac[1] * AUX_INT__p_s_p_s[1] + 1 * one_over_2q * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__d_s_p_s[4];

                    // ( P_100 S_000 | D_002 S_000 )^0 = z * ( P_100 S_000 | P_001 S_000 )^0 + ( P_100 S_000 | S_000 S_000 )^0_{e} - ( D_101 S_000 | P_001 S_000 )^0
                    AUX_INT__p_s_d_s[5] = etfac[2] * AUX_INT__p_s_p_s[2] + 1 * one_over_2q * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__d_s_p_s[8];

                    // ( P_010 S_000 | D_200 S_000 )^0 = x * ( P_010 S_000 | P_100 S_000 )^0 + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( D_110 S_000 | P_100 S_000 )^0
                    AUX_INT__p_s_d_s[6] = etfac[0] * AUX_INT__p_s_p_s[3] + 1 * one_over_2q * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__d_s_p_s[3];

                    // ( P_010 S_000 | D_110 S_000 )^0 = y * ( P_010 S_000 | P_100 S_000 )^0 + ( S_000 S_000 | P_100 S_000 )^0 - ( D_020 S_000 | P_100 S_000 )^0
                    AUX_INT__p_s_d_s[7] = etfac[1] * AUX_INT__p_s_p_s[3] + 1 * one_over_2q * AUX_INT__s_s_p_s[0] - p_over_q * AUX_INT__d_s_p_s[9];

                    // ( P_010 S_000 | D_020 S_000 )^0 = y * ( P_010 S_000 | P_010 S_000 )^0 + ( S_000 S_000 | P_010 S_000 )^0 + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( D_020 S_000 | P_010 S_000 )^0
                    AUX_INT__p_s_d_s[9] = etfac[1] * AUX_INT__p_s_p_s[4] + 1 * one_over_2q * AUX_INT__s_s_p_s[1] + 1 * one_over_2q * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__d_s_p_s[10];

                    // ( P_010 S_000 | D_002 S_000 )^0 = z * ( P_010 S_000 | P_001 S_000 )^0 + ( P_010 S_000 | S_000 S_000 )^0_{e} - ( D_011 S_000 | P_001 S_000 )^0
                    AUX_INT__p_s_d_s[11] = etfac[2] * AUX_INT__p_s_p_s[5] + 1 * one_over_2q * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__d_s_p_s[14];

                    // ( P_001 S_000 | D_200 S_000 )^0 = x * ( P_001 S_000 | P_100 S_000 )^0 + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( D_101 S_000 | P_100 S_000 )^0
                    AUX_INT__p_s_d_s[12] = etfac[0] * AUX_INT__p_s_p_s[6] + 1 * one_over_2q * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__d_s_p_s[6];

                    // ( P_001 S_000 | D_110 S_000 )^0 = y * ( P_001 S_000 | P_100 S_000 )^0 - ( D_011 S_000 | P_100 S_000 )^0
                    AUX_INT__p_s_d_s[13] = etfac[1] * AUX_INT__p_s_p_s[6] - p_over_q * AUX_INT__d_s_p_s[12];

                    // ( P_001 S_000 | D_020 S_000 )^0 = y * ( P_001 S_000 | P_010 S_000 )^0 + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( D_011 S_000 | P_010 S_000 )^0
                    AUX_INT__p_s_d_s[15] = etfac[1] * AUX_INT__p_s_p_s[7] + 1 * one_over_2q * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__d_s_p_s[13];

                    // ( P_001 S_000 | D_002 S_000 )^0 = z * ( P_001 S_000 | P_001 S_000 )^0 + ( S_000 S_000 | P_001 S_000 )^0 + ( P_001 S_000 | S_000 S_000 )^0_{e} - ( D_002 S_000 | P_001 S_000 )^0
                    AUX_INT__p_s_d_s[17] = etfac[2] * AUX_INT__p_s_p_s[8] + 1 * one_over_2q * AUX_INT__s_s_p_s[2] + 1 * one_over_2q * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__d_s_p_s[17];

                    // ( D_200 S_000 | F_300 S_000 )^0_{t} = x * ( D_200 S_000 | D_200 S_000 )^0 + ( P_100 S_000 | D_200 S_000 )^0 + ( D_200 S_000 | P_100 S_000 )^0 - ( F_300 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[0] = etfac[0] * AUX_INT__d_s_d_s[0] + 2 * one_over_2q * AUX_INT__p_s_d_s[0] + 2 * one_over_2q * AUX_INT__d_s_p_s[0] - p_over_q * AUX_INT__f_s_d_s[0];

                    // ( D_200 S_000 | F_210 S_000 )^0_{t} = y * ( D_200 S_000 | D_200 S_000 )^0 - ( F_210 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[1] = etfac[1] * AUX_INT__d_s_d_s[0] - p_over_q * AUX_INT__f_s_d_s[6];

                    // ( D_200 S_000 | F_201 S_000 )^0_{t} = z * ( D_200 S_000 | D_200 S_000 )^0 - ( F_201 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[2] = etfac[2] * AUX_INT__d_s_d_s[0] - p_over_q * AUX_INT__f_s_d_s[12];

                    // ( D_200 S_000 | F_120 S_000 )^0_{t} = x * ( D_200 S_000 | D_020 S_000 )^0 + ( P_100 S_000 | D_020 S_000 )^0 - ( F_300 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[3] = etfac[0] * AUX_INT__d_s_d_s[3] + 2 * one_over_2q * AUX_INT__p_s_d_s[3] - p_over_q * AUX_INT__f_s_d_s[3];

                    // ( D_200 S_000 | F_111 S_000 )^0_{t} = z * ( D_200 S_000 | D_110 S_000 )^0 - ( F_201 S_000 | D_110 S_000 )^0
                    AUX_INT__d_s_f_s[4] = etfac[2] * AUX_INT__d_s_d_s[1] - p_over_q * AUX_INT__f_s_d_s[13];

                    // ( D_200 S_000 | F_102 S_000 )^0_{t} = x * ( D_200 S_000 | D_002 S_000 )^0 + ( P_100 S_000 | D_002 S_000 )^0 - ( F_300 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[5] = etfac[0] * AUX_INT__d_s_d_s[5] + 2 * one_over_2q * AUX_INT__p_s_d_s[5] - p_over_q * AUX_INT__f_s_d_s[5];

                    // ( D_200 S_000 | F_030 S_000 )^0_{t} = y * ( D_200 S_000 | D_020 S_000 )^0 + ( D_200 S_000 | P_010 S_000 )^0 - ( F_210 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[6] = etfac[1] * AUX_INT__d_s_d_s[3] + 2 * one_over_2q * AUX_INT__d_s_p_s[1] - p_over_q * AUX_INT__f_s_d_s[9];

                    // ( D_200 S_000 | F_021 S_000 )^0_{t} = z * ( D_200 S_000 | D_020 S_000 )^0 - ( F_201 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[7] = etfac[2] * AUX_INT__d_s_d_s[3] - p_over_q * AUX_INT__f_s_d_s[15];

                    // ( D_200 S_000 | F_012 S_000 )^0_{t} = y * ( D_200 S_000 | D_002 S_000 )^0 - ( F_210 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[8] = etfac[1] * AUX_INT__d_s_d_s[5] - p_over_q * AUX_INT__f_s_d_s[11];

                    // ( D_200 S_000 | F_003 S_000 )^0_{t} = z * ( D_200 S_000 | D_002 S_000 )^0 + ( D_200 S_000 | P_001 S_000 )^0 - ( F_201 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[9] = etfac[2] * AUX_INT__d_s_d_s[5] + 2 * one_over_2q * AUX_INT__d_s_p_s[2] - p_over_q * AUX_INT__f_s_d_s[17];

                    // ( D_110 S_000 | F_300 S_000 )^0_{t} = x * ( D_110 S_000 | D_200 S_000 )^0 + ( P_010 S_000 | D_200 S_000 )^0 + ( D_110 S_000 | P_100 S_000 )^0 - ( F_210 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[10] = etfac[0] * AUX_INT__d_s_d_s[6] + 1 * one_over_2q * AUX_INT__p_s_d_s[6] + 2 * one_over_2q * AUX_INT__d_s_p_s[3] - p_over_q * AUX_INT__f_s_d_s[6];

                    // ( D_110 S_000 | F_210 S_000 )^0_{t} = y * ( D_110 S_000 | D_200 S_000 )^0 + ( P_100 S_000 | D_200 S_000 )^0 - ( F_120 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[11] = etfac[1] * AUX_INT__d_s_d_s[6] + 1 * one_over_2q * AUX_INT__p_s_d_s[0] - p_over_q * AUX_INT__f_s_d_s[18];

                    // ( D_110 S_000 | F_201 S_000 )^0_{t} = z * ( D_110 S_000 | D_200 S_000 )^0 - ( F_111 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[12] = etfac[2] * AUX_INT__d_s_d_s[6] - p_over_q * AUX_INT__f_s_d_s[24];

                    // ( D_110 S_000 | F_120 S_000 )^0_{t} = x * ( D_110 S_000 | D_020 S_000 )^0 + ( P_010 S_000 | D_020 S_000 )^0 - ( F_210 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[13] = etfac[0] * AUX_INT__d_s_d_s[9] + 1 * one_over_2q * AUX_INT__p_s_d_s[9] - p_over_q * AUX_INT__f_s_d_s[9];

                    // ( D_110 S_000 | F_111 S_000 )^0_{t} = z * ( D_110 S_000 | D_110 S_000 )^0 - ( F_111 S_000 | D_110 S_000 )^0
                    AUX_INT__d_s_f_s[14] = etfac[2] * AUX_INT__d_s_d_s[7] - p_over_q * AUX_INT__f_s_d_s[25];

                    // ( D_110 S_000 | F_102 S_000 )^0_{t} = x * ( D_110 S_000 | D_002 S_000 )^0 + ( P_010 S_000 | D_002 S_000 )^0 - ( F_210 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[15] = etfac[0] * AUX_INT__d_s_d_s[11] + 1 * one_over_2q * AUX_INT__p_s_d_s[11] - p_over_q * AUX_INT__f_s_d_s[11];

                    // ( D_110 S_000 | F_030 S_000 )^0_{t} = y * ( D_110 S_000 | D_020 S_000 )^0 + ( P_100 S_000 | D_020 S_000 )^0 + ( D_110 S_000 | P_010 S_000 )^0 - ( F_120 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[16] = etfac[1] * AUX_INT__d_s_d_s[9] + 1 * one_over_2q * AUX_INT__p_s_d_s[3] + 2 * one_over_2q * AUX_INT__d_s_p_s[4] - p_over_q * AUX_INT__f_s_d_s[21];

                    // ( D_110 S_000 | F_021 S_000 )^0_{t} = z * ( D_110 S_000 | D_020 S_000 )^0 - ( F_111 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[17] = etfac[2] * AUX_INT__d_s_d_s[9] - p_over_q * AUX_INT__f_s_d_s[27];

                    // ( D_110 S_000 | F_012 S_000 )^0_{t} = y * ( D_110 S_000 | D_002 S_000 )^0 + ( P_100 S_000 | D_002 S_000 )^0 - ( F_120 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[18] = etfac[1] * AUX_INT__d_s_d_s[11] + 1 * one_over_2q * AUX_INT__p_s_d_s[5] - p_over_q * AUX_INT__f_s_d_s[23];

                    // ( D_110 S_000 | F_003 S_000 )^0_{t} = z * ( D_110 S_000 | D_002 S_000 )^0 + ( D_110 S_000 | P_001 S_000 )^0 - ( F_111 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[19] = etfac[2] * AUX_INT__d_s_d_s[11] + 2 * one_over_2q * AUX_INT__d_s_p_s[5] - p_over_q * AUX_INT__f_s_d_s[29];

                    // ( D_101 S_000 | F_300 S_000 )^0_{t} = x * ( D_101 S_000 | D_200 S_000 )^0 + ( P_001 S_000 | D_200 S_000 )^0 + ( D_101 S_000 | P_100 S_000 )^0 - ( F_201 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[20] = etfac[0] * AUX_INT__d_s_d_s[12] + 1 * one_over_2q * AUX_INT__p_s_d_s[12] + 2 * one_over_2q * AUX_INT__d_s_p_s[6] - p_over_q * AUX_INT__f_s_d_s[12];

                    // ( D_101 S_000 | F_210 S_000 )^0_{t} = y * ( D_101 S_000 | D_200 S_000 )^0 - ( F_111 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[21] = etfac[1] * AUX_INT__d_s_d_s[12] - p_over_q * AUX_INT__f_s_d_s[24];

                    // ( D_101 S_000 | F_201 S_000 )^0_{t} = z * ( D_101 S_000 | D_200 S_000 )^0 + ( P_100 S_000 | D_200 S_000 )^0 - ( F_102 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[22] = etfac[2] * AUX_INT__d_s_d_s[12] + 1 * one_over_2q * AUX_INT__p_s_d_s[0] - p_over_q * AUX_INT__f_s_d_s[30];

                    // ( D_101 S_000 | F_120 S_000 )^0_{t} = x * ( D_101 S_000 | D_020 S_000 )^0 + ( P_001 S_000 | D_020 S_000 )^0 - ( F_201 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[23] = etfac[0] * AUX_INT__d_s_d_s[15] + 1 * one_over_2q * AUX_INT__p_s_d_s[15] - p_over_q * AUX_INT__f_s_d_s[15];

                    // ( D_101 S_000 | F_111 S_000 )^0_{t} = z * ( D_101 S_000 | D_110 S_000 )^0 + ( P_100 S_000 | D_110 S_000 )^0 - ( F_102 S_000 | D_110 S_000 )^0
                    AUX_INT__d_s_f_s[24] = etfac[2] * AUX_INT__d_s_d_s[13] + 1 * one_over_2q * AUX_INT__p_s_d_s[1] - p_over_q * AUX_INT__f_s_d_s[31];

                    // ( D_101 S_000 | F_102 S_000 )^0_{t} = x * ( D_101 S_000 | D_002 S_000 )^0 + ( P_001 S_000 | D_002 S_000 )^0 - ( F_201 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[25] = etfac[0] * AUX_INT__d_s_d_s[17] + 1 * one_over_2q * AUX_INT__p_s_d_s[17] - p_over_q * AUX_INT__f_s_d_s[17];

                    // ( D_101 S_000 | F_030 S_000 )^0_{t} = y * ( D_101 S_000 | D_020 S_000 )^0 + ( D_101 S_000 | P_010 S_000 )^0 - ( F_111 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[26] = etfac[1] * AUX_INT__d_s_d_s[15] + 2 * one_over_2q * AUX_INT__d_s_p_s[7] - p_over_q * AUX_INT__f_s_d_s[27];

                    // ( D_101 S_000 | F_021 S_000 )^0_{t} = z * ( D_101 S_000 | D_020 S_000 )^0 + ( P_100 S_000 | D_020 S_000 )^0 - ( F_102 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[27] = etfac[2] * AUX_INT__d_s_d_s[15] + 1 * one_over_2q * AUX_INT__p_s_d_s[3] - p_over_q * AUX_INT__f_s_d_s[33];

                    // ( D_101 S_000 | F_012 S_000 )^0_{t} = y * ( D_101 S_000 | D_002 S_000 )^0 - ( F_111 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[28] = etfac[1] * AUX_INT__d_s_d_s[17] - p_over_q * AUX_INT__f_s_d_s[29];

                    // ( D_101 S_000 | F_003 S_000 )^0_{t} = z * ( D_101 S_000 | D_002 S_000 )^0 + ( P_100 S_000 | D_002 S_000 )^0 + ( D_101 S_000 | P_001 S_000 )^0 - ( F_102 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[29] = etfac[2] * AUX_INT__d_s_d_s[17] + 1 * one_over_2q * AUX_INT__p_s_d_s[5] + 2 * one_over_2q * AUX_INT__d_s_p_s[8] - p_over_q * AUX_INT__f_s_d_s[35];

                    // ( D_020 S_000 | F_300 S_000 )^0_{t} = x * ( D_020 S_000 | D_200 S_000 )^0 + ( D_020 S_000 | P_100 S_000 )^0 - ( F_120 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[30] = etfac[0] * AUX_INT__d_s_d_s[18] + 2 * one_over_2q * AUX_INT__d_s_p_s[9] - p_over_q * AUX_INT__f_s_d_s[18];

                    // ( D_020 S_000 | F_210 S_000 )^0_{t} = y * ( D_020 S_000 | D_200 S_000 )^0 + ( P_010 S_000 | D_200 S_000 )^0 - ( F_030 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[31] = etfac[1] * AUX_INT__d_s_d_s[18] + 2 * one_over_2q * AUX_INT__p_s_d_s[6] - p_over_q * AUX_INT__f_s_d_s[36];

                    // ( D_020 S_000 | F_201 S_000 )^0_{t} = z * ( D_020 S_000 | D_200 S_000 )^0 - ( F_021 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[32] = etfac[2] * AUX_INT__d_s_d_s[18] - p_over_q * AUX_INT__f_s_d_s[42];

                    // ( D_020 S_000 | F_120 S_000 )^0_{t} = x * ( D_020 S_000 | D_020 S_000 )^0 - ( F_120 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[33] = etfac[0] * AUX_INT__d_s_d_s[21] - p_over_q * AUX_INT__f_s_d_s[21];

                    // ( D_020 S_000 | F_111 S_000 )^0_{t} = z * ( D_020 S_000 | D_110 S_000 )^0 - ( F_021 S_000 | D_110 S_000 )^0
                    AUX_INT__d_s_f_s[34] = etfac[2] * AUX_INT__d_s_d_s[19] - p_over_q * AUX_INT__f_s_d_s[43];

                    // ( D_020 S_000 | F_102 S_000 )^0_{t} = x * ( D_020 S_000 | D_002 S_000 )^0 - ( F_120 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[35] = etfac[0] * AUX_INT__d_s_d_s[23] - p_over_q * AUX_INT__f_s_d_s[23];

                    // ( D_020 S_000 | F_030 S_000 )^0_{t} = y * ( D_020 S_000 | D_020 S_000 )^0 + ( P_010 S_000 | D_020 S_000 )^0 + ( D_020 S_000 | P_010 S_000 )^0 - ( F_030 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[36] = etfac[1] * AUX_INT__d_s_d_s[21] + 2 * one_over_2q * AUX_INT__p_s_d_s[9] + 2 * one_over_2q * AUX_INT__d_s_p_s[10] - p_over_q * AUX_INT__f_s_d_s[39];

                    // ( D_020 S_000 | F_021 S_000 )^0_{t} = z * ( D_020 S_000 | D_020 S_000 )^0 - ( F_021 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[37] = etfac[2] * AUX_INT__d_s_d_s[21] - p_over_q * AUX_INT__f_s_d_s[45];

                    // ( D_020 S_000 | F_012 S_000 )^0_{t} = y * ( D_020 S_000 | D_002 S_000 )^0 + ( P_010 S_000 | D_002 S_000 )^0 - ( F_030 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[38] = etfac[1] * AUX_INT__d_s_d_s[23] + 2 * one_over_2q * AUX_INT__p_s_d_s[11] - p_over_q * AUX_INT__f_s_d_s[41];

                    // ( D_020 S_000 | F_003 S_000 )^0_{t} = z * ( D_020 S_000 | D_002 S_000 )^0 + ( D_020 S_000 | P_001 S_000 )^0 - ( F_021 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[39] = etfac[2] * AUX_INT__d_s_d_s[23] + 2 * one_over_2q * AUX_INT__d_s_p_s[11] - p_over_q * AUX_INT__f_s_d_s[47];

                    // ( D_011 S_000 | F_300 S_000 )^0_{t} = x * ( D_011 S_000 | D_200 S_000 )^0 + ( D_011 S_000 | P_100 S_000 )^0 - ( F_111 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[40] = etfac[0] * AUX_INT__d_s_d_s[24] + 2 * one_over_2q * AUX_INT__d_s_p_s[12] - p_over_q * AUX_INT__f_s_d_s[24];

                    // ( D_011 S_000 | F_210 S_000 )^0_{t} = y * ( D_011 S_000 | D_200 S_000 )^0 + ( P_001 S_000 | D_200 S_000 )^0 - ( F_021 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[41] = etfac[1] * AUX_INT__d_s_d_s[24] + 1 * one_over_2q * AUX_INT__p_s_d_s[12] - p_over_q * AUX_INT__f_s_d_s[42];

                    // ( D_011 S_000 | F_201 S_000 )^0_{t} = z * ( D_011 S_000 | D_200 S_000 )^0 + ( P_010 S_000 | D_200 S_000 )^0 - ( F_012 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[42] = etfac[2] * AUX_INT__d_s_d_s[24] + 1 * one_over_2q * AUX_INT__p_s_d_s[6] - p_over_q * AUX_INT__f_s_d_s[48];

                    // ( D_011 S_000 | F_120 S_000 )^0_{t} = x * ( D_011 S_000 | D_020 S_000 )^0 - ( F_111 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[43] = etfac[0] * AUX_INT__d_s_d_s[27] - p_over_q * AUX_INT__f_s_d_s[27];

                    // ( D_011 S_000 | F_111 S_000 )^0_{t} = z * ( D_011 S_000 | D_110 S_000 )^0 + ( P_010 S_000 | D_110 S_000 )^0 - ( F_012 S_000 | D_110 S_000 )^0
                    AUX_INT__d_s_f_s[44] = etfac[2] * AUX_INT__d_s_d_s[25] + 1 * one_over_2q * AUX_INT__p_s_d_s[7] - p_over_q * AUX_INT__f_s_d_s[49];

                    // ( D_011 S_000 | F_102 S_000 )^0_{t} = x * ( D_011 S_000 | D_002 S_000 )^0 - ( F_111 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[45] = etfac[0] * AUX_INT__d_s_d_s[29] - p_over_q * AUX_INT__f_s_d_s[29];

                    // ( D_011 S_000 | F_030 S_000 )^0_{t} = y * ( D_011 S_000 | D_020 S_000 )^0 + ( P_001 S_000 | D_020 S_000 )^0 + ( D_011 S_000 | P_010 S_000 )^0 - ( F_021 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[46] = etfac[1] * AUX_INT__d_s_d_s[27] + 1 * one_over_2q * AUX_INT__p_s_d_s[15] + 2 * one_over_2q * AUX_INT__d_s_p_s[13] - p_over_q * AUX_INT__f_s_d_s[45];

                    // ( D_011 S_000 | F_021 S_000 )^0_{t} = z * ( D_011 S_000 | D_020 S_000 )^0 + ( P_010 S_000 | D_020 S_000 )^0 - ( F_012 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[47] = etfac[2] * AUX_INT__d_s_d_s[27] + 1 * one_over_2q * AUX_INT__p_s_d_s[9] - p_over_q * AUX_INT__f_s_d_s[51];

                    // ( D_011 S_000 | F_012 S_000 )^0_{t} = y * ( D_011 S_000 | D_002 S_000 )^0 + ( P_001 S_000 | D_002 S_000 )^0 - ( F_021 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[48] = etfac[1] * AUX_INT__d_s_d_s[29] + 1 * one_over_2q * AUX_INT__p_s_d_s[17] - p_over_q * AUX_INT__f_s_d_s[47];

                    // ( D_011 S_000 | F_003 S_000 )^0_{t} = z * ( D_011 S_000 | D_002 S_000 )^0 + ( P_010 S_000 | D_002 S_000 )^0 + ( D_011 S_000 | P_001 S_000 )^0 - ( F_012 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[49] = etfac[2] * AUX_INT__d_s_d_s[29] + 1 * one_over_2q * AUX_INT__p_s_d_s[11] + 2 * one_over_2q * AUX_INT__d_s_p_s[14] - p_over_q * AUX_INT__f_s_d_s[53];

                    // ( D_002 S_000 | F_300 S_000 )^0_{t} = x * ( D_002 S_000 | D_200 S_000 )^0 + ( D_002 S_000 | P_100 S_000 )^0 - ( F_102 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[50] = etfac[0] * AUX_INT__d_s_d_s[30] + 2 * one_over_2q * AUX_INT__d_s_p_s[15] - p_over_q * AUX_INT__f_s_d_s[30];

                    // ( D_002 S_000 | F_210 S_000 )^0_{t} = y * ( D_002 S_000 | D_200 S_000 )^0 - ( F_012 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[51] = etfac[1] * AUX_INT__d_s_d_s[30] - p_over_q * AUX_INT__f_s_d_s[48];

                    // ( D_002 S_000 | F_201 S_000 )^0_{t} = z * ( D_002 S_000 | D_200 S_000 )^0 + ( P_001 S_000 | D_200 S_000 )^0 - ( F_003 S_000 | D_200 S_000 )^0
                    AUX_INT__d_s_f_s[52] = etfac[2] * AUX_INT__d_s_d_s[30] + 2 * one_over_2q * AUX_INT__p_s_d_s[12] - p_over_q * AUX_INT__f_s_d_s[54];

                    // ( D_002 S_000 | F_120 S_000 )^0_{t} = x * ( D_002 S_000 | D_020 S_000 )^0 - ( F_102 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[53] = etfac[0] * AUX_INT__d_s_d_s[33] - p_over_q * AUX_INT__f_s_d_s[33];

                    // ( D_002 S_000 | F_111 S_000 )^0_{t} = z * ( D_002 S_000 | D_110 S_000 )^0 + ( P_001 S_000 | D_110 S_000 )^0 - ( F_003 S_000 | D_110 S_000 )^0
                    AUX_INT__d_s_f_s[54] = etfac[2] * AUX_INT__d_s_d_s[31] + 2 * one_over_2q * AUX_INT__p_s_d_s[13] - p_over_q * AUX_INT__f_s_d_s[55];

                    // ( D_002 S_000 | F_102 S_000 )^0_{t} = x * ( D_002 S_000 | D_002 S_000 )^0 - ( F_102 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[55] = etfac[0] * AUX_INT__d_s_d_s[35] - p_over_q * AUX_INT__f_s_d_s[35];

                    // ( D_002 S_000 | F_030 S_000 )^0_{t} = y * ( D_002 S_000 | D_020 S_000 )^0 + ( D_002 S_000 | P_010 S_000 )^0 - ( F_012 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[56] = etfac[1] * AUX_INT__d_s_d_s[33] + 2 * one_over_2q * AUX_INT__d_s_p_s[16] - p_over_q * AUX_INT__f_s_d_s[51];

                    // ( D_002 S_000 | F_021 S_000 )^0_{t} = z * ( D_002 S_000 | D_020 S_000 )^0 + ( P_001 S_000 | D_020 S_000 )^0 - ( F_003 S_000 | D_020 S_000 )^0
                    AUX_INT__d_s_f_s[57] = etfac[2] * AUX_INT__d_s_d_s[33] + 2 * one_over_2q * AUX_INT__p_s_d_s[15] - p_over_q * AUX_INT__f_s_d_s[57];

                    // ( D_002 S_000 | F_012 S_000 )^0_{t} = y * ( D_002 S_000 | D_002 S_000 )^0 - ( F_012 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[58] = etfac[1] * AUX_INT__d_s_d_s[35] - p_over_q * AUX_INT__f_s_d_s[53];

                    // ( D_002 S_000 | F_003 S_000 )^0_{t} = z * ( D_002 S_000 | D_002 S_000 )^0 + ( P_001 S_000 | D_002 S_000 )^0 + ( D_002 S_000 | P_001 S_000 )^0 - ( F_003 S_000 | D_002 S_000 )^0
                    AUX_INT__d_s_f_s[59] = etfac[2] * AUX_INT__d_s_d_s[35] + 2 * one_over_2q * AUX_INT__p_s_d_s[17] + 2 * one_over_2q * AUX_INT__d_s_p_s[17] - p_over_q * AUX_INT__f_s_d_s[59];

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

                    // ( H_500 S_000 | P_010 S_000 )^0 = y * ( H_500 S_000 | S_000 S_000 )^0_{e} - ( I_510 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[1] = etfac[1] * AUX_INT__h_s_s_s[0] - p_over_q * AUX_INT__i_s_s_s[1];

                    // ( H_500 S_000 | P_001 S_000 )^0 = z * ( H_500 S_000 | S_000 S_000 )^0_{e} - ( I_501 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[2] = etfac[2] * AUX_INT__h_s_s_s[0] - p_over_q * AUX_INT__i_s_s_s[2];

                    // ( H_410 S_000 | P_001 S_000 )^0 = z * ( H_410 S_000 | S_000 S_000 )^0_{e} - ( I_411 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[5] = etfac[2] * AUX_INT__h_s_s_s[1] - p_over_q * AUX_INT__i_s_s_s[4];

                    // ( H_320 S_000 | P_001 S_000 )^0 = z * ( H_320 S_000 | S_000 S_000 )^0_{e} - ( I_321 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[11] = etfac[2] * AUX_INT__h_s_s_s[3] - p_over_q * AUX_INT__i_s_s_s[7];

                    // ( H_230 S_000 | P_001 S_000 )^0 = z * ( H_230 S_000 | S_000 S_000 )^0_{e} - ( I_231 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[20] = etfac[2] * AUX_INT__h_s_s_s[6] - p_over_q * AUX_INT__i_s_s_s[11];

                    // ( H_140 S_000 | P_001 S_000 )^0 = z * ( H_140 S_000 | S_000 S_000 )^0_{e} - ( I_141 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[32] = etfac[2] * AUX_INT__h_s_s_s[10] - p_over_q * AUX_INT__i_s_s_s[16];

                    // ( H_050 S_000 | P_001 S_000 )^0 = z * ( H_050 S_000 | S_000 S_000 )^0_{e} - ( I_051 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__h_s_p_s[47] = etfac[2] * AUX_INT__h_s_s_s[15] - p_over_q * AUX_INT__i_s_s_s[22];

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

                    // ( F_300 S_000 | F_300 S_000 )^0_{t} = x * ( F_300 S_000 | D_200 S_000 )^0 + ( D_200 S_000 | D_200 S_000 )^0 + ( F_300 S_000 | P_100 S_000 )^0 - ( G_400 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[0] = etfac[0] * AUX_INT__f_s_d_s[0] + 3 * one_over_2q * AUX_INT__d_s_d_s[0] + 2 * one_over_2q * AUX_INT__f_s_p_s[0] - p_over_q * AUX_INT__g_s_d_s[0];

                    // ( F_300 S_000 | F_210 S_000 )^0_{t} = y * ( F_300 S_000 | D_200 S_000 )^0 - ( G_310 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[1] = etfac[1] * AUX_INT__f_s_d_s[0] - p_over_q * AUX_INT__g_s_d_s[6];

                    // ( F_300 S_000 | F_201 S_000 )^0_{t} = z * ( F_300 S_000 | D_200 S_000 )^0 - ( G_301 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[2] = etfac[2] * AUX_INT__f_s_d_s[0] - p_over_q * AUX_INT__g_s_d_s[12];

                    // ( F_300 S_000 | F_120 S_000 )^0_{t} = x * ( F_300 S_000 | D_020 S_000 )^0 + ( D_200 S_000 | D_020 S_000 )^0 - ( G_400 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[3] = etfac[0] * AUX_INT__f_s_d_s[3] + 3 * one_over_2q * AUX_INT__d_s_d_s[3] - p_over_q * AUX_INT__g_s_d_s[3];

                    // ( F_300 S_000 | F_111 S_000 )^0_{t} = z * ( F_300 S_000 | D_110 S_000 )^0 - ( G_301 S_000 | D_110 S_000 )^0
                    AUX_INT__f_s_f_s[4] = etfac[2] * AUX_INT__f_s_d_s[1] - p_over_q * AUX_INT__g_s_d_s[13];

                    // ( F_300 S_000 | F_102 S_000 )^0_{t} = x * ( F_300 S_000 | D_002 S_000 )^0 + ( D_200 S_000 | D_002 S_000 )^0 - ( G_400 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[5] = etfac[0] * AUX_INT__f_s_d_s[5] + 3 * one_over_2q * AUX_INT__d_s_d_s[5] - p_over_q * AUX_INT__g_s_d_s[5];

                    // ( F_300 S_000 | F_030 S_000 )^0_{t} = y * ( F_300 S_000 | D_020 S_000 )^0 + ( F_300 S_000 | P_010 S_000 )^0 - ( G_310 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[6] = etfac[1] * AUX_INT__f_s_d_s[3] + 2 * one_over_2q * AUX_INT__f_s_p_s[1] - p_over_q * AUX_INT__g_s_d_s[9];

                    // ( F_300 S_000 | F_021 S_000 )^0_{t} = z * ( F_300 S_000 | D_020 S_000 )^0 - ( G_301 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[7] = etfac[2] * AUX_INT__f_s_d_s[3] - p_over_q * AUX_INT__g_s_d_s[15];

                    // ( F_300 S_000 | F_012 S_000 )^0_{t} = y * ( F_300 S_000 | D_002 S_000 )^0 - ( G_310 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[8] = etfac[1] * AUX_INT__f_s_d_s[5] - p_over_q * AUX_INT__g_s_d_s[11];

                    // ( F_300 S_000 | F_003 S_000 )^0_{t} = z * ( F_300 S_000 | D_002 S_000 )^0 + ( F_300 S_000 | P_001 S_000 )^0 - ( G_301 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[9] = etfac[2] * AUX_INT__f_s_d_s[5] + 2 * one_over_2q * AUX_INT__f_s_p_s[2] - p_over_q * AUX_INT__g_s_d_s[17];

                    // ( F_210 S_000 | F_300 S_000 )^0_{t} = x * ( F_210 S_000 | D_200 S_000 )^0 + ( D_110 S_000 | D_200 S_000 )^0 + ( F_210 S_000 | P_100 S_000 )^0 - ( G_310 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[10] = etfac[0] * AUX_INT__f_s_d_s[6] + 2 * one_over_2q * AUX_INT__d_s_d_s[6] + 2 * one_over_2q * AUX_INT__f_s_p_s[3] - p_over_q * AUX_INT__g_s_d_s[6];

                    // ( F_210 S_000 | F_210 S_000 )^0_{t} = y * ( F_210 S_000 | D_200 S_000 )^0 + ( D_200 S_000 | D_200 S_000 )^0 - ( G_220 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[11] = etfac[1] * AUX_INT__f_s_d_s[6] + 1 * one_over_2q * AUX_INT__d_s_d_s[0] - p_over_q * AUX_INT__g_s_d_s[18];

                    // ( F_210 S_000 | F_201 S_000 )^0_{t} = z * ( F_210 S_000 | D_200 S_000 )^0 - ( G_211 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[12] = etfac[2] * AUX_INT__f_s_d_s[6] - p_over_q * AUX_INT__g_s_d_s[24];

                    // ( F_210 S_000 | F_120 S_000 )^0_{t} = x * ( F_210 S_000 | D_020 S_000 )^0 + ( D_110 S_000 | D_020 S_000 )^0 - ( G_310 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[13] = etfac[0] * AUX_INT__f_s_d_s[9] + 2 * one_over_2q * AUX_INT__d_s_d_s[9] - p_over_q * AUX_INT__g_s_d_s[9];

                    // ( F_210 S_000 | F_111 S_000 )^0_{t} = z * ( F_210 S_000 | D_110 S_000 )^0 - ( G_211 S_000 | D_110 S_000 )^0
                    AUX_INT__f_s_f_s[14] = etfac[2] * AUX_INT__f_s_d_s[7] - p_over_q * AUX_INT__g_s_d_s[25];

                    // ( F_210 S_000 | F_102 S_000 )^0_{t} = x * ( F_210 S_000 | D_002 S_000 )^0 + ( D_110 S_000 | D_002 S_000 )^0 - ( G_310 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[15] = etfac[0] * AUX_INT__f_s_d_s[11] + 2 * one_over_2q * AUX_INT__d_s_d_s[11] - p_over_q * AUX_INT__g_s_d_s[11];

                    // ( F_210 S_000 | F_030 S_000 )^0_{t} = y * ( F_210 S_000 | D_020 S_000 )^0 + ( D_200 S_000 | D_020 S_000 )^0 + ( F_210 S_000 | P_010 S_000 )^0 - ( G_220 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[16] = etfac[1] * AUX_INT__f_s_d_s[9] + 1 * one_over_2q * AUX_INT__d_s_d_s[3] + 2 * one_over_2q * AUX_INT__f_s_p_s[4] - p_over_q * AUX_INT__g_s_d_s[21];

                    // ( F_210 S_000 | F_021 S_000 )^0_{t} = z * ( F_210 S_000 | D_020 S_000 )^0 - ( G_211 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[17] = etfac[2] * AUX_INT__f_s_d_s[9] - p_over_q * AUX_INT__g_s_d_s[27];

                    // ( F_210 S_000 | F_012 S_000 )^0_{t} = y * ( F_210 S_000 | D_002 S_000 )^0 + ( D_200 S_000 | D_002 S_000 )^0 - ( G_220 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[18] = etfac[1] * AUX_INT__f_s_d_s[11] + 1 * one_over_2q * AUX_INT__d_s_d_s[5] - p_over_q * AUX_INT__g_s_d_s[23];

                    // ( F_210 S_000 | F_003 S_000 )^0_{t} = z * ( F_210 S_000 | D_002 S_000 )^0 + ( F_210 S_000 | P_001 S_000 )^0 - ( G_211 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[19] = etfac[2] * AUX_INT__f_s_d_s[11] + 2 * one_over_2q * AUX_INT__f_s_p_s[5] - p_over_q * AUX_INT__g_s_d_s[29];

                    // ( F_201 S_000 | F_300 S_000 )^0_{t} = x * ( F_201 S_000 | D_200 S_000 )^0 + ( D_101 S_000 | D_200 S_000 )^0 + ( F_201 S_000 | P_100 S_000 )^0 - ( G_301 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[20] = etfac[0] * AUX_INT__f_s_d_s[12] + 2 * one_over_2q * AUX_INT__d_s_d_s[12] + 2 * one_over_2q * AUX_INT__f_s_p_s[6] - p_over_q * AUX_INT__g_s_d_s[12];

                    // ( F_201 S_000 | F_210 S_000 )^0_{t} = y * ( F_201 S_000 | D_200 S_000 )^0 - ( G_211 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[21] = etfac[1] * AUX_INT__f_s_d_s[12] - p_over_q * AUX_INT__g_s_d_s[24];

                    // ( F_201 S_000 | F_201 S_000 )^0_{t} = z * ( F_201 S_000 | D_200 S_000 )^0 + ( D_200 S_000 | D_200 S_000 )^0 - ( G_202 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[22] = etfac[2] * AUX_INT__f_s_d_s[12] + 1 * one_over_2q * AUX_INT__d_s_d_s[0] - p_over_q * AUX_INT__g_s_d_s[30];

                    // ( F_201 S_000 | F_120 S_000 )^0_{t} = x * ( F_201 S_000 | D_020 S_000 )^0 + ( D_101 S_000 | D_020 S_000 )^0 - ( G_301 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[23] = etfac[0] * AUX_INT__f_s_d_s[15] + 2 * one_over_2q * AUX_INT__d_s_d_s[15] - p_over_q * AUX_INT__g_s_d_s[15];

                    // ( F_201 S_000 | F_111 S_000 )^0_{t} = z * ( F_201 S_000 | D_110 S_000 )^0 + ( D_200 S_000 | D_110 S_000 )^0 - ( G_202 S_000 | D_110 S_000 )^0
                    AUX_INT__f_s_f_s[24] = etfac[2] * AUX_INT__f_s_d_s[13] + 1 * one_over_2q * AUX_INT__d_s_d_s[1] - p_over_q * AUX_INT__g_s_d_s[31];

                    // ( F_201 S_000 | F_102 S_000 )^0_{t} = x * ( F_201 S_000 | D_002 S_000 )^0 + ( D_101 S_000 | D_002 S_000 )^0 - ( G_301 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[25] = etfac[0] * AUX_INT__f_s_d_s[17] + 2 * one_over_2q * AUX_INT__d_s_d_s[17] - p_over_q * AUX_INT__g_s_d_s[17];

                    // ( F_201 S_000 | F_030 S_000 )^0_{t} = y * ( F_201 S_000 | D_020 S_000 )^0 + ( F_201 S_000 | P_010 S_000 )^0 - ( G_211 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[26] = etfac[1] * AUX_INT__f_s_d_s[15] + 2 * one_over_2q * AUX_INT__f_s_p_s[7] - p_over_q * AUX_INT__g_s_d_s[27];

                    // ( F_201 S_000 | F_021 S_000 )^0_{t} = z * ( F_201 S_000 | D_020 S_000 )^0 + ( D_200 S_000 | D_020 S_000 )^0 - ( G_202 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[27] = etfac[2] * AUX_INT__f_s_d_s[15] + 1 * one_over_2q * AUX_INT__d_s_d_s[3] - p_over_q * AUX_INT__g_s_d_s[33];

                    // ( F_201 S_000 | F_012 S_000 )^0_{t} = y * ( F_201 S_000 | D_002 S_000 )^0 - ( G_211 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[28] = etfac[1] * AUX_INT__f_s_d_s[17] - p_over_q * AUX_INT__g_s_d_s[29];

                    // ( F_201 S_000 | F_003 S_000 )^0_{t} = z * ( F_201 S_000 | D_002 S_000 )^0 + ( D_200 S_000 | D_002 S_000 )^0 + ( F_201 S_000 | P_001 S_000 )^0 - ( G_202 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[29] = etfac[2] * AUX_INT__f_s_d_s[17] + 1 * one_over_2q * AUX_INT__d_s_d_s[5] + 2 * one_over_2q * AUX_INT__f_s_p_s[8] - p_over_q * AUX_INT__g_s_d_s[35];

                    // ( F_120 S_000 | F_300 S_000 )^0_{t} = x * ( F_120 S_000 | D_200 S_000 )^0 + ( D_020 S_000 | D_200 S_000 )^0 + ( F_120 S_000 | P_100 S_000 )^0 - ( G_220 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[30] = etfac[0] * AUX_INT__f_s_d_s[18] + 1 * one_over_2q * AUX_INT__d_s_d_s[18] + 2 * one_over_2q * AUX_INT__f_s_p_s[9] - p_over_q * AUX_INT__g_s_d_s[18];

                    // ( F_120 S_000 | F_210 S_000 )^0_{t} = y * ( F_120 S_000 | D_200 S_000 )^0 + ( D_110 S_000 | D_200 S_000 )^0 - ( G_130 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[31] = etfac[1] * AUX_INT__f_s_d_s[18] + 2 * one_over_2q * AUX_INT__d_s_d_s[6] - p_over_q * AUX_INT__g_s_d_s[36];

                    // ( F_120 S_000 | F_201 S_000 )^0_{t} = z * ( F_120 S_000 | D_200 S_000 )^0 - ( G_121 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[32] = etfac[2] * AUX_INT__f_s_d_s[18] - p_over_q * AUX_INT__g_s_d_s[42];

                    // ( F_120 S_000 | F_120 S_000 )^0_{t} = x * ( F_120 S_000 | D_020 S_000 )^0 + ( D_020 S_000 | D_020 S_000 )^0 - ( G_220 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[33] = etfac[0] * AUX_INT__f_s_d_s[21] + 1 * one_over_2q * AUX_INT__d_s_d_s[21] - p_over_q * AUX_INT__g_s_d_s[21];

                    // ( F_120 S_000 | F_111 S_000 )^0_{t} = z * ( F_120 S_000 | D_110 S_000 )^0 - ( G_121 S_000 | D_110 S_000 )^0
                    AUX_INT__f_s_f_s[34] = etfac[2] * AUX_INT__f_s_d_s[19] - p_over_q * AUX_INT__g_s_d_s[43];

                    // ( F_120 S_000 | F_102 S_000 )^0_{t} = x * ( F_120 S_000 | D_002 S_000 )^0 + ( D_020 S_000 | D_002 S_000 )^0 - ( G_220 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[35] = etfac[0] * AUX_INT__f_s_d_s[23] + 1 * one_over_2q * AUX_INT__d_s_d_s[23] - p_over_q * AUX_INT__g_s_d_s[23];

                    // ( F_120 S_000 | F_030 S_000 )^0_{t} = y * ( F_120 S_000 | D_020 S_000 )^0 + ( D_110 S_000 | D_020 S_000 )^0 + ( F_120 S_000 | P_010 S_000 )^0 - ( G_130 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[36] = etfac[1] * AUX_INT__f_s_d_s[21] + 2 * one_over_2q * AUX_INT__d_s_d_s[9] + 2 * one_over_2q * AUX_INT__f_s_p_s[10] - p_over_q * AUX_INT__g_s_d_s[39];

                    // ( F_120 S_000 | F_021 S_000 )^0_{t} = z * ( F_120 S_000 | D_020 S_000 )^0 - ( G_121 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[37] = etfac[2] * AUX_INT__f_s_d_s[21] - p_over_q * AUX_INT__g_s_d_s[45];

                    // ( F_120 S_000 | F_012 S_000 )^0_{t} = y * ( F_120 S_000 | D_002 S_000 )^0 + ( D_110 S_000 | D_002 S_000 )^0 - ( G_130 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[38] = etfac[1] * AUX_INT__f_s_d_s[23] + 2 * one_over_2q * AUX_INT__d_s_d_s[11] - p_over_q * AUX_INT__g_s_d_s[41];

                    // ( F_120 S_000 | F_003 S_000 )^0_{t} = z * ( F_120 S_000 | D_002 S_000 )^0 + ( F_120 S_000 | P_001 S_000 )^0 - ( G_121 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[39] = etfac[2] * AUX_INT__f_s_d_s[23] + 2 * one_over_2q * AUX_INT__f_s_p_s[11] - p_over_q * AUX_INT__g_s_d_s[47];

                    // ( F_111 S_000 | F_300 S_000 )^0_{t} = x * ( F_111 S_000 | D_200 S_000 )^0 + ( D_011 S_000 | D_200 S_000 )^0 + ( F_111 S_000 | P_100 S_000 )^0 - ( G_211 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[40] = etfac[0] * AUX_INT__f_s_d_s[24] + 1 * one_over_2q * AUX_INT__d_s_d_s[24] + 2 * one_over_2q * AUX_INT__f_s_p_s[12] - p_over_q * AUX_INT__g_s_d_s[24];

                    // ( F_111 S_000 | F_210 S_000 )^0_{t} = y * ( F_111 S_000 | D_200 S_000 )^0 + ( D_101 S_000 | D_200 S_000 )^0 - ( G_121 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[41] = etfac[1] * AUX_INT__f_s_d_s[24] + 1 * one_over_2q * AUX_INT__d_s_d_s[12] - p_over_q * AUX_INT__g_s_d_s[42];

                    // ( F_111 S_000 | F_201 S_000 )^0_{t} = z * ( F_111 S_000 | D_200 S_000 )^0 + ( D_110 S_000 | D_200 S_000 )^0 - ( G_112 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[42] = etfac[2] * AUX_INT__f_s_d_s[24] + 1 * one_over_2q * AUX_INT__d_s_d_s[6] - p_over_q * AUX_INT__g_s_d_s[48];

                    // ( F_111 S_000 | F_120 S_000 )^0_{t} = x * ( F_111 S_000 | D_020 S_000 )^0 + ( D_011 S_000 | D_020 S_000 )^0 - ( G_211 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[43] = etfac[0] * AUX_INT__f_s_d_s[27] + 1 * one_over_2q * AUX_INT__d_s_d_s[27] - p_over_q * AUX_INT__g_s_d_s[27];

                    // ( F_111 S_000 | F_111 S_000 )^0_{t} = z * ( F_111 S_000 | D_110 S_000 )^0 + ( D_110 S_000 | D_110 S_000 )^0 - ( G_112 S_000 | D_110 S_000 )^0
                    AUX_INT__f_s_f_s[44] = etfac[2] * AUX_INT__f_s_d_s[25] + 1 * one_over_2q * AUX_INT__d_s_d_s[7] - p_over_q * AUX_INT__g_s_d_s[49];

                    // ( F_111 S_000 | F_102 S_000 )^0_{t} = x * ( F_111 S_000 | D_002 S_000 )^0 + ( D_011 S_000 | D_002 S_000 )^0 - ( G_211 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[45] = etfac[0] * AUX_INT__f_s_d_s[29] + 1 * one_over_2q * AUX_INT__d_s_d_s[29] - p_over_q * AUX_INT__g_s_d_s[29];

                    // ( F_111 S_000 | F_030 S_000 )^0_{t} = y * ( F_111 S_000 | D_020 S_000 )^0 + ( D_101 S_000 | D_020 S_000 )^0 + ( F_111 S_000 | P_010 S_000 )^0 - ( G_121 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[46] = etfac[1] * AUX_INT__f_s_d_s[27] + 1 * one_over_2q * AUX_INT__d_s_d_s[15] + 2 * one_over_2q * AUX_INT__f_s_p_s[13] - p_over_q * AUX_INT__g_s_d_s[45];

                    // ( F_111 S_000 | F_021 S_000 )^0_{t} = z * ( F_111 S_000 | D_020 S_000 )^0 + ( D_110 S_000 | D_020 S_000 )^0 - ( G_112 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[47] = etfac[2] * AUX_INT__f_s_d_s[27] + 1 * one_over_2q * AUX_INT__d_s_d_s[9] - p_over_q * AUX_INT__g_s_d_s[51];

                    // ( F_111 S_000 | F_012 S_000 )^0_{t} = y * ( F_111 S_000 | D_002 S_000 )^0 + ( D_101 S_000 | D_002 S_000 )^0 - ( G_121 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[48] = etfac[1] * AUX_INT__f_s_d_s[29] + 1 * one_over_2q * AUX_INT__d_s_d_s[17] - p_over_q * AUX_INT__g_s_d_s[47];

                    // ( F_111 S_000 | F_003 S_000 )^0_{t} = z * ( F_111 S_000 | D_002 S_000 )^0 + ( D_110 S_000 | D_002 S_000 )^0 + ( F_111 S_000 | P_001 S_000 )^0 - ( G_112 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[49] = etfac[2] * AUX_INT__f_s_d_s[29] + 1 * one_over_2q * AUX_INT__d_s_d_s[11] + 2 * one_over_2q * AUX_INT__f_s_p_s[14] - p_over_q * AUX_INT__g_s_d_s[53];

                    // ( F_102 S_000 | F_300 S_000 )^0_{t} = x * ( F_102 S_000 | D_200 S_000 )^0 + ( D_002 S_000 | D_200 S_000 )^0 + ( F_102 S_000 | P_100 S_000 )^0 - ( G_202 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[50] = etfac[0] * AUX_INT__f_s_d_s[30] + 1 * one_over_2q * AUX_INT__d_s_d_s[30] + 2 * one_over_2q * AUX_INT__f_s_p_s[15] - p_over_q * AUX_INT__g_s_d_s[30];

                    // ( F_102 S_000 | F_210 S_000 )^0_{t} = y * ( F_102 S_000 | D_200 S_000 )^0 - ( G_112 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[51] = etfac[1] * AUX_INT__f_s_d_s[30] - p_over_q * AUX_INT__g_s_d_s[48];

                    // ( F_102 S_000 | F_201 S_000 )^0_{t} = z * ( F_102 S_000 | D_200 S_000 )^0 + ( D_101 S_000 | D_200 S_000 )^0 - ( G_103 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[52] = etfac[2] * AUX_INT__f_s_d_s[30] + 2 * one_over_2q * AUX_INT__d_s_d_s[12] - p_over_q * AUX_INT__g_s_d_s[54];

                    // ( F_102 S_000 | F_120 S_000 )^0_{t} = x * ( F_102 S_000 | D_020 S_000 )^0 + ( D_002 S_000 | D_020 S_000 )^0 - ( G_202 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[53] = etfac[0] * AUX_INT__f_s_d_s[33] + 1 * one_over_2q * AUX_INT__d_s_d_s[33] - p_over_q * AUX_INT__g_s_d_s[33];

                    // ( F_102 S_000 | F_111 S_000 )^0_{t} = z * ( F_102 S_000 | D_110 S_000 )^0 + ( D_101 S_000 | D_110 S_000 )^0 - ( G_103 S_000 | D_110 S_000 )^0
                    AUX_INT__f_s_f_s[54] = etfac[2] * AUX_INT__f_s_d_s[31] + 2 * one_over_2q * AUX_INT__d_s_d_s[13] - p_over_q * AUX_INT__g_s_d_s[55];

                    // ( F_102 S_000 | F_102 S_000 )^0_{t} = x * ( F_102 S_000 | D_002 S_000 )^0 + ( D_002 S_000 | D_002 S_000 )^0 - ( G_202 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[55] = etfac[0] * AUX_INT__f_s_d_s[35] + 1 * one_over_2q * AUX_INT__d_s_d_s[35] - p_over_q * AUX_INT__g_s_d_s[35];

                    // ( F_102 S_000 | F_030 S_000 )^0_{t} = y * ( F_102 S_000 | D_020 S_000 )^0 + ( F_102 S_000 | P_010 S_000 )^0 - ( G_112 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[56] = etfac[1] * AUX_INT__f_s_d_s[33] + 2 * one_over_2q * AUX_INT__f_s_p_s[16] - p_over_q * AUX_INT__g_s_d_s[51];

                    // ( F_102 S_000 | F_021 S_000 )^0_{t} = z * ( F_102 S_000 | D_020 S_000 )^0 + ( D_101 S_000 | D_020 S_000 )^0 - ( G_103 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[57] = etfac[2] * AUX_INT__f_s_d_s[33] + 2 * one_over_2q * AUX_INT__d_s_d_s[15] - p_over_q * AUX_INT__g_s_d_s[57];

                    // ( F_102 S_000 | F_012 S_000 )^0_{t} = y * ( F_102 S_000 | D_002 S_000 )^0 - ( G_112 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[58] = etfac[1] * AUX_INT__f_s_d_s[35] - p_over_q * AUX_INT__g_s_d_s[53];

                    // ( F_102 S_000 | F_003 S_000 )^0_{t} = z * ( F_102 S_000 | D_002 S_000 )^0 + ( D_101 S_000 | D_002 S_000 )^0 + ( F_102 S_000 | P_001 S_000 )^0 - ( G_103 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[59] = etfac[2] * AUX_INT__f_s_d_s[35] + 2 * one_over_2q * AUX_INT__d_s_d_s[17] + 2 * one_over_2q * AUX_INT__f_s_p_s[17] - p_over_q * AUX_INT__g_s_d_s[59];

                    // ( F_030 S_000 | F_300 S_000 )^0_{t} = x * ( F_030 S_000 | D_200 S_000 )^0 + ( F_030 S_000 | P_100 S_000 )^0 - ( G_130 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[60] = etfac[0] * AUX_INT__f_s_d_s[36] + 2 * one_over_2q * AUX_INT__f_s_p_s[18] - p_over_q * AUX_INT__g_s_d_s[36];

                    // ( F_030 S_000 | F_210 S_000 )^0_{t} = y * ( F_030 S_000 | D_200 S_000 )^0 + ( D_020 S_000 | D_200 S_000 )^0 - ( G_040 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[61] = etfac[1] * AUX_INT__f_s_d_s[36] + 3 * one_over_2q * AUX_INT__d_s_d_s[18] - p_over_q * AUX_INT__g_s_d_s[60];

                    // ( F_030 S_000 | F_201 S_000 )^0_{t} = z * ( F_030 S_000 | D_200 S_000 )^0 - ( G_031 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[62] = etfac[2] * AUX_INT__f_s_d_s[36] - p_over_q * AUX_INT__g_s_d_s[66];

                    // ( F_030 S_000 | F_120 S_000 )^0_{t} = x * ( F_030 S_000 | D_020 S_000 )^0 - ( G_130 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[63] = etfac[0] * AUX_INT__f_s_d_s[39] - p_over_q * AUX_INT__g_s_d_s[39];

                    // ( F_030 S_000 | F_111 S_000 )^0_{t} = z * ( F_030 S_000 | D_110 S_000 )^0 - ( G_031 S_000 | D_110 S_000 )^0
                    AUX_INT__f_s_f_s[64] = etfac[2] * AUX_INT__f_s_d_s[37] - p_over_q * AUX_INT__g_s_d_s[67];

                    // ( F_030 S_000 | F_102 S_000 )^0_{t} = x * ( F_030 S_000 | D_002 S_000 )^0 - ( G_130 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[65] = etfac[0] * AUX_INT__f_s_d_s[41] - p_over_q * AUX_INT__g_s_d_s[41];

                    // ( F_030 S_000 | F_030 S_000 )^0_{t} = y * ( F_030 S_000 | D_020 S_000 )^0 + ( D_020 S_000 | D_020 S_000 )^0 + ( F_030 S_000 | P_010 S_000 )^0 - ( G_040 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[66] = etfac[1] * AUX_INT__f_s_d_s[39] + 3 * one_over_2q * AUX_INT__d_s_d_s[21] + 2 * one_over_2q * AUX_INT__f_s_p_s[19] - p_over_q * AUX_INT__g_s_d_s[63];

                    // ( F_030 S_000 | F_021 S_000 )^0_{t} = z * ( F_030 S_000 | D_020 S_000 )^0 - ( G_031 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[67] = etfac[2] * AUX_INT__f_s_d_s[39] - p_over_q * AUX_INT__g_s_d_s[69];

                    // ( F_030 S_000 | F_012 S_000 )^0_{t} = y * ( F_030 S_000 | D_002 S_000 )^0 + ( D_020 S_000 | D_002 S_000 )^0 - ( G_040 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[68] = etfac[1] * AUX_INT__f_s_d_s[41] + 3 * one_over_2q * AUX_INT__d_s_d_s[23] - p_over_q * AUX_INT__g_s_d_s[65];

                    // ( F_030 S_000 | F_003 S_000 )^0_{t} = z * ( F_030 S_000 | D_002 S_000 )^0 + ( F_030 S_000 | P_001 S_000 )^0 - ( G_031 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[69] = etfac[2] * AUX_INT__f_s_d_s[41] + 2 * one_over_2q * AUX_INT__f_s_p_s[20] - p_over_q * AUX_INT__g_s_d_s[71];

                    // ( F_021 S_000 | F_300 S_000 )^0_{t} = x * ( F_021 S_000 | D_200 S_000 )^0 + ( F_021 S_000 | P_100 S_000 )^0 - ( G_121 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[70] = etfac[0] * AUX_INT__f_s_d_s[42] + 2 * one_over_2q * AUX_INT__f_s_p_s[21] - p_over_q * AUX_INT__g_s_d_s[42];

                    // ( F_021 S_000 | F_210 S_000 )^0_{t} = y * ( F_021 S_000 | D_200 S_000 )^0 + ( D_011 S_000 | D_200 S_000 )^0 - ( G_031 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[71] = etfac[1] * AUX_INT__f_s_d_s[42] + 2 * one_over_2q * AUX_INT__d_s_d_s[24] - p_over_q * AUX_INT__g_s_d_s[66];

                    // ( F_021 S_000 | F_201 S_000 )^0_{t} = z * ( F_021 S_000 | D_200 S_000 )^0 + ( D_020 S_000 | D_200 S_000 )^0 - ( G_022 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[72] = etfac[2] * AUX_INT__f_s_d_s[42] + 1 * one_over_2q * AUX_INT__d_s_d_s[18] - p_over_q * AUX_INT__g_s_d_s[72];

                    // ( F_021 S_000 | F_120 S_000 )^0_{t} = x * ( F_021 S_000 | D_020 S_000 )^0 - ( G_121 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[73] = etfac[0] * AUX_INT__f_s_d_s[45] - p_over_q * AUX_INT__g_s_d_s[45];

                    // ( F_021 S_000 | F_111 S_000 )^0_{t} = z * ( F_021 S_000 | D_110 S_000 )^0 + ( D_020 S_000 | D_110 S_000 )^0 - ( G_022 S_000 | D_110 S_000 )^0
                    AUX_INT__f_s_f_s[74] = etfac[2] * AUX_INT__f_s_d_s[43] + 1 * one_over_2q * AUX_INT__d_s_d_s[19] - p_over_q * AUX_INT__g_s_d_s[73];

                    // ( F_021 S_000 | F_102 S_000 )^0_{t} = x * ( F_021 S_000 | D_002 S_000 )^0 - ( G_121 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[75] = etfac[0] * AUX_INT__f_s_d_s[47] - p_over_q * AUX_INT__g_s_d_s[47];

                    // ( F_021 S_000 | F_030 S_000 )^0_{t} = y * ( F_021 S_000 | D_020 S_000 )^0 + ( D_011 S_000 | D_020 S_000 )^0 + ( F_021 S_000 | P_010 S_000 )^0 - ( G_031 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[76] = etfac[1] * AUX_INT__f_s_d_s[45] + 2 * one_over_2q * AUX_INT__d_s_d_s[27] + 2 * one_over_2q * AUX_INT__f_s_p_s[22] - p_over_q * AUX_INT__g_s_d_s[69];

                    // ( F_021 S_000 | F_021 S_000 )^0_{t} = z * ( F_021 S_000 | D_020 S_000 )^0 + ( D_020 S_000 | D_020 S_000 )^0 - ( G_022 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[77] = etfac[2] * AUX_INT__f_s_d_s[45] + 1 * one_over_2q * AUX_INT__d_s_d_s[21] - p_over_q * AUX_INT__g_s_d_s[75];

                    // ( F_021 S_000 | F_012 S_000 )^0_{t} = y * ( F_021 S_000 | D_002 S_000 )^0 + ( D_011 S_000 | D_002 S_000 )^0 - ( G_031 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[78] = etfac[1] * AUX_INT__f_s_d_s[47] + 2 * one_over_2q * AUX_INT__d_s_d_s[29] - p_over_q * AUX_INT__g_s_d_s[71];

                    // ( F_021 S_000 | F_003 S_000 )^0_{t} = z * ( F_021 S_000 | D_002 S_000 )^0 + ( D_020 S_000 | D_002 S_000 )^0 + ( F_021 S_000 | P_001 S_000 )^0 - ( G_022 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[79] = etfac[2] * AUX_INT__f_s_d_s[47] + 1 * one_over_2q * AUX_INT__d_s_d_s[23] + 2 * one_over_2q * AUX_INT__f_s_p_s[23] - p_over_q * AUX_INT__g_s_d_s[77];

                    // ( F_012 S_000 | F_300 S_000 )^0_{t} = x * ( F_012 S_000 | D_200 S_000 )^0 + ( F_012 S_000 | P_100 S_000 )^0 - ( G_112 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[80] = etfac[0] * AUX_INT__f_s_d_s[48] + 2 * one_over_2q * AUX_INT__f_s_p_s[24] - p_over_q * AUX_INT__g_s_d_s[48];

                    // ( F_012 S_000 | F_210 S_000 )^0_{t} = y * ( F_012 S_000 | D_200 S_000 )^0 + ( D_002 S_000 | D_200 S_000 )^0 - ( G_022 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[81] = etfac[1] * AUX_INT__f_s_d_s[48] + 1 * one_over_2q * AUX_INT__d_s_d_s[30] - p_over_q * AUX_INT__g_s_d_s[72];

                    // ( F_012 S_000 | F_201 S_000 )^0_{t} = z * ( F_012 S_000 | D_200 S_000 )^0 + ( D_011 S_000 | D_200 S_000 )^0 - ( G_013 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[82] = etfac[2] * AUX_INT__f_s_d_s[48] + 2 * one_over_2q * AUX_INT__d_s_d_s[24] - p_over_q * AUX_INT__g_s_d_s[78];

                    // ( F_012 S_000 | F_120 S_000 )^0_{t} = x * ( F_012 S_000 | D_020 S_000 )^0 - ( G_112 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[83] = etfac[0] * AUX_INT__f_s_d_s[51] - p_over_q * AUX_INT__g_s_d_s[51];

                    // ( F_012 S_000 | F_111 S_000 )^0_{t} = z * ( F_012 S_000 | D_110 S_000 )^0 + ( D_011 S_000 | D_110 S_000 )^0 - ( G_013 S_000 | D_110 S_000 )^0
                    AUX_INT__f_s_f_s[84] = etfac[2] * AUX_INT__f_s_d_s[49] + 2 * one_over_2q * AUX_INT__d_s_d_s[25] - p_over_q * AUX_INT__g_s_d_s[79];

                    // ( F_012 S_000 | F_102 S_000 )^0_{t} = x * ( F_012 S_000 | D_002 S_000 )^0 - ( G_112 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[85] = etfac[0] * AUX_INT__f_s_d_s[53] - p_over_q * AUX_INT__g_s_d_s[53];

                    // ( F_012 S_000 | F_030 S_000 )^0_{t} = y * ( F_012 S_000 | D_020 S_000 )^0 + ( D_002 S_000 | D_020 S_000 )^0 + ( F_012 S_000 | P_010 S_000 )^0 - ( G_022 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[86] = etfac[1] * AUX_INT__f_s_d_s[51] + 1 * one_over_2q * AUX_INT__d_s_d_s[33] + 2 * one_over_2q * AUX_INT__f_s_p_s[25] - p_over_q * AUX_INT__g_s_d_s[75];

                    // ( F_012 S_000 | F_021 S_000 )^0_{t} = z * ( F_012 S_000 | D_020 S_000 )^0 + ( D_011 S_000 | D_020 S_000 )^0 - ( G_013 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[87] = etfac[2] * AUX_INT__f_s_d_s[51] + 2 * one_over_2q * AUX_INT__d_s_d_s[27] - p_over_q * AUX_INT__g_s_d_s[81];

                    // ( F_012 S_000 | F_012 S_000 )^0_{t} = y * ( F_012 S_000 | D_002 S_000 )^0 + ( D_002 S_000 | D_002 S_000 )^0 - ( G_022 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[88] = etfac[1] * AUX_INT__f_s_d_s[53] + 1 * one_over_2q * AUX_INT__d_s_d_s[35] - p_over_q * AUX_INT__g_s_d_s[77];

                    // ( F_012 S_000 | F_003 S_000 )^0_{t} = z * ( F_012 S_000 | D_002 S_000 )^0 + ( D_011 S_000 | D_002 S_000 )^0 + ( F_012 S_000 | P_001 S_000 )^0 - ( G_013 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[89] = etfac[2] * AUX_INT__f_s_d_s[53] + 2 * one_over_2q * AUX_INT__d_s_d_s[29] + 2 * one_over_2q * AUX_INT__f_s_p_s[26] - p_over_q * AUX_INT__g_s_d_s[83];

                    // ( F_003 S_000 | F_300 S_000 )^0_{t} = x * ( F_003 S_000 | D_200 S_000 )^0 + ( F_003 S_000 | P_100 S_000 )^0 - ( G_103 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[90] = etfac[0] * AUX_INT__f_s_d_s[54] + 2 * one_over_2q * AUX_INT__f_s_p_s[27] - p_over_q * AUX_INT__g_s_d_s[54];

                    // ( F_003 S_000 | F_210 S_000 )^0_{t} = y * ( F_003 S_000 | D_200 S_000 )^0 - ( G_013 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[91] = etfac[1] * AUX_INT__f_s_d_s[54] - p_over_q * AUX_INT__g_s_d_s[78];

                    // ( F_003 S_000 | F_201 S_000 )^0_{t} = z * ( F_003 S_000 | D_200 S_000 )^0 + ( D_002 S_000 | D_200 S_000 )^0 - ( G_004 S_000 | D_200 S_000 )^0
                    AUX_INT__f_s_f_s[92] = etfac[2] * AUX_INT__f_s_d_s[54] + 3 * one_over_2q * AUX_INT__d_s_d_s[30] - p_over_q * AUX_INT__g_s_d_s[84];

                    // ( F_003 S_000 | F_120 S_000 )^0_{t} = x * ( F_003 S_000 | D_020 S_000 )^0 - ( G_103 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[93] = etfac[0] * AUX_INT__f_s_d_s[57] - p_over_q * AUX_INT__g_s_d_s[57];

                    // ( F_003 S_000 | F_111 S_000 )^0_{t} = z * ( F_003 S_000 | D_110 S_000 )^0 + ( D_002 S_000 | D_110 S_000 )^0 - ( G_004 S_000 | D_110 S_000 )^0
                    AUX_INT__f_s_f_s[94] = etfac[2] * AUX_INT__f_s_d_s[55] + 3 * one_over_2q * AUX_INT__d_s_d_s[31] - p_over_q * AUX_INT__g_s_d_s[85];

                    // ( F_003 S_000 | F_102 S_000 )^0_{t} = x * ( F_003 S_000 | D_002 S_000 )^0 - ( G_103 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[95] = etfac[0] * AUX_INT__f_s_d_s[59] - p_over_q * AUX_INT__g_s_d_s[59];

                    // ( F_003 S_000 | F_030 S_000 )^0_{t} = y * ( F_003 S_000 | D_020 S_000 )^0 + ( F_003 S_000 | P_010 S_000 )^0 - ( G_013 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[96] = etfac[1] * AUX_INT__f_s_d_s[57] + 2 * one_over_2q * AUX_INT__f_s_p_s[28] - p_over_q * AUX_INT__g_s_d_s[81];

                    // ( F_003 S_000 | F_021 S_000 )^0_{t} = z * ( F_003 S_000 | D_020 S_000 )^0 + ( D_002 S_000 | D_020 S_000 )^0 - ( G_004 S_000 | D_020 S_000 )^0
                    AUX_INT__f_s_f_s[97] = etfac[2] * AUX_INT__f_s_d_s[57] + 3 * one_over_2q * AUX_INT__d_s_d_s[33] - p_over_q * AUX_INT__g_s_d_s[87];

                    // ( F_003 S_000 | F_012 S_000 )^0_{t} = y * ( F_003 S_000 | D_002 S_000 )^0 - ( G_013 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[98] = etfac[1] * AUX_INT__f_s_d_s[59] - p_over_q * AUX_INT__g_s_d_s[83];

                    // ( F_003 S_000 | F_003 S_000 )^0_{t} = z * ( F_003 S_000 | D_002 S_000 )^0 + ( D_002 S_000 | D_002 S_000 )^0 + ( F_003 S_000 | P_001 S_000 )^0 - ( G_004 S_000 | D_002 S_000 )^0
                    AUX_INT__f_s_f_s[99] = etfac[2] * AUX_INT__f_s_d_s[59] + 3 * one_over_2q * AUX_INT__d_s_d_s[35] + 2 * one_over_2q * AUX_INT__f_s_p_s[29] - p_over_q * AUX_INT__g_s_d_s[89];

                    // ( S_000 S_000 | D_200 S_000 )^0 = x * ( S_000 S_000 | P_100 S_000 )^0 + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( P_100 S_000 | P_100 S_000 )^0
                    AUX_INT__s_s_d_s[0] = etfac[0] * AUX_INT__s_s_p_s[0] + 1 * one_over_2q * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__p_s_p_s[0];

                    // ( S_000 S_000 | D_020 S_000 )^0 = y * ( S_000 S_000 | P_010 S_000 )^0 + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( P_010 S_000 | P_010 S_000 )^0
                    AUX_INT__s_s_d_s[3] = etfac[1] * AUX_INT__s_s_p_s[1] + 1 * one_over_2q * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__p_s_p_s[4];

                    // ( S_000 S_000 | D_002 S_000 )^0 = z * ( S_000 S_000 | P_001 S_000 )^0 + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( P_001 S_000 | P_001 S_000 )^0
                    AUX_INT__s_s_d_s[5] = etfac[2] * AUX_INT__s_s_p_s[2] + 1 * one_over_2q * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__p_s_p_s[8];

                    // ( P_100 S_000 | F_300 S_000 )^0 = x * ( P_100 S_000 | D_200 S_000 )^0 + ( S_000 S_000 | D_200 S_000 )^0 + ( P_100 S_000 | P_100 S_000 )^0 - ( D_200 S_000 | D_200 S_000 )^0
                    AUX_INT__p_s_f_s[0] = etfac[0] * AUX_INT__p_s_d_s[0] + 1 * one_over_2q * AUX_INT__s_s_d_s[0] + 2 * one_over_2q * AUX_INT__p_s_p_s[0] - p_over_q * AUX_INT__d_s_d_s[0];

                    // ( P_100 S_000 | F_210 S_000 )^0 = y * ( P_100 S_000 | D_200 S_000 )^0 - ( D_110 S_000 | D_200 S_000 )^0
                    AUX_INT__p_s_f_s[1] = etfac[1] * AUX_INT__p_s_d_s[0] - p_over_q * AUX_INT__d_s_d_s[6];

                    // ( P_100 S_000 | F_201 S_000 )^0 = z * ( P_100 S_000 | D_200 S_000 )^0 - ( D_101 S_000 | D_200 S_000 )^0
                    AUX_INT__p_s_f_s[2] = etfac[2] * AUX_INT__p_s_d_s[0] - p_over_q * AUX_INT__d_s_d_s[12];

                    // ( P_100 S_000 | F_120 S_000 )^0 = x * ( P_100 S_000 | D_020 S_000 )^0 + ( S_000 S_000 | D_020 S_000 )^0 - ( D_200 S_000 | D_020 S_000 )^0
                    AUX_INT__p_s_f_s[3] = etfac[0] * AUX_INT__p_s_d_s[3] + 1 * one_over_2q * AUX_INT__s_s_d_s[3] - p_over_q * AUX_INT__d_s_d_s[3];

                    // ( P_100 S_000 | F_102 S_000 )^0 = x * ( P_100 S_000 | D_002 S_000 )^0 + ( S_000 S_000 | D_002 S_000 )^0 - ( D_200 S_000 | D_002 S_000 )^0
                    AUX_INT__p_s_f_s[5] = etfac[0] * AUX_INT__p_s_d_s[5] + 1 * one_over_2q * AUX_INT__s_s_d_s[5] - p_over_q * AUX_INT__d_s_d_s[5];

                    // ( P_100 S_000 | F_030 S_000 )^0 = y * ( P_100 S_000 | D_020 S_000 )^0 + ( P_100 S_000 | P_010 S_000 )^0 - ( D_110 S_000 | D_020 S_000 )^0
                    AUX_INT__p_s_f_s[6] = etfac[1] * AUX_INT__p_s_d_s[3] + 2 * one_over_2q * AUX_INT__p_s_p_s[1] - p_over_q * AUX_INT__d_s_d_s[9];

                    // ( P_100 S_000 | F_021 S_000 )^0 = z * ( P_100 S_000 | D_020 S_000 )^0 - ( D_101 S_000 | D_020 S_000 )^0
                    AUX_INT__p_s_f_s[7] = etfac[2] * AUX_INT__p_s_d_s[3] - p_over_q * AUX_INT__d_s_d_s[15];

                    // ( P_100 S_000 | F_003 S_000 )^0 = z * ( P_100 S_000 | D_002 S_000 )^0 + ( P_100 S_000 | P_001 S_000 )^0 - ( D_101 S_000 | D_002 S_000 )^0
                    AUX_INT__p_s_f_s[9] = etfac[2] * AUX_INT__p_s_d_s[5] + 2 * one_over_2q * AUX_INT__p_s_p_s[2] - p_over_q * AUX_INT__d_s_d_s[17];

                    // ( P_010 S_000 | F_300 S_000 )^0 = x * ( P_010 S_000 | D_200 S_000 )^0 + ( P_010 S_000 | P_100 S_000 )^0 - ( D_110 S_000 | D_200 S_000 )^0
                    AUX_INT__p_s_f_s[10] = etfac[0] * AUX_INT__p_s_d_s[6] + 2 * one_over_2q * AUX_INT__p_s_p_s[3] - p_over_q * AUX_INT__d_s_d_s[6];

                    // ( P_010 S_000 | F_210 S_000 )^0 = y * ( P_010 S_000 | D_200 S_000 )^0 + ( S_000 S_000 | D_200 S_000 )^0 - ( D_020 S_000 | D_200 S_000 )^0
                    AUX_INT__p_s_f_s[11] = etfac[1] * AUX_INT__p_s_d_s[6] + 1 * one_over_2q * AUX_INT__s_s_d_s[0] - p_over_q * AUX_INT__d_s_d_s[18];

                    // ( P_010 S_000 | F_201 S_000 )^0 = z * ( P_010 S_000 | D_200 S_000 )^0 - ( D_011 S_000 | D_200 S_000 )^0
                    AUX_INT__p_s_f_s[12] = etfac[2] * AUX_INT__p_s_d_s[6] - p_over_q * AUX_INT__d_s_d_s[24];

                    // ( P_010 S_000 | F_120 S_000 )^0 = x * ( P_010 S_000 | D_020 S_000 )^0 - ( D_110 S_000 | D_020 S_000 )^0
                    AUX_INT__p_s_f_s[13] = etfac[0] * AUX_INT__p_s_d_s[9] - p_over_q * AUX_INT__d_s_d_s[9];

                    // ( P_010 S_000 | F_102 S_000 )^0 = x * ( P_010 S_000 | D_002 S_000 )^0 - ( D_110 S_000 | D_002 S_000 )^0
                    AUX_INT__p_s_f_s[15] = etfac[0] * AUX_INT__p_s_d_s[11] - p_over_q * AUX_INT__d_s_d_s[11];

                    // ( P_010 S_000 | F_030 S_000 )^0 = y * ( P_010 S_000 | D_020 S_000 )^0 + ( S_000 S_000 | D_020 S_000 )^0 + ( P_010 S_000 | P_010 S_000 )^0 - ( D_020 S_000 | D_020 S_000 )^0
                    AUX_INT__p_s_f_s[16] = etfac[1] * AUX_INT__p_s_d_s[9] + 1 * one_over_2q * AUX_INT__s_s_d_s[3] + 2 * one_over_2q * AUX_INT__p_s_p_s[4] - p_over_q * AUX_INT__d_s_d_s[21];

                    // ( P_010 S_000 | F_021 S_000 )^0 = z * ( P_010 S_000 | D_020 S_000 )^0 - ( D_011 S_000 | D_020 S_000 )^0
                    AUX_INT__p_s_f_s[17] = etfac[2] * AUX_INT__p_s_d_s[9] - p_over_q * AUX_INT__d_s_d_s[27];

                    // ( P_010 S_000 | F_003 S_000 )^0 = z * ( P_010 S_000 | D_002 S_000 )^0 + ( P_010 S_000 | P_001 S_000 )^0 - ( D_011 S_000 | D_002 S_000 )^0
                    AUX_INT__p_s_f_s[19] = etfac[2] * AUX_INT__p_s_d_s[11] + 2 * one_over_2q * AUX_INT__p_s_p_s[5] - p_over_q * AUX_INT__d_s_d_s[29];

                    // ( P_001 S_000 | F_300 S_000 )^0 = x * ( P_001 S_000 | D_200 S_000 )^0 + ( P_001 S_000 | P_100 S_000 )^0 - ( D_101 S_000 | D_200 S_000 )^0
                    AUX_INT__p_s_f_s[20] = etfac[0] * AUX_INT__p_s_d_s[12] + 2 * one_over_2q * AUX_INT__p_s_p_s[6] - p_over_q * AUX_INT__d_s_d_s[12];

                    // ( P_001 S_000 | F_210 S_000 )^0 = y * ( P_001 S_000 | D_200 S_000 )^0 - ( D_011 S_000 | D_200 S_000 )^0
                    AUX_INT__p_s_f_s[21] = etfac[1] * AUX_INT__p_s_d_s[12] - p_over_q * AUX_INT__d_s_d_s[24];

                    // ( P_001 S_000 | F_201 S_000 )^0 = z * ( P_001 S_000 | D_200 S_000 )^0 + ( S_000 S_000 | D_200 S_000 )^0 - ( D_002 S_000 | D_200 S_000 )^0
                    AUX_INT__p_s_f_s[22] = etfac[2] * AUX_INT__p_s_d_s[12] + 1 * one_over_2q * AUX_INT__s_s_d_s[0] - p_over_q * AUX_INT__d_s_d_s[30];

                    // ( P_001 S_000 | F_120 S_000 )^0 = x * ( P_001 S_000 | D_020 S_000 )^0 - ( D_101 S_000 | D_020 S_000 )^0
                    AUX_INT__p_s_f_s[23] = etfac[0] * AUX_INT__p_s_d_s[15] - p_over_q * AUX_INT__d_s_d_s[15];

                    // ( P_001 S_000 | F_102 S_000 )^0 = x * ( P_001 S_000 | D_002 S_000 )^0 - ( D_101 S_000 | D_002 S_000 )^0
                    AUX_INT__p_s_f_s[25] = etfac[0] * AUX_INT__p_s_d_s[17] - p_over_q * AUX_INT__d_s_d_s[17];

                    // ( P_001 S_000 | F_030 S_000 )^0 = y * ( P_001 S_000 | D_020 S_000 )^0 + ( P_001 S_000 | P_010 S_000 )^0 - ( D_011 S_000 | D_020 S_000 )^0
                    AUX_INT__p_s_f_s[26] = etfac[1] * AUX_INT__p_s_d_s[15] + 2 * one_over_2q * AUX_INT__p_s_p_s[7] - p_over_q * AUX_INT__d_s_d_s[27];

                    // ( P_001 S_000 | F_021 S_000 )^0 = z * ( P_001 S_000 | D_020 S_000 )^0 + ( S_000 S_000 | D_020 S_000 )^0 - ( D_002 S_000 | D_020 S_000 )^0
                    AUX_INT__p_s_f_s[27] = etfac[2] * AUX_INT__p_s_d_s[15] + 1 * one_over_2q * AUX_INT__s_s_d_s[3] - p_over_q * AUX_INT__d_s_d_s[33];

                    // ( P_001 S_000 | F_003 S_000 )^0 = z * ( P_001 S_000 | D_002 S_000 )^0 + ( S_000 S_000 | D_002 S_000 )^0 + ( P_001 S_000 | P_001 S_000 )^0 - ( D_002 S_000 | D_002 S_000 )^0
                    AUX_INT__p_s_f_s[29] = etfac[2] * AUX_INT__p_s_d_s[17] + 1 * one_over_2q * AUX_INT__s_s_d_s[5] + 2 * one_over_2q * AUX_INT__p_s_p_s[8] - p_over_q * AUX_INT__d_s_d_s[35];

                    // ( D_200 S_000 | G_400 S_000 )^0_{t} = x * ( D_200 S_000 | F_300 S_000 )^0 + ( P_100 S_000 | F_300 S_000 )^0 + ( D_200 S_000 | D_200 S_000 )^0 - ( F_300 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[0] = etfac[0] * AUX_INT__d_s_f_s[0] + 2 * one_over_2q * AUX_INT__p_s_f_s[0] + 3 * one_over_2q * AUX_INT__d_s_d_s[0] - p_over_q * AUX_INT__f_s_f_s[0];

                    // ( D_200 S_000 | G_310 S_000 )^0_{t} = y * ( D_200 S_000 | F_300 S_000 )^0 - ( F_210 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[1] = etfac[1] * AUX_INT__d_s_f_s[0] - p_over_q * AUX_INT__f_s_f_s[10];

                    // ( D_200 S_000 | G_301 S_000 )^0_{t} = z * ( D_200 S_000 | F_300 S_000 )^0 - ( F_201 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[2] = etfac[2] * AUX_INT__d_s_f_s[0] - p_over_q * AUX_INT__f_s_f_s[20];

                    // ( D_200 S_000 | G_220 S_000 )^0_{t} = y * ( D_200 S_000 | F_210 S_000 )^0 + ( D_200 S_000 | D_200 S_000 )^0 - ( F_210 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[3] = etfac[1] * AUX_INT__d_s_f_s[1] + 1 * one_over_2q * AUX_INT__d_s_d_s[0] - p_over_q * AUX_INT__f_s_f_s[11];

                    // ( D_200 S_000 | G_211 S_000 )^0_{t} = z * ( D_200 S_000 | F_210 S_000 )^0 - ( F_201 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[4] = etfac[2] * AUX_INT__d_s_f_s[1] - p_over_q * AUX_INT__f_s_f_s[21];

                    // ( D_200 S_000 | G_202 S_000 )^0_{t} = z * ( D_200 S_000 | F_201 S_000 )^0 + ( D_200 S_000 | D_200 S_000 )^0 - ( F_201 S_000 | F_201 S_000 )^0
                    AUX_INT__d_s_g_s[5] = etfac[2] * AUX_INT__d_s_f_s[2] + 1 * one_over_2q * AUX_INT__d_s_d_s[0] - p_over_q * AUX_INT__f_s_f_s[22];

                    // ( D_200 S_000 | G_130 S_000 )^0_{t} = x * ( D_200 S_000 | F_030 S_000 )^0 + ( P_100 S_000 | F_030 S_000 )^0 - ( F_300 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[6] = etfac[0] * AUX_INT__d_s_f_s[6] + 2 * one_over_2q * AUX_INT__p_s_f_s[6] - p_over_q * AUX_INT__f_s_f_s[6];

                    // ( D_200 S_000 | G_121 S_000 )^0_{t} = z * ( D_200 S_000 | F_120 S_000 )^0 - ( F_201 S_000 | F_120 S_000 )^0
                    AUX_INT__d_s_g_s[7] = etfac[2] * AUX_INT__d_s_f_s[3] - p_over_q * AUX_INT__f_s_f_s[23];

                    // ( D_200 S_000 | G_112 S_000 )^0_{t} = y * ( D_200 S_000 | F_102 S_000 )^0 - ( F_210 S_000 | F_102 S_000 )^0
                    AUX_INT__d_s_g_s[8] = etfac[1] * AUX_INT__d_s_f_s[5] - p_over_q * AUX_INT__f_s_f_s[15];

                    // ( D_200 S_000 | G_103 S_000 )^0_{t} = x * ( D_200 S_000 | F_003 S_000 )^0 + ( P_100 S_000 | F_003 S_000 )^0 - ( F_300 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[9] = etfac[0] * AUX_INT__d_s_f_s[9] + 2 * one_over_2q * AUX_INT__p_s_f_s[9] - p_over_q * AUX_INT__f_s_f_s[9];

                    // ( D_200 S_000 | G_040 S_000 )^0_{t} = y * ( D_200 S_000 | F_030 S_000 )^0 + ( D_200 S_000 | D_020 S_000 )^0 - ( F_210 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[10] = etfac[1] * AUX_INT__d_s_f_s[6] + 3 * one_over_2q * AUX_INT__d_s_d_s[3] - p_over_q * AUX_INT__f_s_f_s[16];

                    // ( D_200 S_000 | G_031 S_000 )^0_{t} = z * ( D_200 S_000 | F_030 S_000 )^0 - ( F_201 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[11] = etfac[2] * AUX_INT__d_s_f_s[6] - p_over_q * AUX_INT__f_s_f_s[26];

                    // ( D_200 S_000 | G_022 S_000 )^0_{t} = z * ( D_200 S_000 | F_021 S_000 )^0 + ( D_200 S_000 | D_020 S_000 )^0 - ( F_201 S_000 | F_021 S_000 )^0
                    AUX_INT__d_s_g_s[12] = etfac[2] * AUX_INT__d_s_f_s[7] + 1 * one_over_2q * AUX_INT__d_s_d_s[3] - p_over_q * AUX_INT__f_s_f_s[27];

                    // ( D_200 S_000 | G_013 S_000 )^0_{t} = y * ( D_200 S_000 | F_003 S_000 )^0 - ( F_210 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[13] = etfac[1] * AUX_INT__d_s_f_s[9] - p_over_q * AUX_INT__f_s_f_s[19];

                    // ( D_200 S_000 | G_004 S_000 )^0_{t} = z * ( D_200 S_000 | F_003 S_000 )^0 + ( D_200 S_000 | D_002 S_000 )^0 - ( F_201 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[14] = etfac[2] * AUX_INT__d_s_f_s[9] + 3 * one_over_2q * AUX_INT__d_s_d_s[5] - p_over_q * AUX_INT__f_s_f_s[29];

                    // ( D_110 S_000 | G_400 S_000 )^0_{t} = x * ( D_110 S_000 | F_300 S_000 )^0 + ( P_010 S_000 | F_300 S_000 )^0 + ( D_110 S_000 | D_200 S_000 )^0 - ( F_210 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[15] = etfac[0] * AUX_INT__d_s_f_s[10] + 1 * one_over_2q * AUX_INT__p_s_f_s[10] + 3 * one_over_2q * AUX_INT__d_s_d_s[6] - p_over_q * AUX_INT__f_s_f_s[10];

                    // ( D_110 S_000 | G_310 S_000 )^0_{t} = y * ( D_110 S_000 | F_300 S_000 )^0 + ( P_100 S_000 | F_300 S_000 )^0 - ( F_120 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[16] = etfac[1] * AUX_INT__d_s_f_s[10] + 1 * one_over_2q * AUX_INT__p_s_f_s[0] - p_over_q * AUX_INT__f_s_f_s[30];

                    // ( D_110 S_000 | G_301 S_000 )^0_{t} = z * ( D_110 S_000 | F_300 S_000 )^0 - ( F_111 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[17] = etfac[2] * AUX_INT__d_s_f_s[10] - p_over_q * AUX_INT__f_s_f_s[40];

                    // ( D_110 S_000 | G_220 S_000 )^0_{t} = y * ( D_110 S_000 | F_210 S_000 )^0 + ( P_100 S_000 | F_210 S_000 )^0 + ( D_110 S_000 | D_200 S_000 )^0 - ( F_120 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[18] = etfac[1] * AUX_INT__d_s_f_s[11] + 1 * one_over_2q * AUX_INT__p_s_f_s[1] + 1 * one_over_2q * AUX_INT__d_s_d_s[6] - p_over_q * AUX_INT__f_s_f_s[31];

                    // ( D_110 S_000 | G_211 S_000 )^0_{t} = z * ( D_110 S_000 | F_210 S_000 )^0 - ( F_111 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[19] = etfac[2] * AUX_INT__d_s_f_s[11] - p_over_q * AUX_INT__f_s_f_s[41];

                    // ( D_110 S_000 | G_202 S_000 )^0_{t} = z * ( D_110 S_000 | F_201 S_000 )^0 + ( D_110 S_000 | D_200 S_000 )^0 - ( F_111 S_000 | F_201 S_000 )^0
                    AUX_INT__d_s_g_s[20] = etfac[2] * AUX_INT__d_s_f_s[12] + 1 * one_over_2q * AUX_INT__d_s_d_s[6] - p_over_q * AUX_INT__f_s_f_s[42];

                    // ( D_110 S_000 | G_130 S_000 )^0_{t} = x * ( D_110 S_000 | F_030 S_000 )^0 + ( P_010 S_000 | F_030 S_000 )^0 - ( F_210 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[21] = etfac[0] * AUX_INT__d_s_f_s[16] + 1 * one_over_2q * AUX_INT__p_s_f_s[16] - p_over_q * AUX_INT__f_s_f_s[16];

                    // ( D_110 S_000 | G_121 S_000 )^0_{t} = z * ( D_110 S_000 | F_120 S_000 )^0 - ( F_111 S_000 | F_120 S_000 )^0
                    AUX_INT__d_s_g_s[22] = etfac[2] * AUX_INT__d_s_f_s[13] - p_over_q * AUX_INT__f_s_f_s[43];

                    // ( D_110 S_000 | G_112 S_000 )^0_{t} = y * ( D_110 S_000 | F_102 S_000 )^0 + ( P_100 S_000 | F_102 S_000 )^0 - ( F_120 S_000 | F_102 S_000 )^0
                    AUX_INT__d_s_g_s[23] = etfac[1] * AUX_INT__d_s_f_s[15] + 1 * one_over_2q * AUX_INT__p_s_f_s[5] - p_over_q * AUX_INT__f_s_f_s[35];

                    // ( D_110 S_000 | G_103 S_000 )^0_{t} = x * ( D_110 S_000 | F_003 S_000 )^0 + ( P_010 S_000 | F_003 S_000 )^0 - ( F_210 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[24] = etfac[0] * AUX_INT__d_s_f_s[19] + 1 * one_over_2q * AUX_INT__p_s_f_s[19] - p_over_q * AUX_INT__f_s_f_s[19];

                    // ( D_110 S_000 | G_040 S_000 )^0_{t} = y * ( D_110 S_000 | F_030 S_000 )^0 + ( P_100 S_000 | F_030 S_000 )^0 + ( D_110 S_000 | D_020 S_000 )^0 - ( F_120 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[25] = etfac[1] * AUX_INT__d_s_f_s[16] + 1 * one_over_2q * AUX_INT__p_s_f_s[6] + 3 * one_over_2q * AUX_INT__d_s_d_s[9] - p_over_q * AUX_INT__f_s_f_s[36];

                    // ( D_110 S_000 | G_031 S_000 )^0_{t} = z * ( D_110 S_000 | F_030 S_000 )^0 - ( F_111 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[26] = etfac[2] * AUX_INT__d_s_f_s[16] - p_over_q * AUX_INT__f_s_f_s[46];

                    // ( D_110 S_000 | G_022 S_000 )^0_{t} = z * ( D_110 S_000 | F_021 S_000 )^0 + ( D_110 S_000 | D_020 S_000 )^0 - ( F_111 S_000 | F_021 S_000 )^0
                    AUX_INT__d_s_g_s[27] = etfac[2] * AUX_INT__d_s_f_s[17] + 1 * one_over_2q * AUX_INT__d_s_d_s[9] - p_over_q * AUX_INT__f_s_f_s[47];

                    // ( D_110 S_000 | G_013 S_000 )^0_{t} = y * ( D_110 S_000 | F_003 S_000 )^0 + ( P_100 S_000 | F_003 S_000 )^0 - ( F_120 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[28] = etfac[1] * AUX_INT__d_s_f_s[19] + 1 * one_over_2q * AUX_INT__p_s_f_s[9] - p_over_q * AUX_INT__f_s_f_s[39];

                    // ( D_110 S_000 | G_004 S_000 )^0_{t} = z * ( D_110 S_000 | F_003 S_000 )^0 + ( D_110 S_000 | D_002 S_000 )^0 - ( F_111 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[29] = etfac[2] * AUX_INT__d_s_f_s[19] + 3 * one_over_2q * AUX_INT__d_s_d_s[11] - p_over_q * AUX_INT__f_s_f_s[49];

                    // ( D_101 S_000 | G_400 S_000 )^0_{t} = x * ( D_101 S_000 | F_300 S_000 )^0 + ( P_001 S_000 | F_300 S_000 )^0 + ( D_101 S_000 | D_200 S_000 )^0 - ( F_201 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[30] = etfac[0] * AUX_INT__d_s_f_s[20] + 1 * one_over_2q * AUX_INT__p_s_f_s[20] + 3 * one_over_2q * AUX_INT__d_s_d_s[12] - p_over_q * AUX_INT__f_s_f_s[20];

                    // ( D_101 S_000 | G_310 S_000 )^0_{t} = y * ( D_101 S_000 | F_300 S_000 )^0 - ( F_111 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[31] = etfac[1] * AUX_INT__d_s_f_s[20] - p_over_q * AUX_INT__f_s_f_s[40];

                    // ( D_101 S_000 | G_301 S_000 )^0_{t} = z * ( D_101 S_000 | F_300 S_000 )^0 + ( P_100 S_000 | F_300 S_000 )^0 - ( F_102 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[32] = etfac[2] * AUX_INT__d_s_f_s[20] + 1 * one_over_2q * AUX_INT__p_s_f_s[0] - p_over_q * AUX_INT__f_s_f_s[50];

                    // ( D_101 S_000 | G_220 S_000 )^0_{t} = y * ( D_101 S_000 | F_210 S_000 )^0 + ( D_101 S_000 | D_200 S_000 )^0 - ( F_111 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[33] = etfac[1] * AUX_INT__d_s_f_s[21] + 1 * one_over_2q * AUX_INT__d_s_d_s[12] - p_over_q * AUX_INT__f_s_f_s[41];

                    // ( D_101 S_000 | G_211 S_000 )^0_{t} = z * ( D_101 S_000 | F_210 S_000 )^0 + ( P_100 S_000 | F_210 S_000 )^0 - ( F_102 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[34] = etfac[2] * AUX_INT__d_s_f_s[21] + 1 * one_over_2q * AUX_INT__p_s_f_s[1] - p_over_q * AUX_INT__f_s_f_s[51];

                    // ( D_101 S_000 | G_202 S_000 )^0_{t} = z * ( D_101 S_000 | F_201 S_000 )^0 + ( P_100 S_000 | F_201 S_000 )^0 + ( D_101 S_000 | D_200 S_000 )^0 - ( F_102 S_000 | F_201 S_000 )^0
                    AUX_INT__d_s_g_s[35] = etfac[2] * AUX_INT__d_s_f_s[22] + 1 * one_over_2q * AUX_INT__p_s_f_s[2] + 1 * one_over_2q * AUX_INT__d_s_d_s[12] - p_over_q * AUX_INT__f_s_f_s[52];

                    // ( D_101 S_000 | G_130 S_000 )^0_{t} = x * ( D_101 S_000 | F_030 S_000 )^0 + ( P_001 S_000 | F_030 S_000 )^0 - ( F_201 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[36] = etfac[0] * AUX_INT__d_s_f_s[26] + 1 * one_over_2q * AUX_INT__p_s_f_s[26] - p_over_q * AUX_INT__f_s_f_s[26];

                    // ( D_101 S_000 | G_121 S_000 )^0_{t} = z * ( D_101 S_000 | F_120 S_000 )^0 + ( P_100 S_000 | F_120 S_000 )^0 - ( F_102 S_000 | F_120 S_000 )^0
                    AUX_INT__d_s_g_s[37] = etfac[2] * AUX_INT__d_s_f_s[23] + 1 * one_over_2q * AUX_INT__p_s_f_s[3] - p_over_q * AUX_INT__f_s_f_s[53];

                    // ( D_101 S_000 | G_112 S_000 )^0_{t} = y * ( D_101 S_000 | F_102 S_000 )^0 - ( F_111 S_000 | F_102 S_000 )^0
                    AUX_INT__d_s_g_s[38] = etfac[1] * AUX_INT__d_s_f_s[25] - p_over_q * AUX_INT__f_s_f_s[45];

                    // ( D_101 S_000 | G_103 S_000 )^0_{t} = x * ( D_101 S_000 | F_003 S_000 )^0 + ( P_001 S_000 | F_003 S_000 )^0 - ( F_201 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[39] = etfac[0] * AUX_INT__d_s_f_s[29] + 1 * one_over_2q * AUX_INT__p_s_f_s[29] - p_over_q * AUX_INT__f_s_f_s[29];

                    // ( D_101 S_000 | G_040 S_000 )^0_{t} = y * ( D_101 S_000 | F_030 S_000 )^0 + ( D_101 S_000 | D_020 S_000 )^0 - ( F_111 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[40] = etfac[1] * AUX_INT__d_s_f_s[26] + 3 * one_over_2q * AUX_INT__d_s_d_s[15] - p_over_q * AUX_INT__f_s_f_s[46];

                    // ( D_101 S_000 | G_031 S_000 )^0_{t} = z * ( D_101 S_000 | F_030 S_000 )^0 + ( P_100 S_000 | F_030 S_000 )^0 - ( F_102 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[41] = etfac[2] * AUX_INT__d_s_f_s[26] + 1 * one_over_2q * AUX_INT__p_s_f_s[6] - p_over_q * AUX_INT__f_s_f_s[56];

                    // ( D_101 S_000 | G_022 S_000 )^0_{t} = z * ( D_101 S_000 | F_021 S_000 )^0 + ( P_100 S_000 | F_021 S_000 )^0 + ( D_101 S_000 | D_020 S_000 )^0 - ( F_102 S_000 | F_021 S_000 )^0
                    AUX_INT__d_s_g_s[42] = etfac[2] * AUX_INT__d_s_f_s[27] + 1 * one_over_2q * AUX_INT__p_s_f_s[7] + 1 * one_over_2q * AUX_INT__d_s_d_s[15] - p_over_q * AUX_INT__f_s_f_s[57];

                    // ( D_101 S_000 | G_013 S_000 )^0_{t} = y * ( D_101 S_000 | F_003 S_000 )^0 - ( F_111 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[43] = etfac[1] * AUX_INT__d_s_f_s[29] - p_over_q * AUX_INT__f_s_f_s[49];

                    // ( D_101 S_000 | G_004 S_000 )^0_{t} = z * ( D_101 S_000 | F_003 S_000 )^0 + ( P_100 S_000 | F_003 S_000 )^0 + ( D_101 S_000 | D_002 S_000 )^0 - ( F_102 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[44] = etfac[2] * AUX_INT__d_s_f_s[29] + 1 * one_over_2q * AUX_INT__p_s_f_s[9] + 3 * one_over_2q * AUX_INT__d_s_d_s[17] - p_over_q * AUX_INT__f_s_f_s[59];

                    // ( D_020 S_000 | G_400 S_000 )^0_{t} = x * ( D_020 S_000 | F_300 S_000 )^0 + ( D_020 S_000 | D_200 S_000 )^0 - ( F_120 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[45] = etfac[0] * AUX_INT__d_s_f_s[30] + 3 * one_over_2q * AUX_INT__d_s_d_s[18] - p_over_q * AUX_INT__f_s_f_s[30];

                    // ( D_020 S_000 | G_310 S_000 )^0_{t} = y * ( D_020 S_000 | F_300 S_000 )^0 + ( P_010 S_000 | F_300 S_000 )^0 - ( F_030 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[46] = etfac[1] * AUX_INT__d_s_f_s[30] + 2 * one_over_2q * AUX_INT__p_s_f_s[10] - p_over_q * AUX_INT__f_s_f_s[60];

                    // ( D_020 S_000 | G_301 S_000 )^0_{t} = z * ( D_020 S_000 | F_300 S_000 )^0 - ( F_021 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[47] = etfac[2] * AUX_INT__d_s_f_s[30] - p_over_q * AUX_INT__f_s_f_s[70];

                    // ( D_020 S_000 | G_220 S_000 )^0_{t} = y * ( D_020 S_000 | F_210 S_000 )^0 + ( P_010 S_000 | F_210 S_000 )^0 + ( D_020 S_000 | D_200 S_000 )^0 - ( F_030 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[48] = etfac[1] * AUX_INT__d_s_f_s[31] + 2 * one_over_2q * AUX_INT__p_s_f_s[11] + 1 * one_over_2q * AUX_INT__d_s_d_s[18] - p_over_q * AUX_INT__f_s_f_s[61];

                    // ( D_020 S_000 | G_211 S_000 )^0_{t} = z * ( D_020 S_000 | F_210 S_000 )^0 - ( F_021 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[49] = etfac[2] * AUX_INT__d_s_f_s[31] - p_over_q * AUX_INT__f_s_f_s[71];

                    // ( D_020 S_000 | G_202 S_000 )^0_{t} = z * ( D_020 S_000 | F_201 S_000 )^0 + ( D_020 S_000 | D_200 S_000 )^0 - ( F_021 S_000 | F_201 S_000 )^0
                    AUX_INT__d_s_g_s[50] = etfac[2] * AUX_INT__d_s_f_s[32] + 1 * one_over_2q * AUX_INT__d_s_d_s[18] - p_over_q * AUX_INT__f_s_f_s[72];

                    // ( D_020 S_000 | G_130 S_000 )^0_{t} = x * ( D_020 S_000 | F_030 S_000 )^0 - ( F_120 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[51] = etfac[0] * AUX_INT__d_s_f_s[36] - p_over_q * AUX_INT__f_s_f_s[36];

                    // ( D_020 S_000 | G_121 S_000 )^0_{t} = z * ( D_020 S_000 | F_120 S_000 )^0 - ( F_021 S_000 | F_120 S_000 )^0
                    AUX_INT__d_s_g_s[52] = etfac[2] * AUX_INT__d_s_f_s[33] - p_over_q * AUX_INT__f_s_f_s[73];

                    // ( D_020 S_000 | G_112 S_000 )^0_{t} = y * ( D_020 S_000 | F_102 S_000 )^0 + ( P_010 S_000 | F_102 S_000 )^0 - ( F_030 S_000 | F_102 S_000 )^0
                    AUX_INT__d_s_g_s[53] = etfac[1] * AUX_INT__d_s_f_s[35] + 2 * one_over_2q * AUX_INT__p_s_f_s[15] - p_over_q * AUX_INT__f_s_f_s[65];

                    // ( D_020 S_000 | G_103 S_000 )^0_{t} = x * ( D_020 S_000 | F_003 S_000 )^0 - ( F_120 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[54] = etfac[0] * AUX_INT__d_s_f_s[39] - p_over_q * AUX_INT__f_s_f_s[39];

                    // ( D_020 S_000 | G_040 S_000 )^0_{t} = y * ( D_020 S_000 | F_030 S_000 )^0 + ( P_010 S_000 | F_030 S_000 )^0 + ( D_020 S_000 | D_020 S_000 )^0 - ( F_030 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[55] = etfac[1] * AUX_INT__d_s_f_s[36] + 2 * one_over_2q * AUX_INT__p_s_f_s[16] + 3 * one_over_2q * AUX_INT__d_s_d_s[21] - p_over_q * AUX_INT__f_s_f_s[66];

                    // ( D_020 S_000 | G_031 S_000 )^0_{t} = z * ( D_020 S_000 | F_030 S_000 )^0 - ( F_021 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[56] = etfac[2] * AUX_INT__d_s_f_s[36] - p_over_q * AUX_INT__f_s_f_s[76];

                    // ( D_020 S_000 | G_022 S_000 )^0_{t} = z * ( D_020 S_000 | F_021 S_000 )^0 + ( D_020 S_000 | D_020 S_000 )^0 - ( F_021 S_000 | F_021 S_000 )^0
                    AUX_INT__d_s_g_s[57] = etfac[2] * AUX_INT__d_s_f_s[37] + 1 * one_over_2q * AUX_INT__d_s_d_s[21] - p_over_q * AUX_INT__f_s_f_s[77];

                    // ( D_020 S_000 | G_013 S_000 )^0_{t} = y * ( D_020 S_000 | F_003 S_000 )^0 + ( P_010 S_000 | F_003 S_000 )^0 - ( F_030 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[58] = etfac[1] * AUX_INT__d_s_f_s[39] + 2 * one_over_2q * AUX_INT__p_s_f_s[19] - p_over_q * AUX_INT__f_s_f_s[69];

                    // ( D_020 S_000 | G_004 S_000 )^0_{t} = z * ( D_020 S_000 | F_003 S_000 )^0 + ( D_020 S_000 | D_002 S_000 )^0 - ( F_021 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[59] = etfac[2] * AUX_INT__d_s_f_s[39] + 3 * one_over_2q * AUX_INT__d_s_d_s[23] - p_over_q * AUX_INT__f_s_f_s[79];

                    // ( D_011 S_000 | G_400 S_000 )^0_{t} = x * ( D_011 S_000 | F_300 S_000 )^0 + ( D_011 S_000 | D_200 S_000 )^0 - ( F_111 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[60] = etfac[0] * AUX_INT__d_s_f_s[40] + 3 * one_over_2q * AUX_INT__d_s_d_s[24] - p_over_q * AUX_INT__f_s_f_s[40];

                    // ( D_011 S_000 | G_310 S_000 )^0_{t} = y * ( D_011 S_000 | F_300 S_000 )^0 + ( P_001 S_000 | F_300 S_000 )^0 - ( F_021 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[61] = etfac[1] * AUX_INT__d_s_f_s[40] + 1 * one_over_2q * AUX_INT__p_s_f_s[20] - p_over_q * AUX_INT__f_s_f_s[70];

                    // ( D_011 S_000 | G_301 S_000 )^0_{t} = z * ( D_011 S_000 | F_300 S_000 )^0 + ( P_010 S_000 | F_300 S_000 )^0 - ( F_012 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[62] = etfac[2] * AUX_INT__d_s_f_s[40] + 1 * one_over_2q * AUX_INT__p_s_f_s[10] - p_over_q * AUX_INT__f_s_f_s[80];

                    // ( D_011 S_000 | G_220 S_000 )^0_{t} = y * ( D_011 S_000 | F_210 S_000 )^0 + ( P_001 S_000 | F_210 S_000 )^0 + ( D_011 S_000 | D_200 S_000 )^0 - ( F_021 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[63] = etfac[1] * AUX_INT__d_s_f_s[41] + 1 * one_over_2q * AUX_INT__p_s_f_s[21] + 1 * one_over_2q * AUX_INT__d_s_d_s[24] - p_over_q * AUX_INT__f_s_f_s[71];

                    // ( D_011 S_000 | G_211 S_000 )^0_{t} = z * ( D_011 S_000 | F_210 S_000 )^0 + ( P_010 S_000 | F_210 S_000 )^0 - ( F_012 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[64] = etfac[2] * AUX_INT__d_s_f_s[41] + 1 * one_over_2q * AUX_INT__p_s_f_s[11] - p_over_q * AUX_INT__f_s_f_s[81];

                    // ( D_011 S_000 | G_202 S_000 )^0_{t} = z * ( D_011 S_000 | F_201 S_000 )^0 + ( P_010 S_000 | F_201 S_000 )^0 + ( D_011 S_000 | D_200 S_000 )^0 - ( F_012 S_000 | F_201 S_000 )^0
                    AUX_INT__d_s_g_s[65] = etfac[2] * AUX_INT__d_s_f_s[42] + 1 * one_over_2q * AUX_INT__p_s_f_s[12] + 1 * one_over_2q * AUX_INT__d_s_d_s[24] - p_over_q * AUX_INT__f_s_f_s[82];

                    // ( D_011 S_000 | G_130 S_000 )^0_{t} = x * ( D_011 S_000 | F_030 S_000 )^0 - ( F_111 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[66] = etfac[0] * AUX_INT__d_s_f_s[46] - p_over_q * AUX_INT__f_s_f_s[46];

                    // ( D_011 S_000 | G_121 S_000 )^0_{t} = z * ( D_011 S_000 | F_120 S_000 )^0 + ( P_010 S_000 | F_120 S_000 )^0 - ( F_012 S_000 | F_120 S_000 )^0
                    AUX_INT__d_s_g_s[67] = etfac[2] * AUX_INT__d_s_f_s[43] + 1 * one_over_2q * AUX_INT__p_s_f_s[13] - p_over_q * AUX_INT__f_s_f_s[83];

                    // ( D_011 S_000 | G_112 S_000 )^0_{t} = y * ( D_011 S_000 | F_102 S_000 )^0 + ( P_001 S_000 | F_102 S_000 )^0 - ( F_021 S_000 | F_102 S_000 )^0
                    AUX_INT__d_s_g_s[68] = etfac[1] * AUX_INT__d_s_f_s[45] + 1 * one_over_2q * AUX_INT__p_s_f_s[25] - p_over_q * AUX_INT__f_s_f_s[75];

                    // ( D_011 S_000 | G_103 S_000 )^0_{t} = x * ( D_011 S_000 | F_003 S_000 )^0 - ( F_111 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[69] = etfac[0] * AUX_INT__d_s_f_s[49] - p_over_q * AUX_INT__f_s_f_s[49];

                    // ( D_011 S_000 | G_040 S_000 )^0_{t} = y * ( D_011 S_000 | F_030 S_000 )^0 + ( P_001 S_000 | F_030 S_000 )^0 + ( D_011 S_000 | D_020 S_000 )^0 - ( F_021 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[70] = etfac[1] * AUX_INT__d_s_f_s[46] + 1 * one_over_2q * AUX_INT__p_s_f_s[26] + 3 * one_over_2q * AUX_INT__d_s_d_s[27] - p_over_q * AUX_INT__f_s_f_s[76];

                    // ( D_011 S_000 | G_031 S_000 )^0_{t} = z * ( D_011 S_000 | F_030 S_000 )^0 + ( P_010 S_000 | F_030 S_000 )^0 - ( F_012 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[71] = etfac[2] * AUX_INT__d_s_f_s[46] + 1 * one_over_2q * AUX_INT__p_s_f_s[16] - p_over_q * AUX_INT__f_s_f_s[86];

                    // ( D_011 S_000 | G_022 S_000 )^0_{t} = z * ( D_011 S_000 | F_021 S_000 )^0 + ( P_010 S_000 | F_021 S_000 )^0 + ( D_011 S_000 | D_020 S_000 )^0 - ( F_012 S_000 | F_021 S_000 )^0
                    AUX_INT__d_s_g_s[72] = etfac[2] * AUX_INT__d_s_f_s[47] + 1 * one_over_2q * AUX_INT__p_s_f_s[17] + 1 * one_over_2q * AUX_INT__d_s_d_s[27] - p_over_q * AUX_INT__f_s_f_s[87];

                    // ( D_011 S_000 | G_013 S_000 )^0_{t} = y * ( D_011 S_000 | F_003 S_000 )^0 + ( P_001 S_000 | F_003 S_000 )^0 - ( F_021 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[73] = etfac[1] * AUX_INT__d_s_f_s[49] + 1 * one_over_2q * AUX_INT__p_s_f_s[29] - p_over_q * AUX_INT__f_s_f_s[79];

                    // ( D_011 S_000 | G_004 S_000 )^0_{t} = z * ( D_011 S_000 | F_003 S_000 )^0 + ( P_010 S_000 | F_003 S_000 )^0 + ( D_011 S_000 | D_002 S_000 )^0 - ( F_012 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[74] = etfac[2] * AUX_INT__d_s_f_s[49] + 1 * one_over_2q * AUX_INT__p_s_f_s[19] + 3 * one_over_2q * AUX_INT__d_s_d_s[29] - p_over_q * AUX_INT__f_s_f_s[89];

                    // ( D_002 S_000 | G_400 S_000 )^0_{t} = x * ( D_002 S_000 | F_300 S_000 )^0 + ( D_002 S_000 | D_200 S_000 )^0 - ( F_102 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[75] = etfac[0] * AUX_INT__d_s_f_s[50] + 3 * one_over_2q * AUX_INT__d_s_d_s[30] - p_over_q * AUX_INT__f_s_f_s[50];

                    // ( D_002 S_000 | G_310 S_000 )^0_{t} = y * ( D_002 S_000 | F_300 S_000 )^0 - ( F_012 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[76] = etfac[1] * AUX_INT__d_s_f_s[50] - p_over_q * AUX_INT__f_s_f_s[80];

                    // ( D_002 S_000 | G_301 S_000 )^0_{t} = z * ( D_002 S_000 | F_300 S_000 )^0 + ( P_001 S_000 | F_300 S_000 )^0 - ( F_003 S_000 | F_300 S_000 )^0
                    AUX_INT__d_s_g_s[77] = etfac[2] * AUX_INT__d_s_f_s[50] + 2 * one_over_2q * AUX_INT__p_s_f_s[20] - p_over_q * AUX_INT__f_s_f_s[90];

                    // ( D_002 S_000 | G_220 S_000 )^0_{t} = y * ( D_002 S_000 | F_210 S_000 )^0 + ( D_002 S_000 | D_200 S_000 )^0 - ( F_012 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[78] = etfac[1] * AUX_INT__d_s_f_s[51] + 1 * one_over_2q * AUX_INT__d_s_d_s[30] - p_over_q * AUX_INT__f_s_f_s[81];

                    // ( D_002 S_000 | G_211 S_000 )^0_{t} = z * ( D_002 S_000 | F_210 S_000 )^0 + ( P_001 S_000 | F_210 S_000 )^0 - ( F_003 S_000 | F_210 S_000 )^0
                    AUX_INT__d_s_g_s[79] = etfac[2] * AUX_INT__d_s_f_s[51] + 2 * one_over_2q * AUX_INT__p_s_f_s[21] - p_over_q * AUX_INT__f_s_f_s[91];

                    // ( D_002 S_000 | G_202 S_000 )^0_{t} = z * ( D_002 S_000 | F_201 S_000 )^0 + ( P_001 S_000 | F_201 S_000 )^0 + ( D_002 S_000 | D_200 S_000 )^0 - ( F_003 S_000 | F_201 S_000 )^0
                    AUX_INT__d_s_g_s[80] = etfac[2] * AUX_INT__d_s_f_s[52] + 2 * one_over_2q * AUX_INT__p_s_f_s[22] + 1 * one_over_2q * AUX_INT__d_s_d_s[30] - p_over_q * AUX_INT__f_s_f_s[92];

                    // ( D_002 S_000 | G_130 S_000 )^0_{t} = x * ( D_002 S_000 | F_030 S_000 )^0 - ( F_102 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[81] = etfac[0] * AUX_INT__d_s_f_s[56] - p_over_q * AUX_INT__f_s_f_s[56];

                    // ( D_002 S_000 | G_121 S_000 )^0_{t} = z * ( D_002 S_000 | F_120 S_000 )^0 + ( P_001 S_000 | F_120 S_000 )^0 - ( F_003 S_000 | F_120 S_000 )^0
                    AUX_INT__d_s_g_s[82] = etfac[2] * AUX_INT__d_s_f_s[53] + 2 * one_over_2q * AUX_INT__p_s_f_s[23] - p_over_q * AUX_INT__f_s_f_s[93];

                    // ( D_002 S_000 | G_112 S_000 )^0_{t} = y * ( D_002 S_000 | F_102 S_000 )^0 - ( F_012 S_000 | F_102 S_000 )^0
                    AUX_INT__d_s_g_s[83] = etfac[1] * AUX_INT__d_s_f_s[55] - p_over_q * AUX_INT__f_s_f_s[85];

                    // ( D_002 S_000 | G_103 S_000 )^0_{t} = x * ( D_002 S_000 | F_003 S_000 )^0 - ( F_102 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[84] = etfac[0] * AUX_INT__d_s_f_s[59] - p_over_q * AUX_INT__f_s_f_s[59];

                    // ( D_002 S_000 | G_040 S_000 )^0_{t} = y * ( D_002 S_000 | F_030 S_000 )^0 + ( D_002 S_000 | D_020 S_000 )^0 - ( F_012 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[85] = etfac[1] * AUX_INT__d_s_f_s[56] + 3 * one_over_2q * AUX_INT__d_s_d_s[33] - p_over_q * AUX_INT__f_s_f_s[86];

                    // ( D_002 S_000 | G_031 S_000 )^0_{t} = z * ( D_002 S_000 | F_030 S_000 )^0 + ( P_001 S_000 | F_030 S_000 )^0 - ( F_003 S_000 | F_030 S_000 )^0
                    AUX_INT__d_s_g_s[86] = etfac[2] * AUX_INT__d_s_f_s[56] + 2 * one_over_2q * AUX_INT__p_s_f_s[26] - p_over_q * AUX_INT__f_s_f_s[96];

                    // ( D_002 S_000 | G_022 S_000 )^0_{t} = z * ( D_002 S_000 | F_021 S_000 )^0 + ( P_001 S_000 | F_021 S_000 )^0 + ( D_002 S_000 | D_020 S_000 )^0 - ( F_003 S_000 | F_021 S_000 )^0
                    AUX_INT__d_s_g_s[87] = etfac[2] * AUX_INT__d_s_f_s[57] + 2 * one_over_2q * AUX_INT__p_s_f_s[27] + 1 * one_over_2q * AUX_INT__d_s_d_s[33] - p_over_q * AUX_INT__f_s_f_s[97];

                    // ( D_002 S_000 | G_013 S_000 )^0_{t} = y * ( D_002 S_000 | F_003 S_000 )^0 - ( F_012 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[88] = etfac[1] * AUX_INT__d_s_f_s[59] - p_over_q * AUX_INT__f_s_f_s[89];

                    // ( D_002 S_000 | G_004 S_000 )^0_{t} = z * ( D_002 S_000 | F_003 S_000 )^0 + ( P_001 S_000 | F_003 S_000 )^0 + ( D_002 S_000 | D_002 S_000 )^0 - ( F_003 S_000 | F_003 S_000 )^0
                    AUX_INT__d_s_g_s[89] = etfac[2] * AUX_INT__d_s_f_s[59] + 2 * one_over_2q * AUX_INT__p_s_f_s[29] + 3 * one_over_2q * AUX_INT__d_s_d_s[35] - p_over_q * AUX_INT__f_s_f_s[99];

                    // ( I_600 S_000 | P_100 S_000 )^0 = x * ( I_600 S_000 | S_000 S_000 )^0_{e} + ( H_500 S_000 | S_000 S_000 )^0_{e} - ( J_700 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[0] = etfac[0] * AUX_INT__i_s_s_s[0] + 6 * one_over_2q * AUX_INT__h_s_s_s[0] - p_over_q * AUX_INT__j_s_s_s[0];

                    // ( I_510 S_000 | P_100 S_000 )^0 = x * ( I_510 S_000 | S_000 S_000 )^0_{e} + ( H_410 S_000 | S_000 S_000 )^0_{e} - ( J_610 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[3] = etfac[0] * AUX_INT__i_s_s_s[1] + 5 * one_over_2q * AUX_INT__h_s_s_s[1] - p_over_q * AUX_INT__j_s_s_s[1];

                    // ( I_510 S_000 | P_010 S_000 )^0 = y * ( I_510 S_000 | S_000 S_000 )^0_{e} + ( H_500 S_000 | S_000 S_000 )^0_{e} - ( J_520 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[4] = etfac[1] * AUX_INT__i_s_s_s[1] + 1 * one_over_2q * AUX_INT__h_s_s_s[0] - p_over_q * AUX_INT__j_s_s_s[3];

                    // ( I_501 S_000 | P_100 S_000 )^0 = x * ( I_501 S_000 | S_000 S_000 )^0_{e} + ( H_401 S_000 | S_000 S_000 )^0_{e} - ( J_601 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[6] = etfac[0] * AUX_INT__i_s_s_s[2] + 5 * one_over_2q * AUX_INT__h_s_s_s[2] - p_over_q * AUX_INT__j_s_s_s[2];

                    // ( I_501 S_000 | P_001 S_000 )^0 = z * ( I_501 S_000 | S_000 S_000 )^0_{e} + ( H_500 S_000 | S_000 S_000 )^0_{e} - ( J_502 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[8] = etfac[2] * AUX_INT__i_s_s_s[2] + 1 * one_over_2q * AUX_INT__h_s_s_s[0] - p_over_q * AUX_INT__j_s_s_s[5];

                    // ( I_420 S_000 | P_100 S_000 )^0 = x * ( I_420 S_000 | S_000 S_000 )^0_{e} + ( H_320 S_000 | S_000 S_000 )^0_{e} - ( J_520 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[9] = etfac[0] * AUX_INT__i_s_s_s[3] + 4 * one_over_2q * AUX_INT__h_s_s_s[3] - p_over_q * AUX_INT__j_s_s_s[3];

                    // ( I_420 S_000 | P_010 S_000 )^0 = y * ( I_420 S_000 | S_000 S_000 )^0_{e} + ( H_410 S_000 | S_000 S_000 )^0_{e} - ( J_430 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[10] = etfac[1] * AUX_INT__i_s_s_s[3] + 2 * one_over_2q * AUX_INT__h_s_s_s[1] - p_over_q * AUX_INT__j_s_s_s[6];

                    // ( I_411 S_000 | P_100 S_000 )^0 = x * ( I_411 S_000 | S_000 S_000 )^0_{e} + ( H_311 S_000 | S_000 S_000 )^0_{e} - ( J_511 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[12] = etfac[0] * AUX_INT__i_s_s_s[4] + 4 * one_over_2q * AUX_INT__h_s_s_s[4] - p_over_q * AUX_INT__j_s_s_s[4];

                    // ( I_411 S_000 | P_010 S_000 )^0 = y * ( I_411 S_000 | S_000 S_000 )^0_{e} + ( H_401 S_000 | S_000 S_000 )^0_{e} - ( J_421 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[13] = etfac[1] * AUX_INT__i_s_s_s[4] + 1 * one_over_2q * AUX_INT__h_s_s_s[2] - p_over_q * AUX_INT__j_s_s_s[7];

                    // ( I_411 S_000 | P_001 S_000 )^0 = z * ( I_411 S_000 | S_000 S_000 )^0_{e} + ( H_410 S_000 | S_000 S_000 )^0_{e} - ( J_412 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[14] = etfac[2] * AUX_INT__i_s_s_s[4] + 1 * one_over_2q * AUX_INT__h_s_s_s[1] - p_over_q * AUX_INT__j_s_s_s[8];

                    // ( I_402 S_000 | P_100 S_000 )^0 = x * ( I_402 S_000 | S_000 S_000 )^0_{e} + ( H_302 S_000 | S_000 S_000 )^0_{e} - ( J_502 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[15] = etfac[0] * AUX_INT__i_s_s_s[5] + 4 * one_over_2q * AUX_INT__h_s_s_s[5] - p_over_q * AUX_INT__j_s_s_s[5];

                    // ( I_402 S_000 | P_001 S_000 )^0 = z * ( I_402 S_000 | S_000 S_000 )^0_{e} + ( H_401 S_000 | S_000 S_000 )^0_{e} - ( J_403 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[17] = etfac[2] * AUX_INT__i_s_s_s[5] + 2 * one_over_2q * AUX_INT__h_s_s_s[2] - p_over_q * AUX_INT__j_s_s_s[9];

                    // ( I_330 S_000 | P_100 S_000 )^0 = x * ( I_330 S_000 | S_000 S_000 )^0_{e} + ( H_230 S_000 | S_000 S_000 )^0_{e} - ( J_430 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[18] = etfac[0] * AUX_INT__i_s_s_s[6] + 3 * one_over_2q * AUX_INT__h_s_s_s[6] - p_over_q * AUX_INT__j_s_s_s[6];

                    // ( I_330 S_000 | P_010 S_000 )^0 = y * ( I_330 S_000 | S_000 S_000 )^0_{e} + ( H_320 S_000 | S_000 S_000 )^0_{e} - ( J_340 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[19] = etfac[1] * AUX_INT__i_s_s_s[6] + 3 * one_over_2q * AUX_INT__h_s_s_s[3] - p_over_q * AUX_INT__j_s_s_s[10];

                    // ( I_321 S_000 | P_100 S_000 )^0 = x * ( I_321 S_000 | S_000 S_000 )^0_{e} + ( H_221 S_000 | S_000 S_000 )^0_{e} - ( J_421 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[21] = etfac[0] * AUX_INT__i_s_s_s[7] + 3 * one_over_2q * AUX_INT__h_s_s_s[7] - p_over_q * AUX_INT__j_s_s_s[7];

                    // ( I_321 S_000 | P_010 S_000 )^0 = y * ( I_321 S_000 | S_000 S_000 )^0_{e} + ( H_311 S_000 | S_000 S_000 )^0_{e} - ( J_331 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[22] = etfac[1] * AUX_INT__i_s_s_s[7] + 2 * one_over_2q * AUX_INT__h_s_s_s[4] - p_over_q * AUX_INT__j_s_s_s[11];

                    // ( I_321 S_000 | P_001 S_000 )^0 = z * ( I_321 S_000 | S_000 S_000 )^0_{e} + ( H_320 S_000 | S_000 S_000 )^0_{e} - ( J_322 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[23] = etfac[2] * AUX_INT__i_s_s_s[7] + 1 * one_over_2q * AUX_INT__h_s_s_s[3] - p_over_q * AUX_INT__j_s_s_s[12];

                    // ( I_312 S_000 | P_100 S_000 )^0 = x * ( I_312 S_000 | S_000 S_000 )^0_{e} + ( H_212 S_000 | S_000 S_000 )^0_{e} - ( J_412 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[24] = etfac[0] * AUX_INT__i_s_s_s[8] + 3 * one_over_2q * AUX_INT__h_s_s_s[8] - p_over_q * AUX_INT__j_s_s_s[8];

                    // ( I_312 S_000 | P_010 S_000 )^0 = y * ( I_312 S_000 | S_000 S_000 )^0_{e} + ( H_302 S_000 | S_000 S_000 )^0_{e} - ( J_322 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[25] = etfac[1] * AUX_INT__i_s_s_s[8] + 1 * one_over_2q * AUX_INT__h_s_s_s[5] - p_over_q * AUX_INT__j_s_s_s[12];

                    // ( I_312 S_000 | P_001 S_000 )^0 = z * ( I_312 S_000 | S_000 S_000 )^0_{e} + ( H_311 S_000 | S_000 S_000 )^0_{e} - ( J_313 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[26] = etfac[2] * AUX_INT__i_s_s_s[8] + 2 * one_over_2q * AUX_INT__h_s_s_s[4] - p_over_q * AUX_INT__j_s_s_s[13];

                    // ( I_303 S_000 | P_100 S_000 )^0 = x * ( I_303 S_000 | S_000 S_000 )^0_{e} + ( H_203 S_000 | S_000 S_000 )^0_{e} - ( J_403 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[27] = etfac[0] * AUX_INT__i_s_s_s[9] + 3 * one_over_2q * AUX_INT__h_s_s_s[9] - p_over_q * AUX_INT__j_s_s_s[9];

                    // ( I_303 S_000 | P_001 S_000 )^0 = z * ( I_303 S_000 | S_000 S_000 )^0_{e} + ( H_302 S_000 | S_000 S_000 )^0_{e} - ( J_304 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[29] = etfac[2] * AUX_INT__i_s_s_s[9] + 3 * one_over_2q * AUX_INT__h_s_s_s[5] - p_over_q * AUX_INT__j_s_s_s[14];

                    // ( I_240 S_000 | P_100 S_000 )^0 = x * ( I_240 S_000 | S_000 S_000 )^0_{e} + ( H_140 S_000 | S_000 S_000 )^0_{e} - ( J_340 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[30] = etfac[0] * AUX_INT__i_s_s_s[10] + 2 * one_over_2q * AUX_INT__h_s_s_s[10] - p_over_q * AUX_INT__j_s_s_s[10];

                    // ( I_240 S_000 | P_010 S_000 )^0 = y * ( I_240 S_000 | S_000 S_000 )^0_{e} + ( H_230 S_000 | S_000 S_000 )^0_{e} - ( J_250 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[31] = etfac[1] * AUX_INT__i_s_s_s[10] + 4 * one_over_2q * AUX_INT__h_s_s_s[6] - p_over_q * AUX_INT__j_s_s_s[15];

                    // ( I_231 S_000 | P_100 S_000 )^0 = x * ( I_231 S_000 | S_000 S_000 )^0_{e} + ( H_131 S_000 | S_000 S_000 )^0_{e} - ( J_331 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[33] = etfac[0] * AUX_INT__i_s_s_s[11] + 2 * one_over_2q * AUX_INT__h_s_s_s[11] - p_over_q * AUX_INT__j_s_s_s[11];

                    // ( I_231 S_000 | P_010 S_000 )^0 = y * ( I_231 S_000 | S_000 S_000 )^0_{e} + ( H_221 S_000 | S_000 S_000 )^0_{e} - ( J_241 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[34] = etfac[1] * AUX_INT__i_s_s_s[11] + 3 * one_over_2q * AUX_INT__h_s_s_s[7] - p_over_q * AUX_INT__j_s_s_s[16];

                    // ( I_231 S_000 | P_001 S_000 )^0 = z * ( I_231 S_000 | S_000 S_000 )^0_{e} + ( H_230 S_000 | S_000 S_000 )^0_{e} - ( J_232 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[35] = etfac[2] * AUX_INT__i_s_s_s[11] + 1 * one_over_2q * AUX_INT__h_s_s_s[6] - p_over_q * AUX_INT__j_s_s_s[17];

                    // ( I_222 S_000 | P_100 S_000 )^0 = x * ( I_222 S_000 | S_000 S_000 )^0_{e} + ( H_122 S_000 | S_000 S_000 )^0_{e} - ( J_322 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[36] = etfac[0] * AUX_INT__i_s_s_s[12] + 2 * one_over_2q * AUX_INT__h_s_s_s[12] - p_over_q * AUX_INT__j_s_s_s[12];

                    // ( I_222 S_000 | P_010 S_000 )^0 = y * ( I_222 S_000 | S_000 S_000 )^0_{e} + ( H_212 S_000 | S_000 S_000 )^0_{e} - ( J_232 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[37] = etfac[1] * AUX_INT__i_s_s_s[12] + 2 * one_over_2q * AUX_INT__h_s_s_s[8] - p_over_q * AUX_INT__j_s_s_s[17];

                    // ( I_222 S_000 | P_001 S_000 )^0 = z * ( I_222 S_000 | S_000 S_000 )^0_{e} + ( H_221 S_000 | S_000 S_000 )^0_{e} - ( J_223 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[38] = etfac[2] * AUX_INT__i_s_s_s[12] + 2 * one_over_2q * AUX_INT__h_s_s_s[7] - p_over_q * AUX_INT__j_s_s_s[18];

                    // ( I_213 S_000 | P_100 S_000 )^0 = x * ( I_213 S_000 | S_000 S_000 )^0_{e} + ( H_113 S_000 | S_000 S_000 )^0_{e} - ( J_313 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[39] = etfac[0] * AUX_INT__i_s_s_s[13] + 2 * one_over_2q * AUX_INT__h_s_s_s[13] - p_over_q * AUX_INT__j_s_s_s[13];

                    // ( I_213 S_000 | P_010 S_000 )^0 = y * ( I_213 S_000 | S_000 S_000 )^0_{e} + ( H_203 S_000 | S_000 S_000 )^0_{e} - ( J_223 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[40] = etfac[1] * AUX_INT__i_s_s_s[13] + 1 * one_over_2q * AUX_INT__h_s_s_s[9] - p_over_q * AUX_INT__j_s_s_s[18];

                    // ( I_213 S_000 | P_001 S_000 )^0 = z * ( I_213 S_000 | S_000 S_000 )^0_{e} + ( H_212 S_000 | S_000 S_000 )^0_{e} - ( J_214 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[41] = etfac[2] * AUX_INT__i_s_s_s[13] + 3 * one_over_2q * AUX_INT__h_s_s_s[8] - p_over_q * AUX_INT__j_s_s_s[19];

                    // ( I_204 S_000 | P_100 S_000 )^0 = x * ( I_204 S_000 | S_000 S_000 )^0_{e} + ( H_104 S_000 | S_000 S_000 )^0_{e} - ( J_304 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[42] = etfac[0] * AUX_INT__i_s_s_s[14] + 2 * one_over_2q * AUX_INT__h_s_s_s[14] - p_over_q * AUX_INT__j_s_s_s[14];

                    // ( I_204 S_000 | P_001 S_000 )^0 = z * ( I_204 S_000 | S_000 S_000 )^0_{e} + ( H_203 S_000 | S_000 S_000 )^0_{e} - ( J_205 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[44] = etfac[2] * AUX_INT__i_s_s_s[14] + 4 * one_over_2q * AUX_INT__h_s_s_s[9] - p_over_q * AUX_INT__j_s_s_s[20];

                    // ( I_150 S_000 | P_100 S_000 )^0 = x * ( I_150 S_000 | S_000 S_000 )^0_{e} + ( H_050 S_000 | S_000 S_000 )^0_{e} - ( J_250 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[45] = etfac[0] * AUX_INT__i_s_s_s[15] + 1 * one_over_2q * AUX_INT__h_s_s_s[15] - p_over_q * AUX_INT__j_s_s_s[15];

                    // ( I_150 S_000 | P_010 S_000 )^0 = y * ( I_150 S_000 | S_000 S_000 )^0_{e} + ( H_140 S_000 | S_000 S_000 )^0_{e} - ( J_160 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[46] = etfac[1] * AUX_INT__i_s_s_s[15] + 5 * one_over_2q * AUX_INT__h_s_s_s[10] - p_over_q * AUX_INT__j_s_s_s[21];

                    // ( I_141 S_000 | P_100 S_000 )^0 = x * ( I_141 S_000 | S_000 S_000 )^0_{e} + ( H_041 S_000 | S_000 S_000 )^0_{e} - ( J_241 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[48] = etfac[0] * AUX_INT__i_s_s_s[16] + 1 * one_over_2q * AUX_INT__h_s_s_s[16] - p_over_q * AUX_INT__j_s_s_s[16];

                    // ( I_141 S_000 | P_010 S_000 )^0 = y * ( I_141 S_000 | S_000 S_000 )^0_{e} + ( H_131 S_000 | S_000 S_000 )^0_{e} - ( J_151 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[49] = etfac[1] * AUX_INT__i_s_s_s[16] + 4 * one_over_2q * AUX_INT__h_s_s_s[11] - p_over_q * AUX_INT__j_s_s_s[22];

                    // ( I_141 S_000 | P_001 S_000 )^0 = z * ( I_141 S_000 | S_000 S_000 )^0_{e} + ( H_140 S_000 | S_000 S_000 )^0_{e} - ( J_142 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[50] = etfac[2] * AUX_INT__i_s_s_s[16] + 1 * one_over_2q * AUX_INT__h_s_s_s[10] - p_over_q * AUX_INT__j_s_s_s[23];

                    // ( I_132 S_000 | P_100 S_000 )^0 = x * ( I_132 S_000 | S_000 S_000 )^0_{e} + ( H_032 S_000 | S_000 S_000 )^0_{e} - ( J_232 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[51] = etfac[0] * AUX_INT__i_s_s_s[17] + 1 * one_over_2q * AUX_INT__h_s_s_s[17] - p_over_q * AUX_INT__j_s_s_s[17];

                    // ( I_132 S_000 | P_010 S_000 )^0 = y * ( I_132 S_000 | S_000 S_000 )^0_{e} + ( H_122 S_000 | S_000 S_000 )^0_{e} - ( J_142 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[52] = etfac[1] * AUX_INT__i_s_s_s[17] + 3 * one_over_2q * AUX_INT__h_s_s_s[12] - p_over_q * AUX_INT__j_s_s_s[23];

                    // ( I_132 S_000 | P_001 S_000 )^0 = z * ( I_132 S_000 | S_000 S_000 )^0_{e} + ( H_131 S_000 | S_000 S_000 )^0_{e} - ( J_133 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[53] = etfac[2] * AUX_INT__i_s_s_s[17] + 2 * one_over_2q * AUX_INT__h_s_s_s[11] - p_over_q * AUX_INT__j_s_s_s[24];

                    // ( I_123 S_000 | P_100 S_000 )^0 = x * ( I_123 S_000 | S_000 S_000 )^0_{e} + ( H_023 S_000 | S_000 S_000 )^0_{e} - ( J_223 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[54] = etfac[0] * AUX_INT__i_s_s_s[18] + 1 * one_over_2q * AUX_INT__h_s_s_s[18] - p_over_q * AUX_INT__j_s_s_s[18];

                    // ( I_123 S_000 | P_010 S_000 )^0 = y * ( I_123 S_000 | S_000 S_000 )^0_{e} + ( H_113 S_000 | S_000 S_000 )^0_{e} - ( J_133 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[55] = etfac[1] * AUX_INT__i_s_s_s[18] + 2 * one_over_2q * AUX_INT__h_s_s_s[13] - p_over_q * AUX_INT__j_s_s_s[24];

                    // ( I_123 S_000 | P_001 S_000 )^0 = z * ( I_123 S_000 | S_000 S_000 )^0_{e} + ( H_122 S_000 | S_000 S_000 )^0_{e} - ( J_124 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[56] = etfac[2] * AUX_INT__i_s_s_s[18] + 3 * one_over_2q * AUX_INT__h_s_s_s[12] - p_over_q * AUX_INT__j_s_s_s[25];

                    // ( I_114 S_000 | P_100 S_000 )^0 = x * ( I_114 S_000 | S_000 S_000 )^0_{e} + ( H_014 S_000 | S_000 S_000 )^0_{e} - ( J_214 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[57] = etfac[0] * AUX_INT__i_s_s_s[19] + 1 * one_over_2q * AUX_INT__h_s_s_s[19] - p_over_q * AUX_INT__j_s_s_s[19];

                    // ( I_114 S_000 | P_010 S_000 )^0 = y * ( I_114 S_000 | S_000 S_000 )^0_{e} + ( H_104 S_000 | S_000 S_000 )^0_{e} - ( J_124 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[58] = etfac[1] * AUX_INT__i_s_s_s[19] + 1 * one_over_2q * AUX_INT__h_s_s_s[14] - p_over_q * AUX_INT__j_s_s_s[25];

                    // ( I_114 S_000 | P_001 S_000 )^0 = z * ( I_114 S_000 | S_000 S_000 )^0_{e} + ( H_113 S_000 | S_000 S_000 )^0_{e} - ( J_115 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[59] = etfac[2] * AUX_INT__i_s_s_s[19] + 4 * one_over_2q * AUX_INT__h_s_s_s[13] - p_over_q * AUX_INT__j_s_s_s[26];

                    // ( I_105 S_000 | P_100 S_000 )^0 = x * ( I_105 S_000 | S_000 S_000 )^0_{e} + ( H_005 S_000 | S_000 S_000 )^0_{e} - ( J_205 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[60] = etfac[0] * AUX_INT__i_s_s_s[20] + 1 * one_over_2q * AUX_INT__h_s_s_s[20] - p_over_q * AUX_INT__j_s_s_s[20];

                    // ( I_105 S_000 | P_001 S_000 )^0 = z * ( I_105 S_000 | S_000 S_000 )^0_{e} + ( H_104 S_000 | S_000 S_000 )^0_{e} - ( J_106 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[62] = etfac[2] * AUX_INT__i_s_s_s[20] + 5 * one_over_2q * AUX_INT__h_s_s_s[14] - p_over_q * AUX_INT__j_s_s_s[27];

                    // ( I_060 S_000 | P_010 S_000 )^0 = y * ( I_060 S_000 | S_000 S_000 )^0_{e} + ( H_050 S_000 | S_000 S_000 )^0_{e} - ( J_070 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[64] = etfac[1] * AUX_INT__i_s_s_s[21] + 6 * one_over_2q * AUX_INT__h_s_s_s[15] - p_over_q * AUX_INT__j_s_s_s[28];

                    // ( I_051 S_000 | P_100 S_000 )^0 = x * ( I_051 S_000 | S_000 S_000 )^0_{e} - ( J_151 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[66] = etfac[0] * AUX_INT__i_s_s_s[22] - p_over_q * AUX_INT__j_s_s_s[22];

                    // ( I_051 S_000 | P_010 S_000 )^0 = y * ( I_051 S_000 | S_000 S_000 )^0_{e} + ( H_041 S_000 | S_000 S_000 )^0_{e} - ( J_061 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[67] = etfac[1] * AUX_INT__i_s_s_s[22] + 5 * one_over_2q * AUX_INT__h_s_s_s[16] - p_over_q * AUX_INT__j_s_s_s[29];

                    // ( I_051 S_000 | P_001 S_000 )^0 = z * ( I_051 S_000 | S_000 S_000 )^0_{e} + ( H_050 S_000 | S_000 S_000 )^0_{e} - ( J_052 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[68] = etfac[2] * AUX_INT__i_s_s_s[22] + 1 * one_over_2q * AUX_INT__h_s_s_s[15] - p_over_q * AUX_INT__j_s_s_s[30];

                    // ( I_042 S_000 | P_100 S_000 )^0 = x * ( I_042 S_000 | S_000 S_000 )^0_{e} - ( J_142 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[69] = etfac[0] * AUX_INT__i_s_s_s[23] - p_over_q * AUX_INT__j_s_s_s[23];

                    // ( I_042 S_000 | P_010 S_000 )^0 = y * ( I_042 S_000 | S_000 S_000 )^0_{e} + ( H_032 S_000 | S_000 S_000 )^0_{e} - ( J_052 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[70] = etfac[1] * AUX_INT__i_s_s_s[23] + 4 * one_over_2q * AUX_INT__h_s_s_s[17] - p_over_q * AUX_INT__j_s_s_s[30];

                    // ( I_042 S_000 | P_001 S_000 )^0 = z * ( I_042 S_000 | S_000 S_000 )^0_{e} + ( H_041 S_000 | S_000 S_000 )^0_{e} - ( J_043 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[71] = etfac[2] * AUX_INT__i_s_s_s[23] + 2 * one_over_2q * AUX_INT__h_s_s_s[16] - p_over_q * AUX_INT__j_s_s_s[31];

                    // ( I_033 S_000 | P_100 S_000 )^0 = x * ( I_033 S_000 | S_000 S_000 )^0_{e} - ( J_133 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[72] = etfac[0] * AUX_INT__i_s_s_s[24] - p_over_q * AUX_INT__j_s_s_s[24];

                    // ( I_033 S_000 | P_010 S_000 )^0 = y * ( I_033 S_000 | S_000 S_000 )^0_{e} + ( H_023 S_000 | S_000 S_000 )^0_{e} - ( J_043 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[73] = etfac[1] * AUX_INT__i_s_s_s[24] + 3 * one_over_2q * AUX_INT__h_s_s_s[18] - p_over_q * AUX_INT__j_s_s_s[31];

                    // ( I_033 S_000 | P_001 S_000 )^0 = z * ( I_033 S_000 | S_000 S_000 )^0_{e} + ( H_032 S_000 | S_000 S_000 )^0_{e} - ( J_034 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[74] = etfac[2] * AUX_INT__i_s_s_s[24] + 3 * one_over_2q * AUX_INT__h_s_s_s[17] - p_over_q * AUX_INT__j_s_s_s[32];

                    // ( I_024 S_000 | P_100 S_000 )^0 = x * ( I_024 S_000 | S_000 S_000 )^0_{e} - ( J_124 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[75] = etfac[0] * AUX_INT__i_s_s_s[25] - p_over_q * AUX_INT__j_s_s_s[25];

                    // ( I_024 S_000 | P_010 S_000 )^0 = y * ( I_024 S_000 | S_000 S_000 )^0_{e} + ( H_014 S_000 | S_000 S_000 )^0_{e} - ( J_034 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[76] = etfac[1] * AUX_INT__i_s_s_s[25] + 2 * one_over_2q * AUX_INT__h_s_s_s[19] - p_over_q * AUX_INT__j_s_s_s[32];

                    // ( I_024 S_000 | P_001 S_000 )^0 = z * ( I_024 S_000 | S_000 S_000 )^0_{e} + ( H_023 S_000 | S_000 S_000 )^0_{e} - ( J_025 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[77] = etfac[2] * AUX_INT__i_s_s_s[25] + 4 * one_over_2q * AUX_INT__h_s_s_s[18] - p_over_q * AUX_INT__j_s_s_s[33];

                    // ( I_015 S_000 | P_100 S_000 )^0 = x * ( I_015 S_000 | S_000 S_000 )^0_{e} - ( J_115 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[78] = etfac[0] * AUX_INT__i_s_s_s[26] - p_over_q * AUX_INT__j_s_s_s[26];

                    // ( I_015 S_000 | P_010 S_000 )^0 = y * ( I_015 S_000 | S_000 S_000 )^0_{e} + ( H_005 S_000 | S_000 S_000 )^0_{e} - ( J_025 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[79] = etfac[1] * AUX_INT__i_s_s_s[26] + 1 * one_over_2q * AUX_INT__h_s_s_s[20] - p_over_q * AUX_INT__j_s_s_s[33];

                    // ( I_015 S_000 | P_001 S_000 )^0 = z * ( I_015 S_000 | S_000 S_000 )^0_{e} + ( H_014 S_000 | S_000 S_000 )^0_{e} - ( J_016 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[80] = etfac[2] * AUX_INT__i_s_s_s[26] + 5 * one_over_2q * AUX_INT__h_s_s_s[19] - p_over_q * AUX_INT__j_s_s_s[34];

                    // ( I_006 S_000 | P_001 S_000 )^0 = z * ( I_006 S_000 | S_000 S_000 )^0_{e} + ( H_005 S_000 | S_000 S_000 )^0_{e} - ( J_007 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[83] = etfac[2] * AUX_INT__i_s_s_s[27] + 6 * one_over_2q * AUX_INT__h_s_s_s[20] - p_over_q * AUX_INT__j_s_s_s[35];

                    // ( I_510 S_000 | P_001 S_000 )^0 = z * ( I_510 S_000 | S_000 S_000 )^0_{e} - ( J_511 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[5] = etfac[2] * AUX_INT__i_s_s_s[1] - p_over_q * AUX_INT__j_s_s_s[4];

                    // ( I_501 S_000 | P_010 S_000 )^0 = y * ( I_501 S_000 | S_000 S_000 )^0_{e} - ( J_511 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[7] = etfac[1] * AUX_INT__i_s_s_s[2] - p_over_q * AUX_INT__j_s_s_s[4];

                    // ( I_420 S_000 | P_001 S_000 )^0 = z * ( I_420 S_000 | S_000 S_000 )^0_{e} - ( J_421 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[11] = etfac[2] * AUX_INT__i_s_s_s[3] - p_over_q * AUX_INT__j_s_s_s[7];

                    // ( I_402 S_000 | P_010 S_000 )^0 = y * ( I_402 S_000 | S_000 S_000 )^0_{e} - ( J_412 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[16] = etfac[1] * AUX_INT__i_s_s_s[5] - p_over_q * AUX_INT__j_s_s_s[8];

                    // ( I_330 S_000 | P_001 S_000 )^0 = z * ( I_330 S_000 | S_000 S_000 )^0_{e} - ( J_331 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[20] = etfac[2] * AUX_INT__i_s_s_s[6] - p_over_q * AUX_INT__j_s_s_s[11];

                    // ( I_303 S_000 | P_010 S_000 )^0 = y * ( I_303 S_000 | S_000 S_000 )^0_{e} - ( J_313 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[28] = etfac[1] * AUX_INT__i_s_s_s[9] - p_over_q * AUX_INT__j_s_s_s[13];

                    // ( I_240 S_000 | P_001 S_000 )^0 = z * ( I_240 S_000 | S_000 S_000 )^0_{e} - ( J_241 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[32] = etfac[2] * AUX_INT__i_s_s_s[10] - p_over_q * AUX_INT__j_s_s_s[16];

                    // ( I_204 S_000 | P_010 S_000 )^0 = y * ( I_204 S_000 | S_000 S_000 )^0_{e} - ( J_214 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[43] = etfac[1] * AUX_INT__i_s_s_s[14] - p_over_q * AUX_INT__j_s_s_s[19];

                    // ( I_150 S_000 | P_001 S_000 )^0 = z * ( I_150 S_000 | S_000 S_000 )^0_{e} - ( J_151 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[47] = etfac[2] * AUX_INT__i_s_s_s[15] - p_over_q * AUX_INT__j_s_s_s[22];

                    // ( I_105 S_000 | P_010 S_000 )^0 = y * ( I_105 S_000 | S_000 S_000 )^0_{e} - ( J_115 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[61] = etfac[1] * AUX_INT__i_s_s_s[20] - p_over_q * AUX_INT__j_s_s_s[26];

                    // ( I_060 S_000 | P_100 S_000 )^0 = x * ( I_060 S_000 | S_000 S_000 )^0_{e} - ( J_160 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[63] = etfac[0] * AUX_INT__i_s_s_s[21] - p_over_q * AUX_INT__j_s_s_s[21];

                    // ( I_006 S_000 | P_100 S_000 )^0 = x * ( I_006 S_000 | S_000 S_000 )^0_{e} - ( J_106 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[81] = etfac[0] * AUX_INT__i_s_s_s[27] - p_over_q * AUX_INT__j_s_s_s[27];

                    // ( I_006 S_000 | P_010 S_000 )^0 = y * ( I_006 S_000 | S_000 S_000 )^0_{e} - ( J_016 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__i_s_p_s[82] = etfac[1] * AUX_INT__i_s_s_s[27] - p_over_q * AUX_INT__j_s_s_s[34];

                    // ( H_500 S_000 | D_200 S_000 )^0 = x * ( H_500 S_000 | P_100 S_000 )^0 + ( G_400 S_000 | P_100 S_000 )^0 + ( H_500 S_000 | S_000 S_000 )^0_{e} - ( I_600 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[0] = etfac[0] * AUX_INT__h_s_p_s[0] + 5 * one_over_2q * AUX_INT__g_s_p_s[0] + 1 * one_over_2q * AUX_INT__h_s_s_s[0] - p_over_q * AUX_INT__i_s_p_s[0];

                    // ( H_500 S_000 | D_020 S_000 )^0 = y * ( H_500 S_000 | P_010 S_000 )^0 + ( H_500 S_000 | S_000 S_000 )^0_{e} - ( I_510 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[3] = etfac[1] * AUX_INT__h_s_p_s[1] + 1 * one_over_2q * AUX_INT__h_s_s_s[0] - p_over_q * AUX_INT__i_s_p_s[4];

                    // ( H_500 S_000 | D_002 S_000 )^0 = z * ( H_500 S_000 | P_001 S_000 )^0 + ( H_500 S_000 | S_000 S_000 )^0_{e} - ( I_501 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[5] = etfac[2] * AUX_INT__h_s_p_s[2] + 1 * one_over_2q * AUX_INT__h_s_s_s[0] - p_over_q * AUX_INT__i_s_p_s[8];

                    // ( H_410 S_000 | D_200 S_000 )^0 = x * ( H_410 S_000 | P_100 S_000 )^0 + ( G_310 S_000 | P_100 S_000 )^0 + ( H_410 S_000 | S_000 S_000 )^0_{e} - ( I_510 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[6] = etfac[0] * AUX_INT__h_s_p_s[3] + 4 * one_over_2q * AUX_INT__g_s_p_s[3] + 1 * one_over_2q * AUX_INT__h_s_s_s[1] - p_over_q * AUX_INT__i_s_p_s[3];

                    // ( H_410 S_000 | D_020 S_000 )^0 = y * ( H_410 S_000 | P_010 S_000 )^0 + ( G_400 S_000 | P_010 S_000 )^0 + ( H_410 S_000 | S_000 S_000 )^0_{e} - ( I_420 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[9] = etfac[1] * AUX_INT__h_s_p_s[4] + 1 * one_over_2q * AUX_INT__g_s_p_s[1] + 1 * one_over_2q * AUX_INT__h_s_s_s[1] - p_over_q * AUX_INT__i_s_p_s[10];

                    // ( H_410 S_000 | D_002 S_000 )^0 = z * ( H_410 S_000 | P_001 S_000 )^0 + ( H_410 S_000 | S_000 S_000 )^0_{e} - ( I_411 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[11] = etfac[2] * AUX_INT__h_s_p_s[5] + 1 * one_over_2q * AUX_INT__h_s_s_s[1] - p_over_q * AUX_INT__i_s_p_s[14];

                    // ( H_401 S_000 | D_200 S_000 )^0 = x * ( H_401 S_000 | P_100 S_000 )^0 + ( G_301 S_000 | P_100 S_000 )^0 + ( H_401 S_000 | S_000 S_000 )^0_{e} - ( I_501 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[12] = etfac[0] * AUX_INT__h_s_p_s[6] + 4 * one_over_2q * AUX_INT__g_s_p_s[6] + 1 * one_over_2q * AUX_INT__h_s_s_s[2] - p_over_q * AUX_INT__i_s_p_s[6];

                    // ( H_401 S_000 | D_110 S_000 )^0 = y * ( H_401 S_000 | P_100 S_000 )^0 - ( I_411 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[13] = etfac[1] * AUX_INT__h_s_p_s[6] - p_over_q * AUX_INT__i_s_p_s[12];

                    // ( H_401 S_000 | D_020 S_000 )^0 = y * ( H_401 S_000 | P_010 S_000 )^0 + ( H_401 S_000 | S_000 S_000 )^0_{e} - ( I_411 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[15] = etfac[1] * AUX_INT__h_s_p_s[7] + 1 * one_over_2q * AUX_INT__h_s_s_s[2] - p_over_q * AUX_INT__i_s_p_s[13];

                    // ( H_401 S_000 | D_002 S_000 )^0 = z * ( H_401 S_000 | P_001 S_000 )^0 + ( G_400 S_000 | P_001 S_000 )^0 + ( H_401 S_000 | S_000 S_000 )^0_{e} - ( I_402 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[17] = etfac[2] * AUX_INT__h_s_p_s[8] + 1 * one_over_2q * AUX_INT__g_s_p_s[2] + 1 * one_over_2q * AUX_INT__h_s_s_s[2] - p_over_q * AUX_INT__i_s_p_s[17];

                    // ( H_320 S_000 | D_200 S_000 )^0 = x * ( H_320 S_000 | P_100 S_000 )^0 + ( G_220 S_000 | P_100 S_000 )^0 + ( H_320 S_000 | S_000 S_000 )^0_{e} - ( I_420 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[18] = etfac[0] * AUX_INT__h_s_p_s[9] + 3 * one_over_2q * AUX_INT__g_s_p_s[9] + 1 * one_over_2q * AUX_INT__h_s_s_s[3] - p_over_q * AUX_INT__i_s_p_s[9];

                    // ( H_320 S_000 | D_020 S_000 )^0 = y * ( H_320 S_000 | P_010 S_000 )^0 + ( G_310 S_000 | P_010 S_000 )^0 + ( H_320 S_000 | S_000 S_000 )^0_{e} - ( I_330 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[21] = etfac[1] * AUX_INT__h_s_p_s[10] + 2 * one_over_2q * AUX_INT__g_s_p_s[4] + 1 * one_over_2q * AUX_INT__h_s_s_s[3] - p_over_q * AUX_INT__i_s_p_s[19];

                    // ( H_320 S_000 | D_002 S_000 )^0 = z * ( H_320 S_000 | P_001 S_000 )^0 + ( H_320 S_000 | S_000 S_000 )^0_{e} - ( I_321 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[23] = etfac[2] * AUX_INT__h_s_p_s[11] + 1 * one_over_2q * AUX_INT__h_s_s_s[3] - p_over_q * AUX_INT__i_s_p_s[23];

                    // ( H_311 S_000 | D_200 S_000 )^0 = x * ( H_311 S_000 | P_100 S_000 )^0 + ( G_211 S_000 | P_100 S_000 )^0 + ( H_311 S_000 | S_000 S_000 )^0_{e} - ( I_411 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[24] = etfac[0] * AUX_INT__h_s_p_s[12] + 3 * one_over_2q * AUX_INT__g_s_p_s[12] + 1 * one_over_2q * AUX_INT__h_s_s_s[4] - p_over_q * AUX_INT__i_s_p_s[12];

                    // ( H_311 S_000 | D_110 S_000 )^0 = y * ( H_311 S_000 | P_100 S_000 )^0 + ( G_301 S_000 | P_100 S_000 )^0 - ( I_321 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[25] = etfac[1] * AUX_INT__h_s_p_s[12] + 1 * one_over_2q * AUX_INT__g_s_p_s[6] - p_over_q * AUX_INT__i_s_p_s[21];

                    // ( H_311 S_000 | D_020 S_000 )^0 = y * ( H_311 S_000 | P_010 S_000 )^0 + ( G_301 S_000 | P_010 S_000 )^0 + ( H_311 S_000 | S_000 S_000 )^0_{e} - ( I_321 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[27] = etfac[1] * AUX_INT__h_s_p_s[13] + 1 * one_over_2q * AUX_INT__g_s_p_s[7] + 1 * one_over_2q * AUX_INT__h_s_s_s[4] - p_over_q * AUX_INT__i_s_p_s[22];

                    // ( H_311 S_000 | D_002 S_000 )^0 = z * ( H_311 S_000 | P_001 S_000 )^0 + ( G_310 S_000 | P_001 S_000 )^0 + ( H_311 S_000 | S_000 S_000 )^0_{e} - ( I_312 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[29] = etfac[2] * AUX_INT__h_s_p_s[14] + 1 * one_over_2q * AUX_INT__g_s_p_s[5] + 1 * one_over_2q * AUX_INT__h_s_s_s[4] - p_over_q * AUX_INT__i_s_p_s[26];

                    // ( H_302 S_000 | D_200 S_000 )^0 = x * ( H_302 S_000 | P_100 S_000 )^0 + ( G_202 S_000 | P_100 S_000 )^0 + ( H_302 S_000 | S_000 S_000 )^0_{e} - ( I_402 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[30] = etfac[0] * AUX_INT__h_s_p_s[15] + 3 * one_over_2q * AUX_INT__g_s_p_s[15] + 1 * one_over_2q * AUX_INT__h_s_s_s[5] - p_over_q * AUX_INT__i_s_p_s[15];

                    // ( H_302 S_000 | D_110 S_000 )^0 = y * ( H_302 S_000 | P_100 S_000 )^0 - ( I_312 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[31] = etfac[1] * AUX_INT__h_s_p_s[15] - p_over_q * AUX_INT__i_s_p_s[24];

                    // ( H_302 S_000 | D_020 S_000 )^0 = y * ( H_302 S_000 | P_010 S_000 )^0 + ( H_302 S_000 | S_000 S_000 )^0_{e} - ( I_312 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[33] = etfac[1] * AUX_INT__h_s_p_s[16] + 1 * one_over_2q * AUX_INT__h_s_s_s[5] - p_over_q * AUX_INT__i_s_p_s[25];

                    // ( H_302 S_000 | D_002 S_000 )^0 = z * ( H_302 S_000 | P_001 S_000 )^0 + ( G_301 S_000 | P_001 S_000 )^0 + ( H_302 S_000 | S_000 S_000 )^0_{e} - ( I_303 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[35] = etfac[2] * AUX_INT__h_s_p_s[17] + 2 * one_over_2q * AUX_INT__g_s_p_s[8] + 1 * one_over_2q * AUX_INT__h_s_s_s[5] - p_over_q * AUX_INT__i_s_p_s[29];

                    // ( H_230 S_000 | D_200 S_000 )^0 = x * ( H_230 S_000 | P_100 S_000 )^0 + ( G_130 S_000 | P_100 S_000 )^0 + ( H_230 S_000 | S_000 S_000 )^0_{e} - ( I_330 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[36] = etfac[0] * AUX_INT__h_s_p_s[18] + 2 * one_over_2q * AUX_INT__g_s_p_s[18] + 1 * one_over_2q * AUX_INT__h_s_s_s[6] - p_over_q * AUX_INT__i_s_p_s[18];

                    // ( H_230 S_000 | D_020 S_000 )^0 = y * ( H_230 S_000 | P_010 S_000 )^0 + ( G_220 S_000 | P_010 S_000 )^0 + ( H_230 S_000 | S_000 S_000 )^0_{e} - ( I_240 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[39] = etfac[1] * AUX_INT__h_s_p_s[19] + 3 * one_over_2q * AUX_INT__g_s_p_s[10] + 1 * one_over_2q * AUX_INT__h_s_s_s[6] - p_over_q * AUX_INT__i_s_p_s[31];

                    // ( H_230 S_000 | D_002 S_000 )^0 = z * ( H_230 S_000 | P_001 S_000 )^0 + ( H_230 S_000 | S_000 S_000 )^0_{e} - ( I_231 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[41] = etfac[2] * AUX_INT__h_s_p_s[20] + 1 * one_over_2q * AUX_INT__h_s_s_s[6] - p_over_q * AUX_INT__i_s_p_s[35];

                    // ( H_221 S_000 | D_200 S_000 )^0 = x * ( H_221 S_000 | P_100 S_000 )^0 + ( G_121 S_000 | P_100 S_000 )^0 + ( H_221 S_000 | S_000 S_000 )^0_{e} - ( I_321 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[42] = etfac[0] * AUX_INT__h_s_p_s[21] + 2 * one_over_2q * AUX_INT__g_s_p_s[21] + 1 * one_over_2q * AUX_INT__h_s_s_s[7] - p_over_q * AUX_INT__i_s_p_s[21];

                    // ( H_221 S_000 | D_110 S_000 )^0 = y * ( H_221 S_000 | P_100 S_000 )^0 + ( G_211 S_000 | P_100 S_000 )^0 - ( I_231 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[43] = etfac[1] * AUX_INT__h_s_p_s[21] + 2 * one_over_2q * AUX_INT__g_s_p_s[12] - p_over_q * AUX_INT__i_s_p_s[33];

                    // ( H_221 S_000 | D_020 S_000 )^0 = y * ( H_221 S_000 | P_010 S_000 )^0 + ( G_211 S_000 | P_010 S_000 )^0 + ( H_221 S_000 | S_000 S_000 )^0_{e} - ( I_231 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[45] = etfac[1] * AUX_INT__h_s_p_s[22] + 2 * one_over_2q * AUX_INT__g_s_p_s[13] + 1 * one_over_2q * AUX_INT__h_s_s_s[7] - p_over_q * AUX_INT__i_s_p_s[34];

                    // ( H_221 S_000 | D_002 S_000 )^0 = z * ( H_221 S_000 | P_001 S_000 )^0 + ( G_220 S_000 | P_001 S_000 )^0 + ( H_221 S_000 | S_000 S_000 )^0_{e} - ( I_222 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[47] = etfac[2] * AUX_INT__h_s_p_s[23] + 1 * one_over_2q * AUX_INT__g_s_p_s[11] + 1 * one_over_2q * AUX_INT__h_s_s_s[7] - p_over_q * AUX_INT__i_s_p_s[38];

                    // ( H_212 S_000 | D_200 S_000 )^0 = x * ( H_212 S_000 | P_100 S_000 )^0 + ( G_112 S_000 | P_100 S_000 )^0 + ( H_212 S_000 | S_000 S_000 )^0_{e} - ( I_312 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[48] = etfac[0] * AUX_INT__h_s_p_s[24] + 2 * one_over_2q * AUX_INT__g_s_p_s[24] + 1 * one_over_2q * AUX_INT__h_s_s_s[8] - p_over_q * AUX_INT__i_s_p_s[24];

                    // ( H_212 S_000 | D_110 S_000 )^0 = y * ( H_212 S_000 | P_100 S_000 )^0 + ( G_202 S_000 | P_100 S_000 )^0 - ( I_222 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[49] = etfac[1] * AUX_INT__h_s_p_s[24] + 1 * one_over_2q * AUX_INT__g_s_p_s[15] - p_over_q * AUX_INT__i_s_p_s[36];

                    // ( H_212 S_000 | D_020 S_000 )^0 = y * ( H_212 S_000 | P_010 S_000 )^0 + ( G_202 S_000 | P_010 S_000 )^0 + ( H_212 S_000 | S_000 S_000 )^0_{e} - ( I_222 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[51] = etfac[1] * AUX_INT__h_s_p_s[25] + 1 * one_over_2q * AUX_INT__g_s_p_s[16] + 1 * one_over_2q * AUX_INT__h_s_s_s[8] - p_over_q * AUX_INT__i_s_p_s[37];

                    // ( H_212 S_000 | D_002 S_000 )^0 = z * ( H_212 S_000 | P_001 S_000 )^0 + ( G_211 S_000 | P_001 S_000 )^0 + ( H_212 S_000 | S_000 S_000 )^0_{e} - ( I_213 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[53] = etfac[2] * AUX_INT__h_s_p_s[26] + 2 * one_over_2q * AUX_INT__g_s_p_s[14] + 1 * one_over_2q * AUX_INT__h_s_s_s[8] - p_over_q * AUX_INT__i_s_p_s[41];

                    // ( H_203 S_000 | D_200 S_000 )^0 = x * ( H_203 S_000 | P_100 S_000 )^0 + ( G_103 S_000 | P_100 S_000 )^0 + ( H_203 S_000 | S_000 S_000 )^0_{e} - ( I_303 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[54] = etfac[0] * AUX_INT__h_s_p_s[27] + 2 * one_over_2q * AUX_INT__g_s_p_s[27] + 1 * one_over_2q * AUX_INT__h_s_s_s[9] - p_over_q * AUX_INT__i_s_p_s[27];

                    // ( H_203 S_000 | D_110 S_000 )^0 = y * ( H_203 S_000 | P_100 S_000 )^0 - ( I_213 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[55] = etfac[1] * AUX_INT__h_s_p_s[27] - p_over_q * AUX_INT__i_s_p_s[39];

                    // ( H_203 S_000 | D_020 S_000 )^0 = y * ( H_203 S_000 | P_010 S_000 )^0 + ( H_203 S_000 | S_000 S_000 )^0_{e} - ( I_213 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[57] = etfac[1] * AUX_INT__h_s_p_s[28] + 1 * one_over_2q * AUX_INT__h_s_s_s[9] - p_over_q * AUX_INT__i_s_p_s[40];

                    // ( H_203 S_000 | D_002 S_000 )^0 = z * ( H_203 S_000 | P_001 S_000 )^0 + ( G_202 S_000 | P_001 S_000 )^0 + ( H_203 S_000 | S_000 S_000 )^0_{e} - ( I_204 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[59] = etfac[2] * AUX_INT__h_s_p_s[29] + 3 * one_over_2q * AUX_INT__g_s_p_s[17] + 1 * one_over_2q * AUX_INT__h_s_s_s[9] - p_over_q * AUX_INT__i_s_p_s[44];

                    // ( H_140 S_000 | D_200 S_000 )^0 = x * ( H_140 S_000 | P_100 S_000 )^0 + ( G_040 S_000 | P_100 S_000 )^0 + ( H_140 S_000 | S_000 S_000 )^0_{e} - ( I_240 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[60] = etfac[0] * AUX_INT__h_s_p_s[30] + 1 * one_over_2q * AUX_INT__g_s_p_s[30] + 1 * one_over_2q * AUX_INT__h_s_s_s[10] - p_over_q * AUX_INT__i_s_p_s[30];

                    // ( H_140 S_000 | D_020 S_000 )^0 = y * ( H_140 S_000 | P_010 S_000 )^0 + ( G_130 S_000 | P_010 S_000 )^0 + ( H_140 S_000 | S_000 S_000 )^0_{e} - ( I_150 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[63] = etfac[1] * AUX_INT__h_s_p_s[31] + 4 * one_over_2q * AUX_INT__g_s_p_s[19] + 1 * one_over_2q * AUX_INT__h_s_s_s[10] - p_over_q * AUX_INT__i_s_p_s[46];

                    // ( H_140 S_000 | D_002 S_000 )^0 = z * ( H_140 S_000 | P_001 S_000 )^0 + ( H_140 S_000 | S_000 S_000 )^0_{e} - ( I_141 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[65] = etfac[2] * AUX_INT__h_s_p_s[32] + 1 * one_over_2q * AUX_INT__h_s_s_s[10] - p_over_q * AUX_INT__i_s_p_s[50];

                    // ( H_131 S_000 | D_200 S_000 )^0 = x * ( H_131 S_000 | P_100 S_000 )^0 + ( G_031 S_000 | P_100 S_000 )^0 + ( H_131 S_000 | S_000 S_000 )^0_{e} - ( I_231 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[66] = etfac[0] * AUX_INT__h_s_p_s[33] + 1 * one_over_2q * AUX_INT__g_s_p_s[33] + 1 * one_over_2q * AUX_INT__h_s_s_s[11] - p_over_q * AUX_INT__i_s_p_s[33];

                    // ( H_131 S_000 | D_110 S_000 )^0 = y * ( H_131 S_000 | P_100 S_000 )^0 + ( G_121 S_000 | P_100 S_000 )^0 - ( I_141 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[67] = etfac[1] * AUX_INT__h_s_p_s[33] + 3 * one_over_2q * AUX_INT__g_s_p_s[21] - p_over_q * AUX_INT__i_s_p_s[48];

                    // ( H_131 S_000 | D_020 S_000 )^0 = y * ( H_131 S_000 | P_010 S_000 )^0 + ( G_121 S_000 | P_010 S_000 )^0 + ( H_131 S_000 | S_000 S_000 )^0_{e} - ( I_141 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[69] = etfac[1] * AUX_INT__h_s_p_s[34] + 3 * one_over_2q * AUX_INT__g_s_p_s[22] + 1 * one_over_2q * AUX_INT__h_s_s_s[11] - p_over_q * AUX_INT__i_s_p_s[49];

                    // ( H_131 S_000 | D_002 S_000 )^0 = z * ( H_131 S_000 | P_001 S_000 )^0 + ( G_130 S_000 | P_001 S_000 )^0 + ( H_131 S_000 | S_000 S_000 )^0_{e} - ( I_132 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[71] = etfac[2] * AUX_INT__h_s_p_s[35] + 1 * one_over_2q * AUX_INT__g_s_p_s[20] + 1 * one_over_2q * AUX_INT__h_s_s_s[11] - p_over_q * AUX_INT__i_s_p_s[53];

                    // ( H_122 S_000 | D_200 S_000 )^0 = x * ( H_122 S_000 | P_100 S_000 )^0 + ( G_022 S_000 | P_100 S_000 )^0 + ( H_122 S_000 | S_000 S_000 )^0_{e} - ( I_222 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[72] = etfac[0] * AUX_INT__h_s_p_s[36] + 1 * one_over_2q * AUX_INT__g_s_p_s[36] + 1 * one_over_2q * AUX_INT__h_s_s_s[12] - p_over_q * AUX_INT__i_s_p_s[36];

                    // ( H_122 S_000 | D_110 S_000 )^0 = y * ( H_122 S_000 | P_100 S_000 )^0 + ( G_112 S_000 | P_100 S_000 )^0 - ( I_132 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[73] = etfac[1] * AUX_INT__h_s_p_s[36] + 2 * one_over_2q * AUX_INT__g_s_p_s[24] - p_over_q * AUX_INT__i_s_p_s[51];

                    // ( H_122 S_000 | D_020 S_000 )^0 = y * ( H_122 S_000 | P_010 S_000 )^0 + ( G_112 S_000 | P_010 S_000 )^0 + ( H_122 S_000 | S_000 S_000 )^0_{e} - ( I_132 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[75] = etfac[1] * AUX_INT__h_s_p_s[37] + 2 * one_over_2q * AUX_INT__g_s_p_s[25] + 1 * one_over_2q * AUX_INT__h_s_s_s[12] - p_over_q * AUX_INT__i_s_p_s[52];

                    // ( H_122 S_000 | D_002 S_000 )^0 = z * ( H_122 S_000 | P_001 S_000 )^0 + ( G_121 S_000 | P_001 S_000 )^0 + ( H_122 S_000 | S_000 S_000 )^0_{e} - ( I_123 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[77] = etfac[2] * AUX_INT__h_s_p_s[38] + 2 * one_over_2q * AUX_INT__g_s_p_s[23] + 1 * one_over_2q * AUX_INT__h_s_s_s[12] - p_over_q * AUX_INT__i_s_p_s[56];

                    // ( H_113 S_000 | D_200 S_000 )^0 = x * ( H_113 S_000 | P_100 S_000 )^0 + ( G_013 S_000 | P_100 S_000 )^0 + ( H_113 S_000 | S_000 S_000 )^0_{e} - ( I_213 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[78] = etfac[0] * AUX_INT__h_s_p_s[39] + 1 * one_over_2q * AUX_INT__g_s_p_s[39] + 1 * one_over_2q * AUX_INT__h_s_s_s[13] - p_over_q * AUX_INT__i_s_p_s[39];

                    // ( H_113 S_000 | D_110 S_000 )^0 = y * ( H_113 S_000 | P_100 S_000 )^0 + ( G_103 S_000 | P_100 S_000 )^0 - ( I_123 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[79] = etfac[1] * AUX_INT__h_s_p_s[39] + 1 * one_over_2q * AUX_INT__g_s_p_s[27] - p_over_q * AUX_INT__i_s_p_s[54];

                    // ( H_113 S_000 | D_020 S_000 )^0 = y * ( H_113 S_000 | P_010 S_000 )^0 + ( G_103 S_000 | P_010 S_000 )^0 + ( H_113 S_000 | S_000 S_000 )^0_{e} - ( I_123 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[81] = etfac[1] * AUX_INT__h_s_p_s[40] + 1 * one_over_2q * AUX_INT__g_s_p_s[28] + 1 * one_over_2q * AUX_INT__h_s_s_s[13] - p_over_q * AUX_INT__i_s_p_s[55];

                    // ( H_113 S_000 | D_002 S_000 )^0 = z * ( H_113 S_000 | P_001 S_000 )^0 + ( G_112 S_000 | P_001 S_000 )^0 + ( H_113 S_000 | S_000 S_000 )^0_{e} - ( I_114 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[83] = etfac[2] * AUX_INT__h_s_p_s[41] + 3 * one_over_2q * AUX_INT__g_s_p_s[26] + 1 * one_over_2q * AUX_INT__h_s_s_s[13] - p_over_q * AUX_INT__i_s_p_s[59];

                    // ( H_104 S_000 | D_200 S_000 )^0 = x * ( H_104 S_000 | P_100 S_000 )^0 + ( G_004 S_000 | P_100 S_000 )^0 + ( H_104 S_000 | S_000 S_000 )^0_{e} - ( I_204 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[84] = etfac[0] * AUX_INT__h_s_p_s[42] + 1 * one_over_2q * AUX_INT__g_s_p_s[42] + 1 * one_over_2q * AUX_INT__h_s_s_s[14] - p_over_q * AUX_INT__i_s_p_s[42];

                    // ( H_104 S_000 | D_110 S_000 )^0 = y * ( H_104 S_000 | P_100 S_000 )^0 - ( I_114 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[85] = etfac[1] * AUX_INT__h_s_p_s[42] - p_over_q * AUX_INT__i_s_p_s[57];

                    // ( H_104 S_000 | D_020 S_000 )^0 = y * ( H_104 S_000 | P_010 S_000 )^0 + ( H_104 S_000 | S_000 S_000 )^0_{e} - ( I_114 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[87] = etfac[1] * AUX_INT__h_s_p_s[43] + 1 * one_over_2q * AUX_INT__h_s_s_s[14] - p_over_q * AUX_INT__i_s_p_s[58];

                    // ( H_104 S_000 | D_002 S_000 )^0 = z * ( H_104 S_000 | P_001 S_000 )^0 + ( G_103 S_000 | P_001 S_000 )^0 + ( H_104 S_000 | S_000 S_000 )^0_{e} - ( I_105 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[89] = etfac[2] * AUX_INT__h_s_p_s[44] + 4 * one_over_2q * AUX_INT__g_s_p_s[29] + 1 * one_over_2q * AUX_INT__h_s_s_s[14] - p_over_q * AUX_INT__i_s_p_s[62];

                    // ( H_050 S_000 | D_200 S_000 )^0 = x * ( H_050 S_000 | P_100 S_000 )^0 + ( H_050 S_000 | S_000 S_000 )^0_{e} - ( I_150 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[90] = etfac[0] * AUX_INT__h_s_p_s[45] + 1 * one_over_2q * AUX_INT__h_s_s_s[15] - p_over_q * AUX_INT__i_s_p_s[45];

                    // ( H_050 S_000 | D_020 S_000 )^0 = y * ( H_050 S_000 | P_010 S_000 )^0 + ( G_040 S_000 | P_010 S_000 )^0 + ( H_050 S_000 | S_000 S_000 )^0_{e} - ( I_060 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[93] = etfac[1] * AUX_INT__h_s_p_s[46] + 5 * one_over_2q * AUX_INT__g_s_p_s[31] + 1 * one_over_2q * AUX_INT__h_s_s_s[15] - p_over_q * AUX_INT__i_s_p_s[64];

                    // ( H_050 S_000 | D_002 S_000 )^0 = z * ( H_050 S_000 | P_001 S_000 )^0 + ( H_050 S_000 | S_000 S_000 )^0_{e} - ( I_051 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[95] = etfac[2] * AUX_INT__h_s_p_s[47] + 1 * one_over_2q * AUX_INT__h_s_s_s[15] - p_over_q * AUX_INT__i_s_p_s[68];

                    // ( H_041 S_000 | D_200 S_000 )^0 = x * ( H_041 S_000 | P_100 S_000 )^0 + ( H_041 S_000 | S_000 S_000 )^0_{e} - ( I_141 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[96] = etfac[0] * AUX_INT__h_s_p_s[48] + 1 * one_over_2q * AUX_INT__h_s_s_s[16] - p_over_q * AUX_INT__i_s_p_s[48];

                    // ( H_041 S_000 | D_110 S_000 )^0 = y * ( H_041 S_000 | P_100 S_000 )^0 + ( G_031 S_000 | P_100 S_000 )^0 - ( I_051 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[97] = etfac[1] * AUX_INT__h_s_p_s[48] + 4 * one_over_2q * AUX_INT__g_s_p_s[33] - p_over_q * AUX_INT__i_s_p_s[66];

                    // ( H_041 S_000 | D_020 S_000 )^0 = y * ( H_041 S_000 | P_010 S_000 )^0 + ( G_031 S_000 | P_010 S_000 )^0 + ( H_041 S_000 | S_000 S_000 )^0_{e} - ( I_051 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[99] = etfac[1] * AUX_INT__h_s_p_s[49] + 4 * one_over_2q * AUX_INT__g_s_p_s[34] + 1 * one_over_2q * AUX_INT__h_s_s_s[16] - p_over_q * AUX_INT__i_s_p_s[67];

                    // ( H_041 S_000 | D_002 S_000 )^0 = z * ( H_041 S_000 | P_001 S_000 )^0 + ( G_040 S_000 | P_001 S_000 )^0 + ( H_041 S_000 | S_000 S_000 )^0_{e} - ( I_042 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[101] = etfac[2] * AUX_INT__h_s_p_s[50] + 1 * one_over_2q * AUX_INT__g_s_p_s[32] + 1 * one_over_2q * AUX_INT__h_s_s_s[16] - p_over_q * AUX_INT__i_s_p_s[71];

                    // ( H_032 S_000 | D_200 S_000 )^0 = x * ( H_032 S_000 | P_100 S_000 )^0 + ( H_032 S_000 | S_000 S_000 )^0_{e} - ( I_132 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[102] = etfac[0] * AUX_INT__h_s_p_s[51] + 1 * one_over_2q * AUX_INT__h_s_s_s[17] - p_over_q * AUX_INT__i_s_p_s[51];

                    // ( H_032 S_000 | D_110 S_000 )^0 = y * ( H_032 S_000 | P_100 S_000 )^0 + ( G_022 S_000 | P_100 S_000 )^0 - ( I_042 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[103] = etfac[1] * AUX_INT__h_s_p_s[51] + 3 * one_over_2q * AUX_INT__g_s_p_s[36] - p_over_q * AUX_INT__i_s_p_s[69];

                    // ( H_032 S_000 | D_020 S_000 )^0 = y * ( H_032 S_000 | P_010 S_000 )^0 + ( G_022 S_000 | P_010 S_000 )^0 + ( H_032 S_000 | S_000 S_000 )^0_{e} - ( I_042 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[105] = etfac[1] * AUX_INT__h_s_p_s[52] + 3 * one_over_2q * AUX_INT__g_s_p_s[37] + 1 * one_over_2q * AUX_INT__h_s_s_s[17] - p_over_q * AUX_INT__i_s_p_s[70];

                    // ( H_032 S_000 | D_002 S_000 )^0 = z * ( H_032 S_000 | P_001 S_000 )^0 + ( G_031 S_000 | P_001 S_000 )^0 + ( H_032 S_000 | S_000 S_000 )^0_{e} - ( I_033 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[107] = etfac[2] * AUX_INT__h_s_p_s[53] + 2 * one_over_2q * AUX_INT__g_s_p_s[35] + 1 * one_over_2q * AUX_INT__h_s_s_s[17] - p_over_q * AUX_INT__i_s_p_s[74];

                    // ( H_023 S_000 | D_200 S_000 )^0 = x * ( H_023 S_000 | P_100 S_000 )^0 + ( H_023 S_000 | S_000 S_000 )^0_{e} - ( I_123 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[108] = etfac[0] * AUX_INT__h_s_p_s[54] + 1 * one_over_2q * AUX_INT__h_s_s_s[18] - p_over_q * AUX_INT__i_s_p_s[54];

                    // ( H_023 S_000 | D_110 S_000 )^0 = y * ( H_023 S_000 | P_100 S_000 )^0 + ( G_013 S_000 | P_100 S_000 )^0 - ( I_033 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[109] = etfac[1] * AUX_INT__h_s_p_s[54] + 2 * one_over_2q * AUX_INT__g_s_p_s[39] - p_over_q * AUX_INT__i_s_p_s[72];

                    // ( H_023 S_000 | D_020 S_000 )^0 = y * ( H_023 S_000 | P_010 S_000 )^0 + ( G_013 S_000 | P_010 S_000 )^0 + ( H_023 S_000 | S_000 S_000 )^0_{e} - ( I_033 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[111] = etfac[1] * AUX_INT__h_s_p_s[55] + 2 * one_over_2q * AUX_INT__g_s_p_s[40] + 1 * one_over_2q * AUX_INT__h_s_s_s[18] - p_over_q * AUX_INT__i_s_p_s[73];

                    // ( H_023 S_000 | D_002 S_000 )^0 = z * ( H_023 S_000 | P_001 S_000 )^0 + ( G_022 S_000 | P_001 S_000 )^0 + ( H_023 S_000 | S_000 S_000 )^0_{e} - ( I_024 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[113] = etfac[2] * AUX_INT__h_s_p_s[56] + 3 * one_over_2q * AUX_INT__g_s_p_s[38] + 1 * one_over_2q * AUX_INT__h_s_s_s[18] - p_over_q * AUX_INT__i_s_p_s[77];

                    // ( H_014 S_000 | D_200 S_000 )^0 = x * ( H_014 S_000 | P_100 S_000 )^0 + ( H_014 S_000 | S_000 S_000 )^0_{e} - ( I_114 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[114] = etfac[0] * AUX_INT__h_s_p_s[57] + 1 * one_over_2q * AUX_INT__h_s_s_s[19] - p_over_q * AUX_INT__i_s_p_s[57];

                    // ( H_014 S_000 | D_110 S_000 )^0 = y * ( H_014 S_000 | P_100 S_000 )^0 + ( G_004 S_000 | P_100 S_000 )^0 - ( I_024 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[115] = etfac[1] * AUX_INT__h_s_p_s[57] + 1 * one_over_2q * AUX_INT__g_s_p_s[42] - p_over_q * AUX_INT__i_s_p_s[75];

                    // ( H_014 S_000 | D_020 S_000 )^0 = y * ( H_014 S_000 | P_010 S_000 )^0 + ( G_004 S_000 | P_010 S_000 )^0 + ( H_014 S_000 | S_000 S_000 )^0_{e} - ( I_024 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[117] = etfac[1] * AUX_INT__h_s_p_s[58] + 1 * one_over_2q * AUX_INT__g_s_p_s[43] + 1 * one_over_2q * AUX_INT__h_s_s_s[19] - p_over_q * AUX_INT__i_s_p_s[76];

                    // ( H_014 S_000 | D_002 S_000 )^0 = z * ( H_014 S_000 | P_001 S_000 )^0 + ( G_013 S_000 | P_001 S_000 )^0 + ( H_014 S_000 | S_000 S_000 )^0_{e} - ( I_015 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[119] = etfac[2] * AUX_INT__h_s_p_s[59] + 4 * one_over_2q * AUX_INT__g_s_p_s[41] + 1 * one_over_2q * AUX_INT__h_s_s_s[19] - p_over_q * AUX_INT__i_s_p_s[80];

                    // ( H_005 S_000 | D_200 S_000 )^0 = x * ( H_005 S_000 | P_100 S_000 )^0 + ( H_005 S_000 | S_000 S_000 )^0_{e} - ( I_105 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[120] = etfac[0] * AUX_INT__h_s_p_s[60] + 1 * one_over_2q * AUX_INT__h_s_s_s[20] - p_over_q * AUX_INT__i_s_p_s[60];

                    // ( H_005 S_000 | D_110 S_000 )^0 = y * ( H_005 S_000 | P_100 S_000 )^0 - ( I_015 S_000 | P_100 S_000 )^0
                    AUX_INT__h_s_d_s[121] = etfac[1] * AUX_INT__h_s_p_s[60] - p_over_q * AUX_INT__i_s_p_s[78];

                    // ( H_005 S_000 | D_020 S_000 )^0 = y * ( H_005 S_000 | P_010 S_000 )^0 + ( H_005 S_000 | S_000 S_000 )^0_{e} - ( I_015 S_000 | P_010 S_000 )^0
                    AUX_INT__h_s_d_s[123] = etfac[1] * AUX_INT__h_s_p_s[61] + 1 * one_over_2q * AUX_INT__h_s_s_s[20] - p_over_q * AUX_INT__i_s_p_s[79];

                    // ( H_005 S_000 | D_002 S_000 )^0 = z * ( H_005 S_000 | P_001 S_000 )^0 + ( G_004 S_000 | P_001 S_000 )^0 + ( H_005 S_000 | S_000 S_000 )^0_{e} - ( I_006 S_000 | P_001 S_000 )^0
                    AUX_INT__h_s_d_s[125] = etfac[2] * AUX_INT__h_s_p_s[62] + 5 * one_over_2q * AUX_INT__g_s_p_s[44] + 1 * one_over_2q * AUX_INT__h_s_s_s[20] - p_over_q * AUX_INT__i_s_p_s[83];

                    // ( G_400 S_000 | F_300 S_000 )^0_{t} = x * ( G_400 S_000 | D_200 S_000 )^0 + ( F_300 S_000 | D_200 S_000 )^0 + ( G_400 S_000 | P_100 S_000 )^0 - ( H_500 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[0] = etfac[0] * AUX_INT__g_s_d_s[0] + 4 * one_over_2q * AUX_INT__f_s_d_s[0] + 2 * one_over_2q * AUX_INT__g_s_p_s[0] - p_over_q * AUX_INT__h_s_d_s[0];

                    // ( G_400 S_000 | F_210 S_000 )^0_{t} = y * ( G_400 S_000 | D_200 S_000 )^0 - ( H_410 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[1] = etfac[1] * AUX_INT__g_s_d_s[0] - p_over_q * AUX_INT__h_s_d_s[6];

                    // ( G_400 S_000 | F_201 S_000 )^0_{t} = z * ( G_400 S_000 | D_200 S_000 )^0 - ( H_401 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[2] = etfac[2] * AUX_INT__g_s_d_s[0] - p_over_q * AUX_INT__h_s_d_s[12];

                    // ( G_400 S_000 | F_120 S_000 )^0_{t} = x * ( G_400 S_000 | D_020 S_000 )^0 + ( F_300 S_000 | D_020 S_000 )^0 - ( H_500 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[3] = etfac[0] * AUX_INT__g_s_d_s[3] + 4 * one_over_2q * AUX_INT__f_s_d_s[3] - p_over_q * AUX_INT__h_s_d_s[3];

                    // ( G_400 S_000 | F_111 S_000 )^0_{t} = z * ( G_400 S_000 | D_110 S_000 )^0 - ( H_401 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[4] = etfac[2] * AUX_INT__g_s_d_s[1] - p_over_q * AUX_INT__h_s_d_s[13];

                    // ( G_400 S_000 | F_102 S_000 )^0_{t} = x * ( G_400 S_000 | D_002 S_000 )^0 + ( F_300 S_000 | D_002 S_000 )^0 - ( H_500 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[5] = etfac[0] * AUX_INT__g_s_d_s[5] + 4 * one_over_2q * AUX_INT__f_s_d_s[5] - p_over_q * AUX_INT__h_s_d_s[5];

                    // ( G_400 S_000 | F_030 S_000 )^0_{t} = y * ( G_400 S_000 | D_020 S_000 )^0 + ( G_400 S_000 | P_010 S_000 )^0 - ( H_410 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[6] = etfac[1] * AUX_INT__g_s_d_s[3] + 2 * one_over_2q * AUX_INT__g_s_p_s[1] - p_over_q * AUX_INT__h_s_d_s[9];

                    // ( G_400 S_000 | F_021 S_000 )^0_{t} = z * ( G_400 S_000 | D_020 S_000 )^0 - ( H_401 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[7] = etfac[2] * AUX_INT__g_s_d_s[3] - p_over_q * AUX_INT__h_s_d_s[15];

                    // ( G_400 S_000 | F_012 S_000 )^0_{t} = y * ( G_400 S_000 | D_002 S_000 )^0 - ( H_410 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[8] = etfac[1] * AUX_INT__g_s_d_s[5] - p_over_q * AUX_INT__h_s_d_s[11];

                    // ( G_400 S_000 | F_003 S_000 )^0_{t} = z * ( G_400 S_000 | D_002 S_000 )^0 + ( G_400 S_000 | P_001 S_000 )^0 - ( H_401 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[9] = etfac[2] * AUX_INT__g_s_d_s[5] + 2 * one_over_2q * AUX_INT__g_s_p_s[2] - p_over_q * AUX_INT__h_s_d_s[17];

                    // ( G_310 S_000 | F_300 S_000 )^0_{t} = x * ( G_310 S_000 | D_200 S_000 )^0 + ( F_210 S_000 | D_200 S_000 )^0 + ( G_310 S_000 | P_100 S_000 )^0 - ( H_410 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[10] = etfac[0] * AUX_INT__g_s_d_s[6] + 3 * one_over_2q * AUX_INT__f_s_d_s[6] + 2 * one_over_2q * AUX_INT__g_s_p_s[3] - p_over_q * AUX_INT__h_s_d_s[6];

                    // ( G_310 S_000 | F_210 S_000 )^0_{t} = y * ( G_310 S_000 | D_200 S_000 )^0 + ( F_300 S_000 | D_200 S_000 )^0 - ( H_320 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[11] = etfac[1] * AUX_INT__g_s_d_s[6] + 1 * one_over_2q * AUX_INT__f_s_d_s[0] - p_over_q * AUX_INT__h_s_d_s[18];

                    // ( G_310 S_000 | F_201 S_000 )^0_{t} = z * ( G_310 S_000 | D_200 S_000 )^0 - ( H_311 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[12] = etfac[2] * AUX_INT__g_s_d_s[6] - p_over_q * AUX_INT__h_s_d_s[24];

                    // ( G_310 S_000 | F_120 S_000 )^0_{t} = x * ( G_310 S_000 | D_020 S_000 )^0 + ( F_210 S_000 | D_020 S_000 )^0 - ( H_410 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[13] = etfac[0] * AUX_INT__g_s_d_s[9] + 3 * one_over_2q * AUX_INT__f_s_d_s[9] - p_over_q * AUX_INT__h_s_d_s[9];

                    // ( G_310 S_000 | F_111 S_000 )^0_{t} = z * ( G_310 S_000 | D_110 S_000 )^0 - ( H_311 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[14] = etfac[2] * AUX_INT__g_s_d_s[7] - p_over_q * AUX_INT__h_s_d_s[25];

                    // ( G_310 S_000 | F_102 S_000 )^0_{t} = x * ( G_310 S_000 | D_002 S_000 )^0 + ( F_210 S_000 | D_002 S_000 )^0 - ( H_410 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[15] = etfac[0] * AUX_INT__g_s_d_s[11] + 3 * one_over_2q * AUX_INT__f_s_d_s[11] - p_over_q * AUX_INT__h_s_d_s[11];

                    // ( G_310 S_000 | F_030 S_000 )^0_{t} = y * ( G_310 S_000 | D_020 S_000 )^0 + ( F_300 S_000 | D_020 S_000 )^0 + ( G_310 S_000 | P_010 S_000 )^0 - ( H_320 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[16] = etfac[1] * AUX_INT__g_s_d_s[9] + 1 * one_over_2q * AUX_INT__f_s_d_s[3] + 2 * one_over_2q * AUX_INT__g_s_p_s[4] - p_over_q * AUX_INT__h_s_d_s[21];

                    // ( G_310 S_000 | F_021 S_000 )^0_{t} = z * ( G_310 S_000 | D_020 S_000 )^0 - ( H_311 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[17] = etfac[2] * AUX_INT__g_s_d_s[9] - p_over_q * AUX_INT__h_s_d_s[27];

                    // ( G_310 S_000 | F_012 S_000 )^0_{t} = y * ( G_310 S_000 | D_002 S_000 )^0 + ( F_300 S_000 | D_002 S_000 )^0 - ( H_320 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[18] = etfac[1] * AUX_INT__g_s_d_s[11] + 1 * one_over_2q * AUX_INT__f_s_d_s[5] - p_over_q * AUX_INT__h_s_d_s[23];

                    // ( G_310 S_000 | F_003 S_000 )^0_{t} = z * ( G_310 S_000 | D_002 S_000 )^0 + ( G_310 S_000 | P_001 S_000 )^0 - ( H_311 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[19] = etfac[2] * AUX_INT__g_s_d_s[11] + 2 * one_over_2q * AUX_INT__g_s_p_s[5] - p_over_q * AUX_INT__h_s_d_s[29];

                    // ( G_301 S_000 | F_300 S_000 )^0_{t} = x * ( G_301 S_000 | D_200 S_000 )^0 + ( F_201 S_000 | D_200 S_000 )^0 + ( G_301 S_000 | P_100 S_000 )^0 - ( H_401 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[20] = etfac[0] * AUX_INT__g_s_d_s[12] + 3 * one_over_2q * AUX_INT__f_s_d_s[12] + 2 * one_over_2q * AUX_INT__g_s_p_s[6] - p_over_q * AUX_INT__h_s_d_s[12];

                    // ( G_301 S_000 | F_210 S_000 )^0_{t} = y * ( G_301 S_000 | D_200 S_000 )^0 - ( H_311 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[21] = etfac[1] * AUX_INT__g_s_d_s[12] - p_over_q * AUX_INT__h_s_d_s[24];

                    // ( G_301 S_000 | F_201 S_000 )^0_{t} = z * ( G_301 S_000 | D_200 S_000 )^0 + ( F_300 S_000 | D_200 S_000 )^0 - ( H_302 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[22] = etfac[2] * AUX_INT__g_s_d_s[12] + 1 * one_over_2q * AUX_INT__f_s_d_s[0] - p_over_q * AUX_INT__h_s_d_s[30];

                    // ( G_301 S_000 | F_120 S_000 )^0_{t} = x * ( G_301 S_000 | D_020 S_000 )^0 + ( F_201 S_000 | D_020 S_000 )^0 - ( H_401 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[23] = etfac[0] * AUX_INT__g_s_d_s[15] + 3 * one_over_2q * AUX_INT__f_s_d_s[15] - p_over_q * AUX_INT__h_s_d_s[15];

                    // ( G_301 S_000 | F_111 S_000 )^0_{t} = z * ( G_301 S_000 | D_110 S_000 )^0 + ( F_300 S_000 | D_110 S_000 )^0 - ( H_302 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[24] = etfac[2] * AUX_INT__g_s_d_s[13] + 1 * one_over_2q * AUX_INT__f_s_d_s[1] - p_over_q * AUX_INT__h_s_d_s[31];

                    // ( G_301 S_000 | F_102 S_000 )^0_{t} = x * ( G_301 S_000 | D_002 S_000 )^0 + ( F_201 S_000 | D_002 S_000 )^0 - ( H_401 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[25] = etfac[0] * AUX_INT__g_s_d_s[17] + 3 * one_over_2q * AUX_INT__f_s_d_s[17] - p_over_q * AUX_INT__h_s_d_s[17];

                    // ( G_301 S_000 | F_030 S_000 )^0_{t} = y * ( G_301 S_000 | D_020 S_000 )^0 + ( G_301 S_000 | P_010 S_000 )^0 - ( H_311 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[26] = etfac[1] * AUX_INT__g_s_d_s[15] + 2 * one_over_2q * AUX_INT__g_s_p_s[7] - p_over_q * AUX_INT__h_s_d_s[27];

                    // ( G_301 S_000 | F_021 S_000 )^0_{t} = z * ( G_301 S_000 | D_020 S_000 )^0 + ( F_300 S_000 | D_020 S_000 )^0 - ( H_302 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[27] = etfac[2] * AUX_INT__g_s_d_s[15] + 1 * one_over_2q * AUX_INT__f_s_d_s[3] - p_over_q * AUX_INT__h_s_d_s[33];

                    // ( G_301 S_000 | F_012 S_000 )^0_{t} = y * ( G_301 S_000 | D_002 S_000 )^0 - ( H_311 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[28] = etfac[1] * AUX_INT__g_s_d_s[17] - p_over_q * AUX_INT__h_s_d_s[29];

                    // ( G_301 S_000 | F_003 S_000 )^0_{t} = z * ( G_301 S_000 | D_002 S_000 )^0 + ( F_300 S_000 | D_002 S_000 )^0 + ( G_301 S_000 | P_001 S_000 )^0 - ( H_302 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[29] = etfac[2] * AUX_INT__g_s_d_s[17] + 1 * one_over_2q * AUX_INT__f_s_d_s[5] + 2 * one_over_2q * AUX_INT__g_s_p_s[8] - p_over_q * AUX_INT__h_s_d_s[35];

                    // ( G_220 S_000 | F_300 S_000 )^0_{t} = x * ( G_220 S_000 | D_200 S_000 )^0 + ( F_120 S_000 | D_200 S_000 )^0 + ( G_220 S_000 | P_100 S_000 )^0 - ( H_320 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[30] = etfac[0] * AUX_INT__g_s_d_s[18] + 2 * one_over_2q * AUX_INT__f_s_d_s[18] + 2 * one_over_2q * AUX_INT__g_s_p_s[9] - p_over_q * AUX_INT__h_s_d_s[18];

                    // ( G_220 S_000 | F_210 S_000 )^0_{t} = y * ( G_220 S_000 | D_200 S_000 )^0 + ( F_210 S_000 | D_200 S_000 )^0 - ( H_230 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[31] = etfac[1] * AUX_INT__g_s_d_s[18] + 2 * one_over_2q * AUX_INT__f_s_d_s[6] - p_over_q * AUX_INT__h_s_d_s[36];

                    // ( G_220 S_000 | F_201 S_000 )^0_{t} = z * ( G_220 S_000 | D_200 S_000 )^0 - ( H_221 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[32] = etfac[2] * AUX_INT__g_s_d_s[18] - p_over_q * AUX_INT__h_s_d_s[42];

                    // ( G_220 S_000 | F_120 S_000 )^0_{t} = x * ( G_220 S_000 | D_020 S_000 )^0 + ( F_120 S_000 | D_020 S_000 )^0 - ( H_320 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[33] = etfac[0] * AUX_INT__g_s_d_s[21] + 2 * one_over_2q * AUX_INT__f_s_d_s[21] - p_over_q * AUX_INT__h_s_d_s[21];

                    // ( G_220 S_000 | F_111 S_000 )^0_{t} = z * ( G_220 S_000 | D_110 S_000 )^0 - ( H_221 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[34] = etfac[2] * AUX_INT__g_s_d_s[19] - p_over_q * AUX_INT__h_s_d_s[43];

                    // ( G_220 S_000 | F_102 S_000 )^0_{t} = x * ( G_220 S_000 | D_002 S_000 )^0 + ( F_120 S_000 | D_002 S_000 )^0 - ( H_320 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[35] = etfac[0] * AUX_INT__g_s_d_s[23] + 2 * one_over_2q * AUX_INT__f_s_d_s[23] - p_over_q * AUX_INT__h_s_d_s[23];

                    // ( G_220 S_000 | F_030 S_000 )^0_{t} = y * ( G_220 S_000 | D_020 S_000 )^0 + ( F_210 S_000 | D_020 S_000 )^0 + ( G_220 S_000 | P_010 S_000 )^0 - ( H_230 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[36] = etfac[1] * AUX_INT__g_s_d_s[21] + 2 * one_over_2q * AUX_INT__f_s_d_s[9] + 2 * one_over_2q * AUX_INT__g_s_p_s[10] - p_over_q * AUX_INT__h_s_d_s[39];

                    // ( G_220 S_000 | F_021 S_000 )^0_{t} = z * ( G_220 S_000 | D_020 S_000 )^0 - ( H_221 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[37] = etfac[2] * AUX_INT__g_s_d_s[21] - p_over_q * AUX_INT__h_s_d_s[45];

                    // ( G_220 S_000 | F_012 S_000 )^0_{t} = y * ( G_220 S_000 | D_002 S_000 )^0 + ( F_210 S_000 | D_002 S_000 )^0 - ( H_230 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[38] = etfac[1] * AUX_INT__g_s_d_s[23] + 2 * one_over_2q * AUX_INT__f_s_d_s[11] - p_over_q * AUX_INT__h_s_d_s[41];

                    // ( G_220 S_000 | F_003 S_000 )^0_{t} = z * ( G_220 S_000 | D_002 S_000 )^0 + ( G_220 S_000 | P_001 S_000 )^0 - ( H_221 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[39] = etfac[2] * AUX_INT__g_s_d_s[23] + 2 * one_over_2q * AUX_INT__g_s_p_s[11] - p_over_q * AUX_INT__h_s_d_s[47];

                    // ( G_211 S_000 | F_300 S_000 )^0_{t} = x * ( G_211 S_000 | D_200 S_000 )^0 + ( F_111 S_000 | D_200 S_000 )^0 + ( G_211 S_000 | P_100 S_000 )^0 - ( H_311 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[40] = etfac[0] * AUX_INT__g_s_d_s[24] + 2 * one_over_2q * AUX_INT__f_s_d_s[24] + 2 * one_over_2q * AUX_INT__g_s_p_s[12] - p_over_q * AUX_INT__h_s_d_s[24];

                    // ( G_211 S_000 | F_210 S_000 )^0_{t} = y * ( G_211 S_000 | D_200 S_000 )^0 + ( F_201 S_000 | D_200 S_000 )^0 - ( H_221 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[41] = etfac[1] * AUX_INT__g_s_d_s[24] + 1 * one_over_2q * AUX_INT__f_s_d_s[12] - p_over_q * AUX_INT__h_s_d_s[42];

                    // ( G_211 S_000 | F_201 S_000 )^0_{t} = z * ( G_211 S_000 | D_200 S_000 )^0 + ( F_210 S_000 | D_200 S_000 )^0 - ( H_212 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[42] = etfac[2] * AUX_INT__g_s_d_s[24] + 1 * one_over_2q * AUX_INT__f_s_d_s[6] - p_over_q * AUX_INT__h_s_d_s[48];

                    // ( G_211 S_000 | F_120 S_000 )^0_{t} = x * ( G_211 S_000 | D_020 S_000 )^0 + ( F_111 S_000 | D_020 S_000 )^0 - ( H_311 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[43] = etfac[0] * AUX_INT__g_s_d_s[27] + 2 * one_over_2q * AUX_INT__f_s_d_s[27] - p_over_q * AUX_INT__h_s_d_s[27];

                    // ( G_211 S_000 | F_111 S_000 )^0_{t} = z * ( G_211 S_000 | D_110 S_000 )^0 + ( F_210 S_000 | D_110 S_000 )^0 - ( H_212 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[44] = etfac[2] * AUX_INT__g_s_d_s[25] + 1 * one_over_2q * AUX_INT__f_s_d_s[7] - p_over_q * AUX_INT__h_s_d_s[49];

                    // ( G_211 S_000 | F_102 S_000 )^0_{t} = x * ( G_211 S_000 | D_002 S_000 )^0 + ( F_111 S_000 | D_002 S_000 )^0 - ( H_311 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[45] = etfac[0] * AUX_INT__g_s_d_s[29] + 2 * one_over_2q * AUX_INT__f_s_d_s[29] - p_over_q * AUX_INT__h_s_d_s[29];

                    // ( G_211 S_000 | F_030 S_000 )^0_{t} = y * ( G_211 S_000 | D_020 S_000 )^0 + ( F_201 S_000 | D_020 S_000 )^0 + ( G_211 S_000 | P_010 S_000 )^0 - ( H_221 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[46] = etfac[1] * AUX_INT__g_s_d_s[27] + 1 * one_over_2q * AUX_INT__f_s_d_s[15] + 2 * one_over_2q * AUX_INT__g_s_p_s[13] - p_over_q * AUX_INT__h_s_d_s[45];

                    // ( G_211 S_000 | F_021 S_000 )^0_{t} = z * ( G_211 S_000 | D_020 S_000 )^0 + ( F_210 S_000 | D_020 S_000 )^0 - ( H_212 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[47] = etfac[2] * AUX_INT__g_s_d_s[27] + 1 * one_over_2q * AUX_INT__f_s_d_s[9] - p_over_q * AUX_INT__h_s_d_s[51];

                    // ( G_211 S_000 | F_012 S_000 )^0_{t} = y * ( G_211 S_000 | D_002 S_000 )^0 + ( F_201 S_000 | D_002 S_000 )^0 - ( H_221 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[48] = etfac[1] * AUX_INT__g_s_d_s[29] + 1 * one_over_2q * AUX_INT__f_s_d_s[17] - p_over_q * AUX_INT__h_s_d_s[47];

                    // ( G_211 S_000 | F_003 S_000 )^0_{t} = z * ( G_211 S_000 | D_002 S_000 )^0 + ( F_210 S_000 | D_002 S_000 )^0 + ( G_211 S_000 | P_001 S_000 )^0 - ( H_212 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[49] = etfac[2] * AUX_INT__g_s_d_s[29] + 1 * one_over_2q * AUX_INT__f_s_d_s[11] + 2 * one_over_2q * AUX_INT__g_s_p_s[14] - p_over_q * AUX_INT__h_s_d_s[53];

                    // ( G_202 S_000 | F_300 S_000 )^0_{t} = x * ( G_202 S_000 | D_200 S_000 )^0 + ( F_102 S_000 | D_200 S_000 )^0 + ( G_202 S_000 | P_100 S_000 )^0 - ( H_302 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[50] = etfac[0] * AUX_INT__g_s_d_s[30] + 2 * one_over_2q * AUX_INT__f_s_d_s[30] + 2 * one_over_2q * AUX_INT__g_s_p_s[15] - p_over_q * AUX_INT__h_s_d_s[30];

                    // ( G_202 S_000 | F_210 S_000 )^0_{t} = y * ( G_202 S_000 | D_200 S_000 )^0 - ( H_212 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[51] = etfac[1] * AUX_INT__g_s_d_s[30] - p_over_q * AUX_INT__h_s_d_s[48];

                    // ( G_202 S_000 | F_201 S_000 )^0_{t} = z * ( G_202 S_000 | D_200 S_000 )^0 + ( F_201 S_000 | D_200 S_000 )^0 - ( H_203 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[52] = etfac[2] * AUX_INT__g_s_d_s[30] + 2 * one_over_2q * AUX_INT__f_s_d_s[12] - p_over_q * AUX_INT__h_s_d_s[54];

                    // ( G_202 S_000 | F_120 S_000 )^0_{t} = x * ( G_202 S_000 | D_020 S_000 )^0 + ( F_102 S_000 | D_020 S_000 )^0 - ( H_302 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[53] = etfac[0] * AUX_INT__g_s_d_s[33] + 2 * one_over_2q * AUX_INT__f_s_d_s[33] - p_over_q * AUX_INT__h_s_d_s[33];

                    // ( G_202 S_000 | F_111 S_000 )^0_{t} = z * ( G_202 S_000 | D_110 S_000 )^0 + ( F_201 S_000 | D_110 S_000 )^0 - ( H_203 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[54] = etfac[2] * AUX_INT__g_s_d_s[31] + 2 * one_over_2q * AUX_INT__f_s_d_s[13] - p_over_q * AUX_INT__h_s_d_s[55];

                    // ( G_202 S_000 | F_102 S_000 )^0_{t} = x * ( G_202 S_000 | D_002 S_000 )^0 + ( F_102 S_000 | D_002 S_000 )^0 - ( H_302 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[55] = etfac[0] * AUX_INT__g_s_d_s[35] + 2 * one_over_2q * AUX_INT__f_s_d_s[35] - p_over_q * AUX_INT__h_s_d_s[35];

                    // ( G_202 S_000 | F_030 S_000 )^0_{t} = y * ( G_202 S_000 | D_020 S_000 )^0 + ( G_202 S_000 | P_010 S_000 )^0 - ( H_212 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[56] = etfac[1] * AUX_INT__g_s_d_s[33] + 2 * one_over_2q * AUX_INT__g_s_p_s[16] - p_over_q * AUX_INT__h_s_d_s[51];

                    // ( G_202 S_000 | F_021 S_000 )^0_{t} = z * ( G_202 S_000 | D_020 S_000 )^0 + ( F_201 S_000 | D_020 S_000 )^0 - ( H_203 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[57] = etfac[2] * AUX_INT__g_s_d_s[33] + 2 * one_over_2q * AUX_INT__f_s_d_s[15] - p_over_q * AUX_INT__h_s_d_s[57];

                    // ( G_202 S_000 | F_012 S_000 )^0_{t} = y * ( G_202 S_000 | D_002 S_000 )^0 - ( H_212 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[58] = etfac[1] * AUX_INT__g_s_d_s[35] - p_over_q * AUX_INT__h_s_d_s[53];

                    // ( G_202 S_000 | F_003 S_000 )^0_{t} = z * ( G_202 S_000 | D_002 S_000 )^0 + ( F_201 S_000 | D_002 S_000 )^0 + ( G_202 S_000 | P_001 S_000 )^0 - ( H_203 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[59] = etfac[2] * AUX_INT__g_s_d_s[35] + 2 * one_over_2q * AUX_INT__f_s_d_s[17] + 2 * one_over_2q * AUX_INT__g_s_p_s[17] - p_over_q * AUX_INT__h_s_d_s[59];

                    // ( G_130 S_000 | F_300 S_000 )^0_{t} = x * ( G_130 S_000 | D_200 S_000 )^0 + ( F_030 S_000 | D_200 S_000 )^0 + ( G_130 S_000 | P_100 S_000 )^0 - ( H_230 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[60] = etfac[0] * AUX_INT__g_s_d_s[36] + 1 * one_over_2q * AUX_INT__f_s_d_s[36] + 2 * one_over_2q * AUX_INT__g_s_p_s[18] - p_over_q * AUX_INT__h_s_d_s[36];

                    // ( G_130 S_000 | F_210 S_000 )^0_{t} = y * ( G_130 S_000 | D_200 S_000 )^0 + ( F_120 S_000 | D_200 S_000 )^0 - ( H_140 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[61] = etfac[1] * AUX_INT__g_s_d_s[36] + 3 * one_over_2q * AUX_INT__f_s_d_s[18] - p_over_q * AUX_INT__h_s_d_s[60];

                    // ( G_130 S_000 | F_201 S_000 )^0_{t} = z * ( G_130 S_000 | D_200 S_000 )^0 - ( H_131 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[62] = etfac[2] * AUX_INT__g_s_d_s[36] - p_over_q * AUX_INT__h_s_d_s[66];

                    // ( G_130 S_000 | F_120 S_000 )^0_{t} = x * ( G_130 S_000 | D_020 S_000 )^0 + ( F_030 S_000 | D_020 S_000 )^0 - ( H_230 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[63] = etfac[0] * AUX_INT__g_s_d_s[39] + 1 * one_over_2q * AUX_INT__f_s_d_s[39] - p_over_q * AUX_INT__h_s_d_s[39];

                    // ( G_130 S_000 | F_111 S_000 )^0_{t} = z * ( G_130 S_000 | D_110 S_000 )^0 - ( H_131 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[64] = etfac[2] * AUX_INT__g_s_d_s[37] - p_over_q * AUX_INT__h_s_d_s[67];

                    // ( G_130 S_000 | F_102 S_000 )^0_{t} = x * ( G_130 S_000 | D_002 S_000 )^0 + ( F_030 S_000 | D_002 S_000 )^0 - ( H_230 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[65] = etfac[0] * AUX_INT__g_s_d_s[41] + 1 * one_over_2q * AUX_INT__f_s_d_s[41] - p_over_q * AUX_INT__h_s_d_s[41];

                    // ( G_130 S_000 | F_030 S_000 )^0_{t} = y * ( G_130 S_000 | D_020 S_000 )^0 + ( F_120 S_000 | D_020 S_000 )^0 + ( G_130 S_000 | P_010 S_000 )^0 - ( H_140 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[66] = etfac[1] * AUX_INT__g_s_d_s[39] + 3 * one_over_2q * AUX_INT__f_s_d_s[21] + 2 * one_over_2q * AUX_INT__g_s_p_s[19] - p_over_q * AUX_INT__h_s_d_s[63];

                    // ( G_130 S_000 | F_021 S_000 )^0_{t} = z * ( G_130 S_000 | D_020 S_000 )^0 - ( H_131 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[67] = etfac[2] * AUX_INT__g_s_d_s[39] - p_over_q * AUX_INT__h_s_d_s[69];

                    // ( G_130 S_000 | F_012 S_000 )^0_{t} = y * ( G_130 S_000 | D_002 S_000 )^0 + ( F_120 S_000 | D_002 S_000 )^0 - ( H_140 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[68] = etfac[1] * AUX_INT__g_s_d_s[41] + 3 * one_over_2q * AUX_INT__f_s_d_s[23] - p_over_q * AUX_INT__h_s_d_s[65];

                    // ( G_130 S_000 | F_003 S_000 )^0_{t} = z * ( G_130 S_000 | D_002 S_000 )^0 + ( G_130 S_000 | P_001 S_000 )^0 - ( H_131 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[69] = etfac[2] * AUX_INT__g_s_d_s[41] + 2 * one_over_2q * AUX_INT__g_s_p_s[20] - p_over_q * AUX_INT__h_s_d_s[71];

                    // ( G_121 S_000 | F_300 S_000 )^0_{t} = x * ( G_121 S_000 | D_200 S_000 )^0 + ( F_021 S_000 | D_200 S_000 )^0 + ( G_121 S_000 | P_100 S_000 )^0 - ( H_221 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[70] = etfac[0] * AUX_INT__g_s_d_s[42] + 1 * one_over_2q * AUX_INT__f_s_d_s[42] + 2 * one_over_2q * AUX_INT__g_s_p_s[21] - p_over_q * AUX_INT__h_s_d_s[42];

                    // ( G_121 S_000 | F_210 S_000 )^0_{t} = y * ( G_121 S_000 | D_200 S_000 )^0 + ( F_111 S_000 | D_200 S_000 )^0 - ( H_131 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[71] = etfac[1] * AUX_INT__g_s_d_s[42] + 2 * one_over_2q * AUX_INT__f_s_d_s[24] - p_over_q * AUX_INT__h_s_d_s[66];

                    // ( G_121 S_000 | F_201 S_000 )^0_{t} = z * ( G_121 S_000 | D_200 S_000 )^0 + ( F_120 S_000 | D_200 S_000 )^0 - ( H_122 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[72] = etfac[2] * AUX_INT__g_s_d_s[42] + 1 * one_over_2q * AUX_INT__f_s_d_s[18] - p_over_q * AUX_INT__h_s_d_s[72];

                    // ( G_121 S_000 | F_120 S_000 )^0_{t} = x * ( G_121 S_000 | D_020 S_000 )^0 + ( F_021 S_000 | D_020 S_000 )^0 - ( H_221 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[73] = etfac[0] * AUX_INT__g_s_d_s[45] + 1 * one_over_2q * AUX_INT__f_s_d_s[45] - p_over_q * AUX_INT__h_s_d_s[45];

                    // ( G_121 S_000 | F_111 S_000 )^0_{t} = z * ( G_121 S_000 | D_110 S_000 )^0 + ( F_120 S_000 | D_110 S_000 )^0 - ( H_122 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[74] = etfac[2] * AUX_INT__g_s_d_s[43] + 1 * one_over_2q * AUX_INT__f_s_d_s[19] - p_over_q * AUX_INT__h_s_d_s[73];

                    // ( G_121 S_000 | F_102 S_000 )^0_{t} = x * ( G_121 S_000 | D_002 S_000 )^0 + ( F_021 S_000 | D_002 S_000 )^0 - ( H_221 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[75] = etfac[0] * AUX_INT__g_s_d_s[47] + 1 * one_over_2q * AUX_INT__f_s_d_s[47] - p_over_q * AUX_INT__h_s_d_s[47];

                    // ( G_121 S_000 | F_030 S_000 )^0_{t} = y * ( G_121 S_000 | D_020 S_000 )^0 + ( F_111 S_000 | D_020 S_000 )^0 + ( G_121 S_000 | P_010 S_000 )^0 - ( H_131 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[76] = etfac[1] * AUX_INT__g_s_d_s[45] + 2 * one_over_2q * AUX_INT__f_s_d_s[27] + 2 * one_over_2q * AUX_INT__g_s_p_s[22] - p_over_q * AUX_INT__h_s_d_s[69];

                    // ( G_121 S_000 | F_021 S_000 )^0_{t} = z * ( G_121 S_000 | D_020 S_000 )^0 + ( F_120 S_000 | D_020 S_000 )^0 - ( H_122 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[77] = etfac[2] * AUX_INT__g_s_d_s[45] + 1 * one_over_2q * AUX_INT__f_s_d_s[21] - p_over_q * AUX_INT__h_s_d_s[75];

                    // ( G_121 S_000 | F_012 S_000 )^0_{t} = y * ( G_121 S_000 | D_002 S_000 )^0 + ( F_111 S_000 | D_002 S_000 )^0 - ( H_131 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[78] = etfac[1] * AUX_INT__g_s_d_s[47] + 2 * one_over_2q * AUX_INT__f_s_d_s[29] - p_over_q * AUX_INT__h_s_d_s[71];

                    // ( G_121 S_000 | F_003 S_000 )^0_{t} = z * ( G_121 S_000 | D_002 S_000 )^0 + ( F_120 S_000 | D_002 S_000 )^0 + ( G_121 S_000 | P_001 S_000 )^0 - ( H_122 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[79] = etfac[2] * AUX_INT__g_s_d_s[47] + 1 * one_over_2q * AUX_INT__f_s_d_s[23] + 2 * one_over_2q * AUX_INT__g_s_p_s[23] - p_over_q * AUX_INT__h_s_d_s[77];

                    // ( G_112 S_000 | F_300 S_000 )^0_{t} = x * ( G_112 S_000 | D_200 S_000 )^0 + ( F_012 S_000 | D_200 S_000 )^0 + ( G_112 S_000 | P_100 S_000 )^0 - ( H_212 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[80] = etfac[0] * AUX_INT__g_s_d_s[48] + 1 * one_over_2q * AUX_INT__f_s_d_s[48] + 2 * one_over_2q * AUX_INT__g_s_p_s[24] - p_over_q * AUX_INT__h_s_d_s[48];

                    // ( G_112 S_000 | F_210 S_000 )^0_{t} = y * ( G_112 S_000 | D_200 S_000 )^0 + ( F_102 S_000 | D_200 S_000 )^0 - ( H_122 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[81] = etfac[1] * AUX_INT__g_s_d_s[48] + 1 * one_over_2q * AUX_INT__f_s_d_s[30] - p_over_q * AUX_INT__h_s_d_s[72];

                    // ( G_112 S_000 | F_201 S_000 )^0_{t} = z * ( G_112 S_000 | D_200 S_000 )^0 + ( F_111 S_000 | D_200 S_000 )^0 - ( H_113 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[82] = etfac[2] * AUX_INT__g_s_d_s[48] + 2 * one_over_2q * AUX_INT__f_s_d_s[24] - p_over_q * AUX_INT__h_s_d_s[78];

                    // ( G_112 S_000 | F_120 S_000 )^0_{t} = x * ( G_112 S_000 | D_020 S_000 )^0 + ( F_012 S_000 | D_020 S_000 )^0 - ( H_212 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[83] = etfac[0] * AUX_INT__g_s_d_s[51] + 1 * one_over_2q * AUX_INT__f_s_d_s[51] - p_over_q * AUX_INT__h_s_d_s[51];

                    // ( G_112 S_000 | F_111 S_000 )^0_{t} = z * ( G_112 S_000 | D_110 S_000 )^0 + ( F_111 S_000 | D_110 S_000 )^0 - ( H_113 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[84] = etfac[2] * AUX_INT__g_s_d_s[49] + 2 * one_over_2q * AUX_INT__f_s_d_s[25] - p_over_q * AUX_INT__h_s_d_s[79];

                    // ( G_112 S_000 | F_102 S_000 )^0_{t} = x * ( G_112 S_000 | D_002 S_000 )^0 + ( F_012 S_000 | D_002 S_000 )^0 - ( H_212 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[85] = etfac[0] * AUX_INT__g_s_d_s[53] + 1 * one_over_2q * AUX_INT__f_s_d_s[53] - p_over_q * AUX_INT__h_s_d_s[53];

                    // ( G_112 S_000 | F_030 S_000 )^0_{t} = y * ( G_112 S_000 | D_020 S_000 )^0 + ( F_102 S_000 | D_020 S_000 )^0 + ( G_112 S_000 | P_010 S_000 )^0 - ( H_122 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[86] = etfac[1] * AUX_INT__g_s_d_s[51] + 1 * one_over_2q * AUX_INT__f_s_d_s[33] + 2 * one_over_2q * AUX_INT__g_s_p_s[25] - p_over_q * AUX_INT__h_s_d_s[75];

                    // ( G_112 S_000 | F_021 S_000 )^0_{t} = z * ( G_112 S_000 | D_020 S_000 )^0 + ( F_111 S_000 | D_020 S_000 )^0 - ( H_113 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[87] = etfac[2] * AUX_INT__g_s_d_s[51] + 2 * one_over_2q * AUX_INT__f_s_d_s[27] - p_over_q * AUX_INT__h_s_d_s[81];

                    // ( G_112 S_000 | F_012 S_000 )^0_{t} = y * ( G_112 S_000 | D_002 S_000 )^0 + ( F_102 S_000 | D_002 S_000 )^0 - ( H_122 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[88] = etfac[1] * AUX_INT__g_s_d_s[53] + 1 * one_over_2q * AUX_INT__f_s_d_s[35] - p_over_q * AUX_INT__h_s_d_s[77];

                    // ( G_112 S_000 | F_003 S_000 )^0_{t} = z * ( G_112 S_000 | D_002 S_000 )^0 + ( F_111 S_000 | D_002 S_000 )^0 + ( G_112 S_000 | P_001 S_000 )^0 - ( H_113 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[89] = etfac[2] * AUX_INT__g_s_d_s[53] + 2 * one_over_2q * AUX_INT__f_s_d_s[29] + 2 * one_over_2q * AUX_INT__g_s_p_s[26] - p_over_q * AUX_INT__h_s_d_s[83];

                    // ( G_103 S_000 | F_300 S_000 )^0_{t} = x * ( G_103 S_000 | D_200 S_000 )^0 + ( F_003 S_000 | D_200 S_000 )^0 + ( G_103 S_000 | P_100 S_000 )^0 - ( H_203 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[90] = etfac[0] * AUX_INT__g_s_d_s[54] + 1 * one_over_2q * AUX_INT__f_s_d_s[54] + 2 * one_over_2q * AUX_INT__g_s_p_s[27] - p_over_q * AUX_INT__h_s_d_s[54];

                    // ( G_103 S_000 | F_210 S_000 )^0_{t} = y * ( G_103 S_000 | D_200 S_000 )^0 - ( H_113 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[91] = etfac[1] * AUX_INT__g_s_d_s[54] - p_over_q * AUX_INT__h_s_d_s[78];

                    // ( G_103 S_000 | F_201 S_000 )^0_{t} = z * ( G_103 S_000 | D_200 S_000 )^0 + ( F_102 S_000 | D_200 S_000 )^0 - ( H_104 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[92] = etfac[2] * AUX_INT__g_s_d_s[54] + 3 * one_over_2q * AUX_INT__f_s_d_s[30] - p_over_q * AUX_INT__h_s_d_s[84];

                    // ( G_103 S_000 | F_120 S_000 )^0_{t} = x * ( G_103 S_000 | D_020 S_000 )^0 + ( F_003 S_000 | D_020 S_000 )^0 - ( H_203 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[93] = etfac[0] * AUX_INT__g_s_d_s[57] + 1 * one_over_2q * AUX_INT__f_s_d_s[57] - p_over_q * AUX_INT__h_s_d_s[57];

                    // ( G_103 S_000 | F_111 S_000 )^0_{t} = z * ( G_103 S_000 | D_110 S_000 )^0 + ( F_102 S_000 | D_110 S_000 )^0 - ( H_104 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[94] = etfac[2] * AUX_INT__g_s_d_s[55] + 3 * one_over_2q * AUX_INT__f_s_d_s[31] - p_over_q * AUX_INT__h_s_d_s[85];

                    // ( G_103 S_000 | F_102 S_000 )^0_{t} = x * ( G_103 S_000 | D_002 S_000 )^0 + ( F_003 S_000 | D_002 S_000 )^0 - ( H_203 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[95] = etfac[0] * AUX_INT__g_s_d_s[59] + 1 * one_over_2q * AUX_INT__f_s_d_s[59] - p_over_q * AUX_INT__h_s_d_s[59];

                    // ( G_103 S_000 | F_030 S_000 )^0_{t} = y * ( G_103 S_000 | D_020 S_000 )^0 + ( G_103 S_000 | P_010 S_000 )^0 - ( H_113 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[96] = etfac[1] * AUX_INT__g_s_d_s[57] + 2 * one_over_2q * AUX_INT__g_s_p_s[28] - p_over_q * AUX_INT__h_s_d_s[81];

                    // ( G_103 S_000 | F_021 S_000 )^0_{t} = z * ( G_103 S_000 | D_020 S_000 )^0 + ( F_102 S_000 | D_020 S_000 )^0 - ( H_104 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[97] = etfac[2] * AUX_INT__g_s_d_s[57] + 3 * one_over_2q * AUX_INT__f_s_d_s[33] - p_over_q * AUX_INT__h_s_d_s[87];

                    // ( G_103 S_000 | F_012 S_000 )^0_{t} = y * ( G_103 S_000 | D_002 S_000 )^0 - ( H_113 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[98] = etfac[1] * AUX_INT__g_s_d_s[59] - p_over_q * AUX_INT__h_s_d_s[83];

                    // ( G_103 S_000 | F_003 S_000 )^0_{t} = z * ( G_103 S_000 | D_002 S_000 )^0 + ( F_102 S_000 | D_002 S_000 )^0 + ( G_103 S_000 | P_001 S_000 )^0 - ( H_104 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[99] = etfac[2] * AUX_INT__g_s_d_s[59] + 3 * one_over_2q * AUX_INT__f_s_d_s[35] + 2 * one_over_2q * AUX_INT__g_s_p_s[29] - p_over_q * AUX_INT__h_s_d_s[89];

                    // ( G_040 S_000 | F_300 S_000 )^0_{t} = x * ( G_040 S_000 | D_200 S_000 )^0 + ( G_040 S_000 | P_100 S_000 )^0 - ( H_140 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[100] = etfac[0] * AUX_INT__g_s_d_s[60] + 2 * one_over_2q * AUX_INT__g_s_p_s[30] - p_over_q * AUX_INT__h_s_d_s[60];

                    // ( G_040 S_000 | F_210 S_000 )^0_{t} = y * ( G_040 S_000 | D_200 S_000 )^0 + ( F_030 S_000 | D_200 S_000 )^0 - ( H_050 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[101] = etfac[1] * AUX_INT__g_s_d_s[60] + 4 * one_over_2q * AUX_INT__f_s_d_s[36] - p_over_q * AUX_INT__h_s_d_s[90];

                    // ( G_040 S_000 | F_201 S_000 )^0_{t} = z * ( G_040 S_000 | D_200 S_000 )^0 - ( H_041 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[102] = etfac[2] * AUX_INT__g_s_d_s[60] - p_over_q * AUX_INT__h_s_d_s[96];

                    // ( G_040 S_000 | F_120 S_000 )^0_{t} = x * ( G_040 S_000 | D_020 S_000 )^0 - ( H_140 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[103] = etfac[0] * AUX_INT__g_s_d_s[63] - p_over_q * AUX_INT__h_s_d_s[63];

                    // ( G_040 S_000 | F_111 S_000 )^0_{t} = z * ( G_040 S_000 | D_110 S_000 )^0 - ( H_041 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[104] = etfac[2] * AUX_INT__g_s_d_s[61] - p_over_q * AUX_INT__h_s_d_s[97];

                    // ( G_040 S_000 | F_102 S_000 )^0_{t} = x * ( G_040 S_000 | D_002 S_000 )^0 - ( H_140 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[105] = etfac[0] * AUX_INT__g_s_d_s[65] - p_over_q * AUX_INT__h_s_d_s[65];

                    // ( G_040 S_000 | F_030 S_000 )^0_{t} = y * ( G_040 S_000 | D_020 S_000 )^0 + ( F_030 S_000 | D_020 S_000 )^0 + ( G_040 S_000 | P_010 S_000 )^0 - ( H_050 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[106] = etfac[1] * AUX_INT__g_s_d_s[63] + 4 * one_over_2q * AUX_INT__f_s_d_s[39] + 2 * one_over_2q * AUX_INT__g_s_p_s[31] - p_over_q * AUX_INT__h_s_d_s[93];

                    // ( G_040 S_000 | F_021 S_000 )^0_{t} = z * ( G_040 S_000 | D_020 S_000 )^0 - ( H_041 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[107] = etfac[2] * AUX_INT__g_s_d_s[63] - p_over_q * AUX_INT__h_s_d_s[99];

                    // ( G_040 S_000 | F_012 S_000 )^0_{t} = y * ( G_040 S_000 | D_002 S_000 )^0 + ( F_030 S_000 | D_002 S_000 )^0 - ( H_050 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[108] = etfac[1] * AUX_INT__g_s_d_s[65] + 4 * one_over_2q * AUX_INT__f_s_d_s[41] - p_over_q * AUX_INT__h_s_d_s[95];

                    // ( G_040 S_000 | F_003 S_000 )^0_{t} = z * ( G_040 S_000 | D_002 S_000 )^0 + ( G_040 S_000 | P_001 S_000 )^0 - ( H_041 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[109] = etfac[2] * AUX_INT__g_s_d_s[65] + 2 * one_over_2q * AUX_INT__g_s_p_s[32] - p_over_q * AUX_INT__h_s_d_s[101];

                    // ( G_031 S_000 | F_300 S_000 )^0_{t} = x * ( G_031 S_000 | D_200 S_000 )^0 + ( G_031 S_000 | P_100 S_000 )^0 - ( H_131 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[110] = etfac[0] * AUX_INT__g_s_d_s[66] + 2 * one_over_2q * AUX_INT__g_s_p_s[33] - p_over_q * AUX_INT__h_s_d_s[66];

                    // ( G_031 S_000 | F_210 S_000 )^0_{t} = y * ( G_031 S_000 | D_200 S_000 )^0 + ( F_021 S_000 | D_200 S_000 )^0 - ( H_041 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[111] = etfac[1] * AUX_INT__g_s_d_s[66] + 3 * one_over_2q * AUX_INT__f_s_d_s[42] - p_over_q * AUX_INT__h_s_d_s[96];

                    // ( G_031 S_000 | F_201 S_000 )^0_{t} = z * ( G_031 S_000 | D_200 S_000 )^0 + ( F_030 S_000 | D_200 S_000 )^0 - ( H_032 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[112] = etfac[2] * AUX_INT__g_s_d_s[66] + 1 * one_over_2q * AUX_INT__f_s_d_s[36] - p_over_q * AUX_INT__h_s_d_s[102];

                    // ( G_031 S_000 | F_120 S_000 )^0_{t} = x * ( G_031 S_000 | D_020 S_000 )^0 - ( H_131 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[113] = etfac[0] * AUX_INT__g_s_d_s[69] - p_over_q * AUX_INT__h_s_d_s[69];

                    // ( G_031 S_000 | F_111 S_000 )^0_{t} = z * ( G_031 S_000 | D_110 S_000 )^0 + ( F_030 S_000 | D_110 S_000 )^0 - ( H_032 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[114] = etfac[2] * AUX_INT__g_s_d_s[67] + 1 * one_over_2q * AUX_INT__f_s_d_s[37] - p_over_q * AUX_INT__h_s_d_s[103];

                    // ( G_031 S_000 | F_102 S_000 )^0_{t} = x * ( G_031 S_000 | D_002 S_000 )^0 - ( H_131 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[115] = etfac[0] * AUX_INT__g_s_d_s[71] - p_over_q * AUX_INT__h_s_d_s[71];

                    // ( G_031 S_000 | F_030 S_000 )^0_{t} = y * ( G_031 S_000 | D_020 S_000 )^0 + ( F_021 S_000 | D_020 S_000 )^0 + ( G_031 S_000 | P_010 S_000 )^0 - ( H_041 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[116] = etfac[1] * AUX_INT__g_s_d_s[69] + 3 * one_over_2q * AUX_INT__f_s_d_s[45] + 2 * one_over_2q * AUX_INT__g_s_p_s[34] - p_over_q * AUX_INT__h_s_d_s[99];

                    // ( G_031 S_000 | F_021 S_000 )^0_{t} = z * ( G_031 S_000 | D_020 S_000 )^0 + ( F_030 S_000 | D_020 S_000 )^0 - ( H_032 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[117] = etfac[2] * AUX_INT__g_s_d_s[69] + 1 * one_over_2q * AUX_INT__f_s_d_s[39] - p_over_q * AUX_INT__h_s_d_s[105];

                    // ( G_031 S_000 | F_012 S_000 )^0_{t} = y * ( G_031 S_000 | D_002 S_000 )^0 + ( F_021 S_000 | D_002 S_000 )^0 - ( H_041 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[118] = etfac[1] * AUX_INT__g_s_d_s[71] + 3 * one_over_2q * AUX_INT__f_s_d_s[47] - p_over_q * AUX_INT__h_s_d_s[101];

                    // ( G_031 S_000 | F_003 S_000 )^0_{t} = z * ( G_031 S_000 | D_002 S_000 )^0 + ( F_030 S_000 | D_002 S_000 )^0 + ( G_031 S_000 | P_001 S_000 )^0 - ( H_032 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[119] = etfac[2] * AUX_INT__g_s_d_s[71] + 1 * one_over_2q * AUX_INT__f_s_d_s[41] + 2 * one_over_2q * AUX_INT__g_s_p_s[35] - p_over_q * AUX_INT__h_s_d_s[107];

                    // ( G_022 S_000 | F_300 S_000 )^0_{t} = x * ( G_022 S_000 | D_200 S_000 )^0 + ( G_022 S_000 | P_100 S_000 )^0 - ( H_122 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[120] = etfac[0] * AUX_INT__g_s_d_s[72] + 2 * one_over_2q * AUX_INT__g_s_p_s[36] - p_over_q * AUX_INT__h_s_d_s[72];

                    // ( G_022 S_000 | F_210 S_000 )^0_{t} = y * ( G_022 S_000 | D_200 S_000 )^0 + ( F_012 S_000 | D_200 S_000 )^0 - ( H_032 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[121] = etfac[1] * AUX_INT__g_s_d_s[72] + 2 * one_over_2q * AUX_INT__f_s_d_s[48] - p_over_q * AUX_INT__h_s_d_s[102];

                    // ( G_022 S_000 | F_201 S_000 )^0_{t} = z * ( G_022 S_000 | D_200 S_000 )^0 + ( F_021 S_000 | D_200 S_000 )^0 - ( H_023 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[122] = etfac[2] * AUX_INT__g_s_d_s[72] + 2 * one_over_2q * AUX_INT__f_s_d_s[42] - p_over_q * AUX_INT__h_s_d_s[108];

                    // ( G_022 S_000 | F_120 S_000 )^0_{t} = x * ( G_022 S_000 | D_020 S_000 )^0 - ( H_122 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[123] = etfac[0] * AUX_INT__g_s_d_s[75] - p_over_q * AUX_INT__h_s_d_s[75];

                    // ( G_022 S_000 | F_111 S_000 )^0_{t} = z * ( G_022 S_000 | D_110 S_000 )^0 + ( F_021 S_000 | D_110 S_000 )^0 - ( H_023 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[124] = etfac[2] * AUX_INT__g_s_d_s[73] + 2 * one_over_2q * AUX_INT__f_s_d_s[43] - p_over_q * AUX_INT__h_s_d_s[109];

                    // ( G_022 S_000 | F_102 S_000 )^0_{t} = x * ( G_022 S_000 | D_002 S_000 )^0 - ( H_122 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[125] = etfac[0] * AUX_INT__g_s_d_s[77] - p_over_q * AUX_INT__h_s_d_s[77];

                    // ( G_022 S_000 | F_030 S_000 )^0_{t} = y * ( G_022 S_000 | D_020 S_000 )^0 + ( F_012 S_000 | D_020 S_000 )^0 + ( G_022 S_000 | P_010 S_000 )^0 - ( H_032 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[126] = etfac[1] * AUX_INT__g_s_d_s[75] + 2 * one_over_2q * AUX_INT__f_s_d_s[51] + 2 * one_over_2q * AUX_INT__g_s_p_s[37] - p_over_q * AUX_INT__h_s_d_s[105];

                    // ( G_022 S_000 | F_021 S_000 )^0_{t} = z * ( G_022 S_000 | D_020 S_000 )^0 + ( F_021 S_000 | D_020 S_000 )^0 - ( H_023 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[127] = etfac[2] * AUX_INT__g_s_d_s[75] + 2 * one_over_2q * AUX_INT__f_s_d_s[45] - p_over_q * AUX_INT__h_s_d_s[111];

                    // ( G_022 S_000 | F_012 S_000 )^0_{t} = y * ( G_022 S_000 | D_002 S_000 )^0 + ( F_012 S_000 | D_002 S_000 )^0 - ( H_032 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[128] = etfac[1] * AUX_INT__g_s_d_s[77] + 2 * one_over_2q * AUX_INT__f_s_d_s[53] - p_over_q * AUX_INT__h_s_d_s[107];

                    // ( G_022 S_000 | F_003 S_000 )^0_{t} = z * ( G_022 S_000 | D_002 S_000 )^0 + ( F_021 S_000 | D_002 S_000 )^0 + ( G_022 S_000 | P_001 S_000 )^0 - ( H_023 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[129] = etfac[2] * AUX_INT__g_s_d_s[77] + 2 * one_over_2q * AUX_INT__f_s_d_s[47] + 2 * one_over_2q * AUX_INT__g_s_p_s[38] - p_over_q * AUX_INT__h_s_d_s[113];

                    // ( G_013 S_000 | F_300 S_000 )^0_{t} = x * ( G_013 S_000 | D_200 S_000 )^0 + ( G_013 S_000 | P_100 S_000 )^0 - ( H_113 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[130] = etfac[0] * AUX_INT__g_s_d_s[78] + 2 * one_over_2q * AUX_INT__g_s_p_s[39] - p_over_q * AUX_INT__h_s_d_s[78];

                    // ( G_013 S_000 | F_210 S_000 )^0_{t} = y * ( G_013 S_000 | D_200 S_000 )^0 + ( F_003 S_000 | D_200 S_000 )^0 - ( H_023 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[131] = etfac[1] * AUX_INT__g_s_d_s[78] + 1 * one_over_2q * AUX_INT__f_s_d_s[54] - p_over_q * AUX_INT__h_s_d_s[108];

                    // ( G_013 S_000 | F_201 S_000 )^0_{t} = z * ( G_013 S_000 | D_200 S_000 )^0 + ( F_012 S_000 | D_200 S_000 )^0 - ( H_014 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[132] = etfac[2] * AUX_INT__g_s_d_s[78] + 3 * one_over_2q * AUX_INT__f_s_d_s[48] - p_over_q * AUX_INT__h_s_d_s[114];

                    // ( G_013 S_000 | F_120 S_000 )^0_{t} = x * ( G_013 S_000 | D_020 S_000 )^0 - ( H_113 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[133] = etfac[0] * AUX_INT__g_s_d_s[81] - p_over_q * AUX_INT__h_s_d_s[81];

                    // ( G_013 S_000 | F_111 S_000 )^0_{t} = z * ( G_013 S_000 | D_110 S_000 )^0 + ( F_012 S_000 | D_110 S_000 )^0 - ( H_014 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[134] = etfac[2] * AUX_INT__g_s_d_s[79] + 3 * one_over_2q * AUX_INT__f_s_d_s[49] - p_over_q * AUX_INT__h_s_d_s[115];

                    // ( G_013 S_000 | F_102 S_000 )^0_{t} = x * ( G_013 S_000 | D_002 S_000 )^0 - ( H_113 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[135] = etfac[0] * AUX_INT__g_s_d_s[83] - p_over_q * AUX_INT__h_s_d_s[83];

                    // ( G_013 S_000 | F_030 S_000 )^0_{t} = y * ( G_013 S_000 | D_020 S_000 )^0 + ( F_003 S_000 | D_020 S_000 )^0 + ( G_013 S_000 | P_010 S_000 )^0 - ( H_023 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[136] = etfac[1] * AUX_INT__g_s_d_s[81] + 1 * one_over_2q * AUX_INT__f_s_d_s[57] + 2 * one_over_2q * AUX_INT__g_s_p_s[40] - p_over_q * AUX_INT__h_s_d_s[111];

                    // ( G_013 S_000 | F_021 S_000 )^0_{t} = z * ( G_013 S_000 | D_020 S_000 )^0 + ( F_012 S_000 | D_020 S_000 )^0 - ( H_014 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[137] = etfac[2] * AUX_INT__g_s_d_s[81] + 3 * one_over_2q * AUX_INT__f_s_d_s[51] - p_over_q * AUX_INT__h_s_d_s[117];

                    // ( G_013 S_000 | F_012 S_000 )^0_{t} = y * ( G_013 S_000 | D_002 S_000 )^0 + ( F_003 S_000 | D_002 S_000 )^0 - ( H_023 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[138] = etfac[1] * AUX_INT__g_s_d_s[83] + 1 * one_over_2q * AUX_INT__f_s_d_s[59] - p_over_q * AUX_INT__h_s_d_s[113];

                    // ( G_013 S_000 | F_003 S_000 )^0_{t} = z * ( G_013 S_000 | D_002 S_000 )^0 + ( F_012 S_000 | D_002 S_000 )^0 + ( G_013 S_000 | P_001 S_000 )^0 - ( H_014 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[139] = etfac[2] * AUX_INT__g_s_d_s[83] + 3 * one_over_2q * AUX_INT__f_s_d_s[53] + 2 * one_over_2q * AUX_INT__g_s_p_s[41] - p_over_q * AUX_INT__h_s_d_s[119];

                    // ( G_004 S_000 | F_300 S_000 )^0_{t} = x * ( G_004 S_000 | D_200 S_000 )^0 + ( G_004 S_000 | P_100 S_000 )^0 - ( H_104 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[140] = etfac[0] * AUX_INT__g_s_d_s[84] + 2 * one_over_2q * AUX_INT__g_s_p_s[42] - p_over_q * AUX_INT__h_s_d_s[84];

                    // ( G_004 S_000 | F_210 S_000 )^0_{t} = y * ( G_004 S_000 | D_200 S_000 )^0 - ( H_014 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[141] = etfac[1] * AUX_INT__g_s_d_s[84] - p_over_q * AUX_INT__h_s_d_s[114];

                    // ( G_004 S_000 | F_201 S_000 )^0_{t} = z * ( G_004 S_000 | D_200 S_000 )^0 + ( F_003 S_000 | D_200 S_000 )^0 - ( H_005 S_000 | D_200 S_000 )^0
                    AUX_INT__g_s_f_s[142] = etfac[2] * AUX_INT__g_s_d_s[84] + 4 * one_over_2q * AUX_INT__f_s_d_s[54] - p_over_q * AUX_INT__h_s_d_s[120];

                    // ( G_004 S_000 | F_120 S_000 )^0_{t} = x * ( G_004 S_000 | D_020 S_000 )^0 - ( H_104 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[143] = etfac[0] * AUX_INT__g_s_d_s[87] - p_over_q * AUX_INT__h_s_d_s[87];

                    // ( G_004 S_000 | F_111 S_000 )^0_{t} = z * ( G_004 S_000 | D_110 S_000 )^0 + ( F_003 S_000 | D_110 S_000 )^0 - ( H_005 S_000 | D_110 S_000 )^0
                    AUX_INT__g_s_f_s[144] = etfac[2] * AUX_INT__g_s_d_s[85] + 4 * one_over_2q * AUX_INT__f_s_d_s[55] - p_over_q * AUX_INT__h_s_d_s[121];

                    // ( G_004 S_000 | F_102 S_000 )^0_{t} = x * ( G_004 S_000 | D_002 S_000 )^0 - ( H_104 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[145] = etfac[0] * AUX_INT__g_s_d_s[89] - p_over_q * AUX_INT__h_s_d_s[89];

                    // ( G_004 S_000 | F_030 S_000 )^0_{t} = y * ( G_004 S_000 | D_020 S_000 )^0 + ( G_004 S_000 | P_010 S_000 )^0 - ( H_014 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[146] = etfac[1] * AUX_INT__g_s_d_s[87] + 2 * one_over_2q * AUX_INT__g_s_p_s[43] - p_over_q * AUX_INT__h_s_d_s[117];

                    // ( G_004 S_000 | F_021 S_000 )^0_{t} = z * ( G_004 S_000 | D_020 S_000 )^0 + ( F_003 S_000 | D_020 S_000 )^0 - ( H_005 S_000 | D_020 S_000 )^0
                    AUX_INT__g_s_f_s[147] = etfac[2] * AUX_INT__g_s_d_s[87] + 4 * one_over_2q * AUX_INT__f_s_d_s[57] - p_over_q * AUX_INT__h_s_d_s[123];

                    // ( G_004 S_000 | F_012 S_000 )^0_{t} = y * ( G_004 S_000 | D_002 S_000 )^0 - ( H_014 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[148] = etfac[1] * AUX_INT__g_s_d_s[89] - p_over_q * AUX_INT__h_s_d_s[119];

                    // ( G_004 S_000 | F_003 S_000 )^0_{t} = z * ( G_004 S_000 | D_002 S_000 )^0 + ( F_003 S_000 | D_002 S_000 )^0 + ( G_004 S_000 | P_001 S_000 )^0 - ( H_005 S_000 | D_002 S_000 )^0
                    AUX_INT__g_s_f_s[149] = etfac[2] * AUX_INT__g_s_d_s[89] + 4 * one_over_2q * AUX_INT__f_s_d_s[59] + 2 * one_over_2q * AUX_INT__g_s_p_s[44] - p_over_q * AUX_INT__h_s_d_s[125];

                    // ( F_300 S_000 | G_400 S_000 )^0_{t} = x * ( F_300 S_000 | F_300 S_000 )^0 + ( D_200 S_000 | F_300 S_000 )^0 + ( F_300 S_000 | D_200 S_000 )^0 - ( G_400 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[0] = etfac[0] * AUX_INT__f_s_f_s[0] + 3 * one_over_2q * AUX_INT__d_s_f_s[0] + 3 * one_over_2q * AUX_INT__f_s_d_s[0] - p_over_q * AUX_INT__g_s_f_s[0];

                    // ( F_300 S_000 | G_310 S_000 )^0_{t} = y * ( F_300 S_000 | F_300 S_000 )^0 - ( G_310 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[1] = etfac[1] * AUX_INT__f_s_f_s[0] - p_over_q * AUX_INT__g_s_f_s[10];

                    // ( F_300 S_000 | G_301 S_000 )^0_{t} = z * ( F_300 S_000 | F_300 S_000 )^0 - ( G_301 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[2] = etfac[2] * AUX_INT__f_s_f_s[0] - p_over_q * AUX_INT__g_s_f_s[20];

                    // ( F_300 S_000 | G_220 S_000 )^0_{t} = y * ( F_300 S_000 | F_210 S_000 )^0 + ( F_300 S_000 | D_200 S_000 )^0 - ( G_310 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[3] = etfac[1] * AUX_INT__f_s_f_s[1] + 1 * one_over_2q * AUX_INT__f_s_d_s[0] - p_over_q * AUX_INT__g_s_f_s[11];

                    // ( F_300 S_000 | G_211 S_000 )^0_{t} = z * ( F_300 S_000 | F_210 S_000 )^0 - ( G_301 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[4] = etfac[2] * AUX_INT__f_s_f_s[1] - p_over_q * AUX_INT__g_s_f_s[21];

                    // ( F_300 S_000 | G_202 S_000 )^0_{t} = z * ( F_300 S_000 | F_201 S_000 )^0 + ( F_300 S_000 | D_200 S_000 )^0 - ( G_301 S_000 | F_201 S_000 )^0
                    AUX_INT__f_s_g_s[5] = etfac[2] * AUX_INT__f_s_f_s[2] + 1 * one_over_2q * AUX_INT__f_s_d_s[0] - p_over_q * AUX_INT__g_s_f_s[22];

                    // ( F_300 S_000 | G_130 S_000 )^0_{t} = x * ( F_300 S_000 | F_030 S_000 )^0 + ( D_200 S_000 | F_030 S_000 )^0 - ( G_400 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[6] = etfac[0] * AUX_INT__f_s_f_s[6] + 3 * one_over_2q * AUX_INT__d_s_f_s[6] - p_over_q * AUX_INT__g_s_f_s[6];

                    // ( F_300 S_000 | G_121 S_000 )^0_{t} = z * ( F_300 S_000 | F_120 S_000 )^0 - ( G_301 S_000 | F_120 S_000 )^0
                    AUX_INT__f_s_g_s[7] = etfac[2] * AUX_INT__f_s_f_s[3] - p_over_q * AUX_INT__g_s_f_s[23];

                    // ( F_300 S_000 | G_112 S_000 )^0_{t} = y * ( F_300 S_000 | F_102 S_000 )^0 - ( G_310 S_000 | F_102 S_000 )^0
                    AUX_INT__f_s_g_s[8] = etfac[1] * AUX_INT__f_s_f_s[5] - p_over_q * AUX_INT__g_s_f_s[15];

                    // ( F_300 S_000 | G_103 S_000 )^0_{t} = x * ( F_300 S_000 | F_003 S_000 )^0 + ( D_200 S_000 | F_003 S_000 )^0 - ( G_400 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[9] = etfac[0] * AUX_INT__f_s_f_s[9] + 3 * one_over_2q * AUX_INT__d_s_f_s[9] - p_over_q * AUX_INT__g_s_f_s[9];

                    // ( F_300 S_000 | G_040 S_000 )^0_{t} = y * ( F_300 S_000 | F_030 S_000 )^0 + ( F_300 S_000 | D_020 S_000 )^0 - ( G_310 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[10] = etfac[1] * AUX_INT__f_s_f_s[6] + 3 * one_over_2q * AUX_INT__f_s_d_s[3] - p_over_q * AUX_INT__g_s_f_s[16];

                    // ( F_300 S_000 | G_031 S_000 )^0_{t} = z * ( F_300 S_000 | F_030 S_000 )^0 - ( G_301 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[11] = etfac[2] * AUX_INT__f_s_f_s[6] - p_over_q * AUX_INT__g_s_f_s[26];

                    // ( F_300 S_000 | G_022 S_000 )^0_{t} = z * ( F_300 S_000 | F_021 S_000 )^0 + ( F_300 S_000 | D_020 S_000 )^0 - ( G_301 S_000 | F_021 S_000 )^0
                    AUX_INT__f_s_g_s[12] = etfac[2] * AUX_INT__f_s_f_s[7] + 1 * one_over_2q * AUX_INT__f_s_d_s[3] - p_over_q * AUX_INT__g_s_f_s[27];

                    // ( F_300 S_000 | G_013 S_000 )^0_{t} = y * ( F_300 S_000 | F_003 S_000 )^0 - ( G_310 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[13] = etfac[1] * AUX_INT__f_s_f_s[9] - p_over_q * AUX_INT__g_s_f_s[19];

                    // ( F_300 S_000 | G_004 S_000 )^0_{t} = z * ( F_300 S_000 | F_003 S_000 )^0 + ( F_300 S_000 | D_002 S_000 )^0 - ( G_301 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[14] = etfac[2] * AUX_INT__f_s_f_s[9] + 3 * one_over_2q * AUX_INT__f_s_d_s[5] - p_over_q * AUX_INT__g_s_f_s[29];

                    // ( F_210 S_000 | G_400 S_000 )^0_{t} = x * ( F_210 S_000 | F_300 S_000 )^0 + ( D_110 S_000 | F_300 S_000 )^0 + ( F_210 S_000 | D_200 S_000 )^0 - ( G_310 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[15] = etfac[0] * AUX_INT__f_s_f_s[10] + 2 * one_over_2q * AUX_INT__d_s_f_s[10] + 3 * one_over_2q * AUX_INT__f_s_d_s[6] - p_over_q * AUX_INT__g_s_f_s[10];

                    // ( F_210 S_000 | G_310 S_000 )^0_{t} = y * ( F_210 S_000 | F_300 S_000 )^0 + ( D_200 S_000 | F_300 S_000 )^0 - ( G_220 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[16] = etfac[1] * AUX_INT__f_s_f_s[10] + 1 * one_over_2q * AUX_INT__d_s_f_s[0] - p_over_q * AUX_INT__g_s_f_s[30];

                    // ( F_210 S_000 | G_301 S_000 )^0_{t} = z * ( F_210 S_000 | F_300 S_000 )^0 - ( G_211 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[17] = etfac[2] * AUX_INT__f_s_f_s[10] - p_over_q * AUX_INT__g_s_f_s[40];

                    // ( F_210 S_000 | G_220 S_000 )^0_{t} = y * ( F_210 S_000 | F_210 S_000 )^0 + ( D_200 S_000 | F_210 S_000 )^0 + ( F_210 S_000 | D_200 S_000 )^0 - ( G_220 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[18] = etfac[1] * AUX_INT__f_s_f_s[11] + 1 * one_over_2q * AUX_INT__d_s_f_s[1] + 1 * one_over_2q * AUX_INT__f_s_d_s[6] - p_over_q * AUX_INT__g_s_f_s[31];

                    // ( F_210 S_000 | G_211 S_000 )^0_{t} = z * ( F_210 S_000 | F_210 S_000 )^0 - ( G_211 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[19] = etfac[2] * AUX_INT__f_s_f_s[11] - p_over_q * AUX_INT__g_s_f_s[41];

                    // ( F_210 S_000 | G_202 S_000 )^0_{t} = z * ( F_210 S_000 | F_201 S_000 )^0 + ( F_210 S_000 | D_200 S_000 )^0 - ( G_211 S_000 | F_201 S_000 )^0
                    AUX_INT__f_s_g_s[20] = etfac[2] * AUX_INT__f_s_f_s[12] + 1 * one_over_2q * AUX_INT__f_s_d_s[6] - p_over_q * AUX_INT__g_s_f_s[42];

                    // ( F_210 S_000 | G_130 S_000 )^0_{t} = x * ( F_210 S_000 | F_030 S_000 )^0 + ( D_110 S_000 | F_030 S_000 )^0 - ( G_310 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[21] = etfac[0] * AUX_INT__f_s_f_s[16] + 2 * one_over_2q * AUX_INT__d_s_f_s[16] - p_over_q * AUX_INT__g_s_f_s[16];

                    // ( F_210 S_000 | G_121 S_000 )^0_{t} = z * ( F_210 S_000 | F_120 S_000 )^0 - ( G_211 S_000 | F_120 S_000 )^0
                    AUX_INT__f_s_g_s[22] = etfac[2] * AUX_INT__f_s_f_s[13] - p_over_q * AUX_INT__g_s_f_s[43];

                    // ( F_210 S_000 | G_112 S_000 )^0_{t} = y * ( F_210 S_000 | F_102 S_000 )^0 + ( D_200 S_000 | F_102 S_000 )^0 - ( G_220 S_000 | F_102 S_000 )^0
                    AUX_INT__f_s_g_s[23] = etfac[1] * AUX_INT__f_s_f_s[15] + 1 * one_over_2q * AUX_INT__d_s_f_s[5] - p_over_q * AUX_INT__g_s_f_s[35];

                    // ( F_210 S_000 | G_103 S_000 )^0_{t} = x * ( F_210 S_000 | F_003 S_000 )^0 + ( D_110 S_000 | F_003 S_000 )^0 - ( G_310 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[24] = etfac[0] * AUX_INT__f_s_f_s[19] + 2 * one_over_2q * AUX_INT__d_s_f_s[19] - p_over_q * AUX_INT__g_s_f_s[19];

                    // ( F_210 S_000 | G_040 S_000 )^0_{t} = y * ( F_210 S_000 | F_030 S_000 )^0 + ( D_200 S_000 | F_030 S_000 )^0 + ( F_210 S_000 | D_020 S_000 )^0 - ( G_220 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[25] = etfac[1] * AUX_INT__f_s_f_s[16] + 1 * one_over_2q * AUX_INT__d_s_f_s[6] + 3 * one_over_2q * AUX_INT__f_s_d_s[9] - p_over_q * AUX_INT__g_s_f_s[36];

                    // ( F_210 S_000 | G_031 S_000 )^0_{t} = z * ( F_210 S_000 | F_030 S_000 )^0 - ( G_211 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[26] = etfac[2] * AUX_INT__f_s_f_s[16] - p_over_q * AUX_INT__g_s_f_s[46];

                    // ( F_210 S_000 | G_022 S_000 )^0_{t} = z * ( F_210 S_000 | F_021 S_000 )^0 + ( F_210 S_000 | D_020 S_000 )^0 - ( G_211 S_000 | F_021 S_000 )^0
                    AUX_INT__f_s_g_s[27] = etfac[2] * AUX_INT__f_s_f_s[17] + 1 * one_over_2q * AUX_INT__f_s_d_s[9] - p_over_q * AUX_INT__g_s_f_s[47];

                    // ( F_210 S_000 | G_013 S_000 )^0_{t} = y * ( F_210 S_000 | F_003 S_000 )^0 + ( D_200 S_000 | F_003 S_000 )^0 - ( G_220 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[28] = etfac[1] * AUX_INT__f_s_f_s[19] + 1 * one_over_2q * AUX_INT__d_s_f_s[9] - p_over_q * AUX_INT__g_s_f_s[39];

                    // ( F_210 S_000 | G_004 S_000 )^0_{t} = z * ( F_210 S_000 | F_003 S_000 )^0 + ( F_210 S_000 | D_002 S_000 )^0 - ( G_211 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[29] = etfac[2] * AUX_INT__f_s_f_s[19] + 3 * one_over_2q * AUX_INT__f_s_d_s[11] - p_over_q * AUX_INT__g_s_f_s[49];

                    // ( F_201 S_000 | G_400 S_000 )^0_{t} = x * ( F_201 S_000 | F_300 S_000 )^0 + ( D_101 S_000 | F_300 S_000 )^0 + ( F_201 S_000 | D_200 S_000 )^0 - ( G_301 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[30] = etfac[0] * AUX_INT__f_s_f_s[20] + 2 * one_over_2q * AUX_INT__d_s_f_s[20] + 3 * one_over_2q * AUX_INT__f_s_d_s[12] - p_over_q * AUX_INT__g_s_f_s[20];

                    // ( F_201 S_000 | G_310 S_000 )^0_{t} = y * ( F_201 S_000 | F_300 S_000 )^0 - ( G_211 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[31] = etfac[1] * AUX_INT__f_s_f_s[20] - p_over_q * AUX_INT__g_s_f_s[40];

                    // ( F_201 S_000 | G_301 S_000 )^0_{t} = z * ( F_201 S_000 | F_300 S_000 )^0 + ( D_200 S_000 | F_300 S_000 )^0 - ( G_202 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[32] = etfac[2] * AUX_INT__f_s_f_s[20] + 1 * one_over_2q * AUX_INT__d_s_f_s[0] - p_over_q * AUX_INT__g_s_f_s[50];

                    // ( F_201 S_000 | G_220 S_000 )^0_{t} = y * ( F_201 S_000 | F_210 S_000 )^0 + ( F_201 S_000 | D_200 S_000 )^0 - ( G_211 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[33] = etfac[1] * AUX_INT__f_s_f_s[21] + 1 * one_over_2q * AUX_INT__f_s_d_s[12] - p_over_q * AUX_INT__g_s_f_s[41];

                    // ( F_201 S_000 | G_211 S_000 )^0_{t} = z * ( F_201 S_000 | F_210 S_000 )^0 + ( D_200 S_000 | F_210 S_000 )^0 - ( G_202 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[34] = etfac[2] * AUX_INT__f_s_f_s[21] + 1 * one_over_2q * AUX_INT__d_s_f_s[1] - p_over_q * AUX_INT__g_s_f_s[51];

                    // ( F_201 S_000 | G_202 S_000 )^0_{t} = z * ( F_201 S_000 | F_201 S_000 )^0 + ( D_200 S_000 | F_201 S_000 )^0 + ( F_201 S_000 | D_200 S_000 )^0 - ( G_202 S_000 | F_201 S_000 )^0
                    AUX_INT__f_s_g_s[35] = etfac[2] * AUX_INT__f_s_f_s[22] + 1 * one_over_2q * AUX_INT__d_s_f_s[2] + 1 * one_over_2q * AUX_INT__f_s_d_s[12] - p_over_q * AUX_INT__g_s_f_s[52];

                    // ( F_201 S_000 | G_130 S_000 )^0_{t} = x * ( F_201 S_000 | F_030 S_000 )^0 + ( D_101 S_000 | F_030 S_000 )^0 - ( G_301 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[36] = etfac[0] * AUX_INT__f_s_f_s[26] + 2 * one_over_2q * AUX_INT__d_s_f_s[26] - p_over_q * AUX_INT__g_s_f_s[26];

                    // ( F_201 S_000 | G_121 S_000 )^0_{t} = z * ( F_201 S_000 | F_120 S_000 )^0 + ( D_200 S_000 | F_120 S_000 )^0 - ( G_202 S_000 | F_120 S_000 )^0
                    AUX_INT__f_s_g_s[37] = etfac[2] * AUX_INT__f_s_f_s[23] + 1 * one_over_2q * AUX_INT__d_s_f_s[3] - p_over_q * AUX_INT__g_s_f_s[53];

                    // ( F_201 S_000 | G_112 S_000 )^0_{t} = y * ( F_201 S_000 | F_102 S_000 )^0 - ( G_211 S_000 | F_102 S_000 )^0
                    AUX_INT__f_s_g_s[38] = etfac[1] * AUX_INT__f_s_f_s[25] - p_over_q * AUX_INT__g_s_f_s[45];

                    // ( F_201 S_000 | G_103 S_000 )^0_{t} = x * ( F_201 S_000 | F_003 S_000 )^0 + ( D_101 S_000 | F_003 S_000 )^0 - ( G_301 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[39] = etfac[0] * AUX_INT__f_s_f_s[29] + 2 * one_over_2q * AUX_INT__d_s_f_s[29] - p_over_q * AUX_INT__g_s_f_s[29];

                    // ( F_201 S_000 | G_040 S_000 )^0_{t} = y * ( F_201 S_000 | F_030 S_000 )^0 + ( F_201 S_000 | D_020 S_000 )^0 - ( G_211 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[40] = etfac[1] * AUX_INT__f_s_f_s[26] + 3 * one_over_2q * AUX_INT__f_s_d_s[15] - p_over_q * AUX_INT__g_s_f_s[46];

                    // ( F_201 S_000 | G_031 S_000 )^0_{t} = z * ( F_201 S_000 | F_030 S_000 )^0 + ( D_200 S_000 | F_030 S_000 )^0 - ( G_202 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[41] = etfac[2] * AUX_INT__f_s_f_s[26] + 1 * one_over_2q * AUX_INT__d_s_f_s[6] - p_over_q * AUX_INT__g_s_f_s[56];

                    // ( F_201 S_000 | G_022 S_000 )^0_{t} = z * ( F_201 S_000 | F_021 S_000 )^0 + ( D_200 S_000 | F_021 S_000 )^0 + ( F_201 S_000 | D_020 S_000 )^0 - ( G_202 S_000 | F_021 S_000 )^0
                    AUX_INT__f_s_g_s[42] = etfac[2] * AUX_INT__f_s_f_s[27] + 1 * one_over_2q * AUX_INT__d_s_f_s[7] + 1 * one_over_2q * AUX_INT__f_s_d_s[15] - p_over_q * AUX_INT__g_s_f_s[57];

                    // ( F_201 S_000 | G_013 S_000 )^0_{t} = y * ( F_201 S_000 | F_003 S_000 )^0 - ( G_211 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[43] = etfac[1] * AUX_INT__f_s_f_s[29] - p_over_q * AUX_INT__g_s_f_s[49];

                    // ( F_201 S_000 | G_004 S_000 )^0_{t} = z * ( F_201 S_000 | F_003 S_000 )^0 + ( D_200 S_000 | F_003 S_000 )^0 + ( F_201 S_000 | D_002 S_000 )^0 - ( G_202 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[44] = etfac[2] * AUX_INT__f_s_f_s[29] + 1 * one_over_2q * AUX_INT__d_s_f_s[9] + 3 * one_over_2q * AUX_INT__f_s_d_s[17] - p_over_q * AUX_INT__g_s_f_s[59];

                    // ( F_120 S_000 | G_400 S_000 )^0_{t} = x * ( F_120 S_000 | F_300 S_000 )^0 + ( D_020 S_000 | F_300 S_000 )^0 + ( F_120 S_000 | D_200 S_000 )^0 - ( G_220 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[45] = etfac[0] * AUX_INT__f_s_f_s[30] + 1 * one_over_2q * AUX_INT__d_s_f_s[30] + 3 * one_over_2q * AUX_INT__f_s_d_s[18] - p_over_q * AUX_INT__g_s_f_s[30];

                    // ( F_120 S_000 | G_310 S_000 )^0_{t} = y * ( F_120 S_000 | F_300 S_000 )^0 + ( D_110 S_000 | F_300 S_000 )^0 - ( G_130 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[46] = etfac[1] * AUX_INT__f_s_f_s[30] + 2 * one_over_2q * AUX_INT__d_s_f_s[10] - p_over_q * AUX_INT__g_s_f_s[60];

                    // ( F_120 S_000 | G_301 S_000 )^0_{t} = z * ( F_120 S_000 | F_300 S_000 )^0 - ( G_121 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[47] = etfac[2] * AUX_INT__f_s_f_s[30] - p_over_q * AUX_INT__g_s_f_s[70];

                    // ( F_120 S_000 | G_220 S_000 )^0_{t} = y * ( F_120 S_000 | F_210 S_000 )^0 + ( D_110 S_000 | F_210 S_000 )^0 + ( F_120 S_000 | D_200 S_000 )^0 - ( G_130 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[48] = etfac[1] * AUX_INT__f_s_f_s[31] + 2 * one_over_2q * AUX_INT__d_s_f_s[11] + 1 * one_over_2q * AUX_INT__f_s_d_s[18] - p_over_q * AUX_INT__g_s_f_s[61];

                    // ( F_120 S_000 | G_211 S_000 )^0_{t} = z * ( F_120 S_000 | F_210 S_000 )^0 - ( G_121 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[49] = etfac[2] * AUX_INT__f_s_f_s[31] - p_over_q * AUX_INT__g_s_f_s[71];

                    // ( F_120 S_000 | G_202 S_000 )^0_{t} = z * ( F_120 S_000 | F_201 S_000 )^0 + ( F_120 S_000 | D_200 S_000 )^0 - ( G_121 S_000 | F_201 S_000 )^0
                    AUX_INT__f_s_g_s[50] = etfac[2] * AUX_INT__f_s_f_s[32] + 1 * one_over_2q * AUX_INT__f_s_d_s[18] - p_over_q * AUX_INT__g_s_f_s[72];

                    // ( F_120 S_000 | G_130 S_000 )^0_{t} = x * ( F_120 S_000 | F_030 S_000 )^0 + ( D_020 S_000 | F_030 S_000 )^0 - ( G_220 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[51] = etfac[0] * AUX_INT__f_s_f_s[36] + 1 * one_over_2q * AUX_INT__d_s_f_s[36] - p_over_q * AUX_INT__g_s_f_s[36];

                    // ( F_120 S_000 | G_121 S_000 )^0_{t} = z * ( F_120 S_000 | F_120 S_000 )^0 - ( G_121 S_000 | F_120 S_000 )^0
                    AUX_INT__f_s_g_s[52] = etfac[2] * AUX_INT__f_s_f_s[33] - p_over_q * AUX_INT__g_s_f_s[73];

                    // ( F_120 S_000 | G_112 S_000 )^0_{t} = y * ( F_120 S_000 | F_102 S_000 )^0 + ( D_110 S_000 | F_102 S_000 )^0 - ( G_130 S_000 | F_102 S_000 )^0
                    AUX_INT__f_s_g_s[53] = etfac[1] * AUX_INT__f_s_f_s[35] + 2 * one_over_2q * AUX_INT__d_s_f_s[15] - p_over_q * AUX_INT__g_s_f_s[65];

                    // ( F_120 S_000 | G_103 S_000 )^0_{t} = x * ( F_120 S_000 | F_003 S_000 )^0 + ( D_020 S_000 | F_003 S_000 )^0 - ( G_220 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[54] = etfac[0] * AUX_INT__f_s_f_s[39] + 1 * one_over_2q * AUX_INT__d_s_f_s[39] - p_over_q * AUX_INT__g_s_f_s[39];

                    // ( F_120 S_000 | G_040 S_000 )^0_{t} = y * ( F_120 S_000 | F_030 S_000 )^0 + ( D_110 S_000 | F_030 S_000 )^0 + ( F_120 S_000 | D_020 S_000 )^0 - ( G_130 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[55] = etfac[1] * AUX_INT__f_s_f_s[36] + 2 * one_over_2q * AUX_INT__d_s_f_s[16] + 3 * one_over_2q * AUX_INT__f_s_d_s[21] - p_over_q * AUX_INT__g_s_f_s[66];

                    // ( F_120 S_000 | G_031 S_000 )^0_{t} = z * ( F_120 S_000 | F_030 S_000 )^0 - ( G_121 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[56] = etfac[2] * AUX_INT__f_s_f_s[36] - p_over_q * AUX_INT__g_s_f_s[76];

                    // ( F_120 S_000 | G_022 S_000 )^0_{t} = z * ( F_120 S_000 | F_021 S_000 )^0 + ( F_120 S_000 | D_020 S_000 )^0 - ( G_121 S_000 | F_021 S_000 )^0
                    AUX_INT__f_s_g_s[57] = etfac[2] * AUX_INT__f_s_f_s[37] + 1 * one_over_2q * AUX_INT__f_s_d_s[21] - p_over_q * AUX_INT__g_s_f_s[77];

                    // ( F_120 S_000 | G_013 S_000 )^0_{t} = y * ( F_120 S_000 | F_003 S_000 )^0 + ( D_110 S_000 | F_003 S_000 )^0 - ( G_130 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[58] = etfac[1] * AUX_INT__f_s_f_s[39] + 2 * one_over_2q * AUX_INT__d_s_f_s[19] - p_over_q * AUX_INT__g_s_f_s[69];

                    // ( F_120 S_000 | G_004 S_000 )^0_{t} = z * ( F_120 S_000 | F_003 S_000 )^0 + ( F_120 S_000 | D_002 S_000 )^0 - ( G_121 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[59] = etfac[2] * AUX_INT__f_s_f_s[39] + 3 * one_over_2q * AUX_INT__f_s_d_s[23] - p_over_q * AUX_INT__g_s_f_s[79];

                    // ( F_111 S_000 | G_400 S_000 )^0_{t} = x * ( F_111 S_000 | F_300 S_000 )^0 + ( D_011 S_000 | F_300 S_000 )^0 + ( F_111 S_000 | D_200 S_000 )^0 - ( G_211 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[60] = etfac[0] * AUX_INT__f_s_f_s[40] + 1 * one_over_2q * AUX_INT__d_s_f_s[40] + 3 * one_over_2q * AUX_INT__f_s_d_s[24] - p_over_q * AUX_INT__g_s_f_s[40];

                    // ( F_111 S_000 | G_310 S_000 )^0_{t} = y * ( F_111 S_000 | F_300 S_000 )^0 + ( D_101 S_000 | F_300 S_000 )^0 - ( G_121 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[61] = etfac[1] * AUX_INT__f_s_f_s[40] + 1 * one_over_2q * AUX_INT__d_s_f_s[20] - p_over_q * AUX_INT__g_s_f_s[70];

                    // ( F_111 S_000 | G_301 S_000 )^0_{t} = z * ( F_111 S_000 | F_300 S_000 )^0 + ( D_110 S_000 | F_300 S_000 )^0 - ( G_112 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[62] = etfac[2] * AUX_INT__f_s_f_s[40] + 1 * one_over_2q * AUX_INT__d_s_f_s[10] - p_over_q * AUX_INT__g_s_f_s[80];

                    // ( F_111 S_000 | G_220 S_000 )^0_{t} = y * ( F_111 S_000 | F_210 S_000 )^0 + ( D_101 S_000 | F_210 S_000 )^0 + ( F_111 S_000 | D_200 S_000 )^0 - ( G_121 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[63] = etfac[1] * AUX_INT__f_s_f_s[41] + 1 * one_over_2q * AUX_INT__d_s_f_s[21] + 1 * one_over_2q * AUX_INT__f_s_d_s[24] - p_over_q * AUX_INT__g_s_f_s[71];

                    // ( F_111 S_000 | G_211 S_000 )^0_{t} = z * ( F_111 S_000 | F_210 S_000 )^0 + ( D_110 S_000 | F_210 S_000 )^0 - ( G_112 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[64] = etfac[2] * AUX_INT__f_s_f_s[41] + 1 * one_over_2q * AUX_INT__d_s_f_s[11] - p_over_q * AUX_INT__g_s_f_s[81];

                    // ( F_111 S_000 | G_202 S_000 )^0_{t} = z * ( F_111 S_000 | F_201 S_000 )^0 + ( D_110 S_000 | F_201 S_000 )^0 + ( F_111 S_000 | D_200 S_000 )^0 - ( G_112 S_000 | F_201 S_000 )^0
                    AUX_INT__f_s_g_s[65] = etfac[2] * AUX_INT__f_s_f_s[42] + 1 * one_over_2q * AUX_INT__d_s_f_s[12] + 1 * one_over_2q * AUX_INT__f_s_d_s[24] - p_over_q * AUX_INT__g_s_f_s[82];

                    // ( F_111 S_000 | G_130 S_000 )^0_{t} = x * ( F_111 S_000 | F_030 S_000 )^0 + ( D_011 S_000 | F_030 S_000 )^0 - ( G_211 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[66] = etfac[0] * AUX_INT__f_s_f_s[46] + 1 * one_over_2q * AUX_INT__d_s_f_s[46] - p_over_q * AUX_INT__g_s_f_s[46];

                    // ( F_111 S_000 | G_121 S_000 )^0_{t} = z * ( F_111 S_000 | F_120 S_000 )^0 + ( D_110 S_000 | F_120 S_000 )^0 - ( G_112 S_000 | F_120 S_000 )^0
                    AUX_INT__f_s_g_s[67] = etfac[2] * AUX_INT__f_s_f_s[43] + 1 * one_over_2q * AUX_INT__d_s_f_s[13] - p_over_q * AUX_INT__g_s_f_s[83];

                    // ( F_111 S_000 | G_112 S_000 )^0_{t} = y * ( F_111 S_000 | F_102 S_000 )^0 + ( D_101 S_000 | F_102 S_000 )^0 - ( G_121 S_000 | F_102 S_000 )^0
                    AUX_INT__f_s_g_s[68] = etfac[1] * AUX_INT__f_s_f_s[45] + 1 * one_over_2q * AUX_INT__d_s_f_s[25] - p_over_q * AUX_INT__g_s_f_s[75];

                    // ( F_111 S_000 | G_103 S_000 )^0_{t} = x * ( F_111 S_000 | F_003 S_000 )^0 + ( D_011 S_000 | F_003 S_000 )^0 - ( G_211 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[69] = etfac[0] * AUX_INT__f_s_f_s[49] + 1 * one_over_2q * AUX_INT__d_s_f_s[49] - p_over_q * AUX_INT__g_s_f_s[49];

                    // ( F_111 S_000 | G_040 S_000 )^0_{t} = y * ( F_111 S_000 | F_030 S_000 )^0 + ( D_101 S_000 | F_030 S_000 )^0 + ( F_111 S_000 | D_020 S_000 )^0 - ( G_121 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[70] = etfac[1] * AUX_INT__f_s_f_s[46] + 1 * one_over_2q * AUX_INT__d_s_f_s[26] + 3 * one_over_2q * AUX_INT__f_s_d_s[27] - p_over_q * AUX_INT__g_s_f_s[76];

                    // ( F_111 S_000 | G_031 S_000 )^0_{t} = z * ( F_111 S_000 | F_030 S_000 )^0 + ( D_110 S_000 | F_030 S_000 )^0 - ( G_112 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[71] = etfac[2] * AUX_INT__f_s_f_s[46] + 1 * one_over_2q * AUX_INT__d_s_f_s[16] - p_over_q * AUX_INT__g_s_f_s[86];

                    // ( F_111 S_000 | G_022 S_000 )^0_{t} = z * ( F_111 S_000 | F_021 S_000 )^0 + ( D_110 S_000 | F_021 S_000 )^0 + ( F_111 S_000 | D_020 S_000 )^0 - ( G_112 S_000 | F_021 S_000 )^0
                    AUX_INT__f_s_g_s[72] = etfac[2] * AUX_INT__f_s_f_s[47] + 1 * one_over_2q * AUX_INT__d_s_f_s[17] + 1 * one_over_2q * AUX_INT__f_s_d_s[27] - p_over_q * AUX_INT__g_s_f_s[87];

                    // ( F_111 S_000 | G_013 S_000 )^0_{t} = y * ( F_111 S_000 | F_003 S_000 )^0 + ( D_101 S_000 | F_003 S_000 )^0 - ( G_121 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[73] = etfac[1] * AUX_INT__f_s_f_s[49] + 1 * one_over_2q * AUX_INT__d_s_f_s[29] - p_over_q * AUX_INT__g_s_f_s[79];

                    // ( F_111 S_000 | G_004 S_000 )^0_{t} = z * ( F_111 S_000 | F_003 S_000 )^0 + ( D_110 S_000 | F_003 S_000 )^0 + ( F_111 S_000 | D_002 S_000 )^0 - ( G_112 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[74] = etfac[2] * AUX_INT__f_s_f_s[49] + 1 * one_over_2q * AUX_INT__d_s_f_s[19] + 3 * one_over_2q * AUX_INT__f_s_d_s[29] - p_over_q * AUX_INT__g_s_f_s[89];

                    // ( F_102 S_000 | G_400 S_000 )^0_{t} = x * ( F_102 S_000 | F_300 S_000 )^0 + ( D_002 S_000 | F_300 S_000 )^0 + ( F_102 S_000 | D_200 S_000 )^0 - ( G_202 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[75] = etfac[0] * AUX_INT__f_s_f_s[50] + 1 * one_over_2q * AUX_INT__d_s_f_s[50] + 3 * one_over_2q * AUX_INT__f_s_d_s[30] - p_over_q * AUX_INT__g_s_f_s[50];

                    // ( F_102 S_000 | G_310 S_000 )^0_{t} = y * ( F_102 S_000 | F_300 S_000 )^0 - ( G_112 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[76] = etfac[1] * AUX_INT__f_s_f_s[50] - p_over_q * AUX_INT__g_s_f_s[80];

                    // ( F_102 S_000 | G_301 S_000 )^0_{t} = z * ( F_102 S_000 | F_300 S_000 )^0 + ( D_101 S_000 | F_300 S_000 )^0 - ( G_103 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[77] = etfac[2] * AUX_INT__f_s_f_s[50] + 2 * one_over_2q * AUX_INT__d_s_f_s[20] - p_over_q * AUX_INT__g_s_f_s[90];

                    // ( F_102 S_000 | G_220 S_000 )^0_{t} = y * ( F_102 S_000 | F_210 S_000 )^0 + ( F_102 S_000 | D_200 S_000 )^0 - ( G_112 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[78] = etfac[1] * AUX_INT__f_s_f_s[51] + 1 * one_over_2q * AUX_INT__f_s_d_s[30] - p_over_q * AUX_INT__g_s_f_s[81];

                    // ( F_102 S_000 | G_211 S_000 )^0_{t} = z * ( F_102 S_000 | F_210 S_000 )^0 + ( D_101 S_000 | F_210 S_000 )^0 - ( G_103 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[79] = etfac[2] * AUX_INT__f_s_f_s[51] + 2 * one_over_2q * AUX_INT__d_s_f_s[21] - p_over_q * AUX_INT__g_s_f_s[91];

                    // ( F_102 S_000 | G_202 S_000 )^0_{t} = z * ( F_102 S_000 | F_201 S_000 )^0 + ( D_101 S_000 | F_201 S_000 )^0 + ( F_102 S_000 | D_200 S_000 )^0 - ( G_103 S_000 | F_201 S_000 )^0
                    AUX_INT__f_s_g_s[80] = etfac[2] * AUX_INT__f_s_f_s[52] + 2 * one_over_2q * AUX_INT__d_s_f_s[22] + 1 * one_over_2q * AUX_INT__f_s_d_s[30] - p_over_q * AUX_INT__g_s_f_s[92];

                    // ( F_102 S_000 | G_130 S_000 )^0_{t} = x * ( F_102 S_000 | F_030 S_000 )^0 + ( D_002 S_000 | F_030 S_000 )^0 - ( G_202 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[81] = etfac[0] * AUX_INT__f_s_f_s[56] + 1 * one_over_2q * AUX_INT__d_s_f_s[56] - p_over_q * AUX_INT__g_s_f_s[56];

                    // ( F_102 S_000 | G_121 S_000 )^0_{t} = z * ( F_102 S_000 | F_120 S_000 )^0 + ( D_101 S_000 | F_120 S_000 )^0 - ( G_103 S_000 | F_120 S_000 )^0
                    AUX_INT__f_s_g_s[82] = etfac[2] * AUX_INT__f_s_f_s[53] + 2 * one_over_2q * AUX_INT__d_s_f_s[23] - p_over_q * AUX_INT__g_s_f_s[93];

                    // ( F_102 S_000 | G_112 S_000 )^0_{t} = y * ( F_102 S_000 | F_102 S_000 )^0 - ( G_112 S_000 | F_102 S_000 )^0
                    AUX_INT__f_s_g_s[83] = etfac[1] * AUX_INT__f_s_f_s[55] - p_over_q * AUX_INT__g_s_f_s[85];

                    // ( F_102 S_000 | G_103 S_000 )^0_{t} = x * ( F_102 S_000 | F_003 S_000 )^0 + ( D_002 S_000 | F_003 S_000 )^0 - ( G_202 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[84] = etfac[0] * AUX_INT__f_s_f_s[59] + 1 * one_over_2q * AUX_INT__d_s_f_s[59] - p_over_q * AUX_INT__g_s_f_s[59];

                    // ( F_102 S_000 | G_040 S_000 )^0_{t} = y * ( F_102 S_000 | F_030 S_000 )^0 + ( F_102 S_000 | D_020 S_000 )^0 - ( G_112 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[85] = etfac[1] * AUX_INT__f_s_f_s[56] + 3 * one_over_2q * AUX_INT__f_s_d_s[33] - p_over_q * AUX_INT__g_s_f_s[86];

                    // ( F_102 S_000 | G_031 S_000 )^0_{t} = z * ( F_102 S_000 | F_030 S_000 )^0 + ( D_101 S_000 | F_030 S_000 )^0 - ( G_103 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[86] = etfac[2] * AUX_INT__f_s_f_s[56] + 2 * one_over_2q * AUX_INT__d_s_f_s[26] - p_over_q * AUX_INT__g_s_f_s[96];

                    // ( F_102 S_000 | G_022 S_000 )^0_{t} = z * ( F_102 S_000 | F_021 S_000 )^0 + ( D_101 S_000 | F_021 S_000 )^0 + ( F_102 S_000 | D_020 S_000 )^0 - ( G_103 S_000 | F_021 S_000 )^0
                    AUX_INT__f_s_g_s[87] = etfac[2] * AUX_INT__f_s_f_s[57] + 2 * one_over_2q * AUX_INT__d_s_f_s[27] + 1 * one_over_2q * AUX_INT__f_s_d_s[33] - p_over_q * AUX_INT__g_s_f_s[97];

                    // ( F_102 S_000 | G_013 S_000 )^0_{t} = y * ( F_102 S_000 | F_003 S_000 )^0 - ( G_112 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[88] = etfac[1] * AUX_INT__f_s_f_s[59] - p_over_q * AUX_INT__g_s_f_s[89];

                    // ( F_102 S_000 | G_004 S_000 )^0_{t} = z * ( F_102 S_000 | F_003 S_000 )^0 + ( D_101 S_000 | F_003 S_000 )^0 + ( F_102 S_000 | D_002 S_000 )^0 - ( G_103 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[89] = etfac[2] * AUX_INT__f_s_f_s[59] + 2 * one_over_2q * AUX_INT__d_s_f_s[29] + 3 * one_over_2q * AUX_INT__f_s_d_s[35] - p_over_q * AUX_INT__g_s_f_s[99];

                    // ( F_030 S_000 | G_400 S_000 )^0_{t} = x * ( F_030 S_000 | F_300 S_000 )^0 + ( F_030 S_000 | D_200 S_000 )^0 - ( G_130 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[90] = etfac[0] * AUX_INT__f_s_f_s[60] + 3 * one_over_2q * AUX_INT__f_s_d_s[36] - p_over_q * AUX_INT__g_s_f_s[60];

                    // ( F_030 S_000 | G_310 S_000 )^0_{t} = y * ( F_030 S_000 | F_300 S_000 )^0 + ( D_020 S_000 | F_300 S_000 )^0 - ( G_040 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[91] = etfac[1] * AUX_INT__f_s_f_s[60] + 3 * one_over_2q * AUX_INT__d_s_f_s[30] - p_over_q * AUX_INT__g_s_f_s[100];

                    // ( F_030 S_000 | G_301 S_000 )^0_{t} = z * ( F_030 S_000 | F_300 S_000 )^0 - ( G_031 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[92] = etfac[2] * AUX_INT__f_s_f_s[60] - p_over_q * AUX_INT__g_s_f_s[110];

                    // ( F_030 S_000 | G_220 S_000 )^0_{t} = y * ( F_030 S_000 | F_210 S_000 )^0 + ( D_020 S_000 | F_210 S_000 )^0 + ( F_030 S_000 | D_200 S_000 )^0 - ( G_040 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[93] = etfac[1] * AUX_INT__f_s_f_s[61] + 3 * one_over_2q * AUX_INT__d_s_f_s[31] + 1 * one_over_2q * AUX_INT__f_s_d_s[36] - p_over_q * AUX_INT__g_s_f_s[101];

                    // ( F_030 S_000 | G_211 S_000 )^0_{t} = z * ( F_030 S_000 | F_210 S_000 )^0 - ( G_031 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[94] = etfac[2] * AUX_INT__f_s_f_s[61] - p_over_q * AUX_INT__g_s_f_s[111];

                    // ( F_030 S_000 | G_202 S_000 )^0_{t} = z * ( F_030 S_000 | F_201 S_000 )^0 + ( F_030 S_000 | D_200 S_000 )^0 - ( G_031 S_000 | F_201 S_000 )^0
                    AUX_INT__f_s_g_s[95] = etfac[2] * AUX_INT__f_s_f_s[62] + 1 * one_over_2q * AUX_INT__f_s_d_s[36] - p_over_q * AUX_INT__g_s_f_s[112];

                    // ( F_030 S_000 | G_130 S_000 )^0_{t} = x * ( F_030 S_000 | F_030 S_000 )^0 - ( G_130 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[96] = etfac[0] * AUX_INT__f_s_f_s[66] - p_over_q * AUX_INT__g_s_f_s[66];

                    // ( F_030 S_000 | G_121 S_000 )^0_{t} = z * ( F_030 S_000 | F_120 S_000 )^0 - ( G_031 S_000 | F_120 S_000 )^0
                    AUX_INT__f_s_g_s[97] = etfac[2] * AUX_INT__f_s_f_s[63] - p_over_q * AUX_INT__g_s_f_s[113];

                    // ( F_030 S_000 | G_112 S_000 )^0_{t} = y * ( F_030 S_000 | F_102 S_000 )^0 + ( D_020 S_000 | F_102 S_000 )^0 - ( G_040 S_000 | F_102 S_000 )^0
                    AUX_INT__f_s_g_s[98] = etfac[1] * AUX_INT__f_s_f_s[65] + 3 * one_over_2q * AUX_INT__d_s_f_s[35] - p_over_q * AUX_INT__g_s_f_s[105];

                    // ( F_030 S_000 | G_103 S_000 )^0_{t} = x * ( F_030 S_000 | F_003 S_000 )^0 - ( G_130 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[99] = etfac[0] * AUX_INT__f_s_f_s[69] - p_over_q * AUX_INT__g_s_f_s[69];

                    // ( F_030 S_000 | G_040 S_000 )^0_{t} = y * ( F_030 S_000 | F_030 S_000 )^0 + ( D_020 S_000 | F_030 S_000 )^0 + ( F_030 S_000 | D_020 S_000 )^0 - ( G_040 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[100] = etfac[1] * AUX_INT__f_s_f_s[66] + 3 * one_over_2q * AUX_INT__d_s_f_s[36] + 3 * one_over_2q * AUX_INT__f_s_d_s[39] - p_over_q * AUX_INT__g_s_f_s[106];

                    // ( F_030 S_000 | G_031 S_000 )^0_{t} = z * ( F_030 S_000 | F_030 S_000 )^0 - ( G_031 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[101] = etfac[2] * AUX_INT__f_s_f_s[66] - p_over_q * AUX_INT__g_s_f_s[116];

                    // ( F_030 S_000 | G_022 S_000 )^0_{t} = z * ( F_030 S_000 | F_021 S_000 )^0 + ( F_030 S_000 | D_020 S_000 )^0 - ( G_031 S_000 | F_021 S_000 )^0
                    AUX_INT__f_s_g_s[102] = etfac[2] * AUX_INT__f_s_f_s[67] + 1 * one_over_2q * AUX_INT__f_s_d_s[39] - p_over_q * AUX_INT__g_s_f_s[117];

                    // ( F_030 S_000 | G_013 S_000 )^0_{t} = y * ( F_030 S_000 | F_003 S_000 )^0 + ( D_020 S_000 | F_003 S_000 )^0 - ( G_040 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[103] = etfac[1] * AUX_INT__f_s_f_s[69] + 3 * one_over_2q * AUX_INT__d_s_f_s[39] - p_over_q * AUX_INT__g_s_f_s[109];

                    // ( F_030 S_000 | G_004 S_000 )^0_{t} = z * ( F_030 S_000 | F_003 S_000 )^0 + ( F_030 S_000 | D_002 S_000 )^0 - ( G_031 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[104] = etfac[2] * AUX_INT__f_s_f_s[69] + 3 * one_over_2q * AUX_INT__f_s_d_s[41] - p_over_q * AUX_INT__g_s_f_s[119];

                    // ( F_021 S_000 | G_400 S_000 )^0_{t} = x * ( F_021 S_000 | F_300 S_000 )^0 + ( F_021 S_000 | D_200 S_000 )^0 - ( G_121 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[105] = etfac[0] * AUX_INT__f_s_f_s[70] + 3 * one_over_2q * AUX_INT__f_s_d_s[42] - p_over_q * AUX_INT__g_s_f_s[70];

                    // ( F_021 S_000 | G_310 S_000 )^0_{t} = y * ( F_021 S_000 | F_300 S_000 )^0 + ( D_011 S_000 | F_300 S_000 )^0 - ( G_031 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[106] = etfac[1] * AUX_INT__f_s_f_s[70] + 2 * one_over_2q * AUX_INT__d_s_f_s[40] - p_over_q * AUX_INT__g_s_f_s[110];

                    // ( F_021 S_000 | G_301 S_000 )^0_{t} = z * ( F_021 S_000 | F_300 S_000 )^0 + ( D_020 S_000 | F_300 S_000 )^0 - ( G_022 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[107] = etfac[2] * AUX_INT__f_s_f_s[70] + 1 * one_over_2q * AUX_INT__d_s_f_s[30] - p_over_q * AUX_INT__g_s_f_s[120];

                    // ( F_021 S_000 | G_220 S_000 )^0_{t} = y * ( F_021 S_000 | F_210 S_000 )^0 + ( D_011 S_000 | F_210 S_000 )^0 + ( F_021 S_000 | D_200 S_000 )^0 - ( G_031 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[108] = etfac[1] * AUX_INT__f_s_f_s[71] + 2 * one_over_2q * AUX_INT__d_s_f_s[41] + 1 * one_over_2q * AUX_INT__f_s_d_s[42] - p_over_q * AUX_INT__g_s_f_s[111];

                    // ( F_021 S_000 | G_211 S_000 )^0_{t} = z * ( F_021 S_000 | F_210 S_000 )^0 + ( D_020 S_000 | F_210 S_000 )^0 - ( G_022 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[109] = etfac[2] * AUX_INT__f_s_f_s[71] + 1 * one_over_2q * AUX_INT__d_s_f_s[31] - p_over_q * AUX_INT__g_s_f_s[121];

                    // ( F_021 S_000 | G_202 S_000 )^0_{t} = z * ( F_021 S_000 | F_201 S_000 )^0 + ( D_020 S_000 | F_201 S_000 )^0 + ( F_021 S_000 | D_200 S_000 )^0 - ( G_022 S_000 | F_201 S_000 )^0
                    AUX_INT__f_s_g_s[110] = etfac[2] * AUX_INT__f_s_f_s[72] + 1 * one_over_2q * AUX_INT__d_s_f_s[32] + 1 * one_over_2q * AUX_INT__f_s_d_s[42] - p_over_q * AUX_INT__g_s_f_s[122];

                    // ( F_021 S_000 | G_130 S_000 )^0_{t} = x * ( F_021 S_000 | F_030 S_000 )^0 - ( G_121 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[111] = etfac[0] * AUX_INT__f_s_f_s[76] - p_over_q * AUX_INT__g_s_f_s[76];

                    // ( F_021 S_000 | G_121 S_000 )^0_{t} = z * ( F_021 S_000 | F_120 S_000 )^0 + ( D_020 S_000 | F_120 S_000 )^0 - ( G_022 S_000 | F_120 S_000 )^0
                    AUX_INT__f_s_g_s[112] = etfac[2] * AUX_INT__f_s_f_s[73] + 1 * one_over_2q * AUX_INT__d_s_f_s[33] - p_over_q * AUX_INT__g_s_f_s[123];

                    // ( F_021 S_000 | G_112 S_000 )^0_{t} = y * ( F_021 S_000 | F_102 S_000 )^0 + ( D_011 S_000 | F_102 S_000 )^0 - ( G_031 S_000 | F_102 S_000 )^0
                    AUX_INT__f_s_g_s[113] = etfac[1] * AUX_INT__f_s_f_s[75] + 2 * one_over_2q * AUX_INT__d_s_f_s[45] - p_over_q * AUX_INT__g_s_f_s[115];

                    // ( F_021 S_000 | G_103 S_000 )^0_{t} = x * ( F_021 S_000 | F_003 S_000 )^0 - ( G_121 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[114] = etfac[0] * AUX_INT__f_s_f_s[79] - p_over_q * AUX_INT__g_s_f_s[79];

                    // ( F_021 S_000 | G_040 S_000 )^0_{t} = y * ( F_021 S_000 | F_030 S_000 )^0 + ( D_011 S_000 | F_030 S_000 )^0 + ( F_021 S_000 | D_020 S_000 )^0 - ( G_031 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[115] = etfac[1] * AUX_INT__f_s_f_s[76] + 2 * one_over_2q * AUX_INT__d_s_f_s[46] + 3 * one_over_2q * AUX_INT__f_s_d_s[45] - p_over_q * AUX_INT__g_s_f_s[116];

                    // ( F_021 S_000 | G_031 S_000 )^0_{t} = z * ( F_021 S_000 | F_030 S_000 )^0 + ( D_020 S_000 | F_030 S_000 )^0 - ( G_022 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[116] = etfac[2] * AUX_INT__f_s_f_s[76] + 1 * one_over_2q * AUX_INT__d_s_f_s[36] - p_over_q * AUX_INT__g_s_f_s[126];

                    // ( F_021 S_000 | G_022 S_000 )^0_{t} = z * ( F_021 S_000 | F_021 S_000 )^0 + ( D_020 S_000 | F_021 S_000 )^0 + ( F_021 S_000 | D_020 S_000 )^0 - ( G_022 S_000 | F_021 S_000 )^0
                    AUX_INT__f_s_g_s[117] = etfac[2] * AUX_INT__f_s_f_s[77] + 1 * one_over_2q * AUX_INT__d_s_f_s[37] + 1 * one_over_2q * AUX_INT__f_s_d_s[45] - p_over_q * AUX_INT__g_s_f_s[127];

                    // ( F_021 S_000 | G_013 S_000 )^0_{t} = y * ( F_021 S_000 | F_003 S_000 )^0 + ( D_011 S_000 | F_003 S_000 )^0 - ( G_031 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[118] = etfac[1] * AUX_INT__f_s_f_s[79] + 2 * one_over_2q * AUX_INT__d_s_f_s[49] - p_over_q * AUX_INT__g_s_f_s[119];

                    // ( F_021 S_000 | G_004 S_000 )^0_{t} = z * ( F_021 S_000 | F_003 S_000 )^0 + ( D_020 S_000 | F_003 S_000 )^0 + ( F_021 S_000 | D_002 S_000 )^0 - ( G_022 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[119] = etfac[2] * AUX_INT__f_s_f_s[79] + 1 * one_over_2q * AUX_INT__d_s_f_s[39] + 3 * one_over_2q * AUX_INT__f_s_d_s[47] - p_over_q * AUX_INT__g_s_f_s[129];

                    // ( F_012 S_000 | G_400 S_000 )^0_{t} = x * ( F_012 S_000 | F_300 S_000 )^0 + ( F_012 S_000 | D_200 S_000 )^0 - ( G_112 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[120] = etfac[0] * AUX_INT__f_s_f_s[80] + 3 * one_over_2q * AUX_INT__f_s_d_s[48] - p_over_q * AUX_INT__g_s_f_s[80];

                    // ( F_012 S_000 | G_310 S_000 )^0_{t} = y * ( F_012 S_000 | F_300 S_000 )^0 + ( D_002 S_000 | F_300 S_000 )^0 - ( G_022 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[121] = etfac[1] * AUX_INT__f_s_f_s[80] + 1 * one_over_2q * AUX_INT__d_s_f_s[50] - p_over_q * AUX_INT__g_s_f_s[120];

                    // ( F_012 S_000 | G_301 S_000 )^0_{t} = z * ( F_012 S_000 | F_300 S_000 )^0 + ( D_011 S_000 | F_300 S_000 )^0 - ( G_013 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[122] = etfac[2] * AUX_INT__f_s_f_s[80] + 2 * one_over_2q * AUX_INT__d_s_f_s[40] - p_over_q * AUX_INT__g_s_f_s[130];

                    // ( F_012 S_000 | G_220 S_000 )^0_{t} = y * ( F_012 S_000 | F_210 S_000 )^0 + ( D_002 S_000 | F_210 S_000 )^0 + ( F_012 S_000 | D_200 S_000 )^0 - ( G_022 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[123] = etfac[1] * AUX_INT__f_s_f_s[81] + 1 * one_over_2q * AUX_INT__d_s_f_s[51] + 1 * one_over_2q * AUX_INT__f_s_d_s[48] - p_over_q * AUX_INT__g_s_f_s[121];

                    // ( F_012 S_000 | G_211 S_000 )^0_{t} = z * ( F_012 S_000 | F_210 S_000 )^0 + ( D_011 S_000 | F_210 S_000 )^0 - ( G_013 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[124] = etfac[2] * AUX_INT__f_s_f_s[81] + 2 * one_over_2q * AUX_INT__d_s_f_s[41] - p_over_q * AUX_INT__g_s_f_s[131];

                    // ( F_012 S_000 | G_202 S_000 )^0_{t} = z * ( F_012 S_000 | F_201 S_000 )^0 + ( D_011 S_000 | F_201 S_000 )^0 + ( F_012 S_000 | D_200 S_000 )^0 - ( G_013 S_000 | F_201 S_000 )^0
                    AUX_INT__f_s_g_s[125] = etfac[2] * AUX_INT__f_s_f_s[82] + 2 * one_over_2q * AUX_INT__d_s_f_s[42] + 1 * one_over_2q * AUX_INT__f_s_d_s[48] - p_over_q * AUX_INT__g_s_f_s[132];

                    // ( F_012 S_000 | G_130 S_000 )^0_{t} = x * ( F_012 S_000 | F_030 S_000 )^0 - ( G_112 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[126] = etfac[0] * AUX_INT__f_s_f_s[86] - p_over_q * AUX_INT__g_s_f_s[86];

                    // ( F_012 S_000 | G_121 S_000 )^0_{t} = z * ( F_012 S_000 | F_120 S_000 )^0 + ( D_011 S_000 | F_120 S_000 )^0 - ( G_013 S_000 | F_120 S_000 )^0
                    AUX_INT__f_s_g_s[127] = etfac[2] * AUX_INT__f_s_f_s[83] + 2 * one_over_2q * AUX_INT__d_s_f_s[43] - p_over_q * AUX_INT__g_s_f_s[133];

                    // ( F_012 S_000 | G_112 S_000 )^0_{t} = y * ( F_012 S_000 | F_102 S_000 )^0 + ( D_002 S_000 | F_102 S_000 )^0 - ( G_022 S_000 | F_102 S_000 )^0
                    AUX_INT__f_s_g_s[128] = etfac[1] * AUX_INT__f_s_f_s[85] + 1 * one_over_2q * AUX_INT__d_s_f_s[55] - p_over_q * AUX_INT__g_s_f_s[125];

                    // ( F_012 S_000 | G_103 S_000 )^0_{t} = x * ( F_012 S_000 | F_003 S_000 )^0 - ( G_112 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[129] = etfac[0] * AUX_INT__f_s_f_s[89] - p_over_q * AUX_INT__g_s_f_s[89];

                    // ( F_012 S_000 | G_040 S_000 )^0_{t} = y * ( F_012 S_000 | F_030 S_000 )^0 + ( D_002 S_000 | F_030 S_000 )^0 + ( F_012 S_000 | D_020 S_000 )^0 - ( G_022 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[130] = etfac[1] * AUX_INT__f_s_f_s[86] + 1 * one_over_2q * AUX_INT__d_s_f_s[56] + 3 * one_over_2q * AUX_INT__f_s_d_s[51] - p_over_q * AUX_INT__g_s_f_s[126];

                    // ( F_012 S_000 | G_031 S_000 )^0_{t} = z * ( F_012 S_000 | F_030 S_000 )^0 + ( D_011 S_000 | F_030 S_000 )^0 - ( G_013 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[131] = etfac[2] * AUX_INT__f_s_f_s[86] + 2 * one_over_2q * AUX_INT__d_s_f_s[46] - p_over_q * AUX_INT__g_s_f_s[136];

                    // ( F_012 S_000 | G_022 S_000 )^0_{t} = z * ( F_012 S_000 | F_021 S_000 )^0 + ( D_011 S_000 | F_021 S_000 )^0 + ( F_012 S_000 | D_020 S_000 )^0 - ( G_013 S_000 | F_021 S_000 )^0
                    AUX_INT__f_s_g_s[132] = etfac[2] * AUX_INT__f_s_f_s[87] + 2 * one_over_2q * AUX_INT__d_s_f_s[47] + 1 * one_over_2q * AUX_INT__f_s_d_s[51] - p_over_q * AUX_INT__g_s_f_s[137];

                    // ( F_012 S_000 | G_013 S_000 )^0_{t} = y * ( F_012 S_000 | F_003 S_000 )^0 + ( D_002 S_000 | F_003 S_000 )^0 - ( G_022 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[133] = etfac[1] * AUX_INT__f_s_f_s[89] + 1 * one_over_2q * AUX_INT__d_s_f_s[59] - p_over_q * AUX_INT__g_s_f_s[129];

                    // ( F_012 S_000 | G_004 S_000 )^0_{t} = z * ( F_012 S_000 | F_003 S_000 )^0 + ( D_011 S_000 | F_003 S_000 )^0 + ( F_012 S_000 | D_002 S_000 )^0 - ( G_013 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[134] = etfac[2] * AUX_INT__f_s_f_s[89] + 2 * one_over_2q * AUX_INT__d_s_f_s[49] + 3 * one_over_2q * AUX_INT__f_s_d_s[53] - p_over_q * AUX_INT__g_s_f_s[139];

                    // ( F_003 S_000 | G_400 S_000 )^0_{t} = x * ( F_003 S_000 | F_300 S_000 )^0 + ( F_003 S_000 | D_200 S_000 )^0 - ( G_103 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[135] = etfac[0] * AUX_INT__f_s_f_s[90] + 3 * one_over_2q * AUX_INT__f_s_d_s[54] - p_over_q * AUX_INT__g_s_f_s[90];

                    // ( F_003 S_000 | G_310 S_000 )^0_{t} = y * ( F_003 S_000 | F_300 S_000 )^0 - ( G_013 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[136] = etfac[1] * AUX_INT__f_s_f_s[90] - p_over_q * AUX_INT__g_s_f_s[130];

                    // ( F_003 S_000 | G_301 S_000 )^0_{t} = z * ( F_003 S_000 | F_300 S_000 )^0 + ( D_002 S_000 | F_300 S_000 )^0 - ( G_004 S_000 | F_300 S_000 )^0
                    AUX_INT__f_s_g_s[137] = etfac[2] * AUX_INT__f_s_f_s[90] + 3 * one_over_2q * AUX_INT__d_s_f_s[50] - p_over_q * AUX_INT__g_s_f_s[140];

                    // ( F_003 S_000 | G_220 S_000 )^0_{t} = y * ( F_003 S_000 | F_210 S_000 )^0 + ( F_003 S_000 | D_200 S_000 )^0 - ( G_013 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[138] = etfac[1] * AUX_INT__f_s_f_s[91] + 1 * one_over_2q * AUX_INT__f_s_d_s[54] - p_over_q * AUX_INT__g_s_f_s[131];

                    // ( F_003 S_000 | G_211 S_000 )^0_{t} = z * ( F_003 S_000 | F_210 S_000 )^0 + ( D_002 S_000 | F_210 S_000 )^0 - ( G_004 S_000 | F_210 S_000 )^0
                    AUX_INT__f_s_g_s[139] = etfac[2] * AUX_INT__f_s_f_s[91] + 3 * one_over_2q * AUX_INT__d_s_f_s[51] - p_over_q * AUX_INT__g_s_f_s[141];

                    // ( F_003 S_000 | G_202 S_000 )^0_{t} = z * ( F_003 S_000 | F_201 S_000 )^0 + ( D_002 S_000 | F_201 S_000 )^0 + ( F_003 S_000 | D_200 S_000 )^0 - ( G_004 S_000 | F_201 S_000 )^0
                    AUX_INT__f_s_g_s[140] = etfac[2] * AUX_INT__f_s_f_s[92] + 3 * one_over_2q * AUX_INT__d_s_f_s[52] + 1 * one_over_2q * AUX_INT__f_s_d_s[54] - p_over_q * AUX_INT__g_s_f_s[142];

                    // ( F_003 S_000 | G_130 S_000 )^0_{t} = x * ( F_003 S_000 | F_030 S_000 )^0 - ( G_103 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[141] = etfac[0] * AUX_INT__f_s_f_s[96] - p_over_q * AUX_INT__g_s_f_s[96];

                    // ( F_003 S_000 | G_121 S_000 )^0_{t} = z * ( F_003 S_000 | F_120 S_000 )^0 + ( D_002 S_000 | F_120 S_000 )^0 - ( G_004 S_000 | F_120 S_000 )^0
                    AUX_INT__f_s_g_s[142] = etfac[2] * AUX_INT__f_s_f_s[93] + 3 * one_over_2q * AUX_INT__d_s_f_s[53] - p_over_q * AUX_INT__g_s_f_s[143];

                    // ( F_003 S_000 | G_112 S_000 )^0_{t} = y * ( F_003 S_000 | F_102 S_000 )^0 - ( G_013 S_000 | F_102 S_000 )^0
                    AUX_INT__f_s_g_s[143] = etfac[1] * AUX_INT__f_s_f_s[95] - p_over_q * AUX_INT__g_s_f_s[135];

                    // ( F_003 S_000 | G_103 S_000 )^0_{t} = x * ( F_003 S_000 | F_003 S_000 )^0 - ( G_103 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[144] = etfac[0] * AUX_INT__f_s_f_s[99] - p_over_q * AUX_INT__g_s_f_s[99];

                    // ( F_003 S_000 | G_040 S_000 )^0_{t} = y * ( F_003 S_000 | F_030 S_000 )^0 + ( F_003 S_000 | D_020 S_000 )^0 - ( G_013 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[145] = etfac[1] * AUX_INT__f_s_f_s[96] + 3 * one_over_2q * AUX_INT__f_s_d_s[57] - p_over_q * AUX_INT__g_s_f_s[136];

                    // ( F_003 S_000 | G_031 S_000 )^0_{t} = z * ( F_003 S_000 | F_030 S_000 )^0 + ( D_002 S_000 | F_030 S_000 )^0 - ( G_004 S_000 | F_030 S_000 )^0
                    AUX_INT__f_s_g_s[146] = etfac[2] * AUX_INT__f_s_f_s[96] + 3 * one_over_2q * AUX_INT__d_s_f_s[56] - p_over_q * AUX_INT__g_s_f_s[146];

                    // ( F_003 S_000 | G_022 S_000 )^0_{t} = z * ( F_003 S_000 | F_021 S_000 )^0 + ( D_002 S_000 | F_021 S_000 )^0 + ( F_003 S_000 | D_020 S_000 )^0 - ( G_004 S_000 | F_021 S_000 )^0
                    AUX_INT__f_s_g_s[147] = etfac[2] * AUX_INT__f_s_f_s[97] + 3 * one_over_2q * AUX_INT__d_s_f_s[57] + 1 * one_over_2q * AUX_INT__f_s_d_s[57] - p_over_q * AUX_INT__g_s_f_s[147];

                    // ( F_003 S_000 | G_013 S_000 )^0_{t} = y * ( F_003 S_000 | F_003 S_000 )^0 - ( G_013 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[148] = etfac[1] * AUX_INT__f_s_f_s[99] - p_over_q * AUX_INT__g_s_f_s[139];

                    // ( F_003 S_000 | G_004 S_000 )^0_{t} = z * ( F_003 S_000 | F_003 S_000 )^0 + ( D_002 S_000 | F_003 S_000 )^0 + ( F_003 S_000 | D_002 S_000 )^0 - ( G_004 S_000 | F_003 S_000 )^0
                    AUX_INT__f_s_g_s[149] = etfac[2] * AUX_INT__f_s_f_s[99] + 3 * one_over_2q * AUX_INT__d_s_f_s[59] + 3 * one_over_2q * AUX_INT__f_s_d_s[59] - p_over_q * AUX_INT__g_s_f_s[149];

                    // ( J_700 S_000 | P_100 S_000 )^0 = x * ( J_700 S_000 | S_000 S_000 )^0_{e} + ( I_600 S_000 | S_000 S_000 )^0_{e} - ( K_800 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[0] = etfac[0] * AUX_INT__j_s_s_s[0] + 7 * one_over_2q * AUX_INT__i_s_s_s[0] - p_over_q * AUX_INT__k_s_s_s[0];

                    // ( J_610 S_000 | P_100 S_000 )^0 = x * ( J_610 S_000 | S_000 S_000 )^0_{e} + ( I_510 S_000 | S_000 S_000 )^0_{e} - ( K_710 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[3] = etfac[0] * AUX_INT__j_s_s_s[1] + 6 * one_over_2q * AUX_INT__i_s_s_s[1] - p_over_q * AUX_INT__k_s_s_s[1];

                    // ( J_601 S_000 | P_100 S_000 )^0 = x * ( J_601 S_000 | S_000 S_000 )^0_{e} + ( I_501 S_000 | S_000 S_000 )^0_{e} - ( K_701 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[6] = etfac[0] * AUX_INT__j_s_s_s[2] + 6 * one_over_2q * AUX_INT__i_s_s_s[2] - p_over_q * AUX_INT__k_s_s_s[2];

                    // ( J_520 S_000 | P_100 S_000 )^0 = x * ( J_520 S_000 | S_000 S_000 )^0_{e} + ( I_420 S_000 | S_000 S_000 )^0_{e} - ( K_620 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[9] = etfac[0] * AUX_INT__j_s_s_s[3] + 5 * one_over_2q * AUX_INT__i_s_s_s[3] - p_over_q * AUX_INT__k_s_s_s[3];

                    // ( J_520 S_000 | P_010 S_000 )^0 = y * ( J_520 S_000 | S_000 S_000 )^0_{e} + ( I_510 S_000 | S_000 S_000 )^0_{e} - ( K_530 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[10] = etfac[1] * AUX_INT__j_s_s_s[3] + 2 * one_over_2q * AUX_INT__i_s_s_s[1] - p_over_q * AUX_INT__k_s_s_s[6];

                    // ( J_511 S_000 | P_100 S_000 )^0 = x * ( J_511 S_000 | S_000 S_000 )^0_{e} + ( I_411 S_000 | S_000 S_000 )^0_{e} - ( K_611 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[12] = etfac[0] * AUX_INT__j_s_s_s[4] + 5 * one_over_2q * AUX_INT__i_s_s_s[4] - p_over_q * AUX_INT__k_s_s_s[4];

                    // ( J_511 S_000 | P_010 S_000 )^0 = y * ( J_511 S_000 | S_000 S_000 )^0_{e} + ( I_501 S_000 | S_000 S_000 )^0_{e} - ( K_521 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[13] = etfac[1] * AUX_INT__j_s_s_s[4] + 1 * one_over_2q * AUX_INT__i_s_s_s[2] - p_over_q * AUX_INT__k_s_s_s[7];

                    // ( J_511 S_000 | P_001 S_000 )^0 = z * ( J_511 S_000 | S_000 S_000 )^0_{e} + ( I_510 S_000 | S_000 S_000 )^0_{e} - ( K_512 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[14] = etfac[2] * AUX_INT__j_s_s_s[4] + 1 * one_over_2q * AUX_INT__i_s_s_s[1] - p_over_q * AUX_INT__k_s_s_s[8];

                    // ( J_502 S_000 | P_100 S_000 )^0 = x * ( J_502 S_000 | S_000 S_000 )^0_{e} + ( I_402 S_000 | S_000 S_000 )^0_{e} - ( K_602 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[15] = etfac[0] * AUX_INT__j_s_s_s[5] + 5 * one_over_2q * AUX_INT__i_s_s_s[5] - p_over_q * AUX_INT__k_s_s_s[5];

                    // ( J_502 S_000 | P_001 S_000 )^0 = z * ( J_502 S_000 | S_000 S_000 )^0_{e} + ( I_501 S_000 | S_000 S_000 )^0_{e} - ( K_503 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[17] = etfac[2] * AUX_INT__j_s_s_s[5] + 2 * one_over_2q * AUX_INT__i_s_s_s[2] - p_over_q * AUX_INT__k_s_s_s[9];

                    // ( J_430 S_000 | P_100 S_000 )^0 = x * ( J_430 S_000 | S_000 S_000 )^0_{e} + ( I_330 S_000 | S_000 S_000 )^0_{e} - ( K_530 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[18] = etfac[0] * AUX_INT__j_s_s_s[6] + 4 * one_over_2q * AUX_INT__i_s_s_s[6] - p_over_q * AUX_INT__k_s_s_s[6];

                    // ( J_430 S_000 | P_010 S_000 )^0 = y * ( J_430 S_000 | S_000 S_000 )^0_{e} + ( I_420 S_000 | S_000 S_000 )^0_{e} - ( K_440 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[19] = etfac[1] * AUX_INT__j_s_s_s[6] + 3 * one_over_2q * AUX_INT__i_s_s_s[3] - p_over_q * AUX_INT__k_s_s_s[10];

                    // ( J_421 S_000 | P_100 S_000 )^0 = x * ( J_421 S_000 | S_000 S_000 )^0_{e} + ( I_321 S_000 | S_000 S_000 )^0_{e} - ( K_521 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[21] = etfac[0] * AUX_INT__j_s_s_s[7] + 4 * one_over_2q * AUX_INT__i_s_s_s[7] - p_over_q * AUX_INT__k_s_s_s[7];

                    // ( J_421 S_000 | P_010 S_000 )^0 = y * ( J_421 S_000 | S_000 S_000 )^0_{e} + ( I_411 S_000 | S_000 S_000 )^0_{e} - ( K_431 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[22] = etfac[1] * AUX_INT__j_s_s_s[7] + 2 * one_over_2q * AUX_INT__i_s_s_s[4] - p_over_q * AUX_INT__k_s_s_s[11];

                    // ( J_421 S_000 | P_001 S_000 )^0 = z * ( J_421 S_000 | S_000 S_000 )^0_{e} + ( I_420 S_000 | S_000 S_000 )^0_{e} - ( K_422 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[23] = etfac[2] * AUX_INT__j_s_s_s[7] + 1 * one_over_2q * AUX_INT__i_s_s_s[3] - p_over_q * AUX_INT__k_s_s_s[12];

                    // ( J_412 S_000 | P_100 S_000 )^0 = x * ( J_412 S_000 | S_000 S_000 )^0_{e} + ( I_312 S_000 | S_000 S_000 )^0_{e} - ( K_512 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[24] = etfac[0] * AUX_INT__j_s_s_s[8] + 4 * one_over_2q * AUX_INT__i_s_s_s[8] - p_over_q * AUX_INT__k_s_s_s[8];

                    // ( J_412 S_000 | P_010 S_000 )^0 = y * ( J_412 S_000 | S_000 S_000 )^0_{e} + ( I_402 S_000 | S_000 S_000 )^0_{e} - ( K_422 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[25] = etfac[1] * AUX_INT__j_s_s_s[8] + 1 * one_over_2q * AUX_INT__i_s_s_s[5] - p_over_q * AUX_INT__k_s_s_s[12];

                    // ( J_412 S_000 | P_001 S_000 )^0 = z * ( J_412 S_000 | S_000 S_000 )^0_{e} + ( I_411 S_000 | S_000 S_000 )^0_{e} - ( K_413 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[26] = etfac[2] * AUX_INT__j_s_s_s[8] + 2 * one_over_2q * AUX_INT__i_s_s_s[4] - p_over_q * AUX_INT__k_s_s_s[13];

                    // ( J_403 S_000 | P_100 S_000 )^0 = x * ( J_403 S_000 | S_000 S_000 )^0_{e} + ( I_303 S_000 | S_000 S_000 )^0_{e} - ( K_503 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[27] = etfac[0] * AUX_INT__j_s_s_s[9] + 4 * one_over_2q * AUX_INT__i_s_s_s[9] - p_over_q * AUX_INT__k_s_s_s[9];

                    // ( J_403 S_000 | P_001 S_000 )^0 = z * ( J_403 S_000 | S_000 S_000 )^0_{e} + ( I_402 S_000 | S_000 S_000 )^0_{e} - ( K_404 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[29] = etfac[2] * AUX_INT__j_s_s_s[9] + 3 * one_over_2q * AUX_INT__i_s_s_s[5] - p_over_q * AUX_INT__k_s_s_s[14];

                    // ( J_340 S_000 | P_100 S_000 )^0 = x * ( J_340 S_000 | S_000 S_000 )^0_{e} + ( I_240 S_000 | S_000 S_000 )^0_{e} - ( K_440 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[30] = etfac[0] * AUX_INT__j_s_s_s[10] + 3 * one_over_2q * AUX_INT__i_s_s_s[10] - p_over_q * AUX_INT__k_s_s_s[10];

                    // ( J_340 S_000 | P_010 S_000 )^0 = y * ( J_340 S_000 | S_000 S_000 )^0_{e} + ( I_330 S_000 | S_000 S_000 )^0_{e} - ( K_350 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[31] = etfac[1] * AUX_INT__j_s_s_s[10] + 4 * one_over_2q * AUX_INT__i_s_s_s[6] - p_over_q * AUX_INT__k_s_s_s[15];

                    // ( J_331 S_000 | P_100 S_000 )^0 = x * ( J_331 S_000 | S_000 S_000 )^0_{e} + ( I_231 S_000 | S_000 S_000 )^0_{e} - ( K_431 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[33] = etfac[0] * AUX_INT__j_s_s_s[11] + 3 * one_over_2q * AUX_INT__i_s_s_s[11] - p_over_q * AUX_INT__k_s_s_s[11];

                    // ( J_331 S_000 | P_010 S_000 )^0 = y * ( J_331 S_000 | S_000 S_000 )^0_{e} + ( I_321 S_000 | S_000 S_000 )^0_{e} - ( K_341 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[34] = etfac[1] * AUX_INT__j_s_s_s[11] + 3 * one_over_2q * AUX_INT__i_s_s_s[7] - p_over_q * AUX_INT__k_s_s_s[16];

                    // ( J_331 S_000 | P_001 S_000 )^0 = z * ( J_331 S_000 | S_000 S_000 )^0_{e} + ( I_330 S_000 | S_000 S_000 )^0_{e} - ( K_332 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[35] = etfac[2] * AUX_INT__j_s_s_s[11] + 1 * one_over_2q * AUX_INT__i_s_s_s[6] - p_over_q * AUX_INT__k_s_s_s[17];

                    // ( J_322 S_000 | P_100 S_000 )^0 = x * ( J_322 S_000 | S_000 S_000 )^0_{e} + ( I_222 S_000 | S_000 S_000 )^0_{e} - ( K_422 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[36] = etfac[0] * AUX_INT__j_s_s_s[12] + 3 * one_over_2q * AUX_INT__i_s_s_s[12] - p_over_q * AUX_INT__k_s_s_s[12];

                    // ( J_322 S_000 | P_010 S_000 )^0 = y * ( J_322 S_000 | S_000 S_000 )^0_{e} + ( I_312 S_000 | S_000 S_000 )^0_{e} - ( K_332 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[37] = etfac[1] * AUX_INT__j_s_s_s[12] + 2 * one_over_2q * AUX_INT__i_s_s_s[8] - p_over_q * AUX_INT__k_s_s_s[17];

                    // ( J_322 S_000 | P_001 S_000 )^0 = z * ( J_322 S_000 | S_000 S_000 )^0_{e} + ( I_321 S_000 | S_000 S_000 )^0_{e} - ( K_323 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[38] = etfac[2] * AUX_INT__j_s_s_s[12] + 2 * one_over_2q * AUX_INT__i_s_s_s[7] - p_over_q * AUX_INT__k_s_s_s[18];

                    // ( J_313 S_000 | P_100 S_000 )^0 = x * ( J_313 S_000 | S_000 S_000 )^0_{e} + ( I_213 S_000 | S_000 S_000 )^0_{e} - ( K_413 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[39] = etfac[0] * AUX_INT__j_s_s_s[13] + 3 * one_over_2q * AUX_INT__i_s_s_s[13] - p_over_q * AUX_INT__k_s_s_s[13];

                    // ( J_313 S_000 | P_010 S_000 )^0 = y * ( J_313 S_000 | S_000 S_000 )^0_{e} + ( I_303 S_000 | S_000 S_000 )^0_{e} - ( K_323 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[40] = etfac[1] * AUX_INT__j_s_s_s[13] + 1 * one_over_2q * AUX_INT__i_s_s_s[9] - p_over_q * AUX_INT__k_s_s_s[18];

                    // ( J_313 S_000 | P_001 S_000 )^0 = z * ( J_313 S_000 | S_000 S_000 )^0_{e} + ( I_312 S_000 | S_000 S_000 )^0_{e} - ( K_314 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[41] = etfac[2] * AUX_INT__j_s_s_s[13] + 3 * one_over_2q * AUX_INT__i_s_s_s[8] - p_over_q * AUX_INT__k_s_s_s[19];

                    // ( J_304 S_000 | P_100 S_000 )^0 = x * ( J_304 S_000 | S_000 S_000 )^0_{e} + ( I_204 S_000 | S_000 S_000 )^0_{e} - ( K_404 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[42] = etfac[0] * AUX_INT__j_s_s_s[14] + 3 * one_over_2q * AUX_INT__i_s_s_s[14] - p_over_q * AUX_INT__k_s_s_s[14];

                    // ( J_304 S_000 | P_001 S_000 )^0 = z * ( J_304 S_000 | S_000 S_000 )^0_{e} + ( I_303 S_000 | S_000 S_000 )^0_{e} - ( K_305 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[44] = etfac[2] * AUX_INT__j_s_s_s[14] + 4 * one_over_2q * AUX_INT__i_s_s_s[9] - p_over_q * AUX_INT__k_s_s_s[20];

                    // ( J_250 S_000 | P_100 S_000 )^0 = x * ( J_250 S_000 | S_000 S_000 )^0_{e} + ( I_150 S_000 | S_000 S_000 )^0_{e} - ( K_350 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[45] = etfac[0] * AUX_INT__j_s_s_s[15] + 2 * one_over_2q * AUX_INT__i_s_s_s[15] - p_over_q * AUX_INT__k_s_s_s[15];

                    // ( J_250 S_000 | P_010 S_000 )^0 = y * ( J_250 S_000 | S_000 S_000 )^0_{e} + ( I_240 S_000 | S_000 S_000 )^0_{e} - ( K_260 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[46] = etfac[1] * AUX_INT__j_s_s_s[15] + 5 * one_over_2q * AUX_INT__i_s_s_s[10] - p_over_q * AUX_INT__k_s_s_s[21];

                    // ( J_241 S_000 | P_100 S_000 )^0 = x * ( J_241 S_000 | S_000 S_000 )^0_{e} + ( I_141 S_000 | S_000 S_000 )^0_{e} - ( K_341 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[48] = etfac[0] * AUX_INT__j_s_s_s[16] + 2 * one_over_2q * AUX_INT__i_s_s_s[16] - p_over_q * AUX_INT__k_s_s_s[16];

                    // ( J_241 S_000 | P_010 S_000 )^0 = y * ( J_241 S_000 | S_000 S_000 )^0_{e} + ( I_231 S_000 | S_000 S_000 )^0_{e} - ( K_251 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[49] = etfac[1] * AUX_INT__j_s_s_s[16] + 4 * one_over_2q * AUX_INT__i_s_s_s[11] - p_over_q * AUX_INT__k_s_s_s[22];

                    // ( J_241 S_000 | P_001 S_000 )^0 = z * ( J_241 S_000 | S_000 S_000 )^0_{e} + ( I_240 S_000 | S_000 S_000 )^0_{e} - ( K_242 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[50] = etfac[2] * AUX_INT__j_s_s_s[16] + 1 * one_over_2q * AUX_INT__i_s_s_s[10] - p_over_q * AUX_INT__k_s_s_s[23];

                    // ( J_232 S_000 | P_100 S_000 )^0 = x * ( J_232 S_000 | S_000 S_000 )^0_{e} + ( I_132 S_000 | S_000 S_000 )^0_{e} - ( K_332 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[51] = etfac[0] * AUX_INT__j_s_s_s[17] + 2 * one_over_2q * AUX_INT__i_s_s_s[17] - p_over_q * AUX_INT__k_s_s_s[17];

                    // ( J_232 S_000 | P_010 S_000 )^0 = y * ( J_232 S_000 | S_000 S_000 )^0_{e} + ( I_222 S_000 | S_000 S_000 )^0_{e} - ( K_242 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[52] = etfac[1] * AUX_INT__j_s_s_s[17] + 3 * one_over_2q * AUX_INT__i_s_s_s[12] - p_over_q * AUX_INT__k_s_s_s[23];

                    // ( J_232 S_000 | P_001 S_000 )^0 = z * ( J_232 S_000 | S_000 S_000 )^0_{e} + ( I_231 S_000 | S_000 S_000 )^0_{e} - ( K_233 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[53] = etfac[2] * AUX_INT__j_s_s_s[17] + 2 * one_over_2q * AUX_INT__i_s_s_s[11] - p_over_q * AUX_INT__k_s_s_s[24];

                    // ( J_223 S_000 | P_100 S_000 )^0 = x * ( J_223 S_000 | S_000 S_000 )^0_{e} + ( I_123 S_000 | S_000 S_000 )^0_{e} - ( K_323 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[54] = etfac[0] * AUX_INT__j_s_s_s[18] + 2 * one_over_2q * AUX_INT__i_s_s_s[18] - p_over_q * AUX_INT__k_s_s_s[18];

                    // ( J_223 S_000 | P_010 S_000 )^0 = y * ( J_223 S_000 | S_000 S_000 )^0_{e} + ( I_213 S_000 | S_000 S_000 )^0_{e} - ( K_233 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[55] = etfac[1] * AUX_INT__j_s_s_s[18] + 2 * one_over_2q * AUX_INT__i_s_s_s[13] - p_over_q * AUX_INT__k_s_s_s[24];

                    // ( J_223 S_000 | P_001 S_000 )^0 = z * ( J_223 S_000 | S_000 S_000 )^0_{e} + ( I_222 S_000 | S_000 S_000 )^0_{e} - ( K_224 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[56] = etfac[2] * AUX_INT__j_s_s_s[18] + 3 * one_over_2q * AUX_INT__i_s_s_s[12] - p_over_q * AUX_INT__k_s_s_s[25];

                    // ( J_214 S_000 | P_100 S_000 )^0 = x * ( J_214 S_000 | S_000 S_000 )^0_{e} + ( I_114 S_000 | S_000 S_000 )^0_{e} - ( K_314 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[57] = etfac[0] * AUX_INT__j_s_s_s[19] + 2 * one_over_2q * AUX_INT__i_s_s_s[19] - p_over_q * AUX_INT__k_s_s_s[19];

                    // ( J_214 S_000 | P_010 S_000 )^0 = y * ( J_214 S_000 | S_000 S_000 )^0_{e} + ( I_204 S_000 | S_000 S_000 )^0_{e} - ( K_224 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[58] = etfac[1] * AUX_INT__j_s_s_s[19] + 1 * one_over_2q * AUX_INT__i_s_s_s[14] - p_over_q * AUX_INT__k_s_s_s[25];

                    // ( J_214 S_000 | P_001 S_000 )^0 = z * ( J_214 S_000 | S_000 S_000 )^0_{e} + ( I_213 S_000 | S_000 S_000 )^0_{e} - ( K_215 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[59] = etfac[2] * AUX_INT__j_s_s_s[19] + 4 * one_over_2q * AUX_INT__i_s_s_s[13] - p_over_q * AUX_INT__k_s_s_s[26];

                    // ( J_205 S_000 | P_100 S_000 )^0 = x * ( J_205 S_000 | S_000 S_000 )^0_{e} + ( I_105 S_000 | S_000 S_000 )^0_{e} - ( K_305 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[60] = etfac[0] * AUX_INT__j_s_s_s[20] + 2 * one_over_2q * AUX_INT__i_s_s_s[20] - p_over_q * AUX_INT__k_s_s_s[20];

                    // ( J_205 S_000 | P_001 S_000 )^0 = z * ( J_205 S_000 | S_000 S_000 )^0_{e} + ( I_204 S_000 | S_000 S_000 )^0_{e} - ( K_206 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[62] = etfac[2] * AUX_INT__j_s_s_s[20] + 5 * one_over_2q * AUX_INT__i_s_s_s[14] - p_over_q * AUX_INT__k_s_s_s[27];

                    // ( J_160 S_000 | P_100 S_000 )^0 = x * ( J_160 S_000 | S_000 S_000 )^0_{e} + ( I_060 S_000 | S_000 S_000 )^0_{e} - ( K_260 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[63] = etfac[0] * AUX_INT__j_s_s_s[21] + 1 * one_over_2q * AUX_INT__i_s_s_s[21] - p_over_q * AUX_INT__k_s_s_s[21];

                    // ( J_160 S_000 | P_010 S_000 )^0 = y * ( J_160 S_000 | S_000 S_000 )^0_{e} + ( I_150 S_000 | S_000 S_000 )^0_{e} - ( K_170 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[64] = etfac[1] * AUX_INT__j_s_s_s[21] + 6 * one_over_2q * AUX_INT__i_s_s_s[15] - p_over_q * AUX_INT__k_s_s_s[28];

                    // ( J_151 S_000 | P_100 S_000 )^0 = x * ( J_151 S_000 | S_000 S_000 )^0_{e} + ( I_051 S_000 | S_000 S_000 )^0_{e} - ( K_251 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[66] = etfac[0] * AUX_INT__j_s_s_s[22] + 1 * one_over_2q * AUX_INT__i_s_s_s[22] - p_over_q * AUX_INT__k_s_s_s[22];

                    // ( J_151 S_000 | P_010 S_000 )^0 = y * ( J_151 S_000 | S_000 S_000 )^0_{e} + ( I_141 S_000 | S_000 S_000 )^0_{e} - ( K_161 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[67] = etfac[1] * AUX_INT__j_s_s_s[22] + 5 * one_over_2q * AUX_INT__i_s_s_s[16] - p_over_q * AUX_INT__k_s_s_s[29];

                    // ( J_151 S_000 | P_001 S_000 )^0 = z * ( J_151 S_000 | S_000 S_000 )^0_{e} + ( I_150 S_000 | S_000 S_000 )^0_{e} - ( K_152 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[68] = etfac[2] * AUX_INT__j_s_s_s[22] + 1 * one_over_2q * AUX_INT__i_s_s_s[15] - p_over_q * AUX_INT__k_s_s_s[30];

                    // ( J_142 S_000 | P_100 S_000 )^0 = x * ( J_142 S_000 | S_000 S_000 )^0_{e} + ( I_042 S_000 | S_000 S_000 )^0_{e} - ( K_242 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[69] = etfac[0] * AUX_INT__j_s_s_s[23] + 1 * one_over_2q * AUX_INT__i_s_s_s[23] - p_over_q * AUX_INT__k_s_s_s[23];

                    // ( J_142 S_000 | P_010 S_000 )^0 = y * ( J_142 S_000 | S_000 S_000 )^0_{e} + ( I_132 S_000 | S_000 S_000 )^0_{e} - ( K_152 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[70] = etfac[1] * AUX_INT__j_s_s_s[23] + 4 * one_over_2q * AUX_INT__i_s_s_s[17] - p_over_q * AUX_INT__k_s_s_s[30];

                    // ( J_142 S_000 | P_001 S_000 )^0 = z * ( J_142 S_000 | S_000 S_000 )^0_{e} + ( I_141 S_000 | S_000 S_000 )^0_{e} - ( K_143 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[71] = etfac[2] * AUX_INT__j_s_s_s[23] + 2 * one_over_2q * AUX_INT__i_s_s_s[16] - p_over_q * AUX_INT__k_s_s_s[31];

                    // ( J_133 S_000 | P_100 S_000 )^0 = x * ( J_133 S_000 | S_000 S_000 )^0_{e} + ( I_033 S_000 | S_000 S_000 )^0_{e} - ( K_233 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[72] = etfac[0] * AUX_INT__j_s_s_s[24] + 1 * one_over_2q * AUX_INT__i_s_s_s[24] - p_over_q * AUX_INT__k_s_s_s[24];

                    // ( J_133 S_000 | P_010 S_000 )^0 = y * ( J_133 S_000 | S_000 S_000 )^0_{e} + ( I_123 S_000 | S_000 S_000 )^0_{e} - ( K_143 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[73] = etfac[1] * AUX_INT__j_s_s_s[24] + 3 * one_over_2q * AUX_INT__i_s_s_s[18] - p_over_q * AUX_INT__k_s_s_s[31];

                    // ( J_133 S_000 | P_001 S_000 )^0 = z * ( J_133 S_000 | S_000 S_000 )^0_{e} + ( I_132 S_000 | S_000 S_000 )^0_{e} - ( K_134 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[74] = etfac[2] * AUX_INT__j_s_s_s[24] + 3 * one_over_2q * AUX_INT__i_s_s_s[17] - p_over_q * AUX_INT__k_s_s_s[32];

                    // ( J_124 S_000 | P_100 S_000 )^0 = x * ( J_124 S_000 | S_000 S_000 )^0_{e} + ( I_024 S_000 | S_000 S_000 )^0_{e} - ( K_224 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[75] = etfac[0] * AUX_INT__j_s_s_s[25] + 1 * one_over_2q * AUX_INT__i_s_s_s[25] - p_over_q * AUX_INT__k_s_s_s[25];

                    // ( J_124 S_000 | P_010 S_000 )^0 = y * ( J_124 S_000 | S_000 S_000 )^0_{e} + ( I_114 S_000 | S_000 S_000 )^0_{e} - ( K_134 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[76] = etfac[1] * AUX_INT__j_s_s_s[25] + 2 * one_over_2q * AUX_INT__i_s_s_s[19] - p_over_q * AUX_INT__k_s_s_s[32];

                    // ( J_124 S_000 | P_001 S_000 )^0 = z * ( J_124 S_000 | S_000 S_000 )^0_{e} + ( I_123 S_000 | S_000 S_000 )^0_{e} - ( K_125 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[77] = etfac[2] * AUX_INT__j_s_s_s[25] + 4 * one_over_2q * AUX_INT__i_s_s_s[18] - p_over_q * AUX_INT__k_s_s_s[33];

                    // ( J_115 S_000 | P_100 S_000 )^0 = x * ( J_115 S_000 | S_000 S_000 )^0_{e} + ( I_015 S_000 | S_000 S_000 )^0_{e} - ( K_215 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[78] = etfac[0] * AUX_INT__j_s_s_s[26] + 1 * one_over_2q * AUX_INT__i_s_s_s[26] - p_over_q * AUX_INT__k_s_s_s[26];

                    // ( J_115 S_000 | P_010 S_000 )^0 = y * ( J_115 S_000 | S_000 S_000 )^0_{e} + ( I_105 S_000 | S_000 S_000 )^0_{e} - ( K_125 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[79] = etfac[1] * AUX_INT__j_s_s_s[26] + 1 * one_over_2q * AUX_INT__i_s_s_s[20] - p_over_q * AUX_INT__k_s_s_s[33];

                    // ( J_115 S_000 | P_001 S_000 )^0 = z * ( J_115 S_000 | S_000 S_000 )^0_{e} + ( I_114 S_000 | S_000 S_000 )^0_{e} - ( K_116 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[80] = etfac[2] * AUX_INT__j_s_s_s[26] + 5 * one_over_2q * AUX_INT__i_s_s_s[19] - p_over_q * AUX_INT__k_s_s_s[34];

                    // ( J_106 S_000 | P_100 S_000 )^0 = x * ( J_106 S_000 | S_000 S_000 )^0_{e} + ( I_006 S_000 | S_000 S_000 )^0_{e} - ( K_206 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[81] = etfac[0] * AUX_INT__j_s_s_s[27] + 1 * one_over_2q * AUX_INT__i_s_s_s[27] - p_over_q * AUX_INT__k_s_s_s[27];

                    // ( J_106 S_000 | P_001 S_000 )^0 = z * ( J_106 S_000 | S_000 S_000 )^0_{e} + ( I_105 S_000 | S_000 S_000 )^0_{e} - ( K_107 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[83] = etfac[2] * AUX_INT__j_s_s_s[27] + 6 * one_over_2q * AUX_INT__i_s_s_s[20] - p_over_q * AUX_INT__k_s_s_s[35];

                    // ( J_070 S_000 | P_010 S_000 )^0 = y * ( J_070 S_000 | S_000 S_000 )^0_{e} + ( I_060 S_000 | S_000 S_000 )^0_{e} - ( K_080 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[85] = etfac[1] * AUX_INT__j_s_s_s[28] + 7 * one_over_2q * AUX_INT__i_s_s_s[21] - p_over_q * AUX_INT__k_s_s_s[36];

                    // ( J_061 S_000 | P_010 S_000 )^0 = y * ( J_061 S_000 | S_000 S_000 )^0_{e} + ( I_051 S_000 | S_000 S_000 )^0_{e} - ( K_071 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[88] = etfac[1] * AUX_INT__j_s_s_s[29] + 6 * one_over_2q * AUX_INT__i_s_s_s[22] - p_over_q * AUX_INT__k_s_s_s[37];

                    // ( J_052 S_000 | P_010 S_000 )^0 = y * ( J_052 S_000 | S_000 S_000 )^0_{e} + ( I_042 S_000 | S_000 S_000 )^0_{e} - ( K_062 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[91] = etfac[1] * AUX_INT__j_s_s_s[30] + 5 * one_over_2q * AUX_INT__i_s_s_s[23] - p_over_q * AUX_INT__k_s_s_s[38];

                    // ( J_052 S_000 | P_001 S_000 )^0 = z * ( J_052 S_000 | S_000 S_000 )^0_{e} + ( I_051 S_000 | S_000 S_000 )^0_{e} - ( K_053 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[92] = etfac[2] * AUX_INT__j_s_s_s[30] + 2 * one_over_2q * AUX_INT__i_s_s_s[22] - p_over_q * AUX_INT__k_s_s_s[39];

                    // ( J_043 S_000 | P_010 S_000 )^0 = y * ( J_043 S_000 | S_000 S_000 )^0_{e} + ( I_033 S_000 | S_000 S_000 )^0_{e} - ( K_053 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[94] = etfac[1] * AUX_INT__j_s_s_s[31] + 4 * one_over_2q * AUX_INT__i_s_s_s[24] - p_over_q * AUX_INT__k_s_s_s[39];

                    // ( J_043 S_000 | P_001 S_000 )^0 = z * ( J_043 S_000 | S_000 S_000 )^0_{e} + ( I_042 S_000 | S_000 S_000 )^0_{e} - ( K_044 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[95] = etfac[2] * AUX_INT__j_s_s_s[31] + 3 * one_over_2q * AUX_INT__i_s_s_s[23] - p_over_q * AUX_INT__k_s_s_s[40];

                    // ( J_034 S_000 | P_010 S_000 )^0 = y * ( J_034 S_000 | S_000 S_000 )^0_{e} + ( I_024 S_000 | S_000 S_000 )^0_{e} - ( K_044 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[97] = etfac[1] * AUX_INT__j_s_s_s[32] + 3 * one_over_2q * AUX_INT__i_s_s_s[25] - p_over_q * AUX_INT__k_s_s_s[40];

                    // ( J_034 S_000 | P_001 S_000 )^0 = z * ( J_034 S_000 | S_000 S_000 )^0_{e} + ( I_033 S_000 | S_000 S_000 )^0_{e} - ( K_035 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[98] = etfac[2] * AUX_INT__j_s_s_s[32] + 4 * one_over_2q * AUX_INT__i_s_s_s[24] - p_over_q * AUX_INT__k_s_s_s[41];

                    // ( J_025 S_000 | P_010 S_000 )^0 = y * ( J_025 S_000 | S_000 S_000 )^0_{e} + ( I_015 S_000 | S_000 S_000 )^0_{e} - ( K_035 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[100] = etfac[1] * AUX_INT__j_s_s_s[33] + 2 * one_over_2q * AUX_INT__i_s_s_s[26] - p_over_q * AUX_INT__k_s_s_s[41];

                    // ( J_025 S_000 | P_001 S_000 )^0 = z * ( J_025 S_000 | S_000 S_000 )^0_{e} + ( I_024 S_000 | S_000 S_000 )^0_{e} - ( K_026 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[101] = etfac[2] * AUX_INT__j_s_s_s[33] + 5 * one_over_2q * AUX_INT__i_s_s_s[25] - p_over_q * AUX_INT__k_s_s_s[42];

                    // ( J_016 S_000 | P_010 S_000 )^0 = y * ( J_016 S_000 | S_000 S_000 )^0_{e} + ( I_006 S_000 | S_000 S_000 )^0_{e} - ( K_026 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[103] = etfac[1] * AUX_INT__j_s_s_s[34] + 1 * one_over_2q * AUX_INT__i_s_s_s[27] - p_over_q * AUX_INT__k_s_s_s[42];

                    // ( J_016 S_000 | P_001 S_000 )^0 = z * ( J_016 S_000 | S_000 S_000 )^0_{e} + ( I_015 S_000 | S_000 S_000 )^0_{e} - ( K_017 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[104] = etfac[2] * AUX_INT__j_s_s_s[34] + 6 * one_over_2q * AUX_INT__i_s_s_s[26] - p_over_q * AUX_INT__k_s_s_s[43];

                    // ( J_007 S_000 | P_001 S_000 )^0 = z * ( J_007 S_000 | S_000 S_000 )^0_{e} + ( I_006 S_000 | S_000 S_000 )^0_{e} - ( K_008 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__j_s_p_s[107] = etfac[2] * AUX_INT__j_s_s_s[35] + 7 * one_over_2q * AUX_INT__i_s_s_s[27] - p_over_q * AUX_INT__k_s_s_s[44];

                    // ( I_600 S_000 | D_200 S_000 )^0 = x * ( I_600 S_000 | P_100 S_000 )^0 + ( H_500 S_000 | P_100 S_000 )^0 + ( I_600 S_000 | S_000 S_000 )^0_{e} - ( J_700 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[0] = etfac[0] * AUX_INT__i_s_p_s[0] + 6 * one_over_2q * AUX_INT__h_s_p_s[0] + 1 * one_over_2q * AUX_INT__i_s_s_s[0] - p_over_q * AUX_INT__j_s_p_s[0];

                    // ( I_510 S_000 | D_200 S_000 )^0 = x * ( I_510 S_000 | P_100 S_000 )^0 + ( H_410 S_000 | P_100 S_000 )^0 + ( I_510 S_000 | S_000 S_000 )^0_{e} - ( J_610 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[6] = etfac[0] * AUX_INT__i_s_p_s[3] + 5 * one_over_2q * AUX_INT__h_s_p_s[3] + 1 * one_over_2q * AUX_INT__i_s_s_s[1] - p_over_q * AUX_INT__j_s_p_s[3];

                    // ( I_510 S_000 | D_020 S_000 )^0 = y * ( I_510 S_000 | P_010 S_000 )^0 + ( H_500 S_000 | P_010 S_000 )^0 + ( I_510 S_000 | S_000 S_000 )^0_{e} - ( J_520 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[9] = etfac[1] * AUX_INT__i_s_p_s[4] + 1 * one_over_2q * AUX_INT__h_s_p_s[1] + 1 * one_over_2q * AUX_INT__i_s_s_s[1] - p_over_q * AUX_INT__j_s_p_s[10];

                    // ( I_510 S_000 | D_002 S_000 )^0 = z * ( I_510 S_000 | P_001 S_000 )^0 + ( I_510 S_000 | S_000 S_000 )^0_{e} - ( J_511 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[11] = etfac[2] * AUX_INT__i_s_p_s[5] + 1 * one_over_2q * AUX_INT__i_s_s_s[1] - p_over_q * AUX_INT__j_s_p_s[14];

                    // ( I_501 S_000 | D_200 S_000 )^0 = x * ( I_501 S_000 | P_100 S_000 )^0 + ( H_401 S_000 | P_100 S_000 )^0 + ( I_501 S_000 | S_000 S_000 )^0_{e} - ( J_601 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[12] = etfac[0] * AUX_INT__i_s_p_s[6] + 5 * one_over_2q * AUX_INT__h_s_p_s[6] + 1 * one_over_2q * AUX_INT__i_s_s_s[2] - p_over_q * AUX_INT__j_s_p_s[6];

                    // ( I_501 S_000 | D_020 S_000 )^0 = y * ( I_501 S_000 | P_010 S_000 )^0 + ( I_501 S_000 | S_000 S_000 )^0_{e} - ( J_511 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[15] = etfac[1] * AUX_INT__i_s_p_s[7] + 1 * one_over_2q * AUX_INT__i_s_s_s[2] - p_over_q * AUX_INT__j_s_p_s[13];

                    // ( I_501 S_000 | D_002 S_000 )^0 = z * ( I_501 S_000 | P_001 S_000 )^0 + ( H_500 S_000 | P_001 S_000 )^0 + ( I_501 S_000 | S_000 S_000 )^0_{e} - ( J_502 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[17] = etfac[2] * AUX_INT__i_s_p_s[8] + 1 * one_over_2q * AUX_INT__h_s_p_s[2] + 1 * one_over_2q * AUX_INT__i_s_s_s[2] - p_over_q * AUX_INT__j_s_p_s[17];

                    // ( I_420 S_000 | D_200 S_000 )^0 = x * ( I_420 S_000 | P_100 S_000 )^0 + ( H_320 S_000 | P_100 S_000 )^0 + ( I_420 S_000 | S_000 S_000 )^0_{e} - ( J_520 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[18] = etfac[0] * AUX_INT__i_s_p_s[9] + 4 * one_over_2q * AUX_INT__h_s_p_s[9] + 1 * one_over_2q * AUX_INT__i_s_s_s[3] - p_over_q * AUX_INT__j_s_p_s[9];

                    // ( I_420 S_000 | D_020 S_000 )^0 = y * ( I_420 S_000 | P_010 S_000 )^0 + ( H_410 S_000 | P_010 S_000 )^0 + ( I_420 S_000 | S_000 S_000 )^0_{e} - ( J_430 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[21] = etfac[1] * AUX_INT__i_s_p_s[10] + 2 * one_over_2q * AUX_INT__h_s_p_s[4] + 1 * one_over_2q * AUX_INT__i_s_s_s[3] - p_over_q * AUX_INT__j_s_p_s[19];

                    // ( I_420 S_000 | D_002 S_000 )^0 = z * ( I_420 S_000 | P_001 S_000 )^0 + ( I_420 S_000 | S_000 S_000 )^0_{e} - ( J_421 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[23] = etfac[2] * AUX_INT__i_s_p_s[11] + 1 * one_over_2q * AUX_INT__i_s_s_s[3] - p_over_q * AUX_INT__j_s_p_s[23];

                    // ( I_411 S_000 | D_200 S_000 )^0 = x * ( I_411 S_000 | P_100 S_000 )^0 + ( H_311 S_000 | P_100 S_000 )^0 + ( I_411 S_000 | S_000 S_000 )^0_{e} - ( J_511 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[24] = etfac[0] * AUX_INT__i_s_p_s[12] + 4 * one_over_2q * AUX_INT__h_s_p_s[12] + 1 * one_over_2q * AUX_INT__i_s_s_s[4] - p_over_q * AUX_INT__j_s_p_s[12];

                    // ( I_411 S_000 | D_020 S_000 )^0 = y * ( I_411 S_000 | P_010 S_000 )^0 + ( H_401 S_000 | P_010 S_000 )^0 + ( I_411 S_000 | S_000 S_000 )^0_{e} - ( J_421 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[27] = etfac[1] * AUX_INT__i_s_p_s[13] + 1 * one_over_2q * AUX_INT__h_s_p_s[7] + 1 * one_over_2q * AUX_INT__i_s_s_s[4] - p_over_q * AUX_INT__j_s_p_s[22];

                    // ( I_411 S_000 | D_002 S_000 )^0 = z * ( I_411 S_000 | P_001 S_000 )^0 + ( H_410 S_000 | P_001 S_000 )^0 + ( I_411 S_000 | S_000 S_000 )^0_{e} - ( J_412 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[29] = etfac[2] * AUX_INT__i_s_p_s[14] + 1 * one_over_2q * AUX_INT__h_s_p_s[5] + 1 * one_over_2q * AUX_INT__i_s_s_s[4] - p_over_q * AUX_INT__j_s_p_s[26];

                    // ( I_402 S_000 | D_200 S_000 )^0 = x * ( I_402 S_000 | P_100 S_000 )^0 + ( H_302 S_000 | P_100 S_000 )^0 + ( I_402 S_000 | S_000 S_000 )^0_{e} - ( J_502 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[30] = etfac[0] * AUX_INT__i_s_p_s[15] + 4 * one_over_2q * AUX_INT__h_s_p_s[15] + 1 * one_over_2q * AUX_INT__i_s_s_s[5] - p_over_q * AUX_INT__j_s_p_s[15];

                    // ( I_402 S_000 | D_020 S_000 )^0 = y * ( I_402 S_000 | P_010 S_000 )^0 + ( I_402 S_000 | S_000 S_000 )^0_{e} - ( J_412 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[33] = etfac[1] * AUX_INT__i_s_p_s[16] + 1 * one_over_2q * AUX_INT__i_s_s_s[5] - p_over_q * AUX_INT__j_s_p_s[25];

                    // ( I_402 S_000 | D_002 S_000 )^0 = z * ( I_402 S_000 | P_001 S_000 )^0 + ( H_401 S_000 | P_001 S_000 )^0 + ( I_402 S_000 | S_000 S_000 )^0_{e} - ( J_403 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[35] = etfac[2] * AUX_INT__i_s_p_s[17] + 2 * one_over_2q * AUX_INT__h_s_p_s[8] + 1 * one_over_2q * AUX_INT__i_s_s_s[5] - p_over_q * AUX_INT__j_s_p_s[29];

                    // ( I_330 S_000 | D_200 S_000 )^0 = x * ( I_330 S_000 | P_100 S_000 )^0 + ( H_230 S_000 | P_100 S_000 )^0 + ( I_330 S_000 | S_000 S_000 )^0_{e} - ( J_430 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[36] = etfac[0] * AUX_INT__i_s_p_s[18] + 3 * one_over_2q * AUX_INT__h_s_p_s[18] + 1 * one_over_2q * AUX_INT__i_s_s_s[6] - p_over_q * AUX_INT__j_s_p_s[18];

                    // ( I_330 S_000 | D_020 S_000 )^0 = y * ( I_330 S_000 | P_010 S_000 )^0 + ( H_320 S_000 | P_010 S_000 )^0 + ( I_330 S_000 | S_000 S_000 )^0_{e} - ( J_340 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[39] = etfac[1] * AUX_INT__i_s_p_s[19] + 3 * one_over_2q * AUX_INT__h_s_p_s[10] + 1 * one_over_2q * AUX_INT__i_s_s_s[6] - p_over_q * AUX_INT__j_s_p_s[31];

                    // ( I_330 S_000 | D_002 S_000 )^0 = z * ( I_330 S_000 | P_001 S_000 )^0 + ( I_330 S_000 | S_000 S_000 )^0_{e} - ( J_331 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[41] = etfac[2] * AUX_INT__i_s_p_s[20] + 1 * one_over_2q * AUX_INT__i_s_s_s[6] - p_over_q * AUX_INT__j_s_p_s[35];

                    // ( I_321 S_000 | D_200 S_000 )^0 = x * ( I_321 S_000 | P_100 S_000 )^0 + ( H_221 S_000 | P_100 S_000 )^0 + ( I_321 S_000 | S_000 S_000 )^0_{e} - ( J_421 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[42] = etfac[0] * AUX_INT__i_s_p_s[21] + 3 * one_over_2q * AUX_INT__h_s_p_s[21] + 1 * one_over_2q * AUX_INT__i_s_s_s[7] - p_over_q * AUX_INT__j_s_p_s[21];

                    // ( I_321 S_000 | D_020 S_000 )^0 = y * ( I_321 S_000 | P_010 S_000 )^0 + ( H_311 S_000 | P_010 S_000 )^0 + ( I_321 S_000 | S_000 S_000 )^0_{e} - ( J_331 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[45] = etfac[1] * AUX_INT__i_s_p_s[22] + 2 * one_over_2q * AUX_INT__h_s_p_s[13] + 1 * one_over_2q * AUX_INT__i_s_s_s[7] - p_over_q * AUX_INT__j_s_p_s[34];

                    // ( I_321 S_000 | D_002 S_000 )^0 = z * ( I_321 S_000 | P_001 S_000 )^0 + ( H_320 S_000 | P_001 S_000 )^0 + ( I_321 S_000 | S_000 S_000 )^0_{e} - ( J_322 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[47] = etfac[2] * AUX_INT__i_s_p_s[23] + 1 * one_over_2q * AUX_INT__h_s_p_s[11] + 1 * one_over_2q * AUX_INT__i_s_s_s[7] - p_over_q * AUX_INT__j_s_p_s[38];

                    // ( I_312 S_000 | D_200 S_000 )^0 = x * ( I_312 S_000 | P_100 S_000 )^0 + ( H_212 S_000 | P_100 S_000 )^0 + ( I_312 S_000 | S_000 S_000 )^0_{e} - ( J_412 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[48] = etfac[0] * AUX_INT__i_s_p_s[24] + 3 * one_over_2q * AUX_INT__h_s_p_s[24] + 1 * one_over_2q * AUX_INT__i_s_s_s[8] - p_over_q * AUX_INT__j_s_p_s[24];

                    // ( I_312 S_000 | D_020 S_000 )^0 = y * ( I_312 S_000 | P_010 S_000 )^0 + ( H_302 S_000 | P_010 S_000 )^0 + ( I_312 S_000 | S_000 S_000 )^0_{e} - ( J_322 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[51] = etfac[1] * AUX_INT__i_s_p_s[25] + 1 * one_over_2q * AUX_INT__h_s_p_s[16] + 1 * one_over_2q * AUX_INT__i_s_s_s[8] - p_over_q * AUX_INT__j_s_p_s[37];

                    // ( I_312 S_000 | D_002 S_000 )^0 = z * ( I_312 S_000 | P_001 S_000 )^0 + ( H_311 S_000 | P_001 S_000 )^0 + ( I_312 S_000 | S_000 S_000 )^0_{e} - ( J_313 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[53] = etfac[2] * AUX_INT__i_s_p_s[26] + 2 * one_over_2q * AUX_INT__h_s_p_s[14] + 1 * one_over_2q * AUX_INT__i_s_s_s[8] - p_over_q * AUX_INT__j_s_p_s[41];

                    // ( I_303 S_000 | D_200 S_000 )^0 = x * ( I_303 S_000 | P_100 S_000 )^0 + ( H_203 S_000 | P_100 S_000 )^0 + ( I_303 S_000 | S_000 S_000 )^0_{e} - ( J_403 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[54] = etfac[0] * AUX_INT__i_s_p_s[27] + 3 * one_over_2q * AUX_INT__h_s_p_s[27] + 1 * one_over_2q * AUX_INT__i_s_s_s[9] - p_over_q * AUX_INT__j_s_p_s[27];

                    // ( I_303 S_000 | D_020 S_000 )^0 = y * ( I_303 S_000 | P_010 S_000 )^0 + ( I_303 S_000 | S_000 S_000 )^0_{e} - ( J_313 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[57] = etfac[1] * AUX_INT__i_s_p_s[28] + 1 * one_over_2q * AUX_INT__i_s_s_s[9] - p_over_q * AUX_INT__j_s_p_s[40];

                    // ( I_303 S_000 | D_002 S_000 )^0 = z * ( I_303 S_000 | P_001 S_000 )^0 + ( H_302 S_000 | P_001 S_000 )^0 + ( I_303 S_000 | S_000 S_000 )^0_{e} - ( J_304 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[59] = etfac[2] * AUX_INT__i_s_p_s[29] + 3 * one_over_2q * AUX_INT__h_s_p_s[17] + 1 * one_over_2q * AUX_INT__i_s_s_s[9] - p_over_q * AUX_INT__j_s_p_s[44];

                    // ( I_240 S_000 | D_200 S_000 )^0 = x * ( I_240 S_000 | P_100 S_000 )^0 + ( H_140 S_000 | P_100 S_000 )^0 + ( I_240 S_000 | S_000 S_000 )^0_{e} - ( J_340 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[60] = etfac[0] * AUX_INT__i_s_p_s[30] + 2 * one_over_2q * AUX_INT__h_s_p_s[30] + 1 * one_over_2q * AUX_INT__i_s_s_s[10] - p_over_q * AUX_INT__j_s_p_s[30];

                    // ( I_240 S_000 | D_020 S_000 )^0 = y * ( I_240 S_000 | P_010 S_000 )^0 + ( H_230 S_000 | P_010 S_000 )^0 + ( I_240 S_000 | S_000 S_000 )^0_{e} - ( J_250 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[63] = etfac[1] * AUX_INT__i_s_p_s[31] + 4 * one_over_2q * AUX_INT__h_s_p_s[19] + 1 * one_over_2q * AUX_INT__i_s_s_s[10] - p_over_q * AUX_INT__j_s_p_s[46];

                    // ( I_240 S_000 | D_002 S_000 )^0 = z * ( I_240 S_000 | P_001 S_000 )^0 + ( I_240 S_000 | S_000 S_000 )^0_{e} - ( J_241 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[65] = etfac[2] * AUX_INT__i_s_p_s[32] + 1 * one_over_2q * AUX_INT__i_s_s_s[10] - p_over_q * AUX_INT__j_s_p_s[50];

                    // ( I_231 S_000 | D_200 S_000 )^0 = x * ( I_231 S_000 | P_100 S_000 )^0 + ( H_131 S_000 | P_100 S_000 )^0 + ( I_231 S_000 | S_000 S_000 )^0_{e} - ( J_331 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[66] = etfac[0] * AUX_INT__i_s_p_s[33] + 2 * one_over_2q * AUX_INT__h_s_p_s[33] + 1 * one_over_2q * AUX_INT__i_s_s_s[11] - p_over_q * AUX_INT__j_s_p_s[33];

                    // ( I_231 S_000 | D_020 S_000 )^0 = y * ( I_231 S_000 | P_010 S_000 )^0 + ( H_221 S_000 | P_010 S_000 )^0 + ( I_231 S_000 | S_000 S_000 )^0_{e} - ( J_241 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[69] = etfac[1] * AUX_INT__i_s_p_s[34] + 3 * one_over_2q * AUX_INT__h_s_p_s[22] + 1 * one_over_2q * AUX_INT__i_s_s_s[11] - p_over_q * AUX_INT__j_s_p_s[49];

                    // ( I_231 S_000 | D_002 S_000 )^0 = z * ( I_231 S_000 | P_001 S_000 )^0 + ( H_230 S_000 | P_001 S_000 )^0 + ( I_231 S_000 | S_000 S_000 )^0_{e} - ( J_232 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[71] = etfac[2] * AUX_INT__i_s_p_s[35] + 1 * one_over_2q * AUX_INT__h_s_p_s[20] + 1 * one_over_2q * AUX_INT__i_s_s_s[11] - p_over_q * AUX_INT__j_s_p_s[53];

                    // ( I_222 S_000 | D_200 S_000 )^0 = x * ( I_222 S_000 | P_100 S_000 )^0 + ( H_122 S_000 | P_100 S_000 )^0 + ( I_222 S_000 | S_000 S_000 )^0_{e} - ( J_322 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[72] = etfac[0] * AUX_INT__i_s_p_s[36] + 2 * one_over_2q * AUX_INT__h_s_p_s[36] + 1 * one_over_2q * AUX_INT__i_s_s_s[12] - p_over_q * AUX_INT__j_s_p_s[36];

                    // ( I_222 S_000 | D_020 S_000 )^0 = y * ( I_222 S_000 | P_010 S_000 )^0 + ( H_212 S_000 | P_010 S_000 )^0 + ( I_222 S_000 | S_000 S_000 )^0_{e} - ( J_232 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[75] = etfac[1] * AUX_INT__i_s_p_s[37] + 2 * one_over_2q * AUX_INT__h_s_p_s[25] + 1 * one_over_2q * AUX_INT__i_s_s_s[12] - p_over_q * AUX_INT__j_s_p_s[52];

                    // ( I_222 S_000 | D_002 S_000 )^0 = z * ( I_222 S_000 | P_001 S_000 )^0 + ( H_221 S_000 | P_001 S_000 )^0 + ( I_222 S_000 | S_000 S_000 )^0_{e} - ( J_223 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[77] = etfac[2] * AUX_INT__i_s_p_s[38] + 2 * one_over_2q * AUX_INT__h_s_p_s[23] + 1 * one_over_2q * AUX_INT__i_s_s_s[12] - p_over_q * AUX_INT__j_s_p_s[56];

                    // ( I_213 S_000 | D_200 S_000 )^0 = x * ( I_213 S_000 | P_100 S_000 )^0 + ( H_113 S_000 | P_100 S_000 )^0 + ( I_213 S_000 | S_000 S_000 )^0_{e} - ( J_313 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[78] = etfac[0] * AUX_INT__i_s_p_s[39] + 2 * one_over_2q * AUX_INT__h_s_p_s[39] + 1 * one_over_2q * AUX_INT__i_s_s_s[13] - p_over_q * AUX_INT__j_s_p_s[39];

                    // ( I_213 S_000 | D_020 S_000 )^0 = y * ( I_213 S_000 | P_010 S_000 )^0 + ( H_203 S_000 | P_010 S_000 )^0 + ( I_213 S_000 | S_000 S_000 )^0_{e} - ( J_223 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[81] = etfac[1] * AUX_INT__i_s_p_s[40] + 1 * one_over_2q * AUX_INT__h_s_p_s[28] + 1 * one_over_2q * AUX_INT__i_s_s_s[13] - p_over_q * AUX_INT__j_s_p_s[55];

                    // ( I_213 S_000 | D_002 S_000 )^0 = z * ( I_213 S_000 | P_001 S_000 )^0 + ( H_212 S_000 | P_001 S_000 )^0 + ( I_213 S_000 | S_000 S_000 )^0_{e} - ( J_214 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[83] = etfac[2] * AUX_INT__i_s_p_s[41] + 3 * one_over_2q * AUX_INT__h_s_p_s[26] + 1 * one_over_2q * AUX_INT__i_s_s_s[13] - p_over_q * AUX_INT__j_s_p_s[59];

                    // ( I_204 S_000 | D_200 S_000 )^0 = x * ( I_204 S_000 | P_100 S_000 )^0 + ( H_104 S_000 | P_100 S_000 )^0 + ( I_204 S_000 | S_000 S_000 )^0_{e} - ( J_304 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[84] = etfac[0] * AUX_INT__i_s_p_s[42] + 2 * one_over_2q * AUX_INT__h_s_p_s[42] + 1 * one_over_2q * AUX_INT__i_s_s_s[14] - p_over_q * AUX_INT__j_s_p_s[42];

                    // ( I_204 S_000 | D_020 S_000 )^0 = y * ( I_204 S_000 | P_010 S_000 )^0 + ( I_204 S_000 | S_000 S_000 )^0_{e} - ( J_214 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[87] = etfac[1] * AUX_INT__i_s_p_s[43] + 1 * one_over_2q * AUX_INT__i_s_s_s[14] - p_over_q * AUX_INT__j_s_p_s[58];

                    // ( I_204 S_000 | D_002 S_000 )^0 = z * ( I_204 S_000 | P_001 S_000 )^0 + ( H_203 S_000 | P_001 S_000 )^0 + ( I_204 S_000 | S_000 S_000 )^0_{e} - ( J_205 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[89] = etfac[2] * AUX_INT__i_s_p_s[44] + 4 * one_over_2q * AUX_INT__h_s_p_s[29] + 1 * one_over_2q * AUX_INT__i_s_s_s[14] - p_over_q * AUX_INT__j_s_p_s[62];

                    // ( I_150 S_000 | D_200 S_000 )^0 = x * ( I_150 S_000 | P_100 S_000 )^0 + ( H_050 S_000 | P_100 S_000 )^0 + ( I_150 S_000 | S_000 S_000 )^0_{e} - ( J_250 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[90] = etfac[0] * AUX_INT__i_s_p_s[45] + 1 * one_over_2q * AUX_INT__h_s_p_s[45] + 1 * one_over_2q * AUX_INT__i_s_s_s[15] - p_over_q * AUX_INT__j_s_p_s[45];

                    // ( I_150 S_000 | D_020 S_000 )^0 = y * ( I_150 S_000 | P_010 S_000 )^0 + ( H_140 S_000 | P_010 S_000 )^0 + ( I_150 S_000 | S_000 S_000 )^0_{e} - ( J_160 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[93] = etfac[1] * AUX_INT__i_s_p_s[46] + 5 * one_over_2q * AUX_INT__h_s_p_s[31] + 1 * one_over_2q * AUX_INT__i_s_s_s[15] - p_over_q * AUX_INT__j_s_p_s[64];

                    // ( I_150 S_000 | D_002 S_000 )^0 = z * ( I_150 S_000 | P_001 S_000 )^0 + ( I_150 S_000 | S_000 S_000 )^0_{e} - ( J_151 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[95] = etfac[2] * AUX_INT__i_s_p_s[47] + 1 * one_over_2q * AUX_INT__i_s_s_s[15] - p_over_q * AUX_INT__j_s_p_s[68];

                    // ( I_141 S_000 | D_200 S_000 )^0 = x * ( I_141 S_000 | P_100 S_000 )^0 + ( H_041 S_000 | P_100 S_000 )^0 + ( I_141 S_000 | S_000 S_000 )^0_{e} - ( J_241 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[96] = etfac[0] * AUX_INT__i_s_p_s[48] + 1 * one_over_2q * AUX_INT__h_s_p_s[48] + 1 * one_over_2q * AUX_INT__i_s_s_s[16] - p_over_q * AUX_INT__j_s_p_s[48];

                    // ( I_141 S_000 | D_020 S_000 )^0 = y * ( I_141 S_000 | P_010 S_000 )^0 + ( H_131 S_000 | P_010 S_000 )^0 + ( I_141 S_000 | S_000 S_000 )^0_{e} - ( J_151 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[99] = etfac[1] * AUX_INT__i_s_p_s[49] + 4 * one_over_2q * AUX_INT__h_s_p_s[34] + 1 * one_over_2q * AUX_INT__i_s_s_s[16] - p_over_q * AUX_INT__j_s_p_s[67];

                    // ( I_141 S_000 | D_002 S_000 )^0 = z * ( I_141 S_000 | P_001 S_000 )^0 + ( H_140 S_000 | P_001 S_000 )^0 + ( I_141 S_000 | S_000 S_000 )^0_{e} - ( J_142 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[101] = etfac[2] * AUX_INT__i_s_p_s[50] + 1 * one_over_2q * AUX_INT__h_s_p_s[32] + 1 * one_over_2q * AUX_INT__i_s_s_s[16] - p_over_q * AUX_INT__j_s_p_s[71];

                    // ( I_132 S_000 | D_200 S_000 )^0 = x * ( I_132 S_000 | P_100 S_000 )^0 + ( H_032 S_000 | P_100 S_000 )^0 + ( I_132 S_000 | S_000 S_000 )^0_{e} - ( J_232 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[102] = etfac[0] * AUX_INT__i_s_p_s[51] + 1 * one_over_2q * AUX_INT__h_s_p_s[51] + 1 * one_over_2q * AUX_INT__i_s_s_s[17] - p_over_q * AUX_INT__j_s_p_s[51];

                    // ( I_132 S_000 | D_020 S_000 )^0 = y * ( I_132 S_000 | P_010 S_000 )^0 + ( H_122 S_000 | P_010 S_000 )^0 + ( I_132 S_000 | S_000 S_000 )^0_{e} - ( J_142 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[105] = etfac[1] * AUX_INT__i_s_p_s[52] + 3 * one_over_2q * AUX_INT__h_s_p_s[37] + 1 * one_over_2q * AUX_INT__i_s_s_s[17] - p_over_q * AUX_INT__j_s_p_s[70];

                    // ( I_132 S_000 | D_002 S_000 )^0 = z * ( I_132 S_000 | P_001 S_000 )^0 + ( H_131 S_000 | P_001 S_000 )^0 + ( I_132 S_000 | S_000 S_000 )^0_{e} - ( J_133 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[107] = etfac[2] * AUX_INT__i_s_p_s[53] + 2 * one_over_2q * AUX_INT__h_s_p_s[35] + 1 * one_over_2q * AUX_INT__i_s_s_s[17] - p_over_q * AUX_INT__j_s_p_s[74];

                    // ( I_123 S_000 | D_200 S_000 )^0 = x * ( I_123 S_000 | P_100 S_000 )^0 + ( H_023 S_000 | P_100 S_000 )^0 + ( I_123 S_000 | S_000 S_000 )^0_{e} - ( J_223 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[108] = etfac[0] * AUX_INT__i_s_p_s[54] + 1 * one_over_2q * AUX_INT__h_s_p_s[54] + 1 * one_over_2q * AUX_INT__i_s_s_s[18] - p_over_q * AUX_INT__j_s_p_s[54];

                    // ( I_123 S_000 | D_020 S_000 )^0 = y * ( I_123 S_000 | P_010 S_000 )^0 + ( H_113 S_000 | P_010 S_000 )^0 + ( I_123 S_000 | S_000 S_000 )^0_{e} - ( J_133 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[111] = etfac[1] * AUX_INT__i_s_p_s[55] + 2 * one_over_2q * AUX_INT__h_s_p_s[40] + 1 * one_over_2q * AUX_INT__i_s_s_s[18] - p_over_q * AUX_INT__j_s_p_s[73];

                    // ( I_123 S_000 | D_002 S_000 )^0 = z * ( I_123 S_000 | P_001 S_000 )^0 + ( H_122 S_000 | P_001 S_000 )^0 + ( I_123 S_000 | S_000 S_000 )^0_{e} - ( J_124 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[113] = etfac[2] * AUX_INT__i_s_p_s[56] + 3 * one_over_2q * AUX_INT__h_s_p_s[38] + 1 * one_over_2q * AUX_INT__i_s_s_s[18] - p_over_q * AUX_INT__j_s_p_s[77];

                    // ( I_114 S_000 | D_200 S_000 )^0 = x * ( I_114 S_000 | P_100 S_000 )^0 + ( H_014 S_000 | P_100 S_000 )^0 + ( I_114 S_000 | S_000 S_000 )^0_{e} - ( J_214 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[114] = etfac[0] * AUX_INT__i_s_p_s[57] + 1 * one_over_2q * AUX_INT__h_s_p_s[57] + 1 * one_over_2q * AUX_INT__i_s_s_s[19] - p_over_q * AUX_INT__j_s_p_s[57];

                    // ( I_114 S_000 | D_020 S_000 )^0 = y * ( I_114 S_000 | P_010 S_000 )^0 + ( H_104 S_000 | P_010 S_000 )^0 + ( I_114 S_000 | S_000 S_000 )^0_{e} - ( J_124 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[117] = etfac[1] * AUX_INT__i_s_p_s[58] + 1 * one_over_2q * AUX_INT__h_s_p_s[43] + 1 * one_over_2q * AUX_INT__i_s_s_s[19] - p_over_q * AUX_INT__j_s_p_s[76];

                    // ( I_114 S_000 | D_002 S_000 )^0 = z * ( I_114 S_000 | P_001 S_000 )^0 + ( H_113 S_000 | P_001 S_000 )^0 + ( I_114 S_000 | S_000 S_000 )^0_{e} - ( J_115 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[119] = etfac[2] * AUX_INT__i_s_p_s[59] + 4 * one_over_2q * AUX_INT__h_s_p_s[41] + 1 * one_over_2q * AUX_INT__i_s_s_s[19] - p_over_q * AUX_INT__j_s_p_s[80];

                    // ( I_105 S_000 | D_200 S_000 )^0 = x * ( I_105 S_000 | P_100 S_000 )^0 + ( H_005 S_000 | P_100 S_000 )^0 + ( I_105 S_000 | S_000 S_000 )^0_{e} - ( J_205 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[120] = etfac[0] * AUX_INT__i_s_p_s[60] + 1 * one_over_2q * AUX_INT__h_s_p_s[60] + 1 * one_over_2q * AUX_INT__i_s_s_s[20] - p_over_q * AUX_INT__j_s_p_s[60];

                    // ( I_105 S_000 | D_020 S_000 )^0 = y * ( I_105 S_000 | P_010 S_000 )^0 + ( I_105 S_000 | S_000 S_000 )^0_{e} - ( J_115 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[123] = etfac[1] * AUX_INT__i_s_p_s[61] + 1 * one_over_2q * AUX_INT__i_s_s_s[20] - p_over_q * AUX_INT__j_s_p_s[79];

                    // ( I_105 S_000 | D_002 S_000 )^0 = z * ( I_105 S_000 | P_001 S_000 )^0 + ( H_104 S_000 | P_001 S_000 )^0 + ( I_105 S_000 | S_000 S_000 )^0_{e} - ( J_106 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[125] = etfac[2] * AUX_INT__i_s_p_s[62] + 5 * one_over_2q * AUX_INT__h_s_p_s[44] + 1 * one_over_2q * AUX_INT__i_s_s_s[20] - p_over_q * AUX_INT__j_s_p_s[83];

                    // ( I_060 S_000 | D_200 S_000 )^0 = x * ( I_060 S_000 | P_100 S_000 )^0 + ( I_060 S_000 | S_000 S_000 )^0_{e} - ( J_160 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[126] = etfac[0] * AUX_INT__i_s_p_s[63] + 1 * one_over_2q * AUX_INT__i_s_s_s[21] - p_over_q * AUX_INT__j_s_p_s[63];

                    // ( I_060 S_000 | D_020 S_000 )^0 = y * ( I_060 S_000 | P_010 S_000 )^0 + ( H_050 S_000 | P_010 S_000 )^0 + ( I_060 S_000 | S_000 S_000 )^0_{e} - ( J_070 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[129] = etfac[1] * AUX_INT__i_s_p_s[64] + 6 * one_over_2q * AUX_INT__h_s_p_s[46] + 1 * one_over_2q * AUX_INT__i_s_s_s[21] - p_over_q * AUX_INT__j_s_p_s[85];

                    // ( I_051 S_000 | D_200 S_000 )^0 = x * ( I_051 S_000 | P_100 S_000 )^0 + ( I_051 S_000 | S_000 S_000 )^0_{e} - ( J_151 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[132] = etfac[0] * AUX_INT__i_s_p_s[66] + 1 * one_over_2q * AUX_INT__i_s_s_s[22] - p_over_q * AUX_INT__j_s_p_s[66];

                    // ( I_051 S_000 | D_020 S_000 )^0 = y * ( I_051 S_000 | P_010 S_000 )^0 + ( H_041 S_000 | P_010 S_000 )^0 + ( I_051 S_000 | S_000 S_000 )^0_{e} - ( J_061 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[135] = etfac[1] * AUX_INT__i_s_p_s[67] + 5 * one_over_2q * AUX_INT__h_s_p_s[49] + 1 * one_over_2q * AUX_INT__i_s_s_s[22] - p_over_q * AUX_INT__j_s_p_s[88];

                    // ( I_051 S_000 | D_002 S_000 )^0 = z * ( I_051 S_000 | P_001 S_000 )^0 + ( H_050 S_000 | P_001 S_000 )^0 + ( I_051 S_000 | S_000 S_000 )^0_{e} - ( J_052 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[137] = etfac[2] * AUX_INT__i_s_p_s[68] + 1 * one_over_2q * AUX_INT__h_s_p_s[47] + 1 * one_over_2q * AUX_INT__i_s_s_s[22] - p_over_q * AUX_INT__j_s_p_s[92];

                    // ( I_042 S_000 | D_200 S_000 )^0 = x * ( I_042 S_000 | P_100 S_000 )^0 + ( I_042 S_000 | S_000 S_000 )^0_{e} - ( J_142 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[138] = etfac[0] * AUX_INT__i_s_p_s[69] + 1 * one_over_2q * AUX_INT__i_s_s_s[23] - p_over_q * AUX_INT__j_s_p_s[69];

                    // ( I_042 S_000 | D_020 S_000 )^0 = y * ( I_042 S_000 | P_010 S_000 )^0 + ( H_032 S_000 | P_010 S_000 )^0 + ( I_042 S_000 | S_000 S_000 )^0_{e} - ( J_052 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[141] = etfac[1] * AUX_INT__i_s_p_s[70] + 4 * one_over_2q * AUX_INT__h_s_p_s[52] + 1 * one_over_2q * AUX_INT__i_s_s_s[23] - p_over_q * AUX_INT__j_s_p_s[91];

                    // ( I_042 S_000 | D_002 S_000 )^0 = z * ( I_042 S_000 | P_001 S_000 )^0 + ( H_041 S_000 | P_001 S_000 )^0 + ( I_042 S_000 | S_000 S_000 )^0_{e} - ( J_043 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[143] = etfac[2] * AUX_INT__i_s_p_s[71] + 2 * one_over_2q * AUX_INT__h_s_p_s[50] + 1 * one_over_2q * AUX_INT__i_s_s_s[23] - p_over_q * AUX_INT__j_s_p_s[95];

                    // ( I_033 S_000 | D_200 S_000 )^0 = x * ( I_033 S_000 | P_100 S_000 )^0 + ( I_033 S_000 | S_000 S_000 )^0_{e} - ( J_133 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[144] = etfac[0] * AUX_INT__i_s_p_s[72] + 1 * one_over_2q * AUX_INT__i_s_s_s[24] - p_over_q * AUX_INT__j_s_p_s[72];

                    // ( I_033 S_000 | D_020 S_000 )^0 = y * ( I_033 S_000 | P_010 S_000 )^0 + ( H_023 S_000 | P_010 S_000 )^0 + ( I_033 S_000 | S_000 S_000 )^0_{e} - ( J_043 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[147] = etfac[1] * AUX_INT__i_s_p_s[73] + 3 * one_over_2q * AUX_INT__h_s_p_s[55] + 1 * one_over_2q * AUX_INT__i_s_s_s[24] - p_over_q * AUX_INT__j_s_p_s[94];

                    // ( I_033 S_000 | D_002 S_000 )^0 = z * ( I_033 S_000 | P_001 S_000 )^0 + ( H_032 S_000 | P_001 S_000 )^0 + ( I_033 S_000 | S_000 S_000 )^0_{e} - ( J_034 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[149] = etfac[2] * AUX_INT__i_s_p_s[74] + 3 * one_over_2q * AUX_INT__h_s_p_s[53] + 1 * one_over_2q * AUX_INT__i_s_s_s[24] - p_over_q * AUX_INT__j_s_p_s[98];

                    // ( I_024 S_000 | D_200 S_000 )^0 = x * ( I_024 S_000 | P_100 S_000 )^0 + ( I_024 S_000 | S_000 S_000 )^0_{e} - ( J_124 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[150] = etfac[0] * AUX_INT__i_s_p_s[75] + 1 * one_over_2q * AUX_INT__i_s_s_s[25] - p_over_q * AUX_INT__j_s_p_s[75];

                    // ( I_024 S_000 | D_020 S_000 )^0 = y * ( I_024 S_000 | P_010 S_000 )^0 + ( H_014 S_000 | P_010 S_000 )^0 + ( I_024 S_000 | S_000 S_000 )^0_{e} - ( J_034 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[153] = etfac[1] * AUX_INT__i_s_p_s[76] + 2 * one_over_2q * AUX_INT__h_s_p_s[58] + 1 * one_over_2q * AUX_INT__i_s_s_s[25] - p_over_q * AUX_INT__j_s_p_s[97];

                    // ( I_024 S_000 | D_002 S_000 )^0 = z * ( I_024 S_000 | P_001 S_000 )^0 + ( H_023 S_000 | P_001 S_000 )^0 + ( I_024 S_000 | S_000 S_000 )^0_{e} - ( J_025 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[155] = etfac[2] * AUX_INT__i_s_p_s[77] + 4 * one_over_2q * AUX_INT__h_s_p_s[56] + 1 * one_over_2q * AUX_INT__i_s_s_s[25] - p_over_q * AUX_INT__j_s_p_s[101];

                    // ( I_015 S_000 | D_200 S_000 )^0 = x * ( I_015 S_000 | P_100 S_000 )^0 + ( I_015 S_000 | S_000 S_000 )^0_{e} - ( J_115 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[156] = etfac[0] * AUX_INT__i_s_p_s[78] + 1 * one_over_2q * AUX_INT__i_s_s_s[26] - p_over_q * AUX_INT__j_s_p_s[78];

                    // ( I_015 S_000 | D_020 S_000 )^0 = y * ( I_015 S_000 | P_010 S_000 )^0 + ( H_005 S_000 | P_010 S_000 )^0 + ( I_015 S_000 | S_000 S_000 )^0_{e} - ( J_025 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[159] = etfac[1] * AUX_INT__i_s_p_s[79] + 1 * one_over_2q * AUX_INT__h_s_p_s[61] + 1 * one_over_2q * AUX_INT__i_s_s_s[26] - p_over_q * AUX_INT__j_s_p_s[100];

                    // ( I_015 S_000 | D_002 S_000 )^0 = z * ( I_015 S_000 | P_001 S_000 )^0 + ( H_014 S_000 | P_001 S_000 )^0 + ( I_015 S_000 | S_000 S_000 )^0_{e} - ( J_016 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[161] = etfac[2] * AUX_INT__i_s_p_s[80] + 5 * one_over_2q * AUX_INT__h_s_p_s[59] + 1 * one_over_2q * AUX_INT__i_s_s_s[26] - p_over_q * AUX_INT__j_s_p_s[104];

                    // ( I_006 S_000 | D_200 S_000 )^0 = x * ( I_006 S_000 | P_100 S_000 )^0 + ( I_006 S_000 | S_000 S_000 )^0_{e} - ( J_106 S_000 | P_100 S_000 )^0
                    AUX_INT__i_s_d_s[162] = etfac[0] * AUX_INT__i_s_p_s[81] + 1 * one_over_2q * AUX_INT__i_s_s_s[27] - p_over_q * AUX_INT__j_s_p_s[81];

                    // ( I_006 S_000 | D_020 S_000 )^0 = y * ( I_006 S_000 | P_010 S_000 )^0 + ( I_006 S_000 | S_000 S_000 )^0_{e} - ( J_016 S_000 | P_010 S_000 )^0
                    AUX_INT__i_s_d_s[165] = etfac[1] * AUX_INT__i_s_p_s[82] + 1 * one_over_2q * AUX_INT__i_s_s_s[27] - p_over_q * AUX_INT__j_s_p_s[103];

                    // ( I_006 S_000 | D_002 S_000 )^0 = z * ( I_006 S_000 | P_001 S_000 )^0 + ( H_005 S_000 | P_001 S_000 )^0 + ( I_006 S_000 | S_000 S_000 )^0_{e} - ( J_007 S_000 | P_001 S_000 )^0
                    AUX_INT__i_s_d_s[167] = etfac[2] * AUX_INT__i_s_p_s[83] + 6 * one_over_2q * AUX_INT__h_s_p_s[62] + 1 * one_over_2q * AUX_INT__i_s_s_s[27] - p_over_q * AUX_INT__j_s_p_s[107];

                    // ( H_500 S_000 | F_300 S_000 )^0 = x * ( H_500 S_000 | D_200 S_000 )^0 + ( G_400 S_000 | D_200 S_000 )^0 + ( H_500 S_000 | P_100 S_000 )^0 - ( I_600 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[0] = etfac[0] * AUX_INT__h_s_d_s[0] + 5 * one_over_2q * AUX_INT__g_s_d_s[0] + 2 * one_over_2q * AUX_INT__h_s_p_s[0] - p_over_q * AUX_INT__i_s_d_s[0];

                    // ( H_500 S_000 | F_030 S_000 )^0 = y * ( H_500 S_000 | D_020 S_000 )^0 + ( H_500 S_000 | P_010 S_000 )^0 - ( I_510 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[6] = etfac[1] * AUX_INT__h_s_d_s[3] + 2 * one_over_2q * AUX_INT__h_s_p_s[1] - p_over_q * AUX_INT__i_s_d_s[9];

                    // ( H_500 S_000 | F_003 S_000 )^0 = z * ( H_500 S_000 | D_002 S_000 )^0 + ( H_500 S_000 | P_001 S_000 )^0 - ( I_501 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[9] = etfac[2] * AUX_INT__h_s_d_s[5] + 2 * one_over_2q * AUX_INT__h_s_p_s[2] - p_over_q * AUX_INT__i_s_d_s[17];

                    // ( H_410 S_000 | F_300 S_000 )^0 = x * ( H_410 S_000 | D_200 S_000 )^0 + ( G_310 S_000 | D_200 S_000 )^0 + ( H_410 S_000 | P_100 S_000 )^0 - ( I_510 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[10] = etfac[0] * AUX_INT__h_s_d_s[6] + 4 * one_over_2q * AUX_INT__g_s_d_s[6] + 2 * one_over_2q * AUX_INT__h_s_p_s[3] - p_over_q * AUX_INT__i_s_d_s[6];

                    // ( H_410 S_000 | F_210 S_000 )^0 = y * ( H_410 S_000 | D_200 S_000 )^0 + ( G_400 S_000 | D_200 S_000 )^0 - ( I_420 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[11] = etfac[1] * AUX_INT__h_s_d_s[6] + 1 * one_over_2q * AUX_INT__g_s_d_s[0] - p_over_q * AUX_INT__i_s_d_s[18];

                    // ( H_410 S_000 | F_102 S_000 )^0 = x * ( H_410 S_000 | D_002 S_000 )^0 + ( G_310 S_000 | D_002 S_000 )^0 - ( I_510 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[15] = etfac[0] * AUX_INT__h_s_d_s[11] + 4 * one_over_2q * AUX_INT__g_s_d_s[11] - p_over_q * AUX_INT__i_s_d_s[11];

                    // ( H_410 S_000 | F_030 S_000 )^0 = y * ( H_410 S_000 | D_020 S_000 )^0 + ( G_400 S_000 | D_020 S_000 )^0 + ( H_410 S_000 | P_010 S_000 )^0 - ( I_420 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[16] = etfac[1] * AUX_INT__h_s_d_s[9] + 1 * one_over_2q * AUX_INT__g_s_d_s[3] + 2 * one_over_2q * AUX_INT__h_s_p_s[4] - p_over_q * AUX_INT__i_s_d_s[21];

                    // ( H_410 S_000 | F_003 S_000 )^0 = z * ( H_410 S_000 | D_002 S_000 )^0 + ( H_410 S_000 | P_001 S_000 )^0 - ( I_411 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[19] = etfac[2] * AUX_INT__h_s_d_s[11] + 2 * one_over_2q * AUX_INT__h_s_p_s[5] - p_over_q * AUX_INT__i_s_d_s[29];

                    // ( H_401 S_000 | F_300 S_000 )^0 = x * ( H_401 S_000 | D_200 S_000 )^0 + ( G_301 S_000 | D_200 S_000 )^0 + ( H_401 S_000 | P_100 S_000 )^0 - ( I_501 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[20] = etfac[0] * AUX_INT__h_s_d_s[12] + 4 * one_over_2q * AUX_INT__g_s_d_s[12] + 2 * one_over_2q * AUX_INT__h_s_p_s[6] - p_over_q * AUX_INT__i_s_d_s[12];

                    // ( H_401 S_000 | F_210 S_000 )^0 = y * ( H_401 S_000 | D_200 S_000 )^0 - ( I_411 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[21] = etfac[1] * AUX_INT__h_s_d_s[12] - p_over_q * AUX_INT__i_s_d_s[24];

                    // ( H_401 S_000 | F_201 S_000 )^0 = z * ( H_401 S_000 | D_200 S_000 )^0 + ( G_400 S_000 | D_200 S_000 )^0 - ( I_402 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[22] = etfac[2] * AUX_INT__h_s_d_s[12] + 1 * one_over_2q * AUX_INT__g_s_d_s[0] - p_over_q * AUX_INT__i_s_d_s[30];

                    // ( H_401 S_000 | F_120 S_000 )^0 = x * ( H_401 S_000 | D_020 S_000 )^0 + ( G_301 S_000 | D_020 S_000 )^0 - ( I_501 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[23] = etfac[0] * AUX_INT__h_s_d_s[15] + 4 * one_over_2q * AUX_INT__g_s_d_s[15] - p_over_q * AUX_INT__i_s_d_s[15];

                    // ( H_401 S_000 | F_030 S_000 )^0 = y * ( H_401 S_000 | D_020 S_000 )^0 + ( H_401 S_000 | P_010 S_000 )^0 - ( I_411 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[26] = etfac[1] * AUX_INT__h_s_d_s[15] + 2 * one_over_2q * AUX_INT__h_s_p_s[7] - p_over_q * AUX_INT__i_s_d_s[27];

                    // ( H_401 S_000 | F_021 S_000 )^0 = z * ( H_401 S_000 | D_020 S_000 )^0 + ( G_400 S_000 | D_020 S_000 )^0 - ( I_402 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[27] = etfac[2] * AUX_INT__h_s_d_s[15] + 1 * one_over_2q * AUX_INT__g_s_d_s[3] - p_over_q * AUX_INT__i_s_d_s[33];

                    // ( H_401 S_000 | F_003 S_000 )^0 = z * ( H_401 S_000 | D_002 S_000 )^0 + ( G_400 S_000 | D_002 S_000 )^0 + ( H_401 S_000 | P_001 S_000 )^0 - ( I_402 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[29] = etfac[2] * AUX_INT__h_s_d_s[17] + 1 * one_over_2q * AUX_INT__g_s_d_s[5] + 2 * one_over_2q * AUX_INT__h_s_p_s[8] - p_over_q * AUX_INT__i_s_d_s[35];

                    // ( H_320 S_000 | F_300 S_000 )^0 = x * ( H_320 S_000 | D_200 S_000 )^0 + ( G_220 S_000 | D_200 S_000 )^0 + ( H_320 S_000 | P_100 S_000 )^0 - ( I_420 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[30] = etfac[0] * AUX_INT__h_s_d_s[18] + 3 * one_over_2q * AUX_INT__g_s_d_s[18] + 2 * one_over_2q * AUX_INT__h_s_p_s[9] - p_over_q * AUX_INT__i_s_d_s[18];

                    // ( H_320 S_000 | F_210 S_000 )^0 = y * ( H_320 S_000 | D_200 S_000 )^0 + ( G_310 S_000 | D_200 S_000 )^0 - ( I_330 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[31] = etfac[1] * AUX_INT__h_s_d_s[18] + 2 * one_over_2q * AUX_INT__g_s_d_s[6] - p_over_q * AUX_INT__i_s_d_s[36];

                    // ( H_320 S_000 | F_102 S_000 )^0 = x * ( H_320 S_000 | D_002 S_000 )^0 + ( G_220 S_000 | D_002 S_000 )^0 - ( I_420 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[35] = etfac[0] * AUX_INT__h_s_d_s[23] + 3 * one_over_2q * AUX_INT__g_s_d_s[23] - p_over_q * AUX_INT__i_s_d_s[23];

                    // ( H_320 S_000 | F_030 S_000 )^0 = y * ( H_320 S_000 | D_020 S_000 )^0 + ( G_310 S_000 | D_020 S_000 )^0 + ( H_320 S_000 | P_010 S_000 )^0 - ( I_330 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[36] = etfac[1] * AUX_INT__h_s_d_s[21] + 2 * one_over_2q * AUX_INT__g_s_d_s[9] + 2 * one_over_2q * AUX_INT__h_s_p_s[10] - p_over_q * AUX_INT__i_s_d_s[39];

                    // ( H_320 S_000 | F_003 S_000 )^0 = z * ( H_320 S_000 | D_002 S_000 )^0 + ( H_320 S_000 | P_001 S_000 )^0 - ( I_321 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[39] = etfac[2] * AUX_INT__h_s_d_s[23] + 2 * one_over_2q * AUX_INT__h_s_p_s[11] - p_over_q * AUX_INT__i_s_d_s[47];

                    // ( H_311 S_000 | F_300 S_000 )^0 = x * ( H_311 S_000 | D_200 S_000 )^0 + ( G_211 S_000 | D_200 S_000 )^0 + ( H_311 S_000 | P_100 S_000 )^0 - ( I_411 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[40] = etfac[0] * AUX_INT__h_s_d_s[24] + 3 * one_over_2q * AUX_INT__g_s_d_s[24] + 2 * one_over_2q * AUX_INT__h_s_p_s[12] - p_over_q * AUX_INT__i_s_d_s[24];

                    // ( H_311 S_000 | F_210 S_000 )^0 = y * ( H_311 S_000 | D_200 S_000 )^0 + ( G_301 S_000 | D_200 S_000 )^0 - ( I_321 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[41] = etfac[1] * AUX_INT__h_s_d_s[24] + 1 * one_over_2q * AUX_INT__g_s_d_s[12] - p_over_q * AUX_INT__i_s_d_s[42];

                    // ( H_311 S_000 | F_201 S_000 )^0 = z * ( H_311 S_000 | D_200 S_000 )^0 + ( G_310 S_000 | D_200 S_000 )^0 - ( I_312 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[42] = etfac[2] * AUX_INT__h_s_d_s[24] + 1 * one_over_2q * AUX_INT__g_s_d_s[6] - p_over_q * AUX_INT__i_s_d_s[48];

                    // ( H_311 S_000 | F_120 S_000 )^0 = x * ( H_311 S_000 | D_020 S_000 )^0 + ( G_211 S_000 | D_020 S_000 )^0 - ( I_411 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[43] = etfac[0] * AUX_INT__h_s_d_s[27] + 3 * one_over_2q * AUX_INT__g_s_d_s[27] - p_over_q * AUX_INT__i_s_d_s[27];

                    // ( H_311 S_000 | F_102 S_000 )^0 = x * ( H_311 S_000 | D_002 S_000 )^0 + ( G_211 S_000 | D_002 S_000 )^0 - ( I_411 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[45] = etfac[0] * AUX_INT__h_s_d_s[29] + 3 * one_over_2q * AUX_INT__g_s_d_s[29] - p_over_q * AUX_INT__i_s_d_s[29];

                    // ( H_311 S_000 | F_030 S_000 )^0 = y * ( H_311 S_000 | D_020 S_000 )^0 + ( G_301 S_000 | D_020 S_000 )^0 + ( H_311 S_000 | P_010 S_000 )^0 - ( I_321 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[46] = etfac[1] * AUX_INT__h_s_d_s[27] + 1 * one_over_2q * AUX_INT__g_s_d_s[15] + 2 * one_over_2q * AUX_INT__h_s_p_s[13] - p_over_q * AUX_INT__i_s_d_s[45];

                    // ( H_311 S_000 | F_021 S_000 )^0 = z * ( H_311 S_000 | D_020 S_000 )^0 + ( G_310 S_000 | D_020 S_000 )^0 - ( I_312 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[47] = etfac[2] * AUX_INT__h_s_d_s[27] + 1 * one_over_2q * AUX_INT__g_s_d_s[9] - p_over_q * AUX_INT__i_s_d_s[51];

                    // ( H_311 S_000 | F_003 S_000 )^0 = z * ( H_311 S_000 | D_002 S_000 )^0 + ( G_310 S_000 | D_002 S_000 )^0 + ( H_311 S_000 | P_001 S_000 )^0 - ( I_312 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[49] = etfac[2] * AUX_INT__h_s_d_s[29] + 1 * one_over_2q * AUX_INT__g_s_d_s[11] + 2 * one_over_2q * AUX_INT__h_s_p_s[14] - p_over_q * AUX_INT__i_s_d_s[53];

                    // ( H_302 S_000 | F_300 S_000 )^0 = x * ( H_302 S_000 | D_200 S_000 )^0 + ( G_202 S_000 | D_200 S_000 )^0 + ( H_302 S_000 | P_100 S_000 )^0 - ( I_402 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[50] = etfac[0] * AUX_INT__h_s_d_s[30] + 3 * one_over_2q * AUX_INT__g_s_d_s[30] + 2 * one_over_2q * AUX_INT__h_s_p_s[15] - p_over_q * AUX_INT__i_s_d_s[30];

                    // ( H_302 S_000 | F_210 S_000 )^0 = y * ( H_302 S_000 | D_200 S_000 )^0 - ( I_312 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[51] = etfac[1] * AUX_INT__h_s_d_s[30] - p_over_q * AUX_INT__i_s_d_s[48];

                    // ( H_302 S_000 | F_201 S_000 )^0 = z * ( H_302 S_000 | D_200 S_000 )^0 + ( G_301 S_000 | D_200 S_000 )^0 - ( I_303 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[52] = etfac[2] * AUX_INT__h_s_d_s[30] + 2 * one_over_2q * AUX_INT__g_s_d_s[12] - p_over_q * AUX_INT__i_s_d_s[54];

                    // ( H_302 S_000 | F_120 S_000 )^0 = x * ( H_302 S_000 | D_020 S_000 )^0 + ( G_202 S_000 | D_020 S_000 )^0 - ( I_402 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[53] = etfac[0] * AUX_INT__h_s_d_s[33] + 3 * one_over_2q * AUX_INT__g_s_d_s[33] - p_over_q * AUX_INT__i_s_d_s[33];

                    // ( H_302 S_000 | F_030 S_000 )^0 = y * ( H_302 S_000 | D_020 S_000 )^0 + ( H_302 S_000 | P_010 S_000 )^0 - ( I_312 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[56] = etfac[1] * AUX_INT__h_s_d_s[33] + 2 * one_over_2q * AUX_INT__h_s_p_s[16] - p_over_q * AUX_INT__i_s_d_s[51];

                    // ( H_302 S_000 | F_021 S_000 )^0 = z * ( H_302 S_000 | D_020 S_000 )^0 + ( G_301 S_000 | D_020 S_000 )^0 - ( I_303 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[57] = etfac[2] * AUX_INT__h_s_d_s[33] + 2 * one_over_2q * AUX_INT__g_s_d_s[15] - p_over_q * AUX_INT__i_s_d_s[57];

                    // ( H_302 S_000 | F_003 S_000 )^0 = z * ( H_302 S_000 | D_002 S_000 )^0 + ( G_301 S_000 | D_002 S_000 )^0 + ( H_302 S_000 | P_001 S_000 )^0 - ( I_303 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[59] = etfac[2] * AUX_INT__h_s_d_s[35] + 2 * one_over_2q * AUX_INT__g_s_d_s[17] + 2 * one_over_2q * AUX_INT__h_s_p_s[17] - p_over_q * AUX_INT__i_s_d_s[59];

                    // ( H_230 S_000 | F_300 S_000 )^0 = x * ( H_230 S_000 | D_200 S_000 )^0 + ( G_130 S_000 | D_200 S_000 )^0 + ( H_230 S_000 | P_100 S_000 )^0 - ( I_330 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[60] = etfac[0] * AUX_INT__h_s_d_s[36] + 2 * one_over_2q * AUX_INT__g_s_d_s[36] + 2 * one_over_2q * AUX_INT__h_s_p_s[18] - p_over_q * AUX_INT__i_s_d_s[36];

                    // ( H_230 S_000 | F_210 S_000 )^0 = y * ( H_230 S_000 | D_200 S_000 )^0 + ( G_220 S_000 | D_200 S_000 )^0 - ( I_240 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[61] = etfac[1] * AUX_INT__h_s_d_s[36] + 3 * one_over_2q * AUX_INT__g_s_d_s[18] - p_over_q * AUX_INT__i_s_d_s[60];

                    // ( H_230 S_000 | F_102 S_000 )^0 = x * ( H_230 S_000 | D_002 S_000 )^0 + ( G_130 S_000 | D_002 S_000 )^0 - ( I_330 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[65] = etfac[0] * AUX_INT__h_s_d_s[41] + 2 * one_over_2q * AUX_INT__g_s_d_s[41] - p_over_q * AUX_INT__i_s_d_s[41];

                    // ( H_230 S_000 | F_030 S_000 )^0 = y * ( H_230 S_000 | D_020 S_000 )^0 + ( G_220 S_000 | D_020 S_000 )^0 + ( H_230 S_000 | P_010 S_000 )^0 - ( I_240 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[66] = etfac[1] * AUX_INT__h_s_d_s[39] + 3 * one_over_2q * AUX_INT__g_s_d_s[21] + 2 * one_over_2q * AUX_INT__h_s_p_s[19] - p_over_q * AUX_INT__i_s_d_s[63];

                    // ( H_230 S_000 | F_003 S_000 )^0 = z * ( H_230 S_000 | D_002 S_000 )^0 + ( H_230 S_000 | P_001 S_000 )^0 - ( I_231 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[69] = etfac[2] * AUX_INT__h_s_d_s[41] + 2 * one_over_2q * AUX_INT__h_s_p_s[20] - p_over_q * AUX_INT__i_s_d_s[71];

                    // ( H_221 S_000 | F_300 S_000 )^0 = x * ( H_221 S_000 | D_200 S_000 )^0 + ( G_121 S_000 | D_200 S_000 )^0 + ( H_221 S_000 | P_100 S_000 )^0 - ( I_321 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[70] = etfac[0] * AUX_INT__h_s_d_s[42] + 2 * one_over_2q * AUX_INT__g_s_d_s[42] + 2 * one_over_2q * AUX_INT__h_s_p_s[21] - p_over_q * AUX_INT__i_s_d_s[42];

                    // ( H_221 S_000 | F_210 S_000 )^0 = y * ( H_221 S_000 | D_200 S_000 )^0 + ( G_211 S_000 | D_200 S_000 )^0 - ( I_231 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[71] = etfac[1] * AUX_INT__h_s_d_s[42] + 2 * one_over_2q * AUX_INT__g_s_d_s[24] - p_over_q * AUX_INT__i_s_d_s[66];

                    // ( H_221 S_000 | F_201 S_000 )^0 = z * ( H_221 S_000 | D_200 S_000 )^0 + ( G_220 S_000 | D_200 S_000 )^0 - ( I_222 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[72] = etfac[2] * AUX_INT__h_s_d_s[42] + 1 * one_over_2q * AUX_INT__g_s_d_s[18] - p_over_q * AUX_INT__i_s_d_s[72];

                    // ( H_221 S_000 | F_120 S_000 )^0 = x * ( H_221 S_000 | D_020 S_000 )^0 + ( G_121 S_000 | D_020 S_000 )^0 - ( I_321 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[73] = etfac[0] * AUX_INT__h_s_d_s[45] + 2 * one_over_2q * AUX_INT__g_s_d_s[45] - p_over_q * AUX_INT__i_s_d_s[45];

                    // ( H_221 S_000 | F_102 S_000 )^0 = x * ( H_221 S_000 | D_002 S_000 )^0 + ( G_121 S_000 | D_002 S_000 )^0 - ( I_321 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[75] = etfac[0] * AUX_INT__h_s_d_s[47] + 2 * one_over_2q * AUX_INT__g_s_d_s[47] - p_over_q * AUX_INT__i_s_d_s[47];

                    // ( H_221 S_000 | F_030 S_000 )^0 = y * ( H_221 S_000 | D_020 S_000 )^0 + ( G_211 S_000 | D_020 S_000 )^0 + ( H_221 S_000 | P_010 S_000 )^0 - ( I_231 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[76] = etfac[1] * AUX_INT__h_s_d_s[45] + 2 * one_over_2q * AUX_INT__g_s_d_s[27] + 2 * one_over_2q * AUX_INT__h_s_p_s[22] - p_over_q * AUX_INT__i_s_d_s[69];

                    // ( H_221 S_000 | F_021 S_000 )^0 = z * ( H_221 S_000 | D_020 S_000 )^0 + ( G_220 S_000 | D_020 S_000 )^0 - ( I_222 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[77] = etfac[2] * AUX_INT__h_s_d_s[45] + 1 * one_over_2q * AUX_INT__g_s_d_s[21] - p_over_q * AUX_INT__i_s_d_s[75];

                    // ( H_221 S_000 | F_003 S_000 )^0 = z * ( H_221 S_000 | D_002 S_000 )^0 + ( G_220 S_000 | D_002 S_000 )^0 + ( H_221 S_000 | P_001 S_000 )^0 - ( I_222 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[79] = etfac[2] * AUX_INT__h_s_d_s[47] + 1 * one_over_2q * AUX_INT__g_s_d_s[23] + 2 * one_over_2q * AUX_INT__h_s_p_s[23] - p_over_q * AUX_INT__i_s_d_s[77];

                    // ( H_212 S_000 | F_300 S_000 )^0 = x * ( H_212 S_000 | D_200 S_000 )^0 + ( G_112 S_000 | D_200 S_000 )^0 + ( H_212 S_000 | P_100 S_000 )^0 - ( I_312 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[80] = etfac[0] * AUX_INT__h_s_d_s[48] + 2 * one_over_2q * AUX_INT__g_s_d_s[48] + 2 * one_over_2q * AUX_INT__h_s_p_s[24] - p_over_q * AUX_INT__i_s_d_s[48];

                    // ( H_212 S_000 | F_210 S_000 )^0 = y * ( H_212 S_000 | D_200 S_000 )^0 + ( G_202 S_000 | D_200 S_000 )^0 - ( I_222 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[81] = etfac[1] * AUX_INT__h_s_d_s[48] + 1 * one_over_2q * AUX_INT__g_s_d_s[30] - p_over_q * AUX_INT__i_s_d_s[72];

                    // ( H_212 S_000 | F_201 S_000 )^0 = z * ( H_212 S_000 | D_200 S_000 )^0 + ( G_211 S_000 | D_200 S_000 )^0 - ( I_213 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[82] = etfac[2] * AUX_INT__h_s_d_s[48] + 2 * one_over_2q * AUX_INT__g_s_d_s[24] - p_over_q * AUX_INT__i_s_d_s[78];

                    // ( H_212 S_000 | F_120 S_000 )^0 = x * ( H_212 S_000 | D_020 S_000 )^0 + ( G_112 S_000 | D_020 S_000 )^0 - ( I_312 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[83] = etfac[0] * AUX_INT__h_s_d_s[51] + 2 * one_over_2q * AUX_INT__g_s_d_s[51] - p_over_q * AUX_INT__i_s_d_s[51];

                    // ( H_212 S_000 | F_102 S_000 )^0 = x * ( H_212 S_000 | D_002 S_000 )^0 + ( G_112 S_000 | D_002 S_000 )^0 - ( I_312 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[85] = etfac[0] * AUX_INT__h_s_d_s[53] + 2 * one_over_2q * AUX_INT__g_s_d_s[53] - p_over_q * AUX_INT__i_s_d_s[53];

                    // ( H_212 S_000 | F_030 S_000 )^0 = y * ( H_212 S_000 | D_020 S_000 )^0 + ( G_202 S_000 | D_020 S_000 )^0 + ( H_212 S_000 | P_010 S_000 )^0 - ( I_222 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[86] = etfac[1] * AUX_INT__h_s_d_s[51] + 1 * one_over_2q * AUX_INT__g_s_d_s[33] + 2 * one_over_2q * AUX_INT__h_s_p_s[25] - p_over_q * AUX_INT__i_s_d_s[75];

                    // ( H_212 S_000 | F_021 S_000 )^0 = z * ( H_212 S_000 | D_020 S_000 )^0 + ( G_211 S_000 | D_020 S_000 )^0 - ( I_213 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[87] = etfac[2] * AUX_INT__h_s_d_s[51] + 2 * one_over_2q * AUX_INT__g_s_d_s[27] - p_over_q * AUX_INT__i_s_d_s[81];

                    // ( H_212 S_000 | F_003 S_000 )^0 = z * ( H_212 S_000 | D_002 S_000 )^0 + ( G_211 S_000 | D_002 S_000 )^0 + ( H_212 S_000 | P_001 S_000 )^0 - ( I_213 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[89] = etfac[2] * AUX_INT__h_s_d_s[53] + 2 * one_over_2q * AUX_INT__g_s_d_s[29] + 2 * one_over_2q * AUX_INT__h_s_p_s[26] - p_over_q * AUX_INT__i_s_d_s[83];

                    // ( H_203 S_000 | F_300 S_000 )^0 = x * ( H_203 S_000 | D_200 S_000 )^0 + ( G_103 S_000 | D_200 S_000 )^0 + ( H_203 S_000 | P_100 S_000 )^0 - ( I_303 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[90] = etfac[0] * AUX_INT__h_s_d_s[54] + 2 * one_over_2q * AUX_INT__g_s_d_s[54] + 2 * one_over_2q * AUX_INT__h_s_p_s[27] - p_over_q * AUX_INT__i_s_d_s[54];

                    // ( H_203 S_000 | F_210 S_000 )^0 = y * ( H_203 S_000 | D_200 S_000 )^0 - ( I_213 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[91] = etfac[1] * AUX_INT__h_s_d_s[54] - p_over_q * AUX_INT__i_s_d_s[78];

                    // ( H_203 S_000 | F_201 S_000 )^0 = z * ( H_203 S_000 | D_200 S_000 )^0 + ( G_202 S_000 | D_200 S_000 )^0 - ( I_204 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[92] = etfac[2] * AUX_INT__h_s_d_s[54] + 3 * one_over_2q * AUX_INT__g_s_d_s[30] - p_over_q * AUX_INT__i_s_d_s[84];

                    // ( H_203 S_000 | F_120 S_000 )^0 = x * ( H_203 S_000 | D_020 S_000 )^0 + ( G_103 S_000 | D_020 S_000 )^0 - ( I_303 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[93] = etfac[0] * AUX_INT__h_s_d_s[57] + 2 * one_over_2q * AUX_INT__g_s_d_s[57] - p_over_q * AUX_INT__i_s_d_s[57];

                    // ( H_203 S_000 | F_030 S_000 )^0 = y * ( H_203 S_000 | D_020 S_000 )^0 + ( H_203 S_000 | P_010 S_000 )^0 - ( I_213 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[96] = etfac[1] * AUX_INT__h_s_d_s[57] + 2 * one_over_2q * AUX_INT__h_s_p_s[28] - p_over_q * AUX_INT__i_s_d_s[81];

                    // ( H_203 S_000 | F_021 S_000 )^0 = z * ( H_203 S_000 | D_020 S_000 )^0 + ( G_202 S_000 | D_020 S_000 )^0 - ( I_204 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[97] = etfac[2] * AUX_INT__h_s_d_s[57] + 3 * one_over_2q * AUX_INT__g_s_d_s[33] - p_over_q * AUX_INT__i_s_d_s[87];

                    // ( H_203 S_000 | F_003 S_000 )^0 = z * ( H_203 S_000 | D_002 S_000 )^0 + ( G_202 S_000 | D_002 S_000 )^0 + ( H_203 S_000 | P_001 S_000 )^0 - ( I_204 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[99] = etfac[2] * AUX_INT__h_s_d_s[59] + 3 * one_over_2q * AUX_INT__g_s_d_s[35] + 2 * one_over_2q * AUX_INT__h_s_p_s[29] - p_over_q * AUX_INT__i_s_d_s[89];

                    // ( H_140 S_000 | F_300 S_000 )^0 = x * ( H_140 S_000 | D_200 S_000 )^0 + ( G_040 S_000 | D_200 S_000 )^0 + ( H_140 S_000 | P_100 S_000 )^0 - ( I_240 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[100] = etfac[0] * AUX_INT__h_s_d_s[60] + 1 * one_over_2q * AUX_INT__g_s_d_s[60] + 2 * one_over_2q * AUX_INT__h_s_p_s[30] - p_over_q * AUX_INT__i_s_d_s[60];

                    // ( H_140 S_000 | F_210 S_000 )^0 = y * ( H_140 S_000 | D_200 S_000 )^0 + ( G_130 S_000 | D_200 S_000 )^0 - ( I_150 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[101] = etfac[1] * AUX_INT__h_s_d_s[60] + 4 * one_over_2q * AUX_INT__g_s_d_s[36] - p_over_q * AUX_INT__i_s_d_s[90];

                    // ( H_140 S_000 | F_102 S_000 )^0 = x * ( H_140 S_000 | D_002 S_000 )^0 + ( G_040 S_000 | D_002 S_000 )^0 - ( I_240 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[105] = etfac[0] * AUX_INT__h_s_d_s[65] + 1 * one_over_2q * AUX_INT__g_s_d_s[65] - p_over_q * AUX_INT__i_s_d_s[65];

                    // ( H_140 S_000 | F_030 S_000 )^0 = y * ( H_140 S_000 | D_020 S_000 )^0 + ( G_130 S_000 | D_020 S_000 )^0 + ( H_140 S_000 | P_010 S_000 )^0 - ( I_150 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[106] = etfac[1] * AUX_INT__h_s_d_s[63] + 4 * one_over_2q * AUX_INT__g_s_d_s[39] + 2 * one_over_2q * AUX_INT__h_s_p_s[31] - p_over_q * AUX_INT__i_s_d_s[93];

                    // ( H_140 S_000 | F_003 S_000 )^0 = z * ( H_140 S_000 | D_002 S_000 )^0 + ( H_140 S_000 | P_001 S_000 )^0 - ( I_141 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[109] = etfac[2] * AUX_INT__h_s_d_s[65] + 2 * one_over_2q * AUX_INT__h_s_p_s[32] - p_over_q * AUX_INT__i_s_d_s[101];

                    // ( H_131 S_000 | F_300 S_000 )^0 = x * ( H_131 S_000 | D_200 S_000 )^0 + ( G_031 S_000 | D_200 S_000 )^0 + ( H_131 S_000 | P_100 S_000 )^0 - ( I_231 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[110] = etfac[0] * AUX_INT__h_s_d_s[66] + 1 * one_over_2q * AUX_INT__g_s_d_s[66] + 2 * one_over_2q * AUX_INT__h_s_p_s[33] - p_over_q * AUX_INT__i_s_d_s[66];

                    // ( H_131 S_000 | F_210 S_000 )^0 = y * ( H_131 S_000 | D_200 S_000 )^0 + ( G_121 S_000 | D_200 S_000 )^0 - ( I_141 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[111] = etfac[1] * AUX_INT__h_s_d_s[66] + 3 * one_over_2q * AUX_INT__g_s_d_s[42] - p_over_q * AUX_INT__i_s_d_s[96];

                    // ( H_131 S_000 | F_201 S_000 )^0 = z * ( H_131 S_000 | D_200 S_000 )^0 + ( G_130 S_000 | D_200 S_000 )^0 - ( I_132 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[112] = etfac[2] * AUX_INT__h_s_d_s[66] + 1 * one_over_2q * AUX_INT__g_s_d_s[36] - p_over_q * AUX_INT__i_s_d_s[102];

                    // ( H_131 S_000 | F_120 S_000 )^0 = x * ( H_131 S_000 | D_020 S_000 )^0 + ( G_031 S_000 | D_020 S_000 )^0 - ( I_231 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[113] = etfac[0] * AUX_INT__h_s_d_s[69] + 1 * one_over_2q * AUX_INT__g_s_d_s[69] - p_over_q * AUX_INT__i_s_d_s[69];

                    // ( H_131 S_000 | F_102 S_000 )^0 = x * ( H_131 S_000 | D_002 S_000 )^0 + ( G_031 S_000 | D_002 S_000 )^0 - ( I_231 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[115] = etfac[0] * AUX_INT__h_s_d_s[71] + 1 * one_over_2q * AUX_INT__g_s_d_s[71] - p_over_q * AUX_INT__i_s_d_s[71];

                    // ( H_131 S_000 | F_030 S_000 )^0 = y * ( H_131 S_000 | D_020 S_000 )^0 + ( G_121 S_000 | D_020 S_000 )^0 + ( H_131 S_000 | P_010 S_000 )^0 - ( I_141 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[116] = etfac[1] * AUX_INT__h_s_d_s[69] + 3 * one_over_2q * AUX_INT__g_s_d_s[45] + 2 * one_over_2q * AUX_INT__h_s_p_s[34] - p_over_q * AUX_INT__i_s_d_s[99];

                    // ( H_131 S_000 | F_021 S_000 )^0 = z * ( H_131 S_000 | D_020 S_000 )^0 + ( G_130 S_000 | D_020 S_000 )^0 - ( I_132 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[117] = etfac[2] * AUX_INT__h_s_d_s[69] + 1 * one_over_2q * AUX_INT__g_s_d_s[39] - p_over_q * AUX_INT__i_s_d_s[105];

                    // ( H_131 S_000 | F_003 S_000 )^0 = z * ( H_131 S_000 | D_002 S_000 )^0 + ( G_130 S_000 | D_002 S_000 )^0 + ( H_131 S_000 | P_001 S_000 )^0 - ( I_132 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[119] = etfac[2] * AUX_INT__h_s_d_s[71] + 1 * one_over_2q * AUX_INT__g_s_d_s[41] + 2 * one_over_2q * AUX_INT__h_s_p_s[35] - p_over_q * AUX_INT__i_s_d_s[107];

                    // ( H_122 S_000 | F_300 S_000 )^0 = x * ( H_122 S_000 | D_200 S_000 )^0 + ( G_022 S_000 | D_200 S_000 )^0 + ( H_122 S_000 | P_100 S_000 )^0 - ( I_222 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[120] = etfac[0] * AUX_INT__h_s_d_s[72] + 1 * one_over_2q * AUX_INT__g_s_d_s[72] + 2 * one_over_2q * AUX_INT__h_s_p_s[36] - p_over_q * AUX_INT__i_s_d_s[72];

                    // ( H_122 S_000 | F_210 S_000 )^0 = y * ( H_122 S_000 | D_200 S_000 )^0 + ( G_112 S_000 | D_200 S_000 )^0 - ( I_132 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[121] = etfac[1] * AUX_INT__h_s_d_s[72] + 2 * one_over_2q * AUX_INT__g_s_d_s[48] - p_over_q * AUX_INT__i_s_d_s[102];

                    // ( H_122 S_000 | F_201 S_000 )^0 = z * ( H_122 S_000 | D_200 S_000 )^0 + ( G_121 S_000 | D_200 S_000 )^0 - ( I_123 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[122] = etfac[2] * AUX_INT__h_s_d_s[72] + 2 * one_over_2q * AUX_INT__g_s_d_s[42] - p_over_q * AUX_INT__i_s_d_s[108];

                    // ( H_122 S_000 | F_120 S_000 )^0 = x * ( H_122 S_000 | D_020 S_000 )^0 + ( G_022 S_000 | D_020 S_000 )^0 - ( I_222 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[123] = etfac[0] * AUX_INT__h_s_d_s[75] + 1 * one_over_2q * AUX_INT__g_s_d_s[75] - p_over_q * AUX_INT__i_s_d_s[75];

                    // ( H_122 S_000 | F_102 S_000 )^0 = x * ( H_122 S_000 | D_002 S_000 )^0 + ( G_022 S_000 | D_002 S_000 )^0 - ( I_222 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[125] = etfac[0] * AUX_INT__h_s_d_s[77] + 1 * one_over_2q * AUX_INT__g_s_d_s[77] - p_over_q * AUX_INT__i_s_d_s[77];

                    // ( H_122 S_000 | F_030 S_000 )^0 = y * ( H_122 S_000 | D_020 S_000 )^0 + ( G_112 S_000 | D_020 S_000 )^0 + ( H_122 S_000 | P_010 S_000 )^0 - ( I_132 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[126] = etfac[1] * AUX_INT__h_s_d_s[75] + 2 * one_over_2q * AUX_INT__g_s_d_s[51] + 2 * one_over_2q * AUX_INT__h_s_p_s[37] - p_over_q * AUX_INT__i_s_d_s[105];

                    // ( H_122 S_000 | F_021 S_000 )^0 = z * ( H_122 S_000 | D_020 S_000 )^0 + ( G_121 S_000 | D_020 S_000 )^0 - ( I_123 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[127] = etfac[2] * AUX_INT__h_s_d_s[75] + 2 * one_over_2q * AUX_INT__g_s_d_s[45] - p_over_q * AUX_INT__i_s_d_s[111];

                    // ( H_122 S_000 | F_003 S_000 )^0 = z * ( H_122 S_000 | D_002 S_000 )^0 + ( G_121 S_000 | D_002 S_000 )^0 + ( H_122 S_000 | P_001 S_000 )^0 - ( I_123 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[129] = etfac[2] * AUX_INT__h_s_d_s[77] + 2 * one_over_2q * AUX_INT__g_s_d_s[47] + 2 * one_over_2q * AUX_INT__h_s_p_s[38] - p_over_q * AUX_INT__i_s_d_s[113];

                    // ( H_113 S_000 | F_300 S_000 )^0 = x * ( H_113 S_000 | D_200 S_000 )^0 + ( G_013 S_000 | D_200 S_000 )^0 + ( H_113 S_000 | P_100 S_000 )^0 - ( I_213 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[130] = etfac[0] * AUX_INT__h_s_d_s[78] + 1 * one_over_2q * AUX_INT__g_s_d_s[78] + 2 * one_over_2q * AUX_INT__h_s_p_s[39] - p_over_q * AUX_INT__i_s_d_s[78];

                    // ( H_113 S_000 | F_210 S_000 )^0 = y * ( H_113 S_000 | D_200 S_000 )^0 + ( G_103 S_000 | D_200 S_000 )^0 - ( I_123 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[131] = etfac[1] * AUX_INT__h_s_d_s[78] + 1 * one_over_2q * AUX_INT__g_s_d_s[54] - p_over_q * AUX_INT__i_s_d_s[108];

                    // ( H_113 S_000 | F_201 S_000 )^0 = z * ( H_113 S_000 | D_200 S_000 )^0 + ( G_112 S_000 | D_200 S_000 )^0 - ( I_114 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[132] = etfac[2] * AUX_INT__h_s_d_s[78] + 3 * one_over_2q * AUX_INT__g_s_d_s[48] - p_over_q * AUX_INT__i_s_d_s[114];

                    // ( H_113 S_000 | F_120 S_000 )^0 = x * ( H_113 S_000 | D_020 S_000 )^0 + ( G_013 S_000 | D_020 S_000 )^0 - ( I_213 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[133] = etfac[0] * AUX_INT__h_s_d_s[81] + 1 * one_over_2q * AUX_INT__g_s_d_s[81] - p_over_q * AUX_INT__i_s_d_s[81];

                    // ( H_113 S_000 | F_102 S_000 )^0 = x * ( H_113 S_000 | D_002 S_000 )^0 + ( G_013 S_000 | D_002 S_000 )^0 - ( I_213 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[135] = etfac[0] * AUX_INT__h_s_d_s[83] + 1 * one_over_2q * AUX_INT__g_s_d_s[83] - p_over_q * AUX_INT__i_s_d_s[83];

                    // ( H_113 S_000 | F_030 S_000 )^0 = y * ( H_113 S_000 | D_020 S_000 )^0 + ( G_103 S_000 | D_020 S_000 )^0 + ( H_113 S_000 | P_010 S_000 )^0 - ( I_123 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[136] = etfac[1] * AUX_INT__h_s_d_s[81] + 1 * one_over_2q * AUX_INT__g_s_d_s[57] + 2 * one_over_2q * AUX_INT__h_s_p_s[40] - p_over_q * AUX_INT__i_s_d_s[111];

                    // ( H_113 S_000 | F_021 S_000 )^0 = z * ( H_113 S_000 | D_020 S_000 )^0 + ( G_112 S_000 | D_020 S_000 )^0 - ( I_114 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[137] = etfac[2] * AUX_INT__h_s_d_s[81] + 3 * one_over_2q * AUX_INT__g_s_d_s[51] - p_over_q * AUX_INT__i_s_d_s[117];

                    // ( H_113 S_000 | F_003 S_000 )^0 = z * ( H_113 S_000 | D_002 S_000 )^0 + ( G_112 S_000 | D_002 S_000 )^0 + ( H_113 S_000 | P_001 S_000 )^0 - ( I_114 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[139] = etfac[2] * AUX_INT__h_s_d_s[83] + 3 * one_over_2q * AUX_INT__g_s_d_s[53] + 2 * one_over_2q * AUX_INT__h_s_p_s[41] - p_over_q * AUX_INT__i_s_d_s[119];

                    // ( H_104 S_000 | F_300 S_000 )^0 = x * ( H_104 S_000 | D_200 S_000 )^0 + ( G_004 S_000 | D_200 S_000 )^0 + ( H_104 S_000 | P_100 S_000 )^0 - ( I_204 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[140] = etfac[0] * AUX_INT__h_s_d_s[84] + 1 * one_over_2q * AUX_INT__g_s_d_s[84] + 2 * one_over_2q * AUX_INT__h_s_p_s[42] - p_over_q * AUX_INT__i_s_d_s[84];

                    // ( H_104 S_000 | F_210 S_000 )^0 = y * ( H_104 S_000 | D_200 S_000 )^0 - ( I_114 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[141] = etfac[1] * AUX_INT__h_s_d_s[84] - p_over_q * AUX_INT__i_s_d_s[114];

                    // ( H_104 S_000 | F_201 S_000 )^0 = z * ( H_104 S_000 | D_200 S_000 )^0 + ( G_103 S_000 | D_200 S_000 )^0 - ( I_105 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[142] = etfac[2] * AUX_INT__h_s_d_s[84] + 4 * one_over_2q * AUX_INT__g_s_d_s[54] - p_over_q * AUX_INT__i_s_d_s[120];

                    // ( H_104 S_000 | F_120 S_000 )^0 = x * ( H_104 S_000 | D_020 S_000 )^0 + ( G_004 S_000 | D_020 S_000 )^0 - ( I_204 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[143] = etfac[0] * AUX_INT__h_s_d_s[87] + 1 * one_over_2q * AUX_INT__g_s_d_s[87] - p_over_q * AUX_INT__i_s_d_s[87];

                    // ( H_104 S_000 | F_030 S_000 )^0 = y * ( H_104 S_000 | D_020 S_000 )^0 + ( H_104 S_000 | P_010 S_000 )^0 - ( I_114 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[146] = etfac[1] * AUX_INT__h_s_d_s[87] + 2 * one_over_2q * AUX_INT__h_s_p_s[43] - p_over_q * AUX_INT__i_s_d_s[117];

                    // ( H_104 S_000 | F_021 S_000 )^0 = z * ( H_104 S_000 | D_020 S_000 )^0 + ( G_103 S_000 | D_020 S_000 )^0 - ( I_105 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[147] = etfac[2] * AUX_INT__h_s_d_s[87] + 4 * one_over_2q * AUX_INT__g_s_d_s[57] - p_over_q * AUX_INT__i_s_d_s[123];

                    // ( H_104 S_000 | F_003 S_000 )^0 = z * ( H_104 S_000 | D_002 S_000 )^0 + ( G_103 S_000 | D_002 S_000 )^0 + ( H_104 S_000 | P_001 S_000 )^0 - ( I_105 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[149] = etfac[2] * AUX_INT__h_s_d_s[89] + 4 * one_over_2q * AUX_INT__g_s_d_s[59] + 2 * one_over_2q * AUX_INT__h_s_p_s[44] - p_over_q * AUX_INT__i_s_d_s[125];

                    // ( H_050 S_000 | F_300 S_000 )^0 = x * ( H_050 S_000 | D_200 S_000 )^0 + ( H_050 S_000 | P_100 S_000 )^0 - ( I_150 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[150] = etfac[0] * AUX_INT__h_s_d_s[90] + 2 * one_over_2q * AUX_INT__h_s_p_s[45] - p_over_q * AUX_INT__i_s_d_s[90];

                    // ( H_050 S_000 | F_210 S_000 )^0 = y * ( H_050 S_000 | D_200 S_000 )^0 + ( G_040 S_000 | D_200 S_000 )^0 - ( I_060 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[151] = etfac[1] * AUX_INT__h_s_d_s[90] + 5 * one_over_2q * AUX_INT__g_s_d_s[60] - p_over_q * AUX_INT__i_s_d_s[126];

                    // ( H_050 S_000 | F_102 S_000 )^0 = x * ( H_050 S_000 | D_002 S_000 )^0 - ( I_150 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[155] = etfac[0] * AUX_INT__h_s_d_s[95] - p_over_q * AUX_INT__i_s_d_s[95];

                    // ( H_050 S_000 | F_030 S_000 )^0 = y * ( H_050 S_000 | D_020 S_000 )^0 + ( G_040 S_000 | D_020 S_000 )^0 + ( H_050 S_000 | P_010 S_000 )^0 - ( I_060 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[156] = etfac[1] * AUX_INT__h_s_d_s[93] + 5 * one_over_2q * AUX_INT__g_s_d_s[63] + 2 * one_over_2q * AUX_INT__h_s_p_s[46] - p_over_q * AUX_INT__i_s_d_s[129];

                    // ( H_050 S_000 | F_003 S_000 )^0 = z * ( H_050 S_000 | D_002 S_000 )^0 + ( H_050 S_000 | P_001 S_000 )^0 - ( I_051 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[159] = etfac[2] * AUX_INT__h_s_d_s[95] + 2 * one_over_2q * AUX_INT__h_s_p_s[47] - p_over_q * AUX_INT__i_s_d_s[137];

                    // ( H_041 S_000 | F_300 S_000 )^0 = x * ( H_041 S_000 | D_200 S_000 )^0 + ( H_041 S_000 | P_100 S_000 )^0 - ( I_141 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[160] = etfac[0] * AUX_INT__h_s_d_s[96] + 2 * one_over_2q * AUX_INT__h_s_p_s[48] - p_over_q * AUX_INT__i_s_d_s[96];

                    // ( H_041 S_000 | F_210 S_000 )^0 = y * ( H_041 S_000 | D_200 S_000 )^0 + ( G_031 S_000 | D_200 S_000 )^0 - ( I_051 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[161] = etfac[1] * AUX_INT__h_s_d_s[96] + 4 * one_over_2q * AUX_INT__g_s_d_s[66] - p_over_q * AUX_INT__i_s_d_s[132];

                    // ( H_041 S_000 | F_201 S_000 )^0 = z * ( H_041 S_000 | D_200 S_000 )^0 + ( G_040 S_000 | D_200 S_000 )^0 - ( I_042 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[162] = etfac[2] * AUX_INT__h_s_d_s[96] + 1 * one_over_2q * AUX_INT__g_s_d_s[60] - p_over_q * AUX_INT__i_s_d_s[138];

                    // ( H_041 S_000 | F_120 S_000 )^0 = x * ( H_041 S_000 | D_020 S_000 )^0 - ( I_141 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[163] = etfac[0] * AUX_INT__h_s_d_s[99] - p_over_q * AUX_INT__i_s_d_s[99];

                    // ( H_041 S_000 | F_102 S_000 )^0 = x * ( H_041 S_000 | D_002 S_000 )^0 - ( I_141 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[165] = etfac[0] * AUX_INT__h_s_d_s[101] - p_over_q * AUX_INT__i_s_d_s[101];

                    // ( H_041 S_000 | F_030 S_000 )^0 = y * ( H_041 S_000 | D_020 S_000 )^0 + ( G_031 S_000 | D_020 S_000 )^0 + ( H_041 S_000 | P_010 S_000 )^0 - ( I_051 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[166] = etfac[1] * AUX_INT__h_s_d_s[99] + 4 * one_over_2q * AUX_INT__g_s_d_s[69] + 2 * one_over_2q * AUX_INT__h_s_p_s[49] - p_over_q * AUX_INT__i_s_d_s[135];

                    // ( H_041 S_000 | F_021 S_000 )^0 = z * ( H_041 S_000 | D_020 S_000 )^0 + ( G_040 S_000 | D_020 S_000 )^0 - ( I_042 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[167] = etfac[2] * AUX_INT__h_s_d_s[99] + 1 * one_over_2q * AUX_INT__g_s_d_s[63] - p_over_q * AUX_INT__i_s_d_s[141];

                    // ( H_041 S_000 | F_003 S_000 )^0 = z * ( H_041 S_000 | D_002 S_000 )^0 + ( G_040 S_000 | D_002 S_000 )^0 + ( H_041 S_000 | P_001 S_000 )^0 - ( I_042 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[169] = etfac[2] * AUX_INT__h_s_d_s[101] + 1 * one_over_2q * AUX_INT__g_s_d_s[65] + 2 * one_over_2q * AUX_INT__h_s_p_s[50] - p_over_q * AUX_INT__i_s_d_s[143];

                    // ( H_032 S_000 | F_300 S_000 )^0 = x * ( H_032 S_000 | D_200 S_000 )^0 + ( H_032 S_000 | P_100 S_000 )^0 - ( I_132 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[170] = etfac[0] * AUX_INT__h_s_d_s[102] + 2 * one_over_2q * AUX_INT__h_s_p_s[51] - p_over_q * AUX_INT__i_s_d_s[102];

                    // ( H_032 S_000 | F_210 S_000 )^0 = y * ( H_032 S_000 | D_200 S_000 )^0 + ( G_022 S_000 | D_200 S_000 )^0 - ( I_042 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[171] = etfac[1] * AUX_INT__h_s_d_s[102] + 3 * one_over_2q * AUX_INT__g_s_d_s[72] - p_over_q * AUX_INT__i_s_d_s[138];

                    // ( H_032 S_000 | F_201 S_000 )^0 = z * ( H_032 S_000 | D_200 S_000 )^0 + ( G_031 S_000 | D_200 S_000 )^0 - ( I_033 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[172] = etfac[2] * AUX_INT__h_s_d_s[102] + 2 * one_over_2q * AUX_INT__g_s_d_s[66] - p_over_q * AUX_INT__i_s_d_s[144];

                    // ( H_032 S_000 | F_120 S_000 )^0 = x * ( H_032 S_000 | D_020 S_000 )^0 - ( I_132 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[173] = etfac[0] * AUX_INT__h_s_d_s[105] - p_over_q * AUX_INT__i_s_d_s[105];

                    // ( H_032 S_000 | F_102 S_000 )^0 = x * ( H_032 S_000 | D_002 S_000 )^0 - ( I_132 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[175] = etfac[0] * AUX_INT__h_s_d_s[107] - p_over_q * AUX_INT__i_s_d_s[107];

                    // ( H_032 S_000 | F_030 S_000 )^0 = y * ( H_032 S_000 | D_020 S_000 )^0 + ( G_022 S_000 | D_020 S_000 )^0 + ( H_032 S_000 | P_010 S_000 )^0 - ( I_042 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[176] = etfac[1] * AUX_INT__h_s_d_s[105] + 3 * one_over_2q * AUX_INT__g_s_d_s[75] + 2 * one_over_2q * AUX_INT__h_s_p_s[52] - p_over_q * AUX_INT__i_s_d_s[141];

                    // ( H_032 S_000 | F_021 S_000 )^0 = z * ( H_032 S_000 | D_020 S_000 )^0 + ( G_031 S_000 | D_020 S_000 )^0 - ( I_033 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[177] = etfac[2] * AUX_INT__h_s_d_s[105] + 2 * one_over_2q * AUX_INT__g_s_d_s[69] - p_over_q * AUX_INT__i_s_d_s[147];

                    // ( H_032 S_000 | F_003 S_000 )^0 = z * ( H_032 S_000 | D_002 S_000 )^0 + ( G_031 S_000 | D_002 S_000 )^0 + ( H_032 S_000 | P_001 S_000 )^0 - ( I_033 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[179] = etfac[2] * AUX_INT__h_s_d_s[107] + 2 * one_over_2q * AUX_INT__g_s_d_s[71] + 2 * one_over_2q * AUX_INT__h_s_p_s[53] - p_over_q * AUX_INT__i_s_d_s[149];

                    // ( H_023 S_000 | F_300 S_000 )^0 = x * ( H_023 S_000 | D_200 S_000 )^0 + ( H_023 S_000 | P_100 S_000 )^0 - ( I_123 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[180] = etfac[0] * AUX_INT__h_s_d_s[108] + 2 * one_over_2q * AUX_INT__h_s_p_s[54] - p_over_q * AUX_INT__i_s_d_s[108];

                    // ( H_023 S_000 | F_210 S_000 )^0 = y * ( H_023 S_000 | D_200 S_000 )^0 + ( G_013 S_000 | D_200 S_000 )^0 - ( I_033 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[181] = etfac[1] * AUX_INT__h_s_d_s[108] + 2 * one_over_2q * AUX_INT__g_s_d_s[78] - p_over_q * AUX_INT__i_s_d_s[144];

                    // ( H_023 S_000 | F_201 S_000 )^0 = z * ( H_023 S_000 | D_200 S_000 )^0 + ( G_022 S_000 | D_200 S_000 )^0 - ( I_024 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[182] = etfac[2] * AUX_INT__h_s_d_s[108] + 3 * one_over_2q * AUX_INT__g_s_d_s[72] - p_over_q * AUX_INT__i_s_d_s[150];

                    // ( H_023 S_000 | F_120 S_000 )^0 = x * ( H_023 S_000 | D_020 S_000 )^0 - ( I_123 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[183] = etfac[0] * AUX_INT__h_s_d_s[111] - p_over_q * AUX_INT__i_s_d_s[111];

                    // ( H_023 S_000 | F_102 S_000 )^0 = x * ( H_023 S_000 | D_002 S_000 )^0 - ( I_123 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[185] = etfac[0] * AUX_INT__h_s_d_s[113] - p_over_q * AUX_INT__i_s_d_s[113];

                    // ( H_023 S_000 | F_030 S_000 )^0 = y * ( H_023 S_000 | D_020 S_000 )^0 + ( G_013 S_000 | D_020 S_000 )^0 + ( H_023 S_000 | P_010 S_000 )^0 - ( I_033 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[186] = etfac[1] * AUX_INT__h_s_d_s[111] + 2 * one_over_2q * AUX_INT__g_s_d_s[81] + 2 * one_over_2q * AUX_INT__h_s_p_s[55] - p_over_q * AUX_INT__i_s_d_s[147];

                    // ( H_023 S_000 | F_021 S_000 )^0 = z * ( H_023 S_000 | D_020 S_000 )^0 + ( G_022 S_000 | D_020 S_000 )^0 - ( I_024 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[187] = etfac[2] * AUX_INT__h_s_d_s[111] + 3 * one_over_2q * AUX_INT__g_s_d_s[75] - p_over_q * AUX_INT__i_s_d_s[153];

                    // ( H_023 S_000 | F_003 S_000 )^0 = z * ( H_023 S_000 | D_002 S_000 )^0 + ( G_022 S_000 | D_002 S_000 )^0 + ( H_023 S_000 | P_001 S_000 )^0 - ( I_024 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[189] = etfac[2] * AUX_INT__h_s_d_s[113] + 3 * one_over_2q * AUX_INT__g_s_d_s[77] + 2 * one_over_2q * AUX_INT__h_s_p_s[56] - p_over_q * AUX_INT__i_s_d_s[155];

                    // ( H_014 S_000 | F_300 S_000 )^0 = x * ( H_014 S_000 | D_200 S_000 )^0 + ( H_014 S_000 | P_100 S_000 )^0 - ( I_114 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[190] = etfac[0] * AUX_INT__h_s_d_s[114] + 2 * one_over_2q * AUX_INT__h_s_p_s[57] - p_over_q * AUX_INT__i_s_d_s[114];

                    // ( H_014 S_000 | F_210 S_000 )^0 = y * ( H_014 S_000 | D_200 S_000 )^0 + ( G_004 S_000 | D_200 S_000 )^0 - ( I_024 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[191] = etfac[1] * AUX_INT__h_s_d_s[114] + 1 * one_over_2q * AUX_INT__g_s_d_s[84] - p_over_q * AUX_INT__i_s_d_s[150];

                    // ( H_014 S_000 | F_201 S_000 )^0 = z * ( H_014 S_000 | D_200 S_000 )^0 + ( G_013 S_000 | D_200 S_000 )^0 - ( I_015 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[192] = etfac[2] * AUX_INT__h_s_d_s[114] + 4 * one_over_2q * AUX_INT__g_s_d_s[78] - p_over_q * AUX_INT__i_s_d_s[156];

                    // ( H_014 S_000 | F_120 S_000 )^0 = x * ( H_014 S_000 | D_020 S_000 )^0 - ( I_114 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[193] = etfac[0] * AUX_INT__h_s_d_s[117] - p_over_q * AUX_INT__i_s_d_s[117];

                    // ( H_014 S_000 | F_102 S_000 )^0 = x * ( H_014 S_000 | D_002 S_000 )^0 - ( I_114 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[195] = etfac[0] * AUX_INT__h_s_d_s[119] - p_over_q * AUX_INT__i_s_d_s[119];

                    // ( H_014 S_000 | F_030 S_000 )^0 = y * ( H_014 S_000 | D_020 S_000 )^0 + ( G_004 S_000 | D_020 S_000 )^0 + ( H_014 S_000 | P_010 S_000 )^0 - ( I_024 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[196] = etfac[1] * AUX_INT__h_s_d_s[117] + 1 * one_over_2q * AUX_INT__g_s_d_s[87] + 2 * one_over_2q * AUX_INT__h_s_p_s[58] - p_over_q * AUX_INT__i_s_d_s[153];

                    // ( H_014 S_000 | F_021 S_000 )^0 = z * ( H_014 S_000 | D_020 S_000 )^0 + ( G_013 S_000 | D_020 S_000 )^0 - ( I_015 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[197] = etfac[2] * AUX_INT__h_s_d_s[117] + 4 * one_over_2q * AUX_INT__g_s_d_s[81] - p_over_q * AUX_INT__i_s_d_s[159];

                    // ( H_014 S_000 | F_003 S_000 )^0 = z * ( H_014 S_000 | D_002 S_000 )^0 + ( G_013 S_000 | D_002 S_000 )^0 + ( H_014 S_000 | P_001 S_000 )^0 - ( I_015 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[199] = etfac[2] * AUX_INT__h_s_d_s[119] + 4 * one_over_2q * AUX_INT__g_s_d_s[83] + 2 * one_over_2q * AUX_INT__h_s_p_s[59] - p_over_q * AUX_INT__i_s_d_s[161];

                    // ( H_005 S_000 | F_300 S_000 )^0 = x * ( H_005 S_000 | D_200 S_000 )^0 + ( H_005 S_000 | P_100 S_000 )^0 - ( I_105 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[200] = etfac[0] * AUX_INT__h_s_d_s[120] + 2 * one_over_2q * AUX_INT__h_s_p_s[60] - p_over_q * AUX_INT__i_s_d_s[120];

                    // ( H_005 S_000 | F_210 S_000 )^0 = y * ( H_005 S_000 | D_200 S_000 )^0 - ( I_015 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[201] = etfac[1] * AUX_INT__h_s_d_s[120] - p_over_q * AUX_INT__i_s_d_s[156];

                    // ( H_005 S_000 | F_201 S_000 )^0 = z * ( H_005 S_000 | D_200 S_000 )^0 + ( G_004 S_000 | D_200 S_000 )^0 - ( I_006 S_000 | D_200 S_000 )^0
                    AUX_INT__h_s_f_s[202] = etfac[2] * AUX_INT__h_s_d_s[120] + 5 * one_over_2q * AUX_INT__g_s_d_s[84] - p_over_q * AUX_INT__i_s_d_s[162];

                    // ( H_005 S_000 | F_120 S_000 )^0 = x * ( H_005 S_000 | D_020 S_000 )^0 - ( I_105 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[203] = etfac[0] * AUX_INT__h_s_d_s[123] - p_over_q * AUX_INT__i_s_d_s[123];

                    // ( H_005 S_000 | F_030 S_000 )^0 = y * ( H_005 S_000 | D_020 S_000 )^0 + ( H_005 S_000 | P_010 S_000 )^0 - ( I_015 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[206] = etfac[1] * AUX_INT__h_s_d_s[123] + 2 * one_over_2q * AUX_INT__h_s_p_s[61] - p_over_q * AUX_INT__i_s_d_s[159];

                    // ( H_005 S_000 | F_021 S_000 )^0 = z * ( H_005 S_000 | D_020 S_000 )^0 + ( G_004 S_000 | D_020 S_000 )^0 - ( I_006 S_000 | D_020 S_000 )^0
                    AUX_INT__h_s_f_s[207] = etfac[2] * AUX_INT__h_s_d_s[123] + 5 * one_over_2q * AUX_INT__g_s_d_s[87] - p_over_q * AUX_INT__i_s_d_s[165];

                    // ( H_005 S_000 | F_003 S_000 )^0 = z * ( H_005 S_000 | D_002 S_000 )^0 + ( G_004 S_000 | D_002 S_000 )^0 + ( H_005 S_000 | P_001 S_000 )^0 - ( I_006 S_000 | D_002 S_000 )^0
                    AUX_INT__h_s_f_s[209] = etfac[2] * AUX_INT__h_s_d_s[125] + 5 * one_over_2q * AUX_INT__g_s_d_s[89] + 2 * one_over_2q * AUX_INT__h_s_p_s[62] - p_over_q * AUX_INT__i_s_d_s[167];

                    // ( G_400 S_000 | G_400 S_000 )^0_{t} = x * ( G_400 S_000 | F_300 S_000 )^0 + ( F_300 S_000 | F_300 S_000 )^0 + ( G_400 S_000 | D_200 S_000 )^0 - ( H_500 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[0] = etfac[0] * AUX_INT__g_s_f_s[0] + 4 * one_over_2q * AUX_INT__f_s_f_s[0] + 3 * one_over_2q * AUX_INT__g_s_d_s[0] - p_over_q * AUX_INT__h_s_f_s[0];

                    // ( G_400 S_000 | G_310 S_000 )^0_{t} = y * ( G_400 S_000 | F_300 S_000 )^0 - ( H_410 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[1] = etfac[1] * AUX_INT__g_s_f_s[0] - p_over_q * AUX_INT__h_s_f_s[10];

                    // ( G_400 S_000 | G_301 S_000 )^0_{t} = z * ( G_400 S_000 | F_300 S_000 )^0 - ( H_401 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[2] = etfac[2] * AUX_INT__g_s_f_s[0] - p_over_q * AUX_INT__h_s_f_s[20];

                    // ( G_400 S_000 | G_220 S_000 )^0_{t} = y * ( G_400 S_000 | F_210 S_000 )^0 + ( G_400 S_000 | D_200 S_000 )^0 - ( H_410 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[3] = etfac[1] * AUX_INT__g_s_f_s[1] + 1 * one_over_2q * AUX_INT__g_s_d_s[0] - p_over_q * AUX_INT__h_s_f_s[11];

                    // ( G_400 S_000 | G_211 S_000 )^0_{t} = z * ( G_400 S_000 | F_210 S_000 )^0 - ( H_401 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[4] = etfac[2] * AUX_INT__g_s_f_s[1] - p_over_q * AUX_INT__h_s_f_s[21];

                    // ( G_400 S_000 | G_202 S_000 )^0_{t} = z * ( G_400 S_000 | F_201 S_000 )^0 + ( G_400 S_000 | D_200 S_000 )^0 - ( H_401 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[5] = etfac[2] * AUX_INT__g_s_f_s[2] + 1 * one_over_2q * AUX_INT__g_s_d_s[0] - p_over_q * AUX_INT__h_s_f_s[22];

                    // ( G_400 S_000 | G_130 S_000 )^0_{t} = x * ( G_400 S_000 | F_030 S_000 )^0 + ( F_300 S_000 | F_030 S_000 )^0 - ( H_500 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[6] = etfac[0] * AUX_INT__g_s_f_s[6] + 4 * one_over_2q * AUX_INT__f_s_f_s[6] - p_over_q * AUX_INT__h_s_f_s[6];

                    // ( G_400 S_000 | G_121 S_000 )^0_{t} = z * ( G_400 S_000 | F_120 S_000 )^0 - ( H_401 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[7] = etfac[2] * AUX_INT__g_s_f_s[3] - p_over_q * AUX_INT__h_s_f_s[23];

                    // ( G_400 S_000 | G_112 S_000 )^0_{t} = y * ( G_400 S_000 | F_102 S_000 )^0 - ( H_410 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[8] = etfac[1] * AUX_INT__g_s_f_s[5] - p_over_q * AUX_INT__h_s_f_s[15];

                    // ( G_400 S_000 | G_103 S_000 )^0_{t} = x * ( G_400 S_000 | F_003 S_000 )^0 + ( F_300 S_000 | F_003 S_000 )^0 - ( H_500 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[9] = etfac[0] * AUX_INT__g_s_f_s[9] + 4 * one_over_2q * AUX_INT__f_s_f_s[9] - p_over_q * AUX_INT__h_s_f_s[9];

                    // ( G_400 S_000 | G_040 S_000 )^0_{t} = y * ( G_400 S_000 | F_030 S_000 )^0 + ( G_400 S_000 | D_020 S_000 )^0 - ( H_410 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[10] = etfac[1] * AUX_INT__g_s_f_s[6] + 3 * one_over_2q * AUX_INT__g_s_d_s[3] - p_over_q * AUX_INT__h_s_f_s[16];

                    // ( G_400 S_000 | G_031 S_000 )^0_{t} = z * ( G_400 S_000 | F_030 S_000 )^0 - ( H_401 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[11] = etfac[2] * AUX_INT__g_s_f_s[6] - p_over_q * AUX_INT__h_s_f_s[26];

                    // ( G_400 S_000 | G_022 S_000 )^0_{t} = z * ( G_400 S_000 | F_021 S_000 )^0 + ( G_400 S_000 | D_020 S_000 )^0 - ( H_401 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[12] = etfac[2] * AUX_INT__g_s_f_s[7] + 1 * one_over_2q * AUX_INT__g_s_d_s[3] - p_over_q * AUX_INT__h_s_f_s[27];

                    // ( G_400 S_000 | G_013 S_000 )^0_{t} = y * ( G_400 S_000 | F_003 S_000 )^0 - ( H_410 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[13] = etfac[1] * AUX_INT__g_s_f_s[9] - p_over_q * AUX_INT__h_s_f_s[19];

                    // ( G_400 S_000 | G_004 S_000 )^0_{t} = z * ( G_400 S_000 | F_003 S_000 )^0 + ( G_400 S_000 | D_002 S_000 )^0 - ( H_401 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[14] = etfac[2] * AUX_INT__g_s_f_s[9] + 3 * one_over_2q * AUX_INT__g_s_d_s[5] - p_over_q * AUX_INT__h_s_f_s[29];

                    // ( G_310 S_000 | G_400 S_000 )^0_{t} = x * ( G_310 S_000 | F_300 S_000 )^0 + ( F_210 S_000 | F_300 S_000 )^0 + ( G_310 S_000 | D_200 S_000 )^0 - ( H_410 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[15] = etfac[0] * AUX_INT__g_s_f_s[10] + 3 * one_over_2q * AUX_INT__f_s_f_s[10] + 3 * one_over_2q * AUX_INT__g_s_d_s[6] - p_over_q * AUX_INT__h_s_f_s[10];

                    // ( G_310 S_000 | G_310 S_000 )^0_{t} = y * ( G_310 S_000 | F_300 S_000 )^0 + ( F_300 S_000 | F_300 S_000 )^0 - ( H_320 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[16] = etfac[1] * AUX_INT__g_s_f_s[10] + 1 * one_over_2q * AUX_INT__f_s_f_s[0] - p_over_q * AUX_INT__h_s_f_s[30];

                    // ( G_310 S_000 | G_301 S_000 )^0_{t} = z * ( G_310 S_000 | F_300 S_000 )^0 - ( H_311 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[17] = etfac[2] * AUX_INT__g_s_f_s[10] - p_over_q * AUX_INT__h_s_f_s[40];

                    // ( G_310 S_000 | G_220 S_000 )^0_{t} = y * ( G_310 S_000 | F_210 S_000 )^0 + ( F_300 S_000 | F_210 S_000 )^0 + ( G_310 S_000 | D_200 S_000 )^0 - ( H_320 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[18] = etfac[1] * AUX_INT__g_s_f_s[11] + 1 * one_over_2q * AUX_INT__f_s_f_s[1] + 1 * one_over_2q * AUX_INT__g_s_d_s[6] - p_over_q * AUX_INT__h_s_f_s[31];

                    // ( G_310 S_000 | G_211 S_000 )^0_{t} = z * ( G_310 S_000 | F_210 S_000 )^0 - ( H_311 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[19] = etfac[2] * AUX_INT__g_s_f_s[11] - p_over_q * AUX_INT__h_s_f_s[41];

                    // ( G_310 S_000 | G_202 S_000 )^0_{t} = z * ( G_310 S_000 | F_201 S_000 )^0 + ( G_310 S_000 | D_200 S_000 )^0 - ( H_311 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[20] = etfac[2] * AUX_INT__g_s_f_s[12] + 1 * one_over_2q * AUX_INT__g_s_d_s[6] - p_over_q * AUX_INT__h_s_f_s[42];

                    // ( G_310 S_000 | G_130 S_000 )^0_{t} = x * ( G_310 S_000 | F_030 S_000 )^0 + ( F_210 S_000 | F_030 S_000 )^0 - ( H_410 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[21] = etfac[0] * AUX_INT__g_s_f_s[16] + 3 * one_over_2q * AUX_INT__f_s_f_s[16] - p_over_q * AUX_INT__h_s_f_s[16];

                    // ( G_310 S_000 | G_121 S_000 )^0_{t} = z * ( G_310 S_000 | F_120 S_000 )^0 - ( H_311 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[22] = etfac[2] * AUX_INT__g_s_f_s[13] - p_over_q * AUX_INT__h_s_f_s[43];

                    // ( G_310 S_000 | G_112 S_000 )^0_{t} = y * ( G_310 S_000 | F_102 S_000 )^0 + ( F_300 S_000 | F_102 S_000 )^0 - ( H_320 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[23] = etfac[1] * AUX_INT__g_s_f_s[15] + 1 * one_over_2q * AUX_INT__f_s_f_s[5] - p_over_q * AUX_INT__h_s_f_s[35];

                    // ( G_310 S_000 | G_103 S_000 )^0_{t} = x * ( G_310 S_000 | F_003 S_000 )^0 + ( F_210 S_000 | F_003 S_000 )^0 - ( H_410 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[24] = etfac[0] * AUX_INT__g_s_f_s[19] + 3 * one_over_2q * AUX_INT__f_s_f_s[19] - p_over_q * AUX_INT__h_s_f_s[19];

                    // ( G_310 S_000 | G_040 S_000 )^0_{t} = y * ( G_310 S_000 | F_030 S_000 )^0 + ( F_300 S_000 | F_030 S_000 )^0 + ( G_310 S_000 | D_020 S_000 )^0 - ( H_320 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[25] = etfac[1] * AUX_INT__g_s_f_s[16] + 1 * one_over_2q * AUX_INT__f_s_f_s[6] + 3 * one_over_2q * AUX_INT__g_s_d_s[9] - p_over_q * AUX_INT__h_s_f_s[36];

                    // ( G_310 S_000 | G_031 S_000 )^0_{t} = z * ( G_310 S_000 | F_030 S_000 )^0 - ( H_311 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[26] = etfac[2] * AUX_INT__g_s_f_s[16] - p_over_q * AUX_INT__h_s_f_s[46];

                    // ( G_310 S_000 | G_022 S_000 )^0_{t} = z * ( G_310 S_000 | F_021 S_000 )^0 + ( G_310 S_000 | D_020 S_000 )^0 - ( H_311 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[27] = etfac[2] * AUX_INT__g_s_f_s[17] + 1 * one_over_2q * AUX_INT__g_s_d_s[9] - p_over_q * AUX_INT__h_s_f_s[47];

                    // ( G_310 S_000 | G_013 S_000 )^0_{t} = y * ( G_310 S_000 | F_003 S_000 )^0 + ( F_300 S_000 | F_003 S_000 )^0 - ( H_320 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[28] = etfac[1] * AUX_INT__g_s_f_s[19] + 1 * one_over_2q * AUX_INT__f_s_f_s[9] - p_over_q * AUX_INT__h_s_f_s[39];

                    // ( G_310 S_000 | G_004 S_000 )^0_{t} = z * ( G_310 S_000 | F_003 S_000 )^0 + ( G_310 S_000 | D_002 S_000 )^0 - ( H_311 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[29] = etfac[2] * AUX_INT__g_s_f_s[19] + 3 * one_over_2q * AUX_INT__g_s_d_s[11] - p_over_q * AUX_INT__h_s_f_s[49];

                    // ( G_301 S_000 | G_400 S_000 )^0_{t} = x * ( G_301 S_000 | F_300 S_000 )^0 + ( F_201 S_000 | F_300 S_000 )^0 + ( G_301 S_000 | D_200 S_000 )^0 - ( H_401 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[30] = etfac[0] * AUX_INT__g_s_f_s[20] + 3 * one_over_2q * AUX_INT__f_s_f_s[20] + 3 * one_over_2q * AUX_INT__g_s_d_s[12] - p_over_q * AUX_INT__h_s_f_s[20];

                    // ( G_301 S_000 | G_310 S_000 )^0_{t} = y * ( G_301 S_000 | F_300 S_000 )^0 - ( H_311 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[31] = etfac[1] * AUX_INT__g_s_f_s[20] - p_over_q * AUX_INT__h_s_f_s[40];

                    // ( G_301 S_000 | G_301 S_000 )^0_{t} = z * ( G_301 S_000 | F_300 S_000 )^0 + ( F_300 S_000 | F_300 S_000 )^0 - ( H_302 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[32] = etfac[2] * AUX_INT__g_s_f_s[20] + 1 * one_over_2q * AUX_INT__f_s_f_s[0] - p_over_q * AUX_INT__h_s_f_s[50];

                    // ( G_301 S_000 | G_220 S_000 )^0_{t} = y * ( G_301 S_000 | F_210 S_000 )^0 + ( G_301 S_000 | D_200 S_000 )^0 - ( H_311 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[33] = etfac[1] * AUX_INT__g_s_f_s[21] + 1 * one_over_2q * AUX_INT__g_s_d_s[12] - p_over_q * AUX_INT__h_s_f_s[41];

                    // ( G_301 S_000 | G_211 S_000 )^0_{t} = z * ( G_301 S_000 | F_210 S_000 )^0 + ( F_300 S_000 | F_210 S_000 )^0 - ( H_302 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[34] = etfac[2] * AUX_INT__g_s_f_s[21] + 1 * one_over_2q * AUX_INT__f_s_f_s[1] - p_over_q * AUX_INT__h_s_f_s[51];

                    // ( G_301 S_000 | G_202 S_000 )^0_{t} = z * ( G_301 S_000 | F_201 S_000 )^0 + ( F_300 S_000 | F_201 S_000 )^0 + ( G_301 S_000 | D_200 S_000 )^0 - ( H_302 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[35] = etfac[2] * AUX_INT__g_s_f_s[22] + 1 * one_over_2q * AUX_INT__f_s_f_s[2] + 1 * one_over_2q * AUX_INT__g_s_d_s[12] - p_over_q * AUX_INT__h_s_f_s[52];

                    // ( G_301 S_000 | G_130 S_000 )^0_{t} = x * ( G_301 S_000 | F_030 S_000 )^0 + ( F_201 S_000 | F_030 S_000 )^0 - ( H_401 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[36] = etfac[0] * AUX_INT__g_s_f_s[26] + 3 * one_over_2q * AUX_INT__f_s_f_s[26] - p_over_q * AUX_INT__h_s_f_s[26];

                    // ( G_301 S_000 | G_121 S_000 )^0_{t} = z * ( G_301 S_000 | F_120 S_000 )^0 + ( F_300 S_000 | F_120 S_000 )^0 - ( H_302 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[37] = etfac[2] * AUX_INT__g_s_f_s[23] + 1 * one_over_2q * AUX_INT__f_s_f_s[3] - p_over_q * AUX_INT__h_s_f_s[53];

                    // ( G_301 S_000 | G_112 S_000 )^0_{t} = y * ( G_301 S_000 | F_102 S_000 )^0 - ( H_311 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[38] = etfac[1] * AUX_INT__g_s_f_s[25] - p_over_q * AUX_INT__h_s_f_s[45];

                    // ( G_301 S_000 | G_103 S_000 )^0_{t} = x * ( G_301 S_000 | F_003 S_000 )^0 + ( F_201 S_000 | F_003 S_000 )^0 - ( H_401 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[39] = etfac[0] * AUX_INT__g_s_f_s[29] + 3 * one_over_2q * AUX_INT__f_s_f_s[29] - p_over_q * AUX_INT__h_s_f_s[29];

                    // ( G_301 S_000 | G_040 S_000 )^0_{t} = y * ( G_301 S_000 | F_030 S_000 )^0 + ( G_301 S_000 | D_020 S_000 )^0 - ( H_311 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[40] = etfac[1] * AUX_INT__g_s_f_s[26] + 3 * one_over_2q * AUX_INT__g_s_d_s[15] - p_over_q * AUX_INT__h_s_f_s[46];

                    // ( G_301 S_000 | G_031 S_000 )^0_{t} = z * ( G_301 S_000 | F_030 S_000 )^0 + ( F_300 S_000 | F_030 S_000 )^0 - ( H_302 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[41] = etfac[2] * AUX_INT__g_s_f_s[26] + 1 * one_over_2q * AUX_INT__f_s_f_s[6] - p_over_q * AUX_INT__h_s_f_s[56];

                    // ( G_301 S_000 | G_022 S_000 )^0_{t} = z * ( G_301 S_000 | F_021 S_000 )^0 + ( F_300 S_000 | F_021 S_000 )^0 + ( G_301 S_000 | D_020 S_000 )^0 - ( H_302 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[42] = etfac[2] * AUX_INT__g_s_f_s[27] + 1 * one_over_2q * AUX_INT__f_s_f_s[7] + 1 * one_over_2q * AUX_INT__g_s_d_s[15] - p_over_q * AUX_INT__h_s_f_s[57];

                    // ( G_301 S_000 | G_013 S_000 )^0_{t} = y * ( G_301 S_000 | F_003 S_000 )^0 - ( H_311 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[43] = etfac[1] * AUX_INT__g_s_f_s[29] - p_over_q * AUX_INT__h_s_f_s[49];

                    // ( G_301 S_000 | G_004 S_000 )^0_{t} = z * ( G_301 S_000 | F_003 S_000 )^0 + ( F_300 S_000 | F_003 S_000 )^0 + ( G_301 S_000 | D_002 S_000 )^0 - ( H_302 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[44] = etfac[2] * AUX_INT__g_s_f_s[29] + 1 * one_over_2q * AUX_INT__f_s_f_s[9] + 3 * one_over_2q * AUX_INT__g_s_d_s[17] - p_over_q * AUX_INT__h_s_f_s[59];

                    // ( G_220 S_000 | G_400 S_000 )^0_{t} = x * ( G_220 S_000 | F_300 S_000 )^0 + ( F_120 S_000 | F_300 S_000 )^0 + ( G_220 S_000 | D_200 S_000 )^0 - ( H_320 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[45] = etfac[0] * AUX_INT__g_s_f_s[30] + 2 * one_over_2q * AUX_INT__f_s_f_s[30] + 3 * one_over_2q * AUX_INT__g_s_d_s[18] - p_over_q * AUX_INT__h_s_f_s[30];

                    // ( G_220 S_000 | G_310 S_000 )^0_{t} = y * ( G_220 S_000 | F_300 S_000 )^0 + ( F_210 S_000 | F_300 S_000 )^0 - ( H_230 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[46] = etfac[1] * AUX_INT__g_s_f_s[30] + 2 * one_over_2q * AUX_INT__f_s_f_s[10] - p_over_q * AUX_INT__h_s_f_s[60];

                    // ( G_220 S_000 | G_301 S_000 )^0_{t} = z * ( G_220 S_000 | F_300 S_000 )^0 - ( H_221 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[47] = etfac[2] * AUX_INT__g_s_f_s[30] - p_over_q * AUX_INT__h_s_f_s[70];

                    // ( G_220 S_000 | G_220 S_000 )^0_{t} = y * ( G_220 S_000 | F_210 S_000 )^0 + ( F_210 S_000 | F_210 S_000 )^0 + ( G_220 S_000 | D_200 S_000 )^0 - ( H_230 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[48] = etfac[1] * AUX_INT__g_s_f_s[31] + 2 * one_over_2q * AUX_INT__f_s_f_s[11] + 1 * one_over_2q * AUX_INT__g_s_d_s[18] - p_over_q * AUX_INT__h_s_f_s[61];

                    // ( G_220 S_000 | G_211 S_000 )^0_{t} = z * ( G_220 S_000 | F_210 S_000 )^0 - ( H_221 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[49] = etfac[2] * AUX_INT__g_s_f_s[31] - p_over_q * AUX_INT__h_s_f_s[71];

                    // ( G_220 S_000 | G_202 S_000 )^0_{t} = z * ( G_220 S_000 | F_201 S_000 )^0 + ( G_220 S_000 | D_200 S_000 )^0 - ( H_221 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[50] = etfac[2] * AUX_INT__g_s_f_s[32] + 1 * one_over_2q * AUX_INT__g_s_d_s[18] - p_over_q * AUX_INT__h_s_f_s[72];

                    // ( G_220 S_000 | G_130 S_000 )^0_{t} = x * ( G_220 S_000 | F_030 S_000 )^0 + ( F_120 S_000 | F_030 S_000 )^0 - ( H_320 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[51] = etfac[0] * AUX_INT__g_s_f_s[36] + 2 * one_over_2q * AUX_INT__f_s_f_s[36] - p_over_q * AUX_INT__h_s_f_s[36];

                    // ( G_220 S_000 | G_121 S_000 )^0_{t} = z * ( G_220 S_000 | F_120 S_000 )^0 - ( H_221 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[52] = etfac[2] * AUX_INT__g_s_f_s[33] - p_over_q * AUX_INT__h_s_f_s[73];

                    // ( G_220 S_000 | G_112 S_000 )^0_{t} = y * ( G_220 S_000 | F_102 S_000 )^0 + ( F_210 S_000 | F_102 S_000 )^0 - ( H_230 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[53] = etfac[1] * AUX_INT__g_s_f_s[35] + 2 * one_over_2q * AUX_INT__f_s_f_s[15] - p_over_q * AUX_INT__h_s_f_s[65];

                    // ( G_220 S_000 | G_103 S_000 )^0_{t} = x * ( G_220 S_000 | F_003 S_000 )^0 + ( F_120 S_000 | F_003 S_000 )^0 - ( H_320 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[54] = etfac[0] * AUX_INT__g_s_f_s[39] + 2 * one_over_2q * AUX_INT__f_s_f_s[39] - p_over_q * AUX_INT__h_s_f_s[39];

                    // ( G_220 S_000 | G_040 S_000 )^0_{t} = y * ( G_220 S_000 | F_030 S_000 )^0 + ( F_210 S_000 | F_030 S_000 )^0 + ( G_220 S_000 | D_020 S_000 )^0 - ( H_230 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[55] = etfac[1] * AUX_INT__g_s_f_s[36] + 2 * one_over_2q * AUX_INT__f_s_f_s[16] + 3 * one_over_2q * AUX_INT__g_s_d_s[21] - p_over_q * AUX_INT__h_s_f_s[66];

                    // ( G_220 S_000 | G_031 S_000 )^0_{t} = z * ( G_220 S_000 | F_030 S_000 )^0 - ( H_221 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[56] = etfac[2] * AUX_INT__g_s_f_s[36] - p_over_q * AUX_INT__h_s_f_s[76];

                    // ( G_220 S_000 | G_022 S_000 )^0_{t} = z * ( G_220 S_000 | F_021 S_000 )^0 + ( G_220 S_000 | D_020 S_000 )^0 - ( H_221 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[57] = etfac[2] * AUX_INT__g_s_f_s[37] + 1 * one_over_2q * AUX_INT__g_s_d_s[21] - p_over_q * AUX_INT__h_s_f_s[77];

                    // ( G_220 S_000 | G_013 S_000 )^0_{t} = y * ( G_220 S_000 | F_003 S_000 )^0 + ( F_210 S_000 | F_003 S_000 )^0 - ( H_230 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[58] = etfac[1] * AUX_INT__g_s_f_s[39] + 2 * one_over_2q * AUX_INT__f_s_f_s[19] - p_over_q * AUX_INT__h_s_f_s[69];

                    // ( G_220 S_000 | G_004 S_000 )^0_{t} = z * ( G_220 S_000 | F_003 S_000 )^0 + ( G_220 S_000 | D_002 S_000 )^0 - ( H_221 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[59] = etfac[2] * AUX_INT__g_s_f_s[39] + 3 * one_over_2q * AUX_INT__g_s_d_s[23] - p_over_q * AUX_INT__h_s_f_s[79];

                    // ( G_211 S_000 | G_400 S_000 )^0_{t} = x * ( G_211 S_000 | F_300 S_000 )^0 + ( F_111 S_000 | F_300 S_000 )^0 + ( G_211 S_000 | D_200 S_000 )^0 - ( H_311 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[60] = etfac[0] * AUX_INT__g_s_f_s[40] + 2 * one_over_2q * AUX_INT__f_s_f_s[40] + 3 * one_over_2q * AUX_INT__g_s_d_s[24] - p_over_q * AUX_INT__h_s_f_s[40];

                    // ( G_211 S_000 | G_310 S_000 )^0_{t} = y * ( G_211 S_000 | F_300 S_000 )^0 + ( F_201 S_000 | F_300 S_000 )^0 - ( H_221 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[61] = etfac[1] * AUX_INT__g_s_f_s[40] + 1 * one_over_2q * AUX_INT__f_s_f_s[20] - p_over_q * AUX_INT__h_s_f_s[70];

                    // ( G_211 S_000 | G_301 S_000 )^0_{t} = z * ( G_211 S_000 | F_300 S_000 )^0 + ( F_210 S_000 | F_300 S_000 )^0 - ( H_212 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[62] = etfac[2] * AUX_INT__g_s_f_s[40] + 1 * one_over_2q * AUX_INT__f_s_f_s[10] - p_over_q * AUX_INT__h_s_f_s[80];

                    // ( G_211 S_000 | G_220 S_000 )^0_{t} = y * ( G_211 S_000 | F_210 S_000 )^0 + ( F_201 S_000 | F_210 S_000 )^0 + ( G_211 S_000 | D_200 S_000 )^0 - ( H_221 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[63] = etfac[1] * AUX_INT__g_s_f_s[41] + 1 * one_over_2q * AUX_INT__f_s_f_s[21] + 1 * one_over_2q * AUX_INT__g_s_d_s[24] - p_over_q * AUX_INT__h_s_f_s[71];

                    // ( G_211 S_000 | G_211 S_000 )^0_{t} = z * ( G_211 S_000 | F_210 S_000 )^0 + ( F_210 S_000 | F_210 S_000 )^0 - ( H_212 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[64] = etfac[2] * AUX_INT__g_s_f_s[41] + 1 * one_over_2q * AUX_INT__f_s_f_s[11] - p_over_q * AUX_INT__h_s_f_s[81];

                    // ( G_211 S_000 | G_202 S_000 )^0_{t} = z * ( G_211 S_000 | F_201 S_000 )^0 + ( F_210 S_000 | F_201 S_000 )^0 + ( G_211 S_000 | D_200 S_000 )^0 - ( H_212 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[65] = etfac[2] * AUX_INT__g_s_f_s[42] + 1 * one_over_2q * AUX_INT__f_s_f_s[12] + 1 * one_over_2q * AUX_INT__g_s_d_s[24] - p_over_q * AUX_INT__h_s_f_s[82];

                    // ( G_211 S_000 | G_130 S_000 )^0_{t} = x * ( G_211 S_000 | F_030 S_000 )^0 + ( F_111 S_000 | F_030 S_000 )^0 - ( H_311 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[66] = etfac[0] * AUX_INT__g_s_f_s[46] + 2 * one_over_2q * AUX_INT__f_s_f_s[46] - p_over_q * AUX_INT__h_s_f_s[46];

                    // ( G_211 S_000 | G_121 S_000 )^0_{t} = z * ( G_211 S_000 | F_120 S_000 )^0 + ( F_210 S_000 | F_120 S_000 )^0 - ( H_212 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[67] = etfac[2] * AUX_INT__g_s_f_s[43] + 1 * one_over_2q * AUX_INT__f_s_f_s[13] - p_over_q * AUX_INT__h_s_f_s[83];

                    // ( G_211 S_000 | G_112 S_000 )^0_{t} = y * ( G_211 S_000 | F_102 S_000 )^0 + ( F_201 S_000 | F_102 S_000 )^0 - ( H_221 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[68] = etfac[1] * AUX_INT__g_s_f_s[45] + 1 * one_over_2q * AUX_INT__f_s_f_s[25] - p_over_q * AUX_INT__h_s_f_s[75];

                    // ( G_211 S_000 | G_103 S_000 )^0_{t} = x * ( G_211 S_000 | F_003 S_000 )^0 + ( F_111 S_000 | F_003 S_000 )^0 - ( H_311 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[69] = etfac[0] * AUX_INT__g_s_f_s[49] + 2 * one_over_2q * AUX_INT__f_s_f_s[49] - p_over_q * AUX_INT__h_s_f_s[49];

                    // ( G_211 S_000 | G_040 S_000 )^0_{t} = y * ( G_211 S_000 | F_030 S_000 )^0 + ( F_201 S_000 | F_030 S_000 )^0 + ( G_211 S_000 | D_020 S_000 )^0 - ( H_221 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[70] = etfac[1] * AUX_INT__g_s_f_s[46] + 1 * one_over_2q * AUX_INT__f_s_f_s[26] + 3 * one_over_2q * AUX_INT__g_s_d_s[27] - p_over_q * AUX_INT__h_s_f_s[76];

                    // ( G_211 S_000 | G_031 S_000 )^0_{t} = z * ( G_211 S_000 | F_030 S_000 )^0 + ( F_210 S_000 | F_030 S_000 )^0 - ( H_212 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[71] = etfac[2] * AUX_INT__g_s_f_s[46] + 1 * one_over_2q * AUX_INT__f_s_f_s[16] - p_over_q * AUX_INT__h_s_f_s[86];

                    // ( G_211 S_000 | G_022 S_000 )^0_{t} = z * ( G_211 S_000 | F_021 S_000 )^0 + ( F_210 S_000 | F_021 S_000 )^0 + ( G_211 S_000 | D_020 S_000 )^0 - ( H_212 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[72] = etfac[2] * AUX_INT__g_s_f_s[47] + 1 * one_over_2q * AUX_INT__f_s_f_s[17] + 1 * one_over_2q * AUX_INT__g_s_d_s[27] - p_over_q * AUX_INT__h_s_f_s[87];

                    // ( G_211 S_000 | G_013 S_000 )^0_{t} = y * ( G_211 S_000 | F_003 S_000 )^0 + ( F_201 S_000 | F_003 S_000 )^0 - ( H_221 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[73] = etfac[1] * AUX_INT__g_s_f_s[49] + 1 * one_over_2q * AUX_INT__f_s_f_s[29] - p_over_q * AUX_INT__h_s_f_s[79];

                    // ( G_211 S_000 | G_004 S_000 )^0_{t} = z * ( G_211 S_000 | F_003 S_000 )^0 + ( F_210 S_000 | F_003 S_000 )^0 + ( G_211 S_000 | D_002 S_000 )^0 - ( H_212 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[74] = etfac[2] * AUX_INT__g_s_f_s[49] + 1 * one_over_2q * AUX_INT__f_s_f_s[19] + 3 * one_over_2q * AUX_INT__g_s_d_s[29] - p_over_q * AUX_INT__h_s_f_s[89];

                    // ( G_202 S_000 | G_400 S_000 )^0_{t} = x * ( G_202 S_000 | F_300 S_000 )^0 + ( F_102 S_000 | F_300 S_000 )^0 + ( G_202 S_000 | D_200 S_000 )^0 - ( H_302 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[75] = etfac[0] * AUX_INT__g_s_f_s[50] + 2 * one_over_2q * AUX_INT__f_s_f_s[50] + 3 * one_over_2q * AUX_INT__g_s_d_s[30] - p_over_q * AUX_INT__h_s_f_s[50];

                    // ( G_202 S_000 | G_310 S_000 )^0_{t} = y * ( G_202 S_000 | F_300 S_000 )^0 - ( H_212 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[76] = etfac[1] * AUX_INT__g_s_f_s[50] - p_over_q * AUX_INT__h_s_f_s[80];

                    // ( G_202 S_000 | G_301 S_000 )^0_{t} = z * ( G_202 S_000 | F_300 S_000 )^0 + ( F_201 S_000 | F_300 S_000 )^0 - ( H_203 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[77] = etfac[2] * AUX_INT__g_s_f_s[50] + 2 * one_over_2q * AUX_INT__f_s_f_s[20] - p_over_q * AUX_INT__h_s_f_s[90];

                    // ( G_202 S_000 | G_220 S_000 )^0_{t} = y * ( G_202 S_000 | F_210 S_000 )^0 + ( G_202 S_000 | D_200 S_000 )^0 - ( H_212 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[78] = etfac[1] * AUX_INT__g_s_f_s[51] + 1 * one_over_2q * AUX_INT__g_s_d_s[30] - p_over_q * AUX_INT__h_s_f_s[81];

                    // ( G_202 S_000 | G_211 S_000 )^0_{t} = z * ( G_202 S_000 | F_210 S_000 )^0 + ( F_201 S_000 | F_210 S_000 )^0 - ( H_203 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[79] = etfac[2] * AUX_INT__g_s_f_s[51] + 2 * one_over_2q * AUX_INT__f_s_f_s[21] - p_over_q * AUX_INT__h_s_f_s[91];

                    // ( G_202 S_000 | G_202 S_000 )^0_{t} = z * ( G_202 S_000 | F_201 S_000 )^0 + ( F_201 S_000 | F_201 S_000 )^0 + ( G_202 S_000 | D_200 S_000 )^0 - ( H_203 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[80] = etfac[2] * AUX_INT__g_s_f_s[52] + 2 * one_over_2q * AUX_INT__f_s_f_s[22] + 1 * one_over_2q * AUX_INT__g_s_d_s[30] - p_over_q * AUX_INT__h_s_f_s[92];

                    // ( G_202 S_000 | G_130 S_000 )^0_{t} = x * ( G_202 S_000 | F_030 S_000 )^0 + ( F_102 S_000 | F_030 S_000 )^0 - ( H_302 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[81] = etfac[0] * AUX_INT__g_s_f_s[56] + 2 * one_over_2q * AUX_INT__f_s_f_s[56] - p_over_q * AUX_INT__h_s_f_s[56];

                    // ( G_202 S_000 | G_121 S_000 )^0_{t} = z * ( G_202 S_000 | F_120 S_000 )^0 + ( F_201 S_000 | F_120 S_000 )^0 - ( H_203 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[82] = etfac[2] * AUX_INT__g_s_f_s[53] + 2 * one_over_2q * AUX_INT__f_s_f_s[23] - p_over_q * AUX_INT__h_s_f_s[93];

                    // ( G_202 S_000 | G_112 S_000 )^0_{t} = y * ( G_202 S_000 | F_102 S_000 )^0 - ( H_212 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[83] = etfac[1] * AUX_INT__g_s_f_s[55] - p_over_q * AUX_INT__h_s_f_s[85];

                    // ( G_202 S_000 | G_103 S_000 )^0_{t} = x * ( G_202 S_000 | F_003 S_000 )^0 + ( F_102 S_000 | F_003 S_000 )^0 - ( H_302 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[84] = etfac[0] * AUX_INT__g_s_f_s[59] + 2 * one_over_2q * AUX_INT__f_s_f_s[59] - p_over_q * AUX_INT__h_s_f_s[59];

                    // ( G_202 S_000 | G_040 S_000 )^0_{t} = y * ( G_202 S_000 | F_030 S_000 )^0 + ( G_202 S_000 | D_020 S_000 )^0 - ( H_212 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[85] = etfac[1] * AUX_INT__g_s_f_s[56] + 3 * one_over_2q * AUX_INT__g_s_d_s[33] - p_over_q * AUX_INT__h_s_f_s[86];

                    // ( G_202 S_000 | G_031 S_000 )^0_{t} = z * ( G_202 S_000 | F_030 S_000 )^0 + ( F_201 S_000 | F_030 S_000 )^0 - ( H_203 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[86] = etfac[2] * AUX_INT__g_s_f_s[56] + 2 * one_over_2q * AUX_INT__f_s_f_s[26] - p_over_q * AUX_INT__h_s_f_s[96];

                    // ( G_202 S_000 | G_022 S_000 )^0_{t} = z * ( G_202 S_000 | F_021 S_000 )^0 + ( F_201 S_000 | F_021 S_000 )^0 + ( G_202 S_000 | D_020 S_000 )^0 - ( H_203 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[87] = etfac[2] * AUX_INT__g_s_f_s[57] + 2 * one_over_2q * AUX_INT__f_s_f_s[27] + 1 * one_over_2q * AUX_INT__g_s_d_s[33] - p_over_q * AUX_INT__h_s_f_s[97];

                    // ( G_202 S_000 | G_013 S_000 )^0_{t} = y * ( G_202 S_000 | F_003 S_000 )^0 - ( H_212 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[88] = etfac[1] * AUX_INT__g_s_f_s[59] - p_over_q * AUX_INT__h_s_f_s[89];

                    // ( G_202 S_000 | G_004 S_000 )^0_{t} = z * ( G_202 S_000 | F_003 S_000 )^0 + ( F_201 S_000 | F_003 S_000 )^0 + ( G_202 S_000 | D_002 S_000 )^0 - ( H_203 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[89] = etfac[2] * AUX_INT__g_s_f_s[59] + 2 * one_over_2q * AUX_INT__f_s_f_s[29] + 3 * one_over_2q * AUX_INT__g_s_d_s[35] - p_over_q * AUX_INT__h_s_f_s[99];

                    // ( G_130 S_000 | G_400 S_000 )^0_{t} = x * ( G_130 S_000 | F_300 S_000 )^0 + ( F_030 S_000 | F_300 S_000 )^0 + ( G_130 S_000 | D_200 S_000 )^0 - ( H_230 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[90] = etfac[0] * AUX_INT__g_s_f_s[60] + 1 * one_over_2q * AUX_INT__f_s_f_s[60] + 3 * one_over_2q * AUX_INT__g_s_d_s[36] - p_over_q * AUX_INT__h_s_f_s[60];

                    // ( G_130 S_000 | G_310 S_000 )^0_{t} = y * ( G_130 S_000 | F_300 S_000 )^0 + ( F_120 S_000 | F_300 S_000 )^0 - ( H_140 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[91] = etfac[1] * AUX_INT__g_s_f_s[60] + 3 * one_over_2q * AUX_INT__f_s_f_s[30] - p_over_q * AUX_INT__h_s_f_s[100];

                    // ( G_130 S_000 | G_301 S_000 )^0_{t} = z * ( G_130 S_000 | F_300 S_000 )^0 - ( H_131 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[92] = etfac[2] * AUX_INT__g_s_f_s[60] - p_over_q * AUX_INT__h_s_f_s[110];

                    // ( G_130 S_000 | G_220 S_000 )^0_{t} = y * ( G_130 S_000 | F_210 S_000 )^0 + ( F_120 S_000 | F_210 S_000 )^0 + ( G_130 S_000 | D_200 S_000 )^0 - ( H_140 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[93] = etfac[1] * AUX_INT__g_s_f_s[61] + 3 * one_over_2q * AUX_INT__f_s_f_s[31] + 1 * one_over_2q * AUX_INT__g_s_d_s[36] - p_over_q * AUX_INT__h_s_f_s[101];

                    // ( G_130 S_000 | G_211 S_000 )^0_{t} = z * ( G_130 S_000 | F_210 S_000 )^0 - ( H_131 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[94] = etfac[2] * AUX_INT__g_s_f_s[61] - p_over_q * AUX_INT__h_s_f_s[111];

                    // ( G_130 S_000 | G_202 S_000 )^0_{t} = z * ( G_130 S_000 | F_201 S_000 )^0 + ( G_130 S_000 | D_200 S_000 )^0 - ( H_131 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[95] = etfac[2] * AUX_INT__g_s_f_s[62] + 1 * one_over_2q * AUX_INT__g_s_d_s[36] - p_over_q * AUX_INT__h_s_f_s[112];

                    // ( G_130 S_000 | G_130 S_000 )^0_{t} = x * ( G_130 S_000 | F_030 S_000 )^0 + ( F_030 S_000 | F_030 S_000 )^0 - ( H_230 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[96] = etfac[0] * AUX_INT__g_s_f_s[66] + 1 * one_over_2q * AUX_INT__f_s_f_s[66] - p_over_q * AUX_INT__h_s_f_s[66];

                    // ( G_130 S_000 | G_121 S_000 )^0_{t} = z * ( G_130 S_000 | F_120 S_000 )^0 - ( H_131 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[97] = etfac[2] * AUX_INT__g_s_f_s[63] - p_over_q * AUX_INT__h_s_f_s[113];

                    // ( G_130 S_000 | G_112 S_000 )^0_{t} = y * ( G_130 S_000 | F_102 S_000 )^0 + ( F_120 S_000 | F_102 S_000 )^0 - ( H_140 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[98] = etfac[1] * AUX_INT__g_s_f_s[65] + 3 * one_over_2q * AUX_INT__f_s_f_s[35] - p_over_q * AUX_INT__h_s_f_s[105];

                    // ( G_130 S_000 | G_103 S_000 )^0_{t} = x * ( G_130 S_000 | F_003 S_000 )^0 + ( F_030 S_000 | F_003 S_000 )^0 - ( H_230 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[99] = etfac[0] * AUX_INT__g_s_f_s[69] + 1 * one_over_2q * AUX_INT__f_s_f_s[69] - p_over_q * AUX_INT__h_s_f_s[69];

                    // ( G_130 S_000 | G_040 S_000 )^0_{t} = y * ( G_130 S_000 | F_030 S_000 )^0 + ( F_120 S_000 | F_030 S_000 )^0 + ( G_130 S_000 | D_020 S_000 )^0 - ( H_140 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[100] = etfac[1] * AUX_INT__g_s_f_s[66] + 3 * one_over_2q * AUX_INT__f_s_f_s[36] + 3 * one_over_2q * AUX_INT__g_s_d_s[39] - p_over_q * AUX_INT__h_s_f_s[106];

                    // ( G_130 S_000 | G_031 S_000 )^0_{t} = z * ( G_130 S_000 | F_030 S_000 )^0 - ( H_131 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[101] = etfac[2] * AUX_INT__g_s_f_s[66] - p_over_q * AUX_INT__h_s_f_s[116];

                    // ( G_130 S_000 | G_022 S_000 )^0_{t} = z * ( G_130 S_000 | F_021 S_000 )^0 + ( G_130 S_000 | D_020 S_000 )^0 - ( H_131 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[102] = etfac[2] * AUX_INT__g_s_f_s[67] + 1 * one_over_2q * AUX_INT__g_s_d_s[39] - p_over_q * AUX_INT__h_s_f_s[117];

                    // ( G_130 S_000 | G_013 S_000 )^0_{t} = y * ( G_130 S_000 | F_003 S_000 )^0 + ( F_120 S_000 | F_003 S_000 )^0 - ( H_140 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[103] = etfac[1] * AUX_INT__g_s_f_s[69] + 3 * one_over_2q * AUX_INT__f_s_f_s[39] - p_over_q * AUX_INT__h_s_f_s[109];

                    // ( G_130 S_000 | G_004 S_000 )^0_{t} = z * ( G_130 S_000 | F_003 S_000 )^0 + ( G_130 S_000 | D_002 S_000 )^0 - ( H_131 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[104] = etfac[2] * AUX_INT__g_s_f_s[69] + 3 * one_over_2q * AUX_INT__g_s_d_s[41] - p_over_q * AUX_INT__h_s_f_s[119];

                    // ( G_121 S_000 | G_400 S_000 )^0_{t} = x * ( G_121 S_000 | F_300 S_000 )^0 + ( F_021 S_000 | F_300 S_000 )^0 + ( G_121 S_000 | D_200 S_000 )^0 - ( H_221 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[105] = etfac[0] * AUX_INT__g_s_f_s[70] + 1 * one_over_2q * AUX_INT__f_s_f_s[70] + 3 * one_over_2q * AUX_INT__g_s_d_s[42] - p_over_q * AUX_INT__h_s_f_s[70];

                    // ( G_121 S_000 | G_310 S_000 )^0_{t} = y * ( G_121 S_000 | F_300 S_000 )^0 + ( F_111 S_000 | F_300 S_000 )^0 - ( H_131 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[106] = etfac[1] * AUX_INT__g_s_f_s[70] + 2 * one_over_2q * AUX_INT__f_s_f_s[40] - p_over_q * AUX_INT__h_s_f_s[110];

                    // ( G_121 S_000 | G_301 S_000 )^0_{t} = z * ( G_121 S_000 | F_300 S_000 )^0 + ( F_120 S_000 | F_300 S_000 )^0 - ( H_122 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[107] = etfac[2] * AUX_INT__g_s_f_s[70] + 1 * one_over_2q * AUX_INT__f_s_f_s[30] - p_over_q * AUX_INT__h_s_f_s[120];

                    // ( G_121 S_000 | G_220 S_000 )^0_{t} = y * ( G_121 S_000 | F_210 S_000 )^0 + ( F_111 S_000 | F_210 S_000 )^0 + ( G_121 S_000 | D_200 S_000 )^0 - ( H_131 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[108] = etfac[1] * AUX_INT__g_s_f_s[71] + 2 * one_over_2q * AUX_INT__f_s_f_s[41] + 1 * one_over_2q * AUX_INT__g_s_d_s[42] - p_over_q * AUX_INT__h_s_f_s[111];

                    // ( G_121 S_000 | G_211 S_000 )^0_{t} = z * ( G_121 S_000 | F_210 S_000 )^0 + ( F_120 S_000 | F_210 S_000 )^0 - ( H_122 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[109] = etfac[2] * AUX_INT__g_s_f_s[71] + 1 * one_over_2q * AUX_INT__f_s_f_s[31] - p_over_q * AUX_INT__h_s_f_s[121];

                    // ( G_121 S_000 | G_202 S_000 )^0_{t} = z * ( G_121 S_000 | F_201 S_000 )^0 + ( F_120 S_000 | F_201 S_000 )^0 + ( G_121 S_000 | D_200 S_000 )^0 - ( H_122 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[110] = etfac[2] * AUX_INT__g_s_f_s[72] + 1 * one_over_2q * AUX_INT__f_s_f_s[32] + 1 * one_over_2q * AUX_INT__g_s_d_s[42] - p_over_q * AUX_INT__h_s_f_s[122];

                    // ( G_121 S_000 | G_130 S_000 )^0_{t} = x * ( G_121 S_000 | F_030 S_000 )^0 + ( F_021 S_000 | F_030 S_000 )^0 - ( H_221 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[111] = etfac[0] * AUX_INT__g_s_f_s[76] + 1 * one_over_2q * AUX_INT__f_s_f_s[76] - p_over_q * AUX_INT__h_s_f_s[76];

                    // ( G_121 S_000 | G_121 S_000 )^0_{t} = z * ( G_121 S_000 | F_120 S_000 )^0 + ( F_120 S_000 | F_120 S_000 )^0 - ( H_122 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[112] = etfac[2] * AUX_INT__g_s_f_s[73] + 1 * one_over_2q * AUX_INT__f_s_f_s[33] - p_over_q * AUX_INT__h_s_f_s[123];

                    // ( G_121 S_000 | G_112 S_000 )^0_{t} = y * ( G_121 S_000 | F_102 S_000 )^0 + ( F_111 S_000 | F_102 S_000 )^0 - ( H_131 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[113] = etfac[1] * AUX_INT__g_s_f_s[75] + 2 * one_over_2q * AUX_INT__f_s_f_s[45] - p_over_q * AUX_INT__h_s_f_s[115];

                    // ( G_121 S_000 | G_103 S_000 )^0_{t} = x * ( G_121 S_000 | F_003 S_000 )^0 + ( F_021 S_000 | F_003 S_000 )^0 - ( H_221 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[114] = etfac[0] * AUX_INT__g_s_f_s[79] + 1 * one_over_2q * AUX_INT__f_s_f_s[79] - p_over_q * AUX_INT__h_s_f_s[79];

                    // ( G_121 S_000 | G_040 S_000 )^0_{t} = y * ( G_121 S_000 | F_030 S_000 )^0 + ( F_111 S_000 | F_030 S_000 )^0 + ( G_121 S_000 | D_020 S_000 )^0 - ( H_131 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[115] = etfac[1] * AUX_INT__g_s_f_s[76] + 2 * one_over_2q * AUX_INT__f_s_f_s[46] + 3 * one_over_2q * AUX_INT__g_s_d_s[45] - p_over_q * AUX_INT__h_s_f_s[116];

                    // ( G_121 S_000 | G_031 S_000 )^0_{t} = z * ( G_121 S_000 | F_030 S_000 )^0 + ( F_120 S_000 | F_030 S_000 )^0 - ( H_122 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[116] = etfac[2] * AUX_INT__g_s_f_s[76] + 1 * one_over_2q * AUX_INT__f_s_f_s[36] - p_over_q * AUX_INT__h_s_f_s[126];

                    // ( G_121 S_000 | G_022 S_000 )^0_{t} = z * ( G_121 S_000 | F_021 S_000 )^0 + ( F_120 S_000 | F_021 S_000 )^0 + ( G_121 S_000 | D_020 S_000 )^0 - ( H_122 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[117] = etfac[2] * AUX_INT__g_s_f_s[77] + 1 * one_over_2q * AUX_INT__f_s_f_s[37] + 1 * one_over_2q * AUX_INT__g_s_d_s[45] - p_over_q * AUX_INT__h_s_f_s[127];

                    // ( G_121 S_000 | G_013 S_000 )^0_{t} = y * ( G_121 S_000 | F_003 S_000 )^0 + ( F_111 S_000 | F_003 S_000 )^0 - ( H_131 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[118] = etfac[1] * AUX_INT__g_s_f_s[79] + 2 * one_over_2q * AUX_INT__f_s_f_s[49] - p_over_q * AUX_INT__h_s_f_s[119];

                    // ( G_121 S_000 | G_004 S_000 )^0_{t} = z * ( G_121 S_000 | F_003 S_000 )^0 + ( F_120 S_000 | F_003 S_000 )^0 + ( G_121 S_000 | D_002 S_000 )^0 - ( H_122 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[119] = etfac[2] * AUX_INT__g_s_f_s[79] + 1 * one_over_2q * AUX_INT__f_s_f_s[39] + 3 * one_over_2q * AUX_INT__g_s_d_s[47] - p_over_q * AUX_INT__h_s_f_s[129];

                    // ( G_112 S_000 | G_400 S_000 )^0_{t} = x * ( G_112 S_000 | F_300 S_000 )^0 + ( F_012 S_000 | F_300 S_000 )^0 + ( G_112 S_000 | D_200 S_000 )^0 - ( H_212 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[120] = etfac[0] * AUX_INT__g_s_f_s[80] + 1 * one_over_2q * AUX_INT__f_s_f_s[80] + 3 * one_over_2q * AUX_INT__g_s_d_s[48] - p_over_q * AUX_INT__h_s_f_s[80];

                    // ( G_112 S_000 | G_310 S_000 )^0_{t} = y * ( G_112 S_000 | F_300 S_000 )^0 + ( F_102 S_000 | F_300 S_000 )^0 - ( H_122 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[121] = etfac[1] * AUX_INT__g_s_f_s[80] + 1 * one_over_2q * AUX_INT__f_s_f_s[50] - p_over_q * AUX_INT__h_s_f_s[120];

                    // ( G_112 S_000 | G_301 S_000 )^0_{t} = z * ( G_112 S_000 | F_300 S_000 )^0 + ( F_111 S_000 | F_300 S_000 )^0 - ( H_113 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[122] = etfac[2] * AUX_INT__g_s_f_s[80] + 2 * one_over_2q * AUX_INT__f_s_f_s[40] - p_over_q * AUX_INT__h_s_f_s[130];

                    // ( G_112 S_000 | G_220 S_000 )^0_{t} = y * ( G_112 S_000 | F_210 S_000 )^0 + ( F_102 S_000 | F_210 S_000 )^0 + ( G_112 S_000 | D_200 S_000 )^0 - ( H_122 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[123] = etfac[1] * AUX_INT__g_s_f_s[81] + 1 * one_over_2q * AUX_INT__f_s_f_s[51] + 1 * one_over_2q * AUX_INT__g_s_d_s[48] - p_over_q * AUX_INT__h_s_f_s[121];

                    // ( G_112 S_000 | G_211 S_000 )^0_{t} = z * ( G_112 S_000 | F_210 S_000 )^0 + ( F_111 S_000 | F_210 S_000 )^0 - ( H_113 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[124] = etfac[2] * AUX_INT__g_s_f_s[81] + 2 * one_over_2q * AUX_INT__f_s_f_s[41] - p_over_q * AUX_INT__h_s_f_s[131];

                    // ( G_112 S_000 | G_202 S_000 )^0_{t} = z * ( G_112 S_000 | F_201 S_000 )^0 + ( F_111 S_000 | F_201 S_000 )^0 + ( G_112 S_000 | D_200 S_000 )^0 - ( H_113 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[125] = etfac[2] * AUX_INT__g_s_f_s[82] + 2 * one_over_2q * AUX_INT__f_s_f_s[42] + 1 * one_over_2q * AUX_INT__g_s_d_s[48] - p_over_q * AUX_INT__h_s_f_s[132];

                    // ( G_112 S_000 | G_130 S_000 )^0_{t} = x * ( G_112 S_000 | F_030 S_000 )^0 + ( F_012 S_000 | F_030 S_000 )^0 - ( H_212 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[126] = etfac[0] * AUX_INT__g_s_f_s[86] + 1 * one_over_2q * AUX_INT__f_s_f_s[86] - p_over_q * AUX_INT__h_s_f_s[86];

                    // ( G_112 S_000 | G_121 S_000 )^0_{t} = z * ( G_112 S_000 | F_120 S_000 )^0 + ( F_111 S_000 | F_120 S_000 )^0 - ( H_113 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[127] = etfac[2] * AUX_INT__g_s_f_s[83] + 2 * one_over_2q * AUX_INT__f_s_f_s[43] - p_over_q * AUX_INT__h_s_f_s[133];

                    // ( G_112 S_000 | G_112 S_000 )^0_{t} = y * ( G_112 S_000 | F_102 S_000 )^0 + ( F_102 S_000 | F_102 S_000 )^0 - ( H_122 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[128] = etfac[1] * AUX_INT__g_s_f_s[85] + 1 * one_over_2q * AUX_INT__f_s_f_s[55] - p_over_q * AUX_INT__h_s_f_s[125];

                    // ( G_112 S_000 | G_103 S_000 )^0_{t} = x * ( G_112 S_000 | F_003 S_000 )^0 + ( F_012 S_000 | F_003 S_000 )^0 - ( H_212 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[129] = etfac[0] * AUX_INT__g_s_f_s[89] + 1 * one_over_2q * AUX_INT__f_s_f_s[89] - p_over_q * AUX_INT__h_s_f_s[89];

                    // ( G_112 S_000 | G_040 S_000 )^0_{t} = y * ( G_112 S_000 | F_030 S_000 )^0 + ( F_102 S_000 | F_030 S_000 )^0 + ( G_112 S_000 | D_020 S_000 )^0 - ( H_122 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[130] = etfac[1] * AUX_INT__g_s_f_s[86] + 1 * one_over_2q * AUX_INT__f_s_f_s[56] + 3 * one_over_2q * AUX_INT__g_s_d_s[51] - p_over_q * AUX_INT__h_s_f_s[126];

                    // ( G_112 S_000 | G_031 S_000 )^0_{t} = z * ( G_112 S_000 | F_030 S_000 )^0 + ( F_111 S_000 | F_030 S_000 )^0 - ( H_113 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[131] = etfac[2] * AUX_INT__g_s_f_s[86] + 2 * one_over_2q * AUX_INT__f_s_f_s[46] - p_over_q * AUX_INT__h_s_f_s[136];

                    // ( G_112 S_000 | G_022 S_000 )^0_{t} = z * ( G_112 S_000 | F_021 S_000 )^0 + ( F_111 S_000 | F_021 S_000 )^0 + ( G_112 S_000 | D_020 S_000 )^0 - ( H_113 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[132] = etfac[2] * AUX_INT__g_s_f_s[87] + 2 * one_over_2q * AUX_INT__f_s_f_s[47] + 1 * one_over_2q * AUX_INT__g_s_d_s[51] - p_over_q * AUX_INT__h_s_f_s[137];

                    // ( G_112 S_000 | G_013 S_000 )^0_{t} = y * ( G_112 S_000 | F_003 S_000 )^0 + ( F_102 S_000 | F_003 S_000 )^0 - ( H_122 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[133] = etfac[1] * AUX_INT__g_s_f_s[89] + 1 * one_over_2q * AUX_INT__f_s_f_s[59] - p_over_q * AUX_INT__h_s_f_s[129];

                    // ( G_112 S_000 | G_004 S_000 )^0_{t} = z * ( G_112 S_000 | F_003 S_000 )^0 + ( F_111 S_000 | F_003 S_000 )^0 + ( G_112 S_000 | D_002 S_000 )^0 - ( H_113 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[134] = etfac[2] * AUX_INT__g_s_f_s[89] + 2 * one_over_2q * AUX_INT__f_s_f_s[49] + 3 * one_over_2q * AUX_INT__g_s_d_s[53] - p_over_q * AUX_INT__h_s_f_s[139];

                    // ( G_103 S_000 | G_400 S_000 )^0_{t} = x * ( G_103 S_000 | F_300 S_000 )^0 + ( F_003 S_000 | F_300 S_000 )^0 + ( G_103 S_000 | D_200 S_000 )^0 - ( H_203 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[135] = etfac[0] * AUX_INT__g_s_f_s[90] + 1 * one_over_2q * AUX_INT__f_s_f_s[90] + 3 * one_over_2q * AUX_INT__g_s_d_s[54] - p_over_q * AUX_INT__h_s_f_s[90];

                    // ( G_103 S_000 | G_310 S_000 )^0_{t} = y * ( G_103 S_000 | F_300 S_000 )^0 - ( H_113 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[136] = etfac[1] * AUX_INT__g_s_f_s[90] - p_over_q * AUX_INT__h_s_f_s[130];

                    // ( G_103 S_000 | G_301 S_000 )^0_{t} = z * ( G_103 S_000 | F_300 S_000 )^0 + ( F_102 S_000 | F_300 S_000 )^0 - ( H_104 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[137] = etfac[2] * AUX_INT__g_s_f_s[90] + 3 * one_over_2q * AUX_INT__f_s_f_s[50] - p_over_q * AUX_INT__h_s_f_s[140];

                    // ( G_103 S_000 | G_220 S_000 )^0_{t} = y * ( G_103 S_000 | F_210 S_000 )^0 + ( G_103 S_000 | D_200 S_000 )^0 - ( H_113 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[138] = etfac[1] * AUX_INT__g_s_f_s[91] + 1 * one_over_2q * AUX_INT__g_s_d_s[54] - p_over_q * AUX_INT__h_s_f_s[131];

                    // ( G_103 S_000 | G_211 S_000 )^0_{t} = z * ( G_103 S_000 | F_210 S_000 )^0 + ( F_102 S_000 | F_210 S_000 )^0 - ( H_104 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[139] = etfac[2] * AUX_INT__g_s_f_s[91] + 3 * one_over_2q * AUX_INT__f_s_f_s[51] - p_over_q * AUX_INT__h_s_f_s[141];

                    // ( G_103 S_000 | G_202 S_000 )^0_{t} = z * ( G_103 S_000 | F_201 S_000 )^0 + ( F_102 S_000 | F_201 S_000 )^0 + ( G_103 S_000 | D_200 S_000 )^0 - ( H_104 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[140] = etfac[2] * AUX_INT__g_s_f_s[92] + 3 * one_over_2q * AUX_INT__f_s_f_s[52] + 1 * one_over_2q * AUX_INT__g_s_d_s[54] - p_over_q * AUX_INT__h_s_f_s[142];

                    // ( G_103 S_000 | G_130 S_000 )^0_{t} = x * ( G_103 S_000 | F_030 S_000 )^0 + ( F_003 S_000 | F_030 S_000 )^0 - ( H_203 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[141] = etfac[0] * AUX_INT__g_s_f_s[96] + 1 * one_over_2q * AUX_INT__f_s_f_s[96] - p_over_q * AUX_INT__h_s_f_s[96];

                    // ( G_103 S_000 | G_121 S_000 )^0_{t} = z * ( G_103 S_000 | F_120 S_000 )^0 + ( F_102 S_000 | F_120 S_000 )^0 - ( H_104 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[142] = etfac[2] * AUX_INT__g_s_f_s[93] + 3 * one_over_2q * AUX_INT__f_s_f_s[53] - p_over_q * AUX_INT__h_s_f_s[143];

                    // ( G_103 S_000 | G_112 S_000 )^0_{t} = y * ( G_103 S_000 | F_102 S_000 )^0 - ( H_113 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[143] = etfac[1] * AUX_INT__g_s_f_s[95] - p_over_q * AUX_INT__h_s_f_s[135];

                    // ( G_103 S_000 | G_103 S_000 )^0_{t} = x * ( G_103 S_000 | F_003 S_000 )^0 + ( F_003 S_000 | F_003 S_000 )^0 - ( H_203 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[144] = etfac[0] * AUX_INT__g_s_f_s[99] + 1 * one_over_2q * AUX_INT__f_s_f_s[99] - p_over_q * AUX_INT__h_s_f_s[99];

                    // ( G_103 S_000 | G_040 S_000 )^0_{t} = y * ( G_103 S_000 | F_030 S_000 )^0 + ( G_103 S_000 | D_020 S_000 )^0 - ( H_113 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[145] = etfac[1] * AUX_INT__g_s_f_s[96] + 3 * one_over_2q * AUX_INT__g_s_d_s[57] - p_over_q * AUX_INT__h_s_f_s[136];

                    // ( G_103 S_000 | G_031 S_000 )^0_{t} = z * ( G_103 S_000 | F_030 S_000 )^0 + ( F_102 S_000 | F_030 S_000 )^0 - ( H_104 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[146] = etfac[2] * AUX_INT__g_s_f_s[96] + 3 * one_over_2q * AUX_INT__f_s_f_s[56] - p_over_q * AUX_INT__h_s_f_s[146];

                    // ( G_103 S_000 | G_022 S_000 )^0_{t} = z * ( G_103 S_000 | F_021 S_000 )^0 + ( F_102 S_000 | F_021 S_000 )^0 + ( G_103 S_000 | D_020 S_000 )^0 - ( H_104 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[147] = etfac[2] * AUX_INT__g_s_f_s[97] + 3 * one_over_2q * AUX_INT__f_s_f_s[57] + 1 * one_over_2q * AUX_INT__g_s_d_s[57] - p_over_q * AUX_INT__h_s_f_s[147];

                    // ( G_103 S_000 | G_013 S_000 )^0_{t} = y * ( G_103 S_000 | F_003 S_000 )^0 - ( H_113 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[148] = etfac[1] * AUX_INT__g_s_f_s[99] - p_over_q * AUX_INT__h_s_f_s[139];

                    // ( G_103 S_000 | G_004 S_000 )^0_{t} = z * ( G_103 S_000 | F_003 S_000 )^0 + ( F_102 S_000 | F_003 S_000 )^0 + ( G_103 S_000 | D_002 S_000 )^0 - ( H_104 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[149] = etfac[2] * AUX_INT__g_s_f_s[99] + 3 * one_over_2q * AUX_INT__f_s_f_s[59] + 3 * one_over_2q * AUX_INT__g_s_d_s[59] - p_over_q * AUX_INT__h_s_f_s[149];

                    // ( G_040 S_000 | G_400 S_000 )^0_{t} = x * ( G_040 S_000 | F_300 S_000 )^0 + ( G_040 S_000 | D_200 S_000 )^0 - ( H_140 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[150] = etfac[0] * AUX_INT__g_s_f_s[100] + 3 * one_over_2q * AUX_INT__g_s_d_s[60] - p_over_q * AUX_INT__h_s_f_s[100];

                    // ( G_040 S_000 | G_310 S_000 )^0_{t} = y * ( G_040 S_000 | F_300 S_000 )^0 + ( F_030 S_000 | F_300 S_000 )^0 - ( H_050 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[151] = etfac[1] * AUX_INT__g_s_f_s[100] + 4 * one_over_2q * AUX_INT__f_s_f_s[60] - p_over_q * AUX_INT__h_s_f_s[150];

                    // ( G_040 S_000 | G_301 S_000 )^0_{t} = z * ( G_040 S_000 | F_300 S_000 )^0 - ( H_041 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[152] = etfac[2] * AUX_INT__g_s_f_s[100] - p_over_q * AUX_INT__h_s_f_s[160];

                    // ( G_040 S_000 | G_220 S_000 )^0_{t} = y * ( G_040 S_000 | F_210 S_000 )^0 + ( F_030 S_000 | F_210 S_000 )^0 + ( G_040 S_000 | D_200 S_000 )^0 - ( H_050 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[153] = etfac[1] * AUX_INT__g_s_f_s[101] + 4 * one_over_2q * AUX_INT__f_s_f_s[61] + 1 * one_over_2q * AUX_INT__g_s_d_s[60] - p_over_q * AUX_INT__h_s_f_s[151];

                    // ( G_040 S_000 | G_211 S_000 )^0_{t} = z * ( G_040 S_000 | F_210 S_000 )^0 - ( H_041 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[154] = etfac[2] * AUX_INT__g_s_f_s[101] - p_over_q * AUX_INT__h_s_f_s[161];

                    // ( G_040 S_000 | G_202 S_000 )^0_{t} = z * ( G_040 S_000 | F_201 S_000 )^0 + ( G_040 S_000 | D_200 S_000 )^0 - ( H_041 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[155] = etfac[2] * AUX_INT__g_s_f_s[102] + 1 * one_over_2q * AUX_INT__g_s_d_s[60] - p_over_q * AUX_INT__h_s_f_s[162];

                    // ( G_040 S_000 | G_130 S_000 )^0_{t} = x * ( G_040 S_000 | F_030 S_000 )^0 - ( H_140 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[156] = etfac[0] * AUX_INT__g_s_f_s[106] - p_over_q * AUX_INT__h_s_f_s[106];

                    // ( G_040 S_000 | G_121 S_000 )^0_{t} = z * ( G_040 S_000 | F_120 S_000 )^0 - ( H_041 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[157] = etfac[2] * AUX_INT__g_s_f_s[103] - p_over_q * AUX_INT__h_s_f_s[163];

                    // ( G_040 S_000 | G_112 S_000 )^0_{t} = y * ( G_040 S_000 | F_102 S_000 )^0 + ( F_030 S_000 | F_102 S_000 )^0 - ( H_050 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[158] = etfac[1] * AUX_INT__g_s_f_s[105] + 4 * one_over_2q * AUX_INT__f_s_f_s[65] - p_over_q * AUX_INT__h_s_f_s[155];

                    // ( G_040 S_000 | G_103 S_000 )^0_{t} = x * ( G_040 S_000 | F_003 S_000 )^0 - ( H_140 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[159] = etfac[0] * AUX_INT__g_s_f_s[109] - p_over_q * AUX_INT__h_s_f_s[109];

                    // ( G_040 S_000 | G_040 S_000 )^0_{t} = y * ( G_040 S_000 | F_030 S_000 )^0 + ( F_030 S_000 | F_030 S_000 )^0 + ( G_040 S_000 | D_020 S_000 )^0 - ( H_050 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[160] = etfac[1] * AUX_INT__g_s_f_s[106] + 4 * one_over_2q * AUX_INT__f_s_f_s[66] + 3 * one_over_2q * AUX_INT__g_s_d_s[63] - p_over_q * AUX_INT__h_s_f_s[156];

                    // ( G_040 S_000 | G_031 S_000 )^0_{t} = z * ( G_040 S_000 | F_030 S_000 )^0 - ( H_041 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[161] = etfac[2] * AUX_INT__g_s_f_s[106] - p_over_q * AUX_INT__h_s_f_s[166];

                    // ( G_040 S_000 | G_022 S_000 )^0_{t} = z * ( G_040 S_000 | F_021 S_000 )^0 + ( G_040 S_000 | D_020 S_000 )^0 - ( H_041 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[162] = etfac[2] * AUX_INT__g_s_f_s[107] + 1 * one_over_2q * AUX_INT__g_s_d_s[63] - p_over_q * AUX_INT__h_s_f_s[167];

                    // ( G_040 S_000 | G_013 S_000 )^0_{t} = y * ( G_040 S_000 | F_003 S_000 )^0 + ( F_030 S_000 | F_003 S_000 )^0 - ( H_050 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[163] = etfac[1] * AUX_INT__g_s_f_s[109] + 4 * one_over_2q * AUX_INT__f_s_f_s[69] - p_over_q * AUX_INT__h_s_f_s[159];

                    // ( G_040 S_000 | G_004 S_000 )^0_{t} = z * ( G_040 S_000 | F_003 S_000 )^0 + ( G_040 S_000 | D_002 S_000 )^0 - ( H_041 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[164] = etfac[2] * AUX_INT__g_s_f_s[109] + 3 * one_over_2q * AUX_INT__g_s_d_s[65] - p_over_q * AUX_INT__h_s_f_s[169];

                    // ( G_031 S_000 | G_400 S_000 )^0_{t} = x * ( G_031 S_000 | F_300 S_000 )^0 + ( G_031 S_000 | D_200 S_000 )^0 - ( H_131 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[165] = etfac[0] * AUX_INT__g_s_f_s[110] + 3 * one_over_2q * AUX_INT__g_s_d_s[66] - p_over_q * AUX_INT__h_s_f_s[110];

                    // ( G_031 S_000 | G_310 S_000 )^0_{t} = y * ( G_031 S_000 | F_300 S_000 )^0 + ( F_021 S_000 | F_300 S_000 )^0 - ( H_041 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[166] = etfac[1] * AUX_INT__g_s_f_s[110] + 3 * one_over_2q * AUX_INT__f_s_f_s[70] - p_over_q * AUX_INT__h_s_f_s[160];

                    // ( G_031 S_000 | G_301 S_000 )^0_{t} = z * ( G_031 S_000 | F_300 S_000 )^0 + ( F_030 S_000 | F_300 S_000 )^0 - ( H_032 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[167] = etfac[2] * AUX_INT__g_s_f_s[110] + 1 * one_over_2q * AUX_INT__f_s_f_s[60] - p_over_q * AUX_INT__h_s_f_s[170];

                    // ( G_031 S_000 | G_220 S_000 )^0_{t} = y * ( G_031 S_000 | F_210 S_000 )^0 + ( F_021 S_000 | F_210 S_000 )^0 + ( G_031 S_000 | D_200 S_000 )^0 - ( H_041 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[168] = etfac[1] * AUX_INT__g_s_f_s[111] + 3 * one_over_2q * AUX_INT__f_s_f_s[71] + 1 * one_over_2q * AUX_INT__g_s_d_s[66] - p_over_q * AUX_INT__h_s_f_s[161];

                    // ( G_031 S_000 | G_211 S_000 )^0_{t} = z * ( G_031 S_000 | F_210 S_000 )^0 + ( F_030 S_000 | F_210 S_000 )^0 - ( H_032 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[169] = etfac[2] * AUX_INT__g_s_f_s[111] + 1 * one_over_2q * AUX_INT__f_s_f_s[61] - p_over_q * AUX_INT__h_s_f_s[171];

                    // ( G_031 S_000 | G_202 S_000 )^0_{t} = z * ( G_031 S_000 | F_201 S_000 )^0 + ( F_030 S_000 | F_201 S_000 )^0 + ( G_031 S_000 | D_200 S_000 )^0 - ( H_032 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[170] = etfac[2] * AUX_INT__g_s_f_s[112] + 1 * one_over_2q * AUX_INT__f_s_f_s[62] + 1 * one_over_2q * AUX_INT__g_s_d_s[66] - p_over_q * AUX_INT__h_s_f_s[172];

                    // ( G_031 S_000 | G_130 S_000 )^0_{t} = x * ( G_031 S_000 | F_030 S_000 )^0 - ( H_131 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[171] = etfac[0] * AUX_INT__g_s_f_s[116] - p_over_q * AUX_INT__h_s_f_s[116];

                    // ( G_031 S_000 | G_121 S_000 )^0_{t} = z * ( G_031 S_000 | F_120 S_000 )^0 + ( F_030 S_000 | F_120 S_000 )^0 - ( H_032 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[172] = etfac[2] * AUX_INT__g_s_f_s[113] + 1 * one_over_2q * AUX_INT__f_s_f_s[63] - p_over_q * AUX_INT__h_s_f_s[173];

                    // ( G_031 S_000 | G_112 S_000 )^0_{t} = y * ( G_031 S_000 | F_102 S_000 )^0 + ( F_021 S_000 | F_102 S_000 )^0 - ( H_041 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[173] = etfac[1] * AUX_INT__g_s_f_s[115] + 3 * one_over_2q * AUX_INT__f_s_f_s[75] - p_over_q * AUX_INT__h_s_f_s[165];

                    // ( G_031 S_000 | G_103 S_000 )^0_{t} = x * ( G_031 S_000 | F_003 S_000 )^0 - ( H_131 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[174] = etfac[0] * AUX_INT__g_s_f_s[119] - p_over_q * AUX_INT__h_s_f_s[119];

                    // ( G_031 S_000 | G_040 S_000 )^0_{t} = y * ( G_031 S_000 | F_030 S_000 )^0 + ( F_021 S_000 | F_030 S_000 )^0 + ( G_031 S_000 | D_020 S_000 )^0 - ( H_041 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[175] = etfac[1] * AUX_INT__g_s_f_s[116] + 3 * one_over_2q * AUX_INT__f_s_f_s[76] + 3 * one_over_2q * AUX_INT__g_s_d_s[69] - p_over_q * AUX_INT__h_s_f_s[166];

                    // ( G_031 S_000 | G_031 S_000 )^0_{t} = z * ( G_031 S_000 | F_030 S_000 )^0 + ( F_030 S_000 | F_030 S_000 )^0 - ( H_032 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[176] = etfac[2] * AUX_INT__g_s_f_s[116] + 1 * one_over_2q * AUX_INT__f_s_f_s[66] - p_over_q * AUX_INT__h_s_f_s[176];

                    // ( G_031 S_000 | G_022 S_000 )^0_{t} = z * ( G_031 S_000 | F_021 S_000 )^0 + ( F_030 S_000 | F_021 S_000 )^0 + ( G_031 S_000 | D_020 S_000 )^0 - ( H_032 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[177] = etfac[2] * AUX_INT__g_s_f_s[117] + 1 * one_over_2q * AUX_INT__f_s_f_s[67] + 1 * one_over_2q * AUX_INT__g_s_d_s[69] - p_over_q * AUX_INT__h_s_f_s[177];

                    // ( G_031 S_000 | G_013 S_000 )^0_{t} = y * ( G_031 S_000 | F_003 S_000 )^0 + ( F_021 S_000 | F_003 S_000 )^0 - ( H_041 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[178] = etfac[1] * AUX_INT__g_s_f_s[119] + 3 * one_over_2q * AUX_INT__f_s_f_s[79] - p_over_q * AUX_INT__h_s_f_s[169];

                    // ( G_031 S_000 | G_004 S_000 )^0_{t} = z * ( G_031 S_000 | F_003 S_000 )^0 + ( F_030 S_000 | F_003 S_000 )^0 + ( G_031 S_000 | D_002 S_000 )^0 - ( H_032 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[179] = etfac[2] * AUX_INT__g_s_f_s[119] + 1 * one_over_2q * AUX_INT__f_s_f_s[69] + 3 * one_over_2q * AUX_INT__g_s_d_s[71] - p_over_q * AUX_INT__h_s_f_s[179];

                    // ( G_022 S_000 | G_400 S_000 )^0_{t} = x * ( G_022 S_000 | F_300 S_000 )^0 + ( G_022 S_000 | D_200 S_000 )^0 - ( H_122 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[180] = etfac[0] * AUX_INT__g_s_f_s[120] + 3 * one_over_2q * AUX_INT__g_s_d_s[72] - p_over_q * AUX_INT__h_s_f_s[120];

                    // ( G_022 S_000 | G_310 S_000 )^0_{t} = y * ( G_022 S_000 | F_300 S_000 )^0 + ( F_012 S_000 | F_300 S_000 )^0 - ( H_032 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[181] = etfac[1] * AUX_INT__g_s_f_s[120] + 2 * one_over_2q * AUX_INT__f_s_f_s[80] - p_over_q * AUX_INT__h_s_f_s[170];

                    // ( G_022 S_000 | G_301 S_000 )^0_{t} = z * ( G_022 S_000 | F_300 S_000 )^0 + ( F_021 S_000 | F_300 S_000 )^0 - ( H_023 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[182] = etfac[2] * AUX_INT__g_s_f_s[120] + 2 * one_over_2q * AUX_INT__f_s_f_s[70] - p_over_q * AUX_INT__h_s_f_s[180];

                    // ( G_022 S_000 | G_220 S_000 )^0_{t} = y * ( G_022 S_000 | F_210 S_000 )^0 + ( F_012 S_000 | F_210 S_000 )^0 + ( G_022 S_000 | D_200 S_000 )^0 - ( H_032 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[183] = etfac[1] * AUX_INT__g_s_f_s[121] + 2 * one_over_2q * AUX_INT__f_s_f_s[81] + 1 * one_over_2q * AUX_INT__g_s_d_s[72] - p_over_q * AUX_INT__h_s_f_s[171];

                    // ( G_022 S_000 | G_211 S_000 )^0_{t} = z * ( G_022 S_000 | F_210 S_000 )^0 + ( F_021 S_000 | F_210 S_000 )^0 - ( H_023 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[184] = etfac[2] * AUX_INT__g_s_f_s[121] + 2 * one_over_2q * AUX_INT__f_s_f_s[71] - p_over_q * AUX_INT__h_s_f_s[181];

                    // ( G_022 S_000 | G_202 S_000 )^0_{t} = z * ( G_022 S_000 | F_201 S_000 )^0 + ( F_021 S_000 | F_201 S_000 )^0 + ( G_022 S_000 | D_200 S_000 )^0 - ( H_023 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[185] = etfac[2] * AUX_INT__g_s_f_s[122] + 2 * one_over_2q * AUX_INT__f_s_f_s[72] + 1 * one_over_2q * AUX_INT__g_s_d_s[72] - p_over_q * AUX_INT__h_s_f_s[182];

                    // ( G_022 S_000 | G_130 S_000 )^0_{t} = x * ( G_022 S_000 | F_030 S_000 )^0 - ( H_122 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[186] = etfac[0] * AUX_INT__g_s_f_s[126] - p_over_q * AUX_INT__h_s_f_s[126];

                    // ( G_022 S_000 | G_121 S_000 )^0_{t} = z * ( G_022 S_000 | F_120 S_000 )^0 + ( F_021 S_000 | F_120 S_000 )^0 - ( H_023 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[187] = etfac[2] * AUX_INT__g_s_f_s[123] + 2 * one_over_2q * AUX_INT__f_s_f_s[73] - p_over_q * AUX_INT__h_s_f_s[183];

                    // ( G_022 S_000 | G_112 S_000 )^0_{t} = y * ( G_022 S_000 | F_102 S_000 )^0 + ( F_012 S_000 | F_102 S_000 )^0 - ( H_032 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[188] = etfac[1] * AUX_INT__g_s_f_s[125] + 2 * one_over_2q * AUX_INT__f_s_f_s[85] - p_over_q * AUX_INT__h_s_f_s[175];

                    // ( G_022 S_000 | G_103 S_000 )^0_{t} = x * ( G_022 S_000 | F_003 S_000 )^0 - ( H_122 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[189] = etfac[0] * AUX_INT__g_s_f_s[129] - p_over_q * AUX_INT__h_s_f_s[129];

                    // ( G_022 S_000 | G_040 S_000 )^0_{t} = y * ( G_022 S_000 | F_030 S_000 )^0 + ( F_012 S_000 | F_030 S_000 )^0 + ( G_022 S_000 | D_020 S_000 )^0 - ( H_032 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[190] = etfac[1] * AUX_INT__g_s_f_s[126] + 2 * one_over_2q * AUX_INT__f_s_f_s[86] + 3 * one_over_2q * AUX_INT__g_s_d_s[75] - p_over_q * AUX_INT__h_s_f_s[176];

                    // ( G_022 S_000 | G_031 S_000 )^0_{t} = z * ( G_022 S_000 | F_030 S_000 )^0 + ( F_021 S_000 | F_030 S_000 )^0 - ( H_023 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[191] = etfac[2] * AUX_INT__g_s_f_s[126] + 2 * one_over_2q * AUX_INT__f_s_f_s[76] - p_over_q * AUX_INT__h_s_f_s[186];

                    // ( G_022 S_000 | G_022 S_000 )^0_{t} = z * ( G_022 S_000 | F_021 S_000 )^0 + ( F_021 S_000 | F_021 S_000 )^0 + ( G_022 S_000 | D_020 S_000 )^0 - ( H_023 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[192] = etfac[2] * AUX_INT__g_s_f_s[127] + 2 * one_over_2q * AUX_INT__f_s_f_s[77] + 1 * one_over_2q * AUX_INT__g_s_d_s[75] - p_over_q * AUX_INT__h_s_f_s[187];

                    // ( G_022 S_000 | G_013 S_000 )^0_{t} = y * ( G_022 S_000 | F_003 S_000 )^0 + ( F_012 S_000 | F_003 S_000 )^0 - ( H_032 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[193] = etfac[1] * AUX_INT__g_s_f_s[129] + 2 * one_over_2q * AUX_INT__f_s_f_s[89] - p_over_q * AUX_INT__h_s_f_s[179];

                    // ( G_022 S_000 | G_004 S_000 )^0_{t} = z * ( G_022 S_000 | F_003 S_000 )^0 + ( F_021 S_000 | F_003 S_000 )^0 + ( G_022 S_000 | D_002 S_000 )^0 - ( H_023 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[194] = etfac[2] * AUX_INT__g_s_f_s[129] + 2 * one_over_2q * AUX_INT__f_s_f_s[79] + 3 * one_over_2q * AUX_INT__g_s_d_s[77] - p_over_q * AUX_INT__h_s_f_s[189];

                    // ( G_013 S_000 | G_400 S_000 )^0_{t} = x * ( G_013 S_000 | F_300 S_000 )^0 + ( G_013 S_000 | D_200 S_000 )^0 - ( H_113 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[195] = etfac[0] * AUX_INT__g_s_f_s[130] + 3 * one_over_2q * AUX_INT__g_s_d_s[78] - p_over_q * AUX_INT__h_s_f_s[130];

                    // ( G_013 S_000 | G_310 S_000 )^0_{t} = y * ( G_013 S_000 | F_300 S_000 )^0 + ( F_003 S_000 | F_300 S_000 )^0 - ( H_023 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[196] = etfac[1] * AUX_INT__g_s_f_s[130] + 1 * one_over_2q * AUX_INT__f_s_f_s[90] - p_over_q * AUX_INT__h_s_f_s[180];

                    // ( G_013 S_000 | G_301 S_000 )^0_{t} = z * ( G_013 S_000 | F_300 S_000 )^0 + ( F_012 S_000 | F_300 S_000 )^0 - ( H_014 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[197] = etfac[2] * AUX_INT__g_s_f_s[130] + 3 * one_over_2q * AUX_INT__f_s_f_s[80] - p_over_q * AUX_INT__h_s_f_s[190];

                    // ( G_013 S_000 | G_220 S_000 )^0_{t} = y * ( G_013 S_000 | F_210 S_000 )^0 + ( F_003 S_000 | F_210 S_000 )^0 + ( G_013 S_000 | D_200 S_000 )^0 - ( H_023 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[198] = etfac[1] * AUX_INT__g_s_f_s[131] + 1 * one_over_2q * AUX_INT__f_s_f_s[91] + 1 * one_over_2q * AUX_INT__g_s_d_s[78] - p_over_q * AUX_INT__h_s_f_s[181];

                    // ( G_013 S_000 | G_211 S_000 )^0_{t} = z * ( G_013 S_000 | F_210 S_000 )^0 + ( F_012 S_000 | F_210 S_000 )^0 - ( H_014 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[199] = etfac[2] * AUX_INT__g_s_f_s[131] + 3 * one_over_2q * AUX_INT__f_s_f_s[81] - p_over_q * AUX_INT__h_s_f_s[191];

                    // ( G_013 S_000 | G_202 S_000 )^0_{t} = z * ( G_013 S_000 | F_201 S_000 )^0 + ( F_012 S_000 | F_201 S_000 )^0 + ( G_013 S_000 | D_200 S_000 )^0 - ( H_014 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[200] = etfac[2] * AUX_INT__g_s_f_s[132] + 3 * one_over_2q * AUX_INT__f_s_f_s[82] + 1 * one_over_2q * AUX_INT__g_s_d_s[78] - p_over_q * AUX_INT__h_s_f_s[192];

                    // ( G_013 S_000 | G_130 S_000 )^0_{t} = x * ( G_013 S_000 | F_030 S_000 )^0 - ( H_113 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[201] = etfac[0] * AUX_INT__g_s_f_s[136] - p_over_q * AUX_INT__h_s_f_s[136];

                    // ( G_013 S_000 | G_121 S_000 )^0_{t} = z * ( G_013 S_000 | F_120 S_000 )^0 + ( F_012 S_000 | F_120 S_000 )^0 - ( H_014 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[202] = etfac[2] * AUX_INT__g_s_f_s[133] + 3 * one_over_2q * AUX_INT__f_s_f_s[83] - p_over_q * AUX_INT__h_s_f_s[193];

                    // ( G_013 S_000 | G_112 S_000 )^0_{t} = y * ( G_013 S_000 | F_102 S_000 )^0 + ( F_003 S_000 | F_102 S_000 )^0 - ( H_023 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[203] = etfac[1] * AUX_INT__g_s_f_s[135] + 1 * one_over_2q * AUX_INT__f_s_f_s[95] - p_over_q * AUX_INT__h_s_f_s[185];

                    // ( G_013 S_000 | G_103 S_000 )^0_{t} = x * ( G_013 S_000 | F_003 S_000 )^0 - ( H_113 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[204] = etfac[0] * AUX_INT__g_s_f_s[139] - p_over_q * AUX_INT__h_s_f_s[139];

                    // ( G_013 S_000 | G_040 S_000 )^0_{t} = y * ( G_013 S_000 | F_030 S_000 )^0 + ( F_003 S_000 | F_030 S_000 )^0 + ( G_013 S_000 | D_020 S_000 )^0 - ( H_023 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[205] = etfac[1] * AUX_INT__g_s_f_s[136] + 1 * one_over_2q * AUX_INT__f_s_f_s[96] + 3 * one_over_2q * AUX_INT__g_s_d_s[81] - p_over_q * AUX_INT__h_s_f_s[186];

                    // ( G_013 S_000 | G_031 S_000 )^0_{t} = z * ( G_013 S_000 | F_030 S_000 )^0 + ( F_012 S_000 | F_030 S_000 )^0 - ( H_014 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[206] = etfac[2] * AUX_INT__g_s_f_s[136] + 3 * one_over_2q * AUX_INT__f_s_f_s[86] - p_over_q * AUX_INT__h_s_f_s[196];

                    // ( G_013 S_000 | G_022 S_000 )^0_{t} = z * ( G_013 S_000 | F_021 S_000 )^0 + ( F_012 S_000 | F_021 S_000 )^0 + ( G_013 S_000 | D_020 S_000 )^0 - ( H_014 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[207] = etfac[2] * AUX_INT__g_s_f_s[137] + 3 * one_over_2q * AUX_INT__f_s_f_s[87] + 1 * one_over_2q * AUX_INT__g_s_d_s[81] - p_over_q * AUX_INT__h_s_f_s[197];

                    // ( G_013 S_000 | G_013 S_000 )^0_{t} = y * ( G_013 S_000 | F_003 S_000 )^0 + ( F_003 S_000 | F_003 S_000 )^0 - ( H_023 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[208] = etfac[1] * AUX_INT__g_s_f_s[139] + 1 * one_over_2q * AUX_INT__f_s_f_s[99] - p_over_q * AUX_INT__h_s_f_s[189];

                    // ( G_013 S_000 | G_004 S_000 )^0_{t} = z * ( G_013 S_000 | F_003 S_000 )^0 + ( F_012 S_000 | F_003 S_000 )^0 + ( G_013 S_000 | D_002 S_000 )^0 - ( H_014 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[209] = etfac[2] * AUX_INT__g_s_f_s[139] + 3 * one_over_2q * AUX_INT__f_s_f_s[89] + 3 * one_over_2q * AUX_INT__g_s_d_s[83] - p_over_q * AUX_INT__h_s_f_s[199];

                    // ( G_004 S_000 | G_400 S_000 )^0_{t} = x * ( G_004 S_000 | F_300 S_000 )^0 + ( G_004 S_000 | D_200 S_000 )^0 - ( H_104 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[210] = etfac[0] * AUX_INT__g_s_f_s[140] + 3 * one_over_2q * AUX_INT__g_s_d_s[84] - p_over_q * AUX_INT__h_s_f_s[140];

                    // ( G_004 S_000 | G_310 S_000 )^0_{t} = y * ( G_004 S_000 | F_300 S_000 )^0 - ( H_014 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[211] = etfac[1] * AUX_INT__g_s_f_s[140] - p_over_q * AUX_INT__h_s_f_s[190];

                    // ( G_004 S_000 | G_301 S_000 )^0_{t} = z * ( G_004 S_000 | F_300 S_000 )^0 + ( F_003 S_000 | F_300 S_000 )^0 - ( H_005 S_000 | F_300 S_000 )^0
                    AUX_INT__g_s_g_s[212] = etfac[2] * AUX_INT__g_s_f_s[140] + 4 * one_over_2q * AUX_INT__f_s_f_s[90] - p_over_q * AUX_INT__h_s_f_s[200];

                    // ( G_004 S_000 | G_220 S_000 )^0_{t} = y * ( G_004 S_000 | F_210 S_000 )^0 + ( G_004 S_000 | D_200 S_000 )^0 - ( H_014 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[213] = etfac[1] * AUX_INT__g_s_f_s[141] + 1 * one_over_2q * AUX_INT__g_s_d_s[84] - p_over_q * AUX_INT__h_s_f_s[191];

                    // ( G_004 S_000 | G_211 S_000 )^0_{t} = z * ( G_004 S_000 | F_210 S_000 )^0 + ( F_003 S_000 | F_210 S_000 )^0 - ( H_005 S_000 | F_210 S_000 )^0
                    AUX_INT__g_s_g_s[214] = etfac[2] * AUX_INT__g_s_f_s[141] + 4 * one_over_2q * AUX_INT__f_s_f_s[91] - p_over_q * AUX_INT__h_s_f_s[201];

                    // ( G_004 S_000 | G_202 S_000 )^0_{t} = z * ( G_004 S_000 | F_201 S_000 )^0 + ( F_003 S_000 | F_201 S_000 )^0 + ( G_004 S_000 | D_200 S_000 )^0 - ( H_005 S_000 | F_201 S_000 )^0
                    AUX_INT__g_s_g_s[215] = etfac[2] * AUX_INT__g_s_f_s[142] + 4 * one_over_2q * AUX_INT__f_s_f_s[92] + 1 * one_over_2q * AUX_INT__g_s_d_s[84] - p_over_q * AUX_INT__h_s_f_s[202];

                    // ( G_004 S_000 | G_130 S_000 )^0_{t} = x * ( G_004 S_000 | F_030 S_000 )^0 - ( H_104 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[216] = etfac[0] * AUX_INT__g_s_f_s[146] - p_over_q * AUX_INT__h_s_f_s[146];

                    // ( G_004 S_000 | G_121 S_000 )^0_{t} = z * ( G_004 S_000 | F_120 S_000 )^0 + ( F_003 S_000 | F_120 S_000 )^0 - ( H_005 S_000 | F_120 S_000 )^0
                    AUX_INT__g_s_g_s[217] = etfac[2] * AUX_INT__g_s_f_s[143] + 4 * one_over_2q * AUX_INT__f_s_f_s[93] - p_over_q * AUX_INT__h_s_f_s[203];

                    // ( G_004 S_000 | G_112 S_000 )^0_{t} = y * ( G_004 S_000 | F_102 S_000 )^0 - ( H_014 S_000 | F_102 S_000 )^0
                    AUX_INT__g_s_g_s[218] = etfac[1] * AUX_INT__g_s_f_s[145] - p_over_q * AUX_INT__h_s_f_s[195];

                    // ( G_004 S_000 | G_103 S_000 )^0_{t} = x * ( G_004 S_000 | F_003 S_000 )^0 - ( H_104 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[219] = etfac[0] * AUX_INT__g_s_f_s[149] - p_over_q * AUX_INT__h_s_f_s[149];

                    // ( G_004 S_000 | G_040 S_000 )^0_{t} = y * ( G_004 S_000 | F_030 S_000 )^0 + ( G_004 S_000 | D_020 S_000 )^0 - ( H_014 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[220] = etfac[1] * AUX_INT__g_s_f_s[146] + 3 * one_over_2q * AUX_INT__g_s_d_s[87] - p_over_q * AUX_INT__h_s_f_s[196];

                    // ( G_004 S_000 | G_031 S_000 )^0_{t} = z * ( G_004 S_000 | F_030 S_000 )^0 + ( F_003 S_000 | F_030 S_000 )^0 - ( H_005 S_000 | F_030 S_000 )^0
                    AUX_INT__g_s_g_s[221] = etfac[2] * AUX_INT__g_s_f_s[146] + 4 * one_over_2q * AUX_INT__f_s_f_s[96] - p_over_q * AUX_INT__h_s_f_s[206];

                    // ( G_004 S_000 | G_022 S_000 )^0_{t} = z * ( G_004 S_000 | F_021 S_000 )^0 + ( F_003 S_000 | F_021 S_000 )^0 + ( G_004 S_000 | D_020 S_000 )^0 - ( H_005 S_000 | F_021 S_000 )^0
                    AUX_INT__g_s_g_s[222] = etfac[2] * AUX_INT__g_s_f_s[147] + 4 * one_over_2q * AUX_INT__f_s_f_s[97] + 1 * one_over_2q * AUX_INT__g_s_d_s[87] - p_over_q * AUX_INT__h_s_f_s[207];

                    // ( G_004 S_000 | G_013 S_000 )^0_{t} = y * ( G_004 S_000 | F_003 S_000 )^0 - ( H_014 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[223] = etfac[1] * AUX_INT__g_s_f_s[149] - p_over_q * AUX_INT__h_s_f_s[199];

                    // ( G_004 S_000 | G_004 S_000 )^0_{t} = z * ( G_004 S_000 | F_003 S_000 )^0 + ( F_003 S_000 | F_003 S_000 )^0 + ( G_004 S_000 | D_002 S_000 )^0 - ( H_005 S_000 | F_003 S_000 )^0
                    AUX_INT__g_s_g_s[224] = etfac[2] * AUX_INT__g_s_f_s[149] + 4 * one_over_2q * AUX_INT__f_s_f_s[99] + 3 * one_over_2q * AUX_INT__g_s_d_s[89] - p_over_q * AUX_INT__h_s_f_s[209];


                    // Accumulating in contracted workspace
                    for(int n = 0; n < 36; n++)
                        PRIM_INT__d_s_d_s[n] += AUX_INT__d_s_d_s[n];

                    // Accumulating in contracted workspace
                    for(int n = 0; n < 60; n++)
                        PRIM_INT__d_s_f_s[n] += AUX_INT__d_s_f_s[n];

                    // Accumulating in contracted workspace
                    for(int n = 0; n < 90; n++)
                        PRIM_INT__d_s_g_s[n] += AUX_INT__d_s_g_s[n];

                    // Accumulating in contracted workspace
                    for(int n = 0; n < 60; n++)
                        PRIM_INT__f_s_d_s[n] += AUX_INT__f_s_d_s[n];

                    // Accumulating in contracted workspace
                    for(int n = 0; n < 100; n++)
                        PRIM_INT__f_s_f_s[n] += AUX_INT__f_s_f_s[n];

                    // Accumulating in contracted workspace
                    for(int n = 0; n < 150; n++)
                        PRIM_INT__f_s_g_s[n] += AUX_INT__f_s_g_s[n];

                    // Accumulating in contracted workspace
                    for(int n = 0; n < 90; n++)
                        PRIM_INT__g_s_d_s[n] += AUX_INT__g_s_d_s[n];

                    // Accumulating in contracted workspace
                    for(int n = 0; n < 150; n++)
                        PRIM_INT__g_s_f_s[n] += AUX_INT__g_s_f_s[n];

                    // Accumulating in contracted workspace
                    for(int n = 0; n < 225; n++)
                        PRIM_INT__g_s_g_s[n] += AUX_INT__g_s_g_s[n];

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
        // form INT__d_d_d_s
        for(int iket = 0; iket < 6; ++iket)
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

        // form INT__d_d_f_s
        for(int iket = 0; iket < 10; ++iket)
        {
            // (D_200 P_100| = (F_300 S_000|_{t} + x_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_p_s_s_f = INT__f_s_f_s[abcd * 100 + 0 * 10 + iket] + ( AB_x[abcd] * INT__d_s_f_s[abcd * 60 + 0 * 10 + iket] );

            // (D_200 P_010| = (F_210 S_000|_{t} + y_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_s_p_s_f = INT__f_s_f_s[abcd * 100 + 1 * 10 + iket] + ( AB_y[abcd] * INT__d_s_f_s[abcd * 60 + 0 * 10 + iket] );

            // (D_200 P_001| = (F_201 S_000|_{t} + z_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_s_s_p_f = INT__f_s_f_s[abcd * 100 + 2 * 10 + iket] + ( AB_z[abcd] * INT__d_s_f_s[abcd * 60 + 0 * 10 + iket] );

            // (D_110 P_100| = (F_210 S_000|_{t} + x_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_p_s_s_f = INT__f_s_f_s[abcd * 100 + 1 * 10 + iket] + ( AB_x[abcd] * INT__d_s_f_s[abcd * 60 + 1 * 10 + iket] );

            // (D_110 P_010| = (F_120 S_000|_{t} + y_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_s_p_s_f = INT__f_s_f_s[abcd * 100 + 3 * 10 + iket] + ( AB_y[abcd] * INT__d_s_f_s[abcd * 60 + 1 * 10 + iket] );

            // (D_110 P_001| = (F_111 S_000|_{t} + z_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_s_s_p_f = INT__f_s_f_s[abcd * 100 + 4 * 10 + iket] + ( AB_z[abcd] * INT__d_s_f_s[abcd * 60 + 1 * 10 + iket] );

            // (D_101 P_100| = (F_201 S_000|_{t} + x_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_p_s_s_f = INT__f_s_f_s[abcd * 100 + 2 * 10 + iket] + ( AB_x[abcd] * INT__d_s_f_s[abcd * 60 + 2 * 10 + iket] );

            // (D_101 P_010| = (F_111 S_000|_{t} + y_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_s_p_s_f = INT__f_s_f_s[abcd * 100 + 4 * 10 + iket] + ( AB_y[abcd] * INT__d_s_f_s[abcd * 60 + 2 * 10 + iket] );

            // (D_101 P_001| = (F_102 S_000|_{t} + z_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_s_s_p_f = INT__f_s_f_s[abcd * 100 + 5 * 10 + iket] + ( AB_z[abcd] * INT__d_s_f_s[abcd * 60 + 2 * 10 + iket] );

            // (D_020 P_100| = (F_120 S_000|_{t} + x_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_p_s_s_f = INT__f_s_f_s[abcd * 100 + 3 * 10 + iket] + ( AB_x[abcd] * INT__d_s_f_s[abcd * 60 + 3 * 10 + iket] );

            // (D_020 P_010| = (F_030 S_000|_{t} + y_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_s_p_s_f = INT__f_s_f_s[abcd * 100 + 6 * 10 + iket] + ( AB_y[abcd] * INT__d_s_f_s[abcd * 60 + 3 * 10 + iket] );

            // (D_020 P_001| = (F_021 S_000|_{t} + z_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_s_s_p_f = INT__f_s_f_s[abcd * 100 + 7 * 10 + iket] + ( AB_z[abcd] * INT__d_s_f_s[abcd * 60 + 3 * 10 + iket] );

            // (D_011 P_100| = (F_111 S_000|_{t} + x_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_p_s_s_f = INT__f_s_f_s[abcd * 100 + 4 * 10 + iket] + ( AB_x[abcd] * INT__d_s_f_s[abcd * 60 + 4 * 10 + iket] );

            // (D_011 P_010| = (F_021 S_000|_{t} + y_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_s_p_s_f = INT__f_s_f_s[abcd * 100 + 7 * 10 + iket] + ( AB_y[abcd] * INT__d_s_f_s[abcd * 60 + 4 * 10 + iket] );

            // (D_011 P_001| = (F_012 S_000|_{t} + z_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_s_s_p_f = INT__f_s_f_s[abcd * 100 + 8 * 10 + iket] + ( AB_z[abcd] * INT__d_s_f_s[abcd * 60 + 4 * 10 + iket] );

            // (D_002 P_100| = (F_102 S_000|_{t} + x_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_p_s_s_f = INT__f_s_f_s[abcd * 100 + 5 * 10 + iket] + ( AB_x[abcd] * INT__d_s_f_s[abcd * 60 + 5 * 10 + iket] );

            // (D_002 P_010| = (F_012 S_000|_{t} + y_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_s_p_s_f = INT__f_s_f_s[abcd * 100 + 8 * 10 + iket] + ( AB_y[abcd] * INT__d_s_f_s[abcd * 60 + 5 * 10 + iket] );

            // (D_002 P_001| = (F_003 S_000|_{t} + z_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_s_s_p_f = INT__f_s_f_s[abcd * 100 + 9 * 10 + iket] + ( AB_z[abcd] * INT__d_s_f_s[abcd * 60 + 5 * 10 + iket] );

            // (F_300 P_100| = (G_400 S_000|_{t} + x_ab * (F_300 S_000|_{t}
            const double Q_f_s_s_p_s_s_f = INT__g_s_f_s[abcd * 150 + 0 * 10 + iket] + ( AB_x[abcd] * INT__f_s_f_s[abcd * 100 + 0 * 10 + iket] );

            // (F_210 P_100| = (G_310 S_000|_{t} + x_ab * (F_210 S_000|_{t}
            const double Q_d_p_s_p_s_s_f = INT__g_s_f_s[abcd * 150 + 1 * 10 + iket] + ( AB_x[abcd] * INT__f_s_f_s[abcd * 100 + 1 * 10 + iket] );

            // (F_210 P_010| = (G_220 S_000|_{t} + y_ab * (F_210 S_000|_{t}
            const double Q_d_p_s_s_p_s_f = INT__g_s_f_s[abcd * 150 + 3 * 10 + iket] + ( AB_y[abcd] * INT__f_s_f_s[abcd * 100 + 1 * 10 + iket] );

            // (F_201 P_100| = (G_301 S_000|_{t} + x_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_p_s_s_f = INT__g_s_f_s[abcd * 150 + 2 * 10 + iket] + ( AB_x[abcd] * INT__f_s_f_s[abcd * 100 + 2 * 10 + iket] );

            // (F_201 P_010| = (G_211 S_000|_{t} + y_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_s_p_s_f = INT__g_s_f_s[abcd * 150 + 4 * 10 + iket] + ( AB_y[abcd] * INT__f_s_f_s[abcd * 100 + 2 * 10 + iket] );

            // (F_201 P_001| = (G_202 S_000|_{t} + z_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_s_s_p_f = INT__g_s_f_s[abcd * 150 + 5 * 10 + iket] + ( AB_z[abcd] * INT__f_s_f_s[abcd * 100 + 2 * 10 + iket] );

            // (F_120 P_100| = (G_220 S_000|_{t} + x_ab * (F_120 S_000|_{t}
            const double Q_p_d_s_p_s_s_f = INT__g_s_f_s[abcd * 150 + 3 * 10 + iket] + ( AB_x[abcd] * INT__f_s_f_s[abcd * 100 + 3 * 10 + iket] );

            // (F_120 P_010| = (G_130 S_000|_{t} + y_ab * (F_120 S_000|_{t}
            const double Q_p_d_s_s_p_s_f = INT__g_s_f_s[abcd * 150 + 6 * 10 + iket] + ( AB_y[abcd] * INT__f_s_f_s[abcd * 100 + 3 * 10 + iket] );

            // (F_111 P_100| = (G_211 S_000|_{t} + x_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_p_s_s_f = INT__g_s_f_s[abcd * 150 + 4 * 10 + iket] + ( AB_x[abcd] * INT__f_s_f_s[abcd * 100 + 4 * 10 + iket] );

            // (F_111 P_010| = (G_121 S_000|_{t} + y_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_s_p_s_f = INT__g_s_f_s[abcd * 150 + 7 * 10 + iket] + ( AB_y[abcd] * INT__f_s_f_s[abcd * 100 + 4 * 10 + iket] );

            // (F_111 P_001| = (G_112 S_000|_{t} + z_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_s_s_p_f = INT__g_s_f_s[abcd * 150 + 8 * 10 + iket] + ( AB_z[abcd] * INT__f_s_f_s[abcd * 100 + 4 * 10 + iket] );

            // (F_102 P_100| = (G_202 S_000|_{t} + x_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_p_s_s_f = INT__g_s_f_s[abcd * 150 + 5 * 10 + iket] + ( AB_x[abcd] * INT__f_s_f_s[abcd * 100 + 5 * 10 + iket] );

            // (F_102 P_010| = (G_112 S_000|_{t} + y_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_s_p_s_f = INT__g_s_f_s[abcd * 150 + 8 * 10 + iket] + ( AB_y[abcd] * INT__f_s_f_s[abcd * 100 + 5 * 10 + iket] );

            // (F_102 P_001| = (G_103 S_000|_{t} + z_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_s_s_p_f = INT__g_s_f_s[abcd * 150 + 9 * 10 + iket] + ( AB_z[abcd] * INT__f_s_f_s[abcd * 100 + 5 * 10 + iket] );

            // (F_030 P_100| = (G_130 S_000|_{t} + x_ab * (F_030 S_000|_{t}
            const double Q_s_f_s_p_s_s_f = INT__g_s_f_s[abcd * 150 + 6 * 10 + iket] + ( AB_x[abcd] * INT__f_s_f_s[abcd * 100 + 6 * 10 + iket] );

            // (F_030 P_010| = (G_040 S_000|_{t} + y_ab * (F_030 S_000|_{t}
            const double Q_s_f_s_s_p_s_f = INT__g_s_f_s[abcd * 150 + 10 * 10 + iket] + ( AB_y[abcd] * INT__f_s_f_s[abcd * 100 + 6 * 10 + iket] );

            // (F_021 P_100| = (G_121 S_000|_{t} + x_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_p_s_s_f = INT__g_s_f_s[abcd * 150 + 7 * 10 + iket] + ( AB_x[abcd] * INT__f_s_f_s[abcd * 100 + 7 * 10 + iket] );

            // (F_021 P_010| = (G_031 S_000|_{t} + y_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_s_p_s_f = INT__g_s_f_s[abcd * 150 + 11 * 10 + iket] + ( AB_y[abcd] * INT__f_s_f_s[abcd * 100 + 7 * 10 + iket] );

            // (F_021 P_001| = (G_022 S_000|_{t} + z_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_s_s_p_f = INT__g_s_f_s[abcd * 150 + 12 * 10 + iket] + ( AB_z[abcd] * INT__f_s_f_s[abcd * 100 + 7 * 10 + iket] );

            // (F_012 P_100| = (G_112 S_000|_{t} + x_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_p_s_s_f = INT__g_s_f_s[abcd * 150 + 8 * 10 + iket] + ( AB_x[abcd] * INT__f_s_f_s[abcd * 100 + 8 * 10 + iket] );

            // (F_012 P_010| = (G_022 S_000|_{t} + y_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_s_p_s_f = INT__g_s_f_s[abcd * 150 + 12 * 10 + iket] + ( AB_y[abcd] * INT__f_s_f_s[abcd * 100 + 8 * 10 + iket] );

            // (F_012 P_001| = (G_013 S_000|_{t} + z_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_s_s_p_f = INT__g_s_f_s[abcd * 150 + 13 * 10 + iket] + ( AB_z[abcd] * INT__f_s_f_s[abcd * 100 + 8 * 10 + iket] );

            // (F_003 P_100| = (G_103 S_000|_{t} + x_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_p_s_s_f = INT__g_s_f_s[abcd * 150 + 9 * 10 + iket] + ( AB_x[abcd] * INT__f_s_f_s[abcd * 100 + 9 * 10 + iket] );

            // (F_003 P_010| = (G_013 S_000|_{t} + y_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_s_p_s_f = INT__g_s_f_s[abcd * 150 + 13 * 10 + iket] + ( AB_y[abcd] * INT__f_s_f_s[abcd * 100 + 9 * 10 + iket] );

            // (F_003 P_001| = (G_004 S_000|_{t} + z_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_s_s_p_f = INT__g_s_f_s[abcd * 150 + 14 * 10 + iket] + ( AB_z[abcd] * INT__f_s_f_s[abcd * 100 + 9 * 10 + iket] );

            // (D_200 D_200|_{i} = (F_300 P_100| + x_ab * (D_200 P_100|
            INT__d_d_f_s[abcd * 360 + 0 * 10 + iket] = Q_f_s_s_p_s_s_f + ( AB_x[abcd] * Q_d_s_s_p_s_s_f );

            // (D_200 D_110|_{i} = (F_210 P_100| + y_ab * (D_200 P_100|
            INT__d_d_f_s[abcd * 360 + 1 * 10 + iket] = Q_d_p_s_p_s_s_f + ( AB_y[abcd] * Q_d_s_s_p_s_s_f );

            // (D_200 D_101|_{i} = (F_201 P_100| + z_ab * (D_200 P_100|
            INT__d_d_f_s[abcd * 360 + 2 * 10 + iket] = Q_d_s_p_p_s_s_f + ( AB_z[abcd] * Q_d_s_s_p_s_s_f );

            // (D_200 D_020|_{i} = (F_210 P_010| + y_ab * (D_200 P_010|
            INT__d_d_f_s[abcd * 360 + 3 * 10 + iket] = Q_d_p_s_s_p_s_f + ( AB_y[abcd] * Q_d_s_s_s_p_s_f );

            // (D_200 D_011|_{i} = (F_201 P_010| + z_ab * (D_200 P_010|
            INT__d_d_f_s[abcd * 360 + 4 * 10 + iket] = Q_d_s_p_s_p_s_f + ( AB_z[abcd] * Q_d_s_s_s_p_s_f );

            // (D_200 D_002|_{i} = (F_201 P_001| + z_ab * (D_200 P_001|
            INT__d_d_f_s[abcd * 360 + 5 * 10 + iket] = Q_d_s_p_s_s_p_f + ( AB_z[abcd] * Q_d_s_s_s_s_p_f );

            // (D_110 D_200|_{i} = (F_210 P_100| + x_ab * (D_110 P_100|
            INT__d_d_f_s[abcd * 360 + 6 * 10 + iket] = Q_d_p_s_p_s_s_f + ( AB_x[abcd] * Q_p_p_s_p_s_s_f );

            // (D_110 D_110|_{i} = (F_120 P_100| + y_ab * (D_110 P_100|
            INT__d_d_f_s[abcd * 360 + 7 * 10 + iket] = Q_p_d_s_p_s_s_f + ( AB_y[abcd] * Q_p_p_s_p_s_s_f );

            // (D_110 D_101|_{i} = (F_111 P_100| + z_ab * (D_110 P_100|
            INT__d_d_f_s[abcd * 360 + 8 * 10 + iket] = Q_p_p_p_p_s_s_f + ( AB_z[abcd] * Q_p_p_s_p_s_s_f );

            // (D_110 D_020|_{i} = (F_120 P_010| + y_ab * (D_110 P_010|
            INT__d_d_f_s[abcd * 360 + 9 * 10 + iket] = Q_p_d_s_s_p_s_f + ( AB_y[abcd] * Q_p_p_s_s_p_s_f );

            // (D_110 D_011|_{i} = (F_111 P_010| + z_ab * (D_110 P_010|
            INT__d_d_f_s[abcd * 360 + 10 * 10 + iket] = Q_p_p_p_s_p_s_f + ( AB_z[abcd] * Q_p_p_s_s_p_s_f );

            // (D_110 D_002|_{i} = (F_111 P_001| + z_ab * (D_110 P_001|
            INT__d_d_f_s[abcd * 360 + 11 * 10 + iket] = Q_p_p_p_s_s_p_f + ( AB_z[abcd] * Q_p_p_s_s_s_p_f );

            // (D_101 D_200|_{i} = (F_201 P_100| + x_ab * (D_101 P_100|
            INT__d_d_f_s[abcd * 360 + 12 * 10 + iket] = Q_d_s_p_p_s_s_f + ( AB_x[abcd] * Q_p_s_p_p_s_s_f );

            // (D_101 D_110|_{i} = (F_111 P_100| + y_ab * (D_101 P_100|
            INT__d_d_f_s[abcd * 360 + 13 * 10 + iket] = Q_p_p_p_p_s_s_f + ( AB_y[abcd] * Q_p_s_p_p_s_s_f );

            // (D_101 D_101|_{i} = (F_102 P_100| + z_ab * (D_101 P_100|
            INT__d_d_f_s[abcd * 360 + 14 * 10 + iket] = Q_p_s_d_p_s_s_f + ( AB_z[abcd] * Q_p_s_p_p_s_s_f );

            // (D_101 D_020|_{i} = (F_111 P_010| + y_ab * (D_101 P_010|
            INT__d_d_f_s[abcd * 360 + 15 * 10 + iket] = Q_p_p_p_s_p_s_f + ( AB_y[abcd] * Q_p_s_p_s_p_s_f );

            // (D_101 D_011|_{i} = (F_102 P_010| + z_ab * (D_101 P_010|
            INT__d_d_f_s[abcd * 360 + 16 * 10 + iket] = Q_p_s_d_s_p_s_f + ( AB_z[abcd] * Q_p_s_p_s_p_s_f );

            // (D_101 D_002|_{i} = (F_102 P_001| + z_ab * (D_101 P_001|
            INT__d_d_f_s[abcd * 360 + 17 * 10 + iket] = Q_p_s_d_s_s_p_f + ( AB_z[abcd] * Q_p_s_p_s_s_p_f );

            // (D_020 D_200|_{i} = (F_120 P_100| + x_ab * (D_020 P_100|
            INT__d_d_f_s[abcd * 360 + 18 * 10 + iket] = Q_p_d_s_p_s_s_f + ( AB_x[abcd] * Q_s_d_s_p_s_s_f );

            // (D_020 D_110|_{i} = (F_030 P_100| + y_ab * (D_020 P_100|
            INT__d_d_f_s[abcd * 360 + 19 * 10 + iket] = Q_s_f_s_p_s_s_f + ( AB_y[abcd] * Q_s_d_s_p_s_s_f );

            // (D_020 D_101|_{i} = (F_021 P_100| + z_ab * (D_020 P_100|
            INT__d_d_f_s[abcd * 360 + 20 * 10 + iket] = Q_s_d_p_p_s_s_f + ( AB_z[abcd] * Q_s_d_s_p_s_s_f );

            // (D_020 D_020|_{i} = (F_030 P_010| + y_ab * (D_020 P_010|
            INT__d_d_f_s[abcd * 360 + 21 * 10 + iket] = Q_s_f_s_s_p_s_f + ( AB_y[abcd] * Q_s_d_s_s_p_s_f );

            // (D_020 D_011|_{i} = (F_021 P_010| + z_ab * (D_020 P_010|
            INT__d_d_f_s[abcd * 360 + 22 * 10 + iket] = Q_s_d_p_s_p_s_f + ( AB_z[abcd] * Q_s_d_s_s_p_s_f );

            // (D_020 D_002|_{i} = (F_021 P_001| + z_ab * (D_020 P_001|
            INT__d_d_f_s[abcd * 360 + 23 * 10 + iket] = Q_s_d_p_s_s_p_f + ( AB_z[abcd] * Q_s_d_s_s_s_p_f );

            // (D_011 D_200|_{i} = (F_111 P_100| + x_ab * (D_011 P_100|
            INT__d_d_f_s[abcd * 360 + 24 * 10 + iket] = Q_p_p_p_p_s_s_f + ( AB_x[abcd] * Q_s_p_p_p_s_s_f );

            // (D_011 D_110|_{i} = (F_021 P_100| + y_ab * (D_011 P_100|
            INT__d_d_f_s[abcd * 360 + 25 * 10 + iket] = Q_s_d_p_p_s_s_f + ( AB_y[abcd] * Q_s_p_p_p_s_s_f );

            // (D_011 D_101|_{i} = (F_012 P_100| + z_ab * (D_011 P_100|
            INT__d_d_f_s[abcd * 360 + 26 * 10 + iket] = Q_s_p_d_p_s_s_f + ( AB_z[abcd] * Q_s_p_p_p_s_s_f );

            // (D_011 D_020|_{i} = (F_021 P_010| + y_ab * (D_011 P_010|
            INT__d_d_f_s[abcd * 360 + 27 * 10 + iket] = Q_s_d_p_s_p_s_f + ( AB_y[abcd] * Q_s_p_p_s_p_s_f );

            // (D_011 D_011|_{i} = (F_012 P_010| + z_ab * (D_011 P_010|
            INT__d_d_f_s[abcd * 360 + 28 * 10 + iket] = Q_s_p_d_s_p_s_f + ( AB_z[abcd] * Q_s_p_p_s_p_s_f );

            // (D_011 D_002|_{i} = (F_012 P_001| + z_ab * (D_011 P_001|
            INT__d_d_f_s[abcd * 360 + 29 * 10 + iket] = Q_s_p_d_s_s_p_f + ( AB_z[abcd] * Q_s_p_p_s_s_p_f );

            // (D_002 D_200|_{i} = (F_102 P_100| + x_ab * (D_002 P_100|
            INT__d_d_f_s[abcd * 360 + 30 * 10 + iket] = Q_p_s_d_p_s_s_f + ( AB_x[abcd] * Q_s_s_d_p_s_s_f );

            // (D_002 D_110|_{i} = (F_012 P_100| + y_ab * (D_002 P_100|
            INT__d_d_f_s[abcd * 360 + 31 * 10 + iket] = Q_s_p_d_p_s_s_f + ( AB_y[abcd] * Q_s_s_d_p_s_s_f );

            // (D_002 D_101|_{i} = (F_003 P_100| + z_ab * (D_002 P_100|
            INT__d_d_f_s[abcd * 360 + 32 * 10 + iket] = Q_s_s_f_p_s_s_f + ( AB_z[abcd] * Q_s_s_d_p_s_s_f );

            // (D_002 D_020|_{i} = (F_012 P_010| + y_ab * (D_002 P_010|
            INT__d_d_f_s[abcd * 360 + 33 * 10 + iket] = Q_s_p_d_s_p_s_f + ( AB_y[abcd] * Q_s_s_d_s_p_s_f );

            // (D_002 D_011|_{i} = (F_003 P_010| + z_ab * (D_002 P_010|
            INT__d_d_f_s[abcd * 360 + 34 * 10 + iket] = Q_s_s_f_s_p_s_f + ( AB_z[abcd] * Q_s_s_d_s_p_s_f );

            // (D_002 D_002|_{i} = (F_003 P_001| + z_ab * (D_002 P_001|
            INT__d_d_f_s[abcd * 360 + 35 * 10 + iket] = Q_s_s_f_s_s_p_f + ( AB_z[abcd] * Q_s_s_d_s_s_p_f );

        }

        // form INT__d_d_g_s
        for(int iket = 0; iket < 15; ++iket)
        {
            // (D_200 P_100| = (F_300 S_000|_{t} + x_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_p_s_s_g = INT__f_s_g_s[abcd * 150 + 0 * 15 + iket] + ( AB_x[abcd] * INT__d_s_g_s[abcd * 90 + 0 * 15 + iket] );

            // (D_200 P_010| = (F_210 S_000|_{t} + y_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_s_p_s_g = INT__f_s_g_s[abcd * 150 + 1 * 15 + iket] + ( AB_y[abcd] * INT__d_s_g_s[abcd * 90 + 0 * 15 + iket] );

            // (D_200 P_001| = (F_201 S_000|_{t} + z_ab * (D_200 S_000|_{t}
            const double Q_d_s_s_s_s_p_g = INT__f_s_g_s[abcd * 150 + 2 * 15 + iket] + ( AB_z[abcd] * INT__d_s_g_s[abcd * 90 + 0 * 15 + iket] );

            // (D_110 P_100| = (F_210 S_000|_{t} + x_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_p_s_s_g = INT__f_s_g_s[abcd * 150 + 1 * 15 + iket] + ( AB_x[abcd] * INT__d_s_g_s[abcd * 90 + 1 * 15 + iket] );

            // (D_110 P_010| = (F_120 S_000|_{t} + y_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_s_p_s_g = INT__f_s_g_s[abcd * 150 + 3 * 15 + iket] + ( AB_y[abcd] * INT__d_s_g_s[abcd * 90 + 1 * 15 + iket] );

            // (D_110 P_001| = (F_111 S_000|_{t} + z_ab * (D_110 S_000|_{t}
            const double Q_p_p_s_s_s_p_g = INT__f_s_g_s[abcd * 150 + 4 * 15 + iket] + ( AB_z[abcd] * INT__d_s_g_s[abcd * 90 + 1 * 15 + iket] );

            // (D_101 P_100| = (F_201 S_000|_{t} + x_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_p_s_s_g = INT__f_s_g_s[abcd * 150 + 2 * 15 + iket] + ( AB_x[abcd] * INT__d_s_g_s[abcd * 90 + 2 * 15 + iket] );

            // (D_101 P_010| = (F_111 S_000|_{t} + y_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_s_p_s_g = INT__f_s_g_s[abcd * 150 + 4 * 15 + iket] + ( AB_y[abcd] * INT__d_s_g_s[abcd * 90 + 2 * 15 + iket] );

            // (D_101 P_001| = (F_102 S_000|_{t} + z_ab * (D_101 S_000|_{t}
            const double Q_p_s_p_s_s_p_g = INT__f_s_g_s[abcd * 150 + 5 * 15 + iket] + ( AB_z[abcd] * INT__d_s_g_s[abcd * 90 + 2 * 15 + iket] );

            // (D_020 P_100| = (F_120 S_000|_{t} + x_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_p_s_s_g = INT__f_s_g_s[abcd * 150 + 3 * 15 + iket] + ( AB_x[abcd] * INT__d_s_g_s[abcd * 90 + 3 * 15 + iket] );

            // (D_020 P_010| = (F_030 S_000|_{t} + y_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_s_p_s_g = INT__f_s_g_s[abcd * 150 + 6 * 15 + iket] + ( AB_y[abcd] * INT__d_s_g_s[abcd * 90 + 3 * 15 + iket] );

            // (D_020 P_001| = (F_021 S_000|_{t} + z_ab * (D_020 S_000|_{t}
            const double Q_s_d_s_s_s_p_g = INT__f_s_g_s[abcd * 150 + 7 * 15 + iket] + ( AB_z[abcd] * INT__d_s_g_s[abcd * 90 + 3 * 15 + iket] );

            // (D_011 P_100| = (F_111 S_000|_{t} + x_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_p_s_s_g = INT__f_s_g_s[abcd * 150 + 4 * 15 + iket] + ( AB_x[abcd] * INT__d_s_g_s[abcd * 90 + 4 * 15 + iket] );

            // (D_011 P_010| = (F_021 S_000|_{t} + y_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_s_p_s_g = INT__f_s_g_s[abcd * 150 + 7 * 15 + iket] + ( AB_y[abcd] * INT__d_s_g_s[abcd * 90 + 4 * 15 + iket] );

            // (D_011 P_001| = (F_012 S_000|_{t} + z_ab * (D_011 S_000|_{t}
            const double Q_s_p_p_s_s_p_g = INT__f_s_g_s[abcd * 150 + 8 * 15 + iket] + ( AB_z[abcd] * INT__d_s_g_s[abcd * 90 + 4 * 15 + iket] );

            // (D_002 P_100| = (F_102 S_000|_{t} + x_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_p_s_s_g = INT__f_s_g_s[abcd * 150 + 5 * 15 + iket] + ( AB_x[abcd] * INT__d_s_g_s[abcd * 90 + 5 * 15 + iket] );

            // (D_002 P_010| = (F_012 S_000|_{t} + y_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_s_p_s_g = INT__f_s_g_s[abcd * 150 + 8 * 15 + iket] + ( AB_y[abcd] * INT__d_s_g_s[abcd * 90 + 5 * 15 + iket] );

            // (D_002 P_001| = (F_003 S_000|_{t} + z_ab * (D_002 S_000|_{t}
            const double Q_s_s_d_s_s_p_g = INT__f_s_g_s[abcd * 150 + 9 * 15 + iket] + ( AB_z[abcd] * INT__d_s_g_s[abcd * 90 + 5 * 15 + iket] );

            // (F_300 P_100| = (G_400 S_000|_{t} + x_ab * (F_300 S_000|_{t}
            const double Q_f_s_s_p_s_s_g = INT__g_s_g_s[abcd * 225 + 0 * 15 + iket] + ( AB_x[abcd] * INT__f_s_g_s[abcd * 150 + 0 * 15 + iket] );

            // (F_210 P_100| = (G_310 S_000|_{t} + x_ab * (F_210 S_000|_{t}
            const double Q_d_p_s_p_s_s_g = INT__g_s_g_s[abcd * 225 + 1 * 15 + iket] + ( AB_x[abcd] * INT__f_s_g_s[abcd * 150 + 1 * 15 + iket] );

            // (F_210 P_010| = (G_220 S_000|_{t} + y_ab * (F_210 S_000|_{t}
            const double Q_d_p_s_s_p_s_g = INT__g_s_g_s[abcd * 225 + 3 * 15 + iket] + ( AB_y[abcd] * INT__f_s_g_s[abcd * 150 + 1 * 15 + iket] );

            // (F_201 P_100| = (G_301 S_000|_{t} + x_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_p_s_s_g = INT__g_s_g_s[abcd * 225 + 2 * 15 + iket] + ( AB_x[abcd] * INT__f_s_g_s[abcd * 150 + 2 * 15 + iket] );

            // (F_201 P_010| = (G_211 S_000|_{t} + y_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_s_p_s_g = INT__g_s_g_s[abcd * 225 + 4 * 15 + iket] + ( AB_y[abcd] * INT__f_s_g_s[abcd * 150 + 2 * 15 + iket] );

            // (F_201 P_001| = (G_202 S_000|_{t} + z_ab * (F_201 S_000|_{t}
            const double Q_d_s_p_s_s_p_g = INT__g_s_g_s[abcd * 225 + 5 * 15 + iket] + ( AB_z[abcd] * INT__f_s_g_s[abcd * 150 + 2 * 15 + iket] );

            // (F_120 P_100| = (G_220 S_000|_{t} + x_ab * (F_120 S_000|_{t}
            const double Q_p_d_s_p_s_s_g = INT__g_s_g_s[abcd * 225 + 3 * 15 + iket] + ( AB_x[abcd] * INT__f_s_g_s[abcd * 150 + 3 * 15 + iket] );

            // (F_120 P_010| = (G_130 S_000|_{t} + y_ab * (F_120 S_000|_{t}
            const double Q_p_d_s_s_p_s_g = INT__g_s_g_s[abcd * 225 + 6 * 15 + iket] + ( AB_y[abcd] * INT__f_s_g_s[abcd * 150 + 3 * 15 + iket] );

            // (F_111 P_100| = (G_211 S_000|_{t} + x_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_p_s_s_g = INT__g_s_g_s[abcd * 225 + 4 * 15 + iket] + ( AB_x[abcd] * INT__f_s_g_s[abcd * 150 + 4 * 15 + iket] );

            // (F_111 P_010| = (G_121 S_000|_{t} + y_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_s_p_s_g = INT__g_s_g_s[abcd * 225 + 7 * 15 + iket] + ( AB_y[abcd] * INT__f_s_g_s[abcd * 150 + 4 * 15 + iket] );

            // (F_111 P_001| = (G_112 S_000|_{t} + z_ab * (F_111 S_000|_{t}
            const double Q_p_p_p_s_s_p_g = INT__g_s_g_s[abcd * 225 + 8 * 15 + iket] + ( AB_z[abcd] * INT__f_s_g_s[abcd * 150 + 4 * 15 + iket] );

            // (F_102 P_100| = (G_202 S_000|_{t} + x_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_p_s_s_g = INT__g_s_g_s[abcd * 225 + 5 * 15 + iket] + ( AB_x[abcd] * INT__f_s_g_s[abcd * 150 + 5 * 15 + iket] );

            // (F_102 P_010| = (G_112 S_000|_{t} + y_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_s_p_s_g = INT__g_s_g_s[abcd * 225 + 8 * 15 + iket] + ( AB_y[abcd] * INT__f_s_g_s[abcd * 150 + 5 * 15 + iket] );

            // (F_102 P_001| = (G_103 S_000|_{t} + z_ab * (F_102 S_000|_{t}
            const double Q_p_s_d_s_s_p_g = INT__g_s_g_s[abcd * 225 + 9 * 15 + iket] + ( AB_z[abcd] * INT__f_s_g_s[abcd * 150 + 5 * 15 + iket] );

            // (F_030 P_100| = (G_130 S_000|_{t} + x_ab * (F_030 S_000|_{t}
            const double Q_s_f_s_p_s_s_g = INT__g_s_g_s[abcd * 225 + 6 * 15 + iket] + ( AB_x[abcd] * INT__f_s_g_s[abcd * 150 + 6 * 15 + iket] );

            // (F_030 P_010| = (G_040 S_000|_{t} + y_ab * (F_030 S_000|_{t}
            const double Q_s_f_s_s_p_s_g = INT__g_s_g_s[abcd * 225 + 10 * 15 + iket] + ( AB_y[abcd] * INT__f_s_g_s[abcd * 150 + 6 * 15 + iket] );

            // (F_021 P_100| = (G_121 S_000|_{t} + x_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_p_s_s_g = INT__g_s_g_s[abcd * 225 + 7 * 15 + iket] + ( AB_x[abcd] * INT__f_s_g_s[abcd * 150 + 7 * 15 + iket] );

            // (F_021 P_010| = (G_031 S_000|_{t} + y_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_s_p_s_g = INT__g_s_g_s[abcd * 225 + 11 * 15 + iket] + ( AB_y[abcd] * INT__f_s_g_s[abcd * 150 + 7 * 15 + iket] );

            // (F_021 P_001| = (G_022 S_000|_{t} + z_ab * (F_021 S_000|_{t}
            const double Q_s_d_p_s_s_p_g = INT__g_s_g_s[abcd * 225 + 12 * 15 + iket] + ( AB_z[abcd] * INT__f_s_g_s[abcd * 150 + 7 * 15 + iket] );

            // (F_012 P_100| = (G_112 S_000|_{t} + x_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_p_s_s_g = INT__g_s_g_s[abcd * 225 + 8 * 15 + iket] + ( AB_x[abcd] * INT__f_s_g_s[abcd * 150 + 8 * 15 + iket] );

            // (F_012 P_010| = (G_022 S_000|_{t} + y_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_s_p_s_g = INT__g_s_g_s[abcd * 225 + 12 * 15 + iket] + ( AB_y[abcd] * INT__f_s_g_s[abcd * 150 + 8 * 15 + iket] );

            // (F_012 P_001| = (G_013 S_000|_{t} + z_ab * (F_012 S_000|_{t}
            const double Q_s_p_d_s_s_p_g = INT__g_s_g_s[abcd * 225 + 13 * 15 + iket] + ( AB_z[abcd] * INT__f_s_g_s[abcd * 150 + 8 * 15 + iket] );

            // (F_003 P_100| = (G_103 S_000|_{t} + x_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_p_s_s_g = INT__g_s_g_s[abcd * 225 + 9 * 15 + iket] + ( AB_x[abcd] * INT__f_s_g_s[abcd * 150 + 9 * 15 + iket] );

            // (F_003 P_010| = (G_013 S_000|_{t} + y_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_s_p_s_g = INT__g_s_g_s[abcd * 225 + 13 * 15 + iket] + ( AB_y[abcd] * INT__f_s_g_s[abcd * 150 + 9 * 15 + iket] );

            // (F_003 P_001| = (G_004 S_000|_{t} + z_ab * (F_003 S_000|_{t}
            const double Q_s_s_f_s_s_p_g = INT__g_s_g_s[abcd * 225 + 14 * 15 + iket] + ( AB_z[abcd] * INT__f_s_g_s[abcd * 150 + 9 * 15 + iket] );

            // (D_200 D_200|_{i} = (F_300 P_100| + x_ab * (D_200 P_100|
            INT__d_d_g_s[abcd * 540 + 0 * 15 + iket] = Q_f_s_s_p_s_s_g + ( AB_x[abcd] * Q_d_s_s_p_s_s_g );

            // (D_200 D_110|_{i} = (F_210 P_100| + y_ab * (D_200 P_100|
            INT__d_d_g_s[abcd * 540 + 1 * 15 + iket] = Q_d_p_s_p_s_s_g + ( AB_y[abcd] * Q_d_s_s_p_s_s_g );

            // (D_200 D_101|_{i} = (F_201 P_100| + z_ab * (D_200 P_100|
            INT__d_d_g_s[abcd * 540 + 2 * 15 + iket] = Q_d_s_p_p_s_s_g + ( AB_z[abcd] * Q_d_s_s_p_s_s_g );

            // (D_200 D_020|_{i} = (F_210 P_010| + y_ab * (D_200 P_010|
            INT__d_d_g_s[abcd * 540 + 3 * 15 + iket] = Q_d_p_s_s_p_s_g + ( AB_y[abcd] * Q_d_s_s_s_p_s_g );

            // (D_200 D_011|_{i} = (F_201 P_010| + z_ab * (D_200 P_010|
            INT__d_d_g_s[abcd * 540 + 4 * 15 + iket] = Q_d_s_p_s_p_s_g + ( AB_z[abcd] * Q_d_s_s_s_p_s_g );

            // (D_200 D_002|_{i} = (F_201 P_001| + z_ab * (D_200 P_001|
            INT__d_d_g_s[abcd * 540 + 5 * 15 + iket] = Q_d_s_p_s_s_p_g + ( AB_z[abcd] * Q_d_s_s_s_s_p_g );

            // (D_110 D_200|_{i} = (F_210 P_100| + x_ab * (D_110 P_100|
            INT__d_d_g_s[abcd * 540 + 6 * 15 + iket] = Q_d_p_s_p_s_s_g + ( AB_x[abcd] * Q_p_p_s_p_s_s_g );

            // (D_110 D_110|_{i} = (F_120 P_100| + y_ab * (D_110 P_100|
            INT__d_d_g_s[abcd * 540 + 7 * 15 + iket] = Q_p_d_s_p_s_s_g + ( AB_y[abcd] * Q_p_p_s_p_s_s_g );

            // (D_110 D_101|_{i} = (F_111 P_100| + z_ab * (D_110 P_100|
            INT__d_d_g_s[abcd * 540 + 8 * 15 + iket] = Q_p_p_p_p_s_s_g + ( AB_z[abcd] * Q_p_p_s_p_s_s_g );

            // (D_110 D_020|_{i} = (F_120 P_010| + y_ab * (D_110 P_010|
            INT__d_d_g_s[abcd * 540 + 9 * 15 + iket] = Q_p_d_s_s_p_s_g + ( AB_y[abcd] * Q_p_p_s_s_p_s_g );

            // (D_110 D_011|_{i} = (F_111 P_010| + z_ab * (D_110 P_010|
            INT__d_d_g_s[abcd * 540 + 10 * 15 + iket] = Q_p_p_p_s_p_s_g + ( AB_z[abcd] * Q_p_p_s_s_p_s_g );

            // (D_110 D_002|_{i} = (F_111 P_001| + z_ab * (D_110 P_001|
            INT__d_d_g_s[abcd * 540 + 11 * 15 + iket] = Q_p_p_p_s_s_p_g + ( AB_z[abcd] * Q_p_p_s_s_s_p_g );

            // (D_101 D_200|_{i} = (F_201 P_100| + x_ab * (D_101 P_100|
            INT__d_d_g_s[abcd * 540 + 12 * 15 + iket] = Q_d_s_p_p_s_s_g + ( AB_x[abcd] * Q_p_s_p_p_s_s_g );

            // (D_101 D_110|_{i} = (F_111 P_100| + y_ab * (D_101 P_100|
            INT__d_d_g_s[abcd * 540 + 13 * 15 + iket] = Q_p_p_p_p_s_s_g + ( AB_y[abcd] * Q_p_s_p_p_s_s_g );

            // (D_101 D_101|_{i} = (F_102 P_100| + z_ab * (D_101 P_100|
            INT__d_d_g_s[abcd * 540 + 14 * 15 + iket] = Q_p_s_d_p_s_s_g + ( AB_z[abcd] * Q_p_s_p_p_s_s_g );

            // (D_101 D_020|_{i} = (F_111 P_010| + y_ab * (D_101 P_010|
            INT__d_d_g_s[abcd * 540 + 15 * 15 + iket] = Q_p_p_p_s_p_s_g + ( AB_y[abcd] * Q_p_s_p_s_p_s_g );

            // (D_101 D_011|_{i} = (F_102 P_010| + z_ab * (D_101 P_010|
            INT__d_d_g_s[abcd * 540 + 16 * 15 + iket] = Q_p_s_d_s_p_s_g + ( AB_z[abcd] * Q_p_s_p_s_p_s_g );

            // (D_101 D_002|_{i} = (F_102 P_001| + z_ab * (D_101 P_001|
            INT__d_d_g_s[abcd * 540 + 17 * 15 + iket] = Q_p_s_d_s_s_p_g + ( AB_z[abcd] * Q_p_s_p_s_s_p_g );

            // (D_020 D_200|_{i} = (F_120 P_100| + x_ab * (D_020 P_100|
            INT__d_d_g_s[abcd * 540 + 18 * 15 + iket] = Q_p_d_s_p_s_s_g + ( AB_x[abcd] * Q_s_d_s_p_s_s_g );

            // (D_020 D_110|_{i} = (F_030 P_100| + y_ab * (D_020 P_100|
            INT__d_d_g_s[abcd * 540 + 19 * 15 + iket] = Q_s_f_s_p_s_s_g + ( AB_y[abcd] * Q_s_d_s_p_s_s_g );

            // (D_020 D_101|_{i} = (F_021 P_100| + z_ab * (D_020 P_100|
            INT__d_d_g_s[abcd * 540 + 20 * 15 + iket] = Q_s_d_p_p_s_s_g + ( AB_z[abcd] * Q_s_d_s_p_s_s_g );

            // (D_020 D_020|_{i} = (F_030 P_010| + y_ab * (D_020 P_010|
            INT__d_d_g_s[abcd * 540 + 21 * 15 + iket] = Q_s_f_s_s_p_s_g + ( AB_y[abcd] * Q_s_d_s_s_p_s_g );

            // (D_020 D_011|_{i} = (F_021 P_010| + z_ab * (D_020 P_010|
            INT__d_d_g_s[abcd * 540 + 22 * 15 + iket] = Q_s_d_p_s_p_s_g + ( AB_z[abcd] * Q_s_d_s_s_p_s_g );

            // (D_020 D_002|_{i} = (F_021 P_001| + z_ab * (D_020 P_001|
            INT__d_d_g_s[abcd * 540 + 23 * 15 + iket] = Q_s_d_p_s_s_p_g + ( AB_z[abcd] * Q_s_d_s_s_s_p_g );

            // (D_011 D_200|_{i} = (F_111 P_100| + x_ab * (D_011 P_100|
            INT__d_d_g_s[abcd * 540 + 24 * 15 + iket] = Q_p_p_p_p_s_s_g + ( AB_x[abcd] * Q_s_p_p_p_s_s_g );

            // (D_011 D_110|_{i} = (F_021 P_100| + y_ab * (D_011 P_100|
            INT__d_d_g_s[abcd * 540 + 25 * 15 + iket] = Q_s_d_p_p_s_s_g + ( AB_y[abcd] * Q_s_p_p_p_s_s_g );

            // (D_011 D_101|_{i} = (F_012 P_100| + z_ab * (D_011 P_100|
            INT__d_d_g_s[abcd * 540 + 26 * 15 + iket] = Q_s_p_d_p_s_s_g + ( AB_z[abcd] * Q_s_p_p_p_s_s_g );

            // (D_011 D_020|_{i} = (F_021 P_010| + y_ab * (D_011 P_010|
            INT__d_d_g_s[abcd * 540 + 27 * 15 + iket] = Q_s_d_p_s_p_s_g + ( AB_y[abcd] * Q_s_p_p_s_p_s_g );

            // (D_011 D_011|_{i} = (F_012 P_010| + z_ab * (D_011 P_010|
            INT__d_d_g_s[abcd * 540 + 28 * 15 + iket] = Q_s_p_d_s_p_s_g + ( AB_z[abcd] * Q_s_p_p_s_p_s_g );

            // (D_011 D_002|_{i} = (F_012 P_001| + z_ab * (D_011 P_001|
            INT__d_d_g_s[abcd * 540 + 29 * 15 + iket] = Q_s_p_d_s_s_p_g + ( AB_z[abcd] * Q_s_p_p_s_s_p_g );

            // (D_002 D_200|_{i} = (F_102 P_100| + x_ab * (D_002 P_100|
            INT__d_d_g_s[abcd * 540 + 30 * 15 + iket] = Q_p_s_d_p_s_s_g + ( AB_x[abcd] * Q_s_s_d_p_s_s_g );

            // (D_002 D_110|_{i} = (F_012 P_100| + y_ab * (D_002 P_100|
            INT__d_d_g_s[abcd * 540 + 31 * 15 + iket] = Q_s_p_d_p_s_s_g + ( AB_y[abcd] * Q_s_s_d_p_s_s_g );

            // (D_002 D_101|_{i} = (F_003 P_100| + z_ab * (D_002 P_100|
            INT__d_d_g_s[abcd * 540 + 32 * 15 + iket] = Q_s_s_f_p_s_s_g + ( AB_z[abcd] * Q_s_s_d_p_s_s_g );

            // (D_002 D_020|_{i} = (F_012 P_010| + y_ab * (D_002 P_010|
            INT__d_d_g_s[abcd * 540 + 33 * 15 + iket] = Q_s_p_d_s_p_s_g + ( AB_y[abcd] * Q_s_s_d_s_p_s_g );

            // (D_002 D_011|_{i} = (F_003 P_010| + z_ab * (D_002 P_010|
            INT__d_d_g_s[abcd * 540 + 34 * 15 + iket] = Q_s_s_f_s_p_s_g + ( AB_z[abcd] * Q_s_s_d_s_p_s_g );

            // (D_002 D_002|_{i} = (F_003 P_001| + z_ab * (D_002 P_001|
            INT__d_d_g_s[abcd * 540 + 35 * 15 + iket] = Q_s_s_f_s_s_p_g + ( AB_z[abcd] * Q_s_s_d_s_s_p_g );

        }


    }


    //////////////////////////////////////////////
    // Contracted integrals: Horizontal recurrance
    // Ket part
    // Steps: 79
    //////////////////////////////////////////////

    for(abcd = 0; abcd < nshell1234; ++abcd)
    {
        for(int ibra = 0; ibra < 36; ++ibra)
        {
            // |D_200 P_100) = |F_300 S_000)_{t} + x_cd * |D_200 S_000)_{t}
            const double Q_d_s_s_p_s_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 0] + ( CD_x[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 0] );

            // |D_200 P_010) = |F_210 S_000)_{t} + y_cd * |D_200 S_000)_{t}
            const double Q_d_s_s_s_p_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 1] + ( CD_y[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 0] );

            // |D_200 P_001) = |F_201 S_000)_{t} + z_cd * |D_200 S_000)_{t}
            const double Q_d_s_s_s_s_p_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 2] + ( CD_z[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 0] );

            // |D_110 P_100) = |F_210 S_000)_{t} + x_cd * |D_110 S_000)_{t}
            const double Q_p_p_s_p_s_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 1] + ( CD_x[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 1] );

            // |D_110 P_010) = |F_120 S_000)_{t} + y_cd * |D_110 S_000)_{t}
            const double Q_p_p_s_s_p_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 3] + ( CD_y[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 1] );

            // |D_110 P_001) = |F_111 S_000)_{t} + z_cd * |D_110 S_000)_{t}
            const double Q_p_p_s_s_s_p_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 4] + ( CD_z[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 1] );

            // |D_101 P_100) = |F_201 S_000)_{t} + x_cd * |D_101 S_000)_{t}
            const double Q_p_s_p_p_s_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 2] + ( CD_x[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 2] );

            // |D_101 P_010) = |F_111 S_000)_{t} + y_cd * |D_101 S_000)_{t}
            const double Q_p_s_p_s_p_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 4] + ( CD_y[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 2] );

            // |D_101 P_001) = |F_102 S_000)_{t} + z_cd * |D_101 S_000)_{t}
            const double Q_p_s_p_s_s_p_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 5] + ( CD_z[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 2] );

            // |D_020 P_100) = |F_120 S_000)_{t} + x_cd * |D_020 S_000)_{t}
            const double Q_s_d_s_p_s_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 3] + ( CD_x[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 3] );

            // |D_020 P_010) = |F_030 S_000)_{t} + y_cd * |D_020 S_000)_{t}
            const double Q_s_d_s_s_p_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 6] + ( CD_y[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 3] );

            // |D_020 P_001) = |F_021 S_000)_{t} + z_cd * |D_020 S_000)_{t}
            const double Q_s_d_s_s_s_p_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 7] + ( CD_z[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 3] );

            // |D_011 P_100) = |F_111 S_000)_{t} + x_cd * |D_011 S_000)_{t}
            const double Q_s_p_p_p_s_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 4] + ( CD_x[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 4] );

            // |D_011 P_010) = |F_021 S_000)_{t} + y_cd * |D_011 S_000)_{t}
            const double Q_s_p_p_s_p_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 7] + ( CD_y[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 4] );

            // |D_011 P_001) = |F_012 S_000)_{t} + z_cd * |D_011 S_000)_{t}
            const double Q_s_p_p_s_s_p_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 8] + ( CD_z[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 4] );

            // |D_002 P_100) = |F_102 S_000)_{t} + x_cd * |D_002 S_000)_{t}
            const double Q_s_s_d_p_s_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 5] + ( CD_x[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 5] );

            // |D_002 P_010) = |F_012 S_000)_{t} + y_cd * |D_002 S_000)_{t}
            const double Q_s_s_d_s_p_s_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 8] + ( CD_y[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 5] );

            // |D_002 P_001) = |F_003 S_000)_{t} + z_cd * |D_002 S_000)_{t}
            const double Q_s_s_d_s_s_p_d_d = INT__d_d_f_s[abcd * 360 + ibra * 10 + 9] + ( CD_z[abcd] * INT__d_d_d_s[abcd * 216 + ibra * 6 + 5] );

            // |F_300 P_100) = |G_400 S_000)_{t} + x_cd * |F_300 S_000)_{t}
            const double Q_f_s_s_p_s_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 0] + ( CD_x[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 0] );

            // |F_210 P_100) = |G_310 S_000)_{t} + x_cd * |F_210 S_000)_{t}
            const double Q_d_p_s_p_s_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 1] + ( CD_x[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 1] );

            // |F_210 P_010) = |G_220 S_000)_{t} + y_cd * |F_210 S_000)_{t}
            const double Q_d_p_s_s_p_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 3] + ( CD_y[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 1] );

            // |F_201 P_100) = |G_301 S_000)_{t} + x_cd * |F_201 S_000)_{t}
            const double Q_d_s_p_p_s_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 2] + ( CD_x[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 2] );

            // |F_201 P_010) = |G_211 S_000)_{t} + y_cd * |F_201 S_000)_{t}
            const double Q_d_s_p_s_p_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 4] + ( CD_y[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 2] );

            // |F_201 P_001) = |G_202 S_000)_{t} + z_cd * |F_201 S_000)_{t}
            const double Q_d_s_p_s_s_p_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 5] + ( CD_z[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 2] );

            // |F_120 P_100) = |G_220 S_000)_{t} + x_cd * |F_120 S_000)_{t}
            const double Q_p_d_s_p_s_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 3] + ( CD_x[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 3] );

            // |F_120 P_010) = |G_130 S_000)_{t} + y_cd * |F_120 S_000)_{t}
            const double Q_p_d_s_s_p_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 6] + ( CD_y[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 3] );

            // |F_111 P_100) = |G_211 S_000)_{t} + x_cd * |F_111 S_000)_{t}
            const double Q_p_p_p_p_s_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 4] + ( CD_x[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 4] );

            // |F_111 P_010) = |G_121 S_000)_{t} + y_cd * |F_111 S_000)_{t}
            const double Q_p_p_p_s_p_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 7] + ( CD_y[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 4] );

            // |F_111 P_001) = |G_112 S_000)_{t} + z_cd * |F_111 S_000)_{t}
            const double Q_p_p_p_s_s_p_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 8] + ( CD_z[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 4] );

            // |F_102 P_100) = |G_202 S_000)_{t} + x_cd * |F_102 S_000)_{t}
            const double Q_p_s_d_p_s_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 5] + ( CD_x[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 5] );

            // |F_102 P_010) = |G_112 S_000)_{t} + y_cd * |F_102 S_000)_{t}
            const double Q_p_s_d_s_p_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 8] + ( CD_y[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 5] );

            // |F_102 P_001) = |G_103 S_000)_{t} + z_cd * |F_102 S_000)_{t}
            const double Q_p_s_d_s_s_p_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 9] + ( CD_z[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 5] );

            // |F_030 P_100) = |G_130 S_000)_{t} + x_cd * |F_030 S_000)_{t}
            const double Q_s_f_s_p_s_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 6] + ( CD_x[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 6] );

            // |F_030 P_010) = |G_040 S_000)_{t} + y_cd * |F_030 S_000)_{t}
            const double Q_s_f_s_s_p_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 10] + ( CD_y[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 6] );

            // |F_021 P_100) = |G_121 S_000)_{t} + x_cd * |F_021 S_000)_{t}
            const double Q_s_d_p_p_s_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 7] + ( CD_x[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 7] );

            // |F_021 P_010) = |G_031 S_000)_{t} + y_cd * |F_021 S_000)_{t}
            const double Q_s_d_p_s_p_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 11] + ( CD_y[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 7] );

            // |F_021 P_001) = |G_022 S_000)_{t} + z_cd * |F_021 S_000)_{t}
            const double Q_s_d_p_s_s_p_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 12] + ( CD_z[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 7] );

            // |F_012 P_100) = |G_112 S_000)_{t} + x_cd * |F_012 S_000)_{t}
            const double Q_s_p_d_p_s_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 8] + ( CD_x[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 8] );

            // |F_012 P_010) = |G_022 S_000)_{t} + y_cd * |F_012 S_000)_{t}
            const double Q_s_p_d_s_p_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 12] + ( CD_y[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 8] );

            // |F_012 P_001) = |G_013 S_000)_{t} + z_cd * |F_012 S_000)_{t}
            const double Q_s_p_d_s_s_p_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 13] + ( CD_z[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 8] );

            // |F_003 P_100) = |G_103 S_000)_{t} + x_cd * |F_003 S_000)_{t}
            const double Q_s_s_f_p_s_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 9] + ( CD_x[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 9] );

            // |F_003 P_010) = |G_013 S_000)_{t} + y_cd * |F_003 S_000)_{t}
            const double Q_s_s_f_s_p_s_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 13] + ( CD_y[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 9] );

            // |F_003 P_001) = |G_004 S_000)_{t} + z_cd * |F_003 S_000)_{t}
            const double Q_s_s_f_s_s_p_d_d = INT__d_d_g_s[abcd * 540 + ibra * 15 + 14] + ( CD_z[abcd] * INT__d_d_f_s[abcd * 360 + ibra * 10 + 9] );

            // |D_200 D_200)_{i} = |F_300 P_100) + x_cd * |D_200 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 0] = Q_f_s_s_p_s_s_d_d + ( CD_x[abcd] * Q_d_s_s_p_s_s_d_d );

            // |D_200 D_110)_{i} = |F_210 P_100) + y_cd * |D_200 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 1] = Q_d_p_s_p_s_s_d_d + ( CD_y[abcd] * Q_d_s_s_p_s_s_d_d );

            // |D_200 D_101)_{i} = |F_201 P_100) + z_cd * |D_200 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 2] = Q_d_s_p_p_s_s_d_d + ( CD_z[abcd] * Q_d_s_s_p_s_s_d_d );

            // |D_200 D_020)_{i} = |F_210 P_010) + y_cd * |D_200 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 3] = Q_d_p_s_s_p_s_d_d + ( CD_y[abcd] * Q_d_s_s_s_p_s_d_d );

            // |D_200 D_011)_{i} = |F_201 P_010) + z_cd * |D_200 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 4] = Q_d_s_p_s_p_s_d_d + ( CD_z[abcd] * Q_d_s_s_s_p_s_d_d );

            // |D_200 D_002)_{i} = |F_201 P_001) + z_cd * |D_200 P_001)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 5] = Q_d_s_p_s_s_p_d_d + ( CD_z[abcd] * Q_d_s_s_s_s_p_d_d );

            // |D_110 D_200)_{i} = |F_210 P_100) + x_cd * |D_110 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 6] = Q_d_p_s_p_s_s_d_d + ( CD_x[abcd] * Q_p_p_s_p_s_s_d_d );

            // |D_110 D_110)_{i} = |F_120 P_100) + y_cd * |D_110 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 7] = Q_p_d_s_p_s_s_d_d + ( CD_y[abcd] * Q_p_p_s_p_s_s_d_d );

            // |D_110 D_101)_{i} = |F_111 P_100) + z_cd * |D_110 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 8] = Q_p_p_p_p_s_s_d_d + ( CD_z[abcd] * Q_p_p_s_p_s_s_d_d );

            // |D_110 D_020)_{i} = |F_120 P_010) + y_cd * |D_110 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 9] = Q_p_d_s_s_p_s_d_d + ( CD_y[abcd] * Q_p_p_s_s_p_s_d_d );

            // |D_110 D_011)_{i} = |F_111 P_010) + z_cd * |D_110 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 10] = Q_p_p_p_s_p_s_d_d + ( CD_z[abcd] * Q_p_p_s_s_p_s_d_d );

            // |D_110 D_002)_{i} = |F_111 P_001) + z_cd * |D_110 P_001)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 11] = Q_p_p_p_s_s_p_d_d + ( CD_z[abcd] * Q_p_p_s_s_s_p_d_d );

            // |D_101 D_200)_{i} = |F_201 P_100) + x_cd * |D_101 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 12] = Q_d_s_p_p_s_s_d_d + ( CD_x[abcd] * Q_p_s_p_p_s_s_d_d );

            // |D_101 D_110)_{i} = |F_111 P_100) + y_cd * |D_101 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 13] = Q_p_p_p_p_s_s_d_d + ( CD_y[abcd] * Q_p_s_p_p_s_s_d_d );

            // |D_101 D_101)_{i} = |F_102 P_100) + z_cd * |D_101 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 14] = Q_p_s_d_p_s_s_d_d + ( CD_z[abcd] * Q_p_s_p_p_s_s_d_d );

            // |D_101 D_020)_{i} = |F_111 P_010) + y_cd * |D_101 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 15] = Q_p_p_p_s_p_s_d_d + ( CD_y[abcd] * Q_p_s_p_s_p_s_d_d );

            // |D_101 D_011)_{i} = |F_102 P_010) + z_cd * |D_101 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 16] = Q_p_s_d_s_p_s_d_d + ( CD_z[abcd] * Q_p_s_p_s_p_s_d_d );

            // |D_101 D_002)_{i} = |F_102 P_001) + z_cd * |D_101 P_001)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 17] = Q_p_s_d_s_s_p_d_d + ( CD_z[abcd] * Q_p_s_p_s_s_p_d_d );

            // |D_020 D_200)_{i} = |F_120 P_100) + x_cd * |D_020 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 18] = Q_p_d_s_p_s_s_d_d + ( CD_x[abcd] * Q_s_d_s_p_s_s_d_d );

            // |D_020 D_110)_{i} = |F_030 P_100) + y_cd * |D_020 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 19] = Q_s_f_s_p_s_s_d_d + ( CD_y[abcd] * Q_s_d_s_p_s_s_d_d );

            // |D_020 D_101)_{i} = |F_021 P_100) + z_cd * |D_020 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 20] = Q_s_d_p_p_s_s_d_d + ( CD_z[abcd] * Q_s_d_s_p_s_s_d_d );

            // |D_020 D_020)_{i} = |F_030 P_010) + y_cd * |D_020 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 21] = Q_s_f_s_s_p_s_d_d + ( CD_y[abcd] * Q_s_d_s_s_p_s_d_d );

            // |D_020 D_011)_{i} = |F_021 P_010) + z_cd * |D_020 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 22] = Q_s_d_p_s_p_s_d_d + ( CD_z[abcd] * Q_s_d_s_s_p_s_d_d );

            // |D_020 D_002)_{i} = |F_021 P_001) + z_cd * |D_020 P_001)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 23] = Q_s_d_p_s_s_p_d_d + ( CD_z[abcd] * Q_s_d_s_s_s_p_d_d );

            // |D_011 D_200)_{i} = |F_111 P_100) + x_cd * |D_011 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 24] = Q_p_p_p_p_s_s_d_d + ( CD_x[abcd] * Q_s_p_p_p_s_s_d_d );

            // |D_011 D_110)_{i} = |F_021 P_100) + y_cd * |D_011 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 25] = Q_s_d_p_p_s_s_d_d + ( CD_y[abcd] * Q_s_p_p_p_s_s_d_d );

            // |D_011 D_101)_{i} = |F_012 P_100) + z_cd * |D_011 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 26] = Q_s_p_d_p_s_s_d_d + ( CD_z[abcd] * Q_s_p_p_p_s_s_d_d );

            // |D_011 D_020)_{i} = |F_021 P_010) + y_cd * |D_011 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 27] = Q_s_d_p_s_p_s_d_d + ( CD_y[abcd] * Q_s_p_p_s_p_s_d_d );

            // |D_011 D_011)_{i} = |F_012 P_010) + z_cd * |D_011 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 28] = Q_s_p_d_s_p_s_d_d + ( CD_z[abcd] * Q_s_p_p_s_p_s_d_d );

            // |D_011 D_002)_{i} = |F_012 P_001) + z_cd * |D_011 P_001)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 29] = Q_s_p_d_s_s_p_d_d + ( CD_z[abcd] * Q_s_p_p_s_s_p_d_d );

            // |D_002 D_200)_{i} = |F_102 P_100) + x_cd * |D_002 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 30] = Q_p_s_d_p_s_s_d_d + ( CD_x[abcd] * Q_s_s_d_p_s_s_d_d );

            // |D_002 D_110)_{i} = |F_012 P_100) + y_cd * |D_002 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 31] = Q_s_p_d_p_s_s_d_d + ( CD_y[abcd] * Q_s_s_d_p_s_s_d_d );

            // |D_002 D_101)_{i} = |F_003 P_100) + z_cd * |D_002 P_100)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 32] = Q_s_s_f_p_s_s_d_d + ( CD_z[abcd] * Q_s_s_d_p_s_s_d_d );

            // |D_002 D_020)_{i} = |F_012 P_010) + y_cd * |D_002 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 33] = Q_s_p_d_s_p_s_d_d + ( CD_y[abcd] * Q_s_s_d_s_p_s_d_d );

            // |D_002 D_011)_{i} = |F_003 P_010) + z_cd * |D_002 P_010)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 34] = Q_s_s_f_s_p_s_d_d + ( CD_z[abcd] * Q_s_s_d_s_p_s_d_d );

            // |D_002 D_002)_{i} = |F_003 P_001) + z_cd * |D_002 P_001)
            INT__d_d_d_d[abcd * 1296 + ibra * 36 + 35] = Q_s_s_f_s_s_p_d_d + ( CD_z[abcd] * Q_s_s_d_s_s_p_d_d );

        }
    }


    // Free contracted work space
    free(contwork);

    return nshell1234;
}

