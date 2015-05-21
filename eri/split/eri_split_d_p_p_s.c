#include <string.h>
#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"
#include "boys/boys_split.h"


int eri_split_d_p_p_s(struct multishell_pair const P,
                      struct multishell_pair const Q,
                      double * const restrict INT__d_p_p_s)
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

    ASSUME_ALIGN(INT__d_p_p_s);

    const int nshell1234 = P.nshell12 * Q.nshell12;


    // Holds AB_{xyz} and CD_{xyz} in a flattened fashion for later
    double AB_x[nshell1234];  double CD_x[nshell1234];
    double AB_y[nshell1234];  double CD_y[nshell1234];
    double AB_z[nshell1234];  double CD_z[nshell1234];

    int ab, cd, abcd;
    int i, j;

    // Workspace for contracted integrals
    double * const contwork = malloc(nshell1234 * 384);
    memset(contwork, 0, nshell1234 * 384);

    // partition workspace into shells
    double * const INT__d_s_p_s = contwork + (nshell1234 * 0);
    double * const INT__f_s_p_s = contwork + (nshell1234 * 18);


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
            double * const restrict PRIM_INT__f_s_p_s = INT__f_s_p_s + (abcd * 30);

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
                    double AUX_INT__s_s_s_s[5 * 1];

                    // AM = 1: Needed from this AM: 3
                    double AUX_INT__p_s_s_s[4 * 3];

                    // AM = 2: Needed from this AM: 6
                    double AUX_INT__d_s_s_s[3 * 6];

                    // AM = 3: Needed from this AM: 10
                    double AUX_INT__f_s_s_s[2 * 10];

                    // AM = 4: Needed from this AM: 15
                    double AUX_INT__g_s_s_s[1 * 15];



                    // Holds temporary integrals for electron transfer
                    double AUX_INT__d_s_p_s[18];
                    double AUX_INT__f_s_p_s[30];


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
                    // Maximum v value: 4
                    //////////////////////////////////////////////
                    // The paremeter to the boys function
                    const double F_x = R2 * alpha;


                    Boys_F_split(AUX_INT__s_s_s_s, 4, F_x);
                    AUX_INT__s_s_s_s[0] *= allprefac;
                    AUX_INT__s_s_s_s[1] *= allprefac;
                    AUX_INT__s_s_s_s[2] *= allprefac;
                    AUX_INT__s_s_s_s[3] *= allprefac;
                    AUX_INT__s_s_s_s[4] *= allprefac;

                    //////////////////////////////////////////////
                    // Primitive integrals: Vertical recurrance
                    //////////////////////////////////////////////

                    // Forming AUX_INT__p_s_s_s[4 * 3];
                    // Needed from this AM:
                    //    P_100
                    //    P_010
                    //    P_001
                    for(int m = 0; m < 4; m++)  // loop over orders of boys function
                    {
                        //P_100 : STEP: x
                        AUX_INT__p_s_s_s[m * 3 + 0] = P_PA_x * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_x * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_010 : STEP: y
                        AUX_INT__p_s_s_s[m * 3 + 1] = P_PA_y * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_y * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_001 : STEP: z
                        AUX_INT__p_s_s_s[m * 3 + 2] = P_PA_z * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_z * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                    }


                    // Forming AUX_INT__d_s_s_s[3 * 6];
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


                    // Forming AUX_INT__f_s_s_s[2 * 10];
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


                    // Forming AUX_INT__g_s_s_s[1 * 15];
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


                    // Accumulating in contracted workspace
                    for(int n = 0; n < 18; n++)
                        PRIM_INT__d_s_p_s[n] += AUX_INT__d_s_p_s[n];

                    // Accumulating in contracted workspace
                    for(int n = 0; n < 30; n++)
                        PRIM_INT__f_s_p_s[n] += AUX_INT__f_s_p_s[n];

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
        // form INT__d_p_p_s
        for(int iket = 0; iket < 3; ++iket)
        {
            // (D_200 P_100|_{i} = (F_300 S_000|_{t} + x_ab * (D_200 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 0 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 0 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 0 * 3 + iket] );

            // (D_200 P_010|_{i} = (F_210 S_000|_{t} + y_ab * (D_200 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 1 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 1 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 0 * 3 + iket] );

            // (D_200 P_001|_{i} = (F_201 S_000|_{t} + z_ab * (D_200 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 2 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 2 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 0 * 3 + iket] );

            // (D_110 P_100|_{i} = (F_210 S_000|_{t} + x_ab * (D_110 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 3 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 1 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 1 * 3 + iket] );

            // (D_110 P_010|_{i} = (F_120 S_000|_{t} + y_ab * (D_110 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 4 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 3 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 1 * 3 + iket] );

            // (D_110 P_001|_{i} = (F_111 S_000|_{t} + z_ab * (D_110 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 5 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 4 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 1 * 3 + iket] );

            // (D_101 P_100|_{i} = (F_201 S_000|_{t} + x_ab * (D_101 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 6 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 2 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 2 * 3 + iket] );

            // (D_101 P_010|_{i} = (F_111 S_000|_{t} + y_ab * (D_101 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 7 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 4 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 2 * 3 + iket] );

            // (D_101 P_001|_{i} = (F_102 S_000|_{t} + z_ab * (D_101 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 8 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 5 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 2 * 3 + iket] );

            // (D_020 P_100|_{i} = (F_120 S_000|_{t} + x_ab * (D_020 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 9 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 3 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 3 * 3 + iket] );

            // (D_020 P_010|_{i} = (F_030 S_000|_{t} + y_ab * (D_020 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 10 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 6 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 3 * 3 + iket] );

            // (D_020 P_001|_{i} = (F_021 S_000|_{t} + z_ab * (D_020 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 11 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 7 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 3 * 3 + iket] );

            // (D_011 P_100|_{i} = (F_111 S_000|_{t} + x_ab * (D_011 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 12 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 4 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 4 * 3 + iket] );

            // (D_011 P_010|_{i} = (F_021 S_000|_{t} + y_ab * (D_011 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 13 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 7 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 4 * 3 + iket] );

            // (D_011 P_001|_{i} = (F_012 S_000|_{t} + z_ab * (D_011 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 14 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 8 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 4 * 3 + iket] );

            // (D_002 P_100|_{i} = (F_102 S_000|_{t} + x_ab * (D_002 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 15 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 5 * 3 + iket] + ( AB_x[abcd] * INT__d_s_p_s[abcd * 18 + 5 * 3 + iket] );

            // (D_002 P_010|_{i} = (F_012 S_000|_{t} + y_ab * (D_002 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 16 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 8 * 3 + iket] + ( AB_y[abcd] * INT__d_s_p_s[abcd * 18 + 5 * 3 + iket] );

            // (D_002 P_001|_{i} = (F_003 S_000|_{t} + z_ab * (D_002 S_000|_{t}
            INT__d_p_p_s[abcd * 54 + 17 * 3 + iket] = INT__f_s_p_s[abcd * 30 + 9 * 3 + iket] + ( AB_z[abcd] * INT__d_s_p_s[abcd * 18 + 5 * 3 + iket] );

        }


    }




    // Free contracted work space
    free(contwork);

    return nshell1234;
}

