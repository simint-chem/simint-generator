#include <string.h>
#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"
#include "boys/boys_split.h"


int eri_split_d_s_p_s(struct multishell_pair const P,
                      struct multishell_pair const Q,
                      double * const restrict S_2_0_1_0)
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

    memset(S_2_0_1_0, 0, nshell1234*18*sizeof(double));

    int ab, cd, abcd;
    int i, j;

    // Workspace for contracted integrals
    double * const contwork = malloc(nshell1234 * 0);
    memset(contwork, 0, nshell1234 * 0);

    // partition workspace into shells


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
        double * const restrict PRIM_S_2_0_1_0 = S_2_0_1_0 + (abcd * 18);

            const int cdstart = Q.primstart[cd];
            const int cdend = Q.primend[cd];

            // this should have been set/aligned in fill_multishell_pair or something else
            ASSUME(cdstart%SIMD_ALIGN_DBL == 0);

            for(i = abstart; i < abend; ++i)
            {
                for(j = cdstart; j < cdend; ++j)
                {

                    // Holds the auxiliary integrals ( i 0 | 0 0 )^m in the primitive basis
                    // with m as the slowest index
                    // AM = 0: Needed from this AM: 1
                    double AUX_S_0_0_0_0[4 * 1];

                    // AM = 1: Needed from this AM: 3
                    double AUX_S_1_0_0_0[3 * 3];

                    // AM = 2: Needed from this AM: 6
                    double AUX_S_2_0_0_0[2 * 6];

                    // AM = 3: Needed from this AM: 10
                    double AUX_S_3_0_0_0[1 * 10];



                    // Holds temporary integrals for electron transfer
                    double AUX_S_2_0_1_0[18];


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
                    // Maximum v value: 3
                    //////////////////////////////////////////////
                    // The paremeter to the boys function
                    const double F_x = R2 * alpha;


                    Boys_F_split(AUX_S_0_0_0_0, 3, F_x);
                    AUX_S_0_0_0_0[0] *= allprefac;
                    AUX_S_0_0_0_0[1] *= allprefac;
                    AUX_S_0_0_0_0[2] *= allprefac;
                    AUX_S_0_0_0_0[3] *= allprefac;

                    //////////////////////////////////////////////
                    // Primitive integrals: Vertical recurrance
                    //////////////////////////////////////////////

                    // Forming AUX_S_1_0_0_0[3 * 3];
                    // Needed from this AM:
                    //    P_100
                    //    P_010
                    //    P_001
                    for(int m = 0; m < 3; m++)  // loop over orders of boys function
                    {
                        //P_100 : STEP: x
                        AUX_S_1_0_0_0[m * 3 + 0] = P.PA_x[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_x * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                        //P_010 : STEP: y
                        AUX_S_1_0_0_0[m * 3 + 1] = P.PA_y[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_y * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                        //P_001 : STEP: z
                        AUX_S_1_0_0_0[m * 3 + 2] = P.PA_z[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_z * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                    }


                    // Forming AUX_S_2_0_0_0[2 * 6];
                    // Needed from this AM:
                    //    D_200
                    //    D_110
                    //    D_101
                    //    D_020
                    //    D_011
                    //    D_002
                    for(int m = 0; m < 2; m++)  // loop over orders of boys function
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


                    // Forming AUX_S_3_0_0_0[1 * 10];
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
                    for(int m = 0; m < 1; m++)  // loop over orders of boys function
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


                    // Accumulating in contracted workspace
                    for(int i = 0; i < 18; i++)
                        PRIM_S_2_0_1_0[i] += AUX_S_2_0_1_0[i];

                 }
            }
        }
    }


    //////////////////////////////////////////////
    // Contracted integrals: Horizontal recurrance
    // Bra part
    // Steps: 0
    //////////////////////////////////////////////

    // Nothing to do.....


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

