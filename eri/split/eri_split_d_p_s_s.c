#include <string.h>
#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"
#include "boys/boys_split.h"


int eri_split_d_p_s_s(struct multishell_pair const P,
                      struct multishell_pair const Q,
                      double * const restrict INT__d_p_s_s)
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
    double * const contwork = malloc(nshell1234 * 128);
    memset(contwork, 0, nshell1234 * 128);

    // partition workspace into shells
    double * const INT__d_s_s_s = contwork + (nshell1234 * 0);
    double * const INT__f_s_s_s = contwork + (nshell1234 * 6);


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
            double * const restrict PRIM_INT__d_s_s_s = INT__d_s_s_s + (abcd * 6);
            double * const restrict PRIM_INT__f_s_s_s = INT__f_s_s_s + (abcd * 10);

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
                    double AUX_INT__s_s_s_s[4 * 1];

                    // AM = 1: Needed from this AM: 3
                    double AUX_INT__p_s_s_s[3 * 3];

                    // AM = 2: Needed from this AM: 6
                    double AUX_INT__d_s_s_s[2 * 6];

                    // AM = 3: Needed from this AM: 10
                    double AUX_INT__f_s_s_s[1 * 10];



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
                    // Maximum v value: 3
                    //////////////////////////////////////////////
                    // The paremeter to the boys function
                    const double F_x = R2 * alpha;


                    Boys_F_split(AUX_INT__s_s_s_s, 3, F_x);
                    AUX_INT__s_s_s_s[0] *= allprefac;
                    AUX_INT__s_s_s_s[1] *= allprefac;
                    AUX_INT__s_s_s_s[2] *= allprefac;
                    AUX_INT__s_s_s_s[3] *= allprefac;

                    //////////////////////////////////////////////
                    // Primitive integrals: Vertical recurrance
                    //////////////////////////////////////////////

                    // Forming AUX_INT__p_s_s_s[3 * 3];
                    // Needed from this AM:
                    //    P_100
                    //    P_010
                    //    P_001
                    for(int m = 0; m < 3; m++)  // loop over orders of boys function
                    {
                        //P_100 : STEP: x
                        AUX_INT__p_s_s_s[m * 3 + 0] = P.PA_x[i] * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_x * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_010 : STEP: y
                        AUX_INT__p_s_s_s[m * 3 + 1] = P.PA_y[i] * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_y * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_001 : STEP: z
                        AUX_INT__p_s_s_s[m * 3 + 2] = P.PA_z[i] * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_z * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                    }


                    // Forming AUX_INT__d_s_s_s[2 * 6];
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
                        AUX_INT__d_s_s_s[m * 6 + 0] = P.PA_x[i] * AUX_INT__p_s_s_s[m * 3 + 0] - a_over_p * PQ_x * AUX_INT__p_s_s_s[(m+1) * 3 + 0]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                        //D_110 : STEP: y
                        AUX_INT__d_s_s_s[m * 6 + 1] = P.PA_y[i] * AUX_INT__p_s_s_s[m * 3 + 0] - a_over_p * PQ_y * AUX_INT__p_s_s_s[(m+1) * 3 + 0];

                        //D_101 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 2] = P.PA_z[i] * AUX_INT__p_s_s_s[m * 3 + 0] - a_over_p * PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 0];

                        //D_020 : STEP: y
                        AUX_INT__d_s_s_s[m * 6 + 3] = P.PA_y[i] * AUX_INT__p_s_s_s[m * 3 + 1] - a_over_p * PQ_y * AUX_INT__p_s_s_s[(m+1) * 3 + 1]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                        //D_011 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 4] = P.PA_z[i] * AUX_INT__p_s_s_s[m * 3 + 1] - a_over_p * PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 1];

                        //D_002 : STEP: z
                        AUX_INT__d_s_s_s[m * 6 + 5] = P.PA_z[i] * AUX_INT__p_s_s_s[m * 3 + 2] - a_over_p * PQ_z * AUX_INT__p_s_s_s[(m+1) * 3 + 2]
                                      + 1 * one_over_2p * ( AUX_INT__s_s_s_s[m * 1 +  0] - a_over_p * AUX_INT__s_s_s_s[(m+1) * 1 + 0] );

                    }

                    // Accumulating in contracted workspace
                    for(int i = 0; i < 6; i++)
                        PRIM_INT__d_s_s_s[i] += AUX_INT__d_s_s_s[i];


                    // Forming AUX_INT__f_s_s_s[1 * 10];
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
                        AUX_INT__f_s_s_s[m * 10 + 0] = P.PA_x[i] * AUX_INT__d_s_s_s[m * 6 + 0] - a_over_p * PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 0]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  0] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 0] );

                        //F_210 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 1] = P.PA_y[i] * AUX_INT__d_s_s_s[m * 6 + 0] - a_over_p * PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 0];

                        //F_201 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 2] = P.PA_z[i] * AUX_INT__d_s_s_s[m * 6 + 0] - a_over_p * PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 0];

                        //F_120 : STEP: x
                        AUX_INT__f_s_s_s[m * 10 + 3] = P.PA_x[i] * AUX_INT__d_s_s_s[m * 6 + 3] - a_over_p * PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 3];

                        //F_111 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 4] = P.PA_z[i] * AUX_INT__d_s_s_s[m * 6 + 1] - a_over_p * PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 1];

                        //F_102 : STEP: x
                        AUX_INT__f_s_s_s[m * 10 + 5] = P.PA_x[i] * AUX_INT__d_s_s_s[m * 6 + 5] - a_over_p * PQ_x * AUX_INT__d_s_s_s[(m+1) * 6 + 5];

                        //F_030 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 6] = P.PA_y[i] * AUX_INT__d_s_s_s[m * 6 + 3] - a_over_p * PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 3]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  1] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 1] );

                        //F_021 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 7] = P.PA_z[i] * AUX_INT__d_s_s_s[m * 6 + 3] - a_over_p * PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 3];

                        //F_012 : STEP: y
                        AUX_INT__f_s_s_s[m * 10 + 8] = P.PA_y[i] * AUX_INT__d_s_s_s[m * 6 + 5] - a_over_p * PQ_y * AUX_INT__d_s_s_s[(m+1) * 6 + 5];

                        //F_003 : STEP: z
                        AUX_INT__f_s_s_s[m * 10 + 9] = P.PA_z[i] * AUX_INT__d_s_s_s[m * 6 + 5] - a_over_p * PQ_z * AUX_INT__d_s_s_s[(m+1) * 6 + 5]
                                      + 2 * one_over_2p * ( AUX_INT__p_s_s_s[m * 3 +  2] - a_over_p * AUX_INT__p_s_s_s[(m+1) * 3 + 2] );

                    }

                    // Accumulating in contracted workspace
                    for(int i = 0; i < 10; i++)
                        PRIM_INT__f_s_s_s[i] += AUX_INT__f_s_s_s[i];




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
    // Steps: 18
    //////////////////////////////////////////////

    for(abcd = 0; abcd < nshell1234; ++abcd)
    {
        // form INT__d_p_s_s
        for(int iket = 0; iket < 1; ++iket)
        {
            // (D_200 P_100|_{i} = (F_300 S_000|_{t} + x_ab * (D_200 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 0 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 0 * 1 + iket] + ( AB_x[abcd] * INT__d_s_s_s[abcd * 6 + 0 * 1 + iket] );

            // (D_200 P_010|_{i} = (F_210 S_000|_{t} + y_ab * (D_200 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 1 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 1 * 1 + iket] + ( AB_y[abcd] * INT__d_s_s_s[abcd * 6 + 0 * 1 + iket] );

            // (D_200 P_001|_{i} = (F_201 S_000|_{t} + z_ab * (D_200 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 2 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 2 * 1 + iket] + ( AB_z[abcd] * INT__d_s_s_s[abcd * 6 + 0 * 1 + iket] );

            // (D_110 P_100|_{i} = (F_210 S_000|_{t} + x_ab * (D_110 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 3 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 1 * 1 + iket] + ( AB_x[abcd] * INT__d_s_s_s[abcd * 6 + 1 * 1 + iket] );

            // (D_110 P_010|_{i} = (F_120 S_000|_{t} + y_ab * (D_110 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 4 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 3 * 1 + iket] + ( AB_y[abcd] * INT__d_s_s_s[abcd * 6 + 1 * 1 + iket] );

            // (D_110 P_001|_{i} = (F_111 S_000|_{t} + z_ab * (D_110 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 5 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 4 * 1 + iket] + ( AB_z[abcd] * INT__d_s_s_s[abcd * 6 + 1 * 1 + iket] );

            // (D_101 P_100|_{i} = (F_201 S_000|_{t} + x_ab * (D_101 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 6 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 2 * 1 + iket] + ( AB_x[abcd] * INT__d_s_s_s[abcd * 6 + 2 * 1 + iket] );

            // (D_101 P_010|_{i} = (F_111 S_000|_{t} + y_ab * (D_101 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 7 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 4 * 1 + iket] + ( AB_y[abcd] * INT__d_s_s_s[abcd * 6 + 2 * 1 + iket] );

            // (D_101 P_001|_{i} = (F_102 S_000|_{t} + z_ab * (D_101 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 8 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 5 * 1 + iket] + ( AB_z[abcd] * INT__d_s_s_s[abcd * 6 + 2 * 1 + iket] );

            // (D_020 P_100|_{i} = (F_120 S_000|_{t} + x_ab * (D_020 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 9 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 3 * 1 + iket] + ( AB_x[abcd] * INT__d_s_s_s[abcd * 6 + 3 * 1 + iket] );

            // (D_020 P_010|_{i} = (F_030 S_000|_{t} + y_ab * (D_020 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 10 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 6 * 1 + iket] + ( AB_y[abcd] * INT__d_s_s_s[abcd * 6 + 3 * 1 + iket] );

            // (D_020 P_001|_{i} = (F_021 S_000|_{t} + z_ab * (D_020 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 11 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 7 * 1 + iket] + ( AB_z[abcd] * INT__d_s_s_s[abcd * 6 + 3 * 1 + iket] );

            // (D_011 P_100|_{i} = (F_111 S_000|_{t} + x_ab * (D_011 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 12 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 4 * 1 + iket] + ( AB_x[abcd] * INT__d_s_s_s[abcd * 6 + 4 * 1 + iket] );

            // (D_011 P_010|_{i} = (F_021 S_000|_{t} + y_ab * (D_011 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 13 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 7 * 1 + iket] + ( AB_y[abcd] * INT__d_s_s_s[abcd * 6 + 4 * 1 + iket] );

            // (D_011 P_001|_{i} = (F_012 S_000|_{t} + z_ab * (D_011 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 14 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 8 * 1 + iket] + ( AB_z[abcd] * INT__d_s_s_s[abcd * 6 + 4 * 1 + iket] );

            // (D_002 P_100|_{i} = (F_102 S_000|_{t} + x_ab * (D_002 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 15 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 5 * 1 + iket] + ( AB_x[abcd] * INT__d_s_s_s[abcd * 6 + 5 * 1 + iket] );

            // (D_002 P_010|_{i} = (F_012 S_000|_{t} + y_ab * (D_002 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 16 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 8 * 1 + iket] + ( AB_y[abcd] * INT__d_s_s_s[abcd * 6 + 5 * 1 + iket] );

            // (D_002 P_001|_{i} = (F_003 S_000|_{t} + z_ab * (D_002 S_000|_{t}
            INT__d_p_s_s[abcd * 18 + 17 * 1 + iket] = INT__f_s_s_s[abcd * 10 + 9 * 1 + iket] + ( AB_z[abcd] * INT__d_s_s_s[abcd * 6 + 5 * 1 + iket] );

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

