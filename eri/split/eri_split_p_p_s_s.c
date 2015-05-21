#include <string.h>
#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"
#include "boys/boys_split.h"


int eri_split_p_p_s_s(struct multishell_pair const P,
                      struct multishell_pair const Q,
                      double * const restrict INT__p_p_s_s)
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
    double * const contwork = malloc(nshell1234 * 72);
    memset(contwork, 0, nshell1234 * 72);

    // partition workspace into shells
    double * const INT__p_s_s_s = contwork + (nshell1234 * 0);
    double * const INT__d_s_s_s = contwork + (nshell1234 * 3);


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
            double * const restrict PRIM_INT__p_s_s_s = INT__p_s_s_s + (abcd * 3);
            double * const restrict PRIM_INT__d_s_s_s = INT__d_s_s_s + (abcd * 6);

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
                    double AUX_INT__s_s_s_s[3 * 1];

                    // AM = 1: Needed from this AM: 3
                    double AUX_INT__p_s_s_s[2 * 3];

                    // AM = 2: Needed from this AM: 6
                    double AUX_INT__d_s_s_s[1 * 6];



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
                    // Maximum v value: 2
                    //////////////////////////////////////////////
                    // The paremeter to the boys function
                    const double F_x = R2 * alpha;


                    Boys_F_split(AUX_INT__s_s_s_s, 2, F_x);
                    AUX_INT__s_s_s_s[0] *= allprefac;
                    AUX_INT__s_s_s_s[1] *= allprefac;
                    AUX_INT__s_s_s_s[2] *= allprefac;

                    //////////////////////////////////////////////
                    // Primitive integrals: Vertical recurrance
                    //////////////////////////////////////////////

                    // Forming AUX_INT__p_s_s_s[2 * 3];
                    // Needed from this AM:
                    //    P_100
                    //    P_010
                    //    P_001
                    for(int m = 0; m < 2; m++)  // loop over orders of boys function
                    {
                        //P_100 : STEP: x
                        AUX_INT__p_s_s_s[m * 3 + 0] = P.PA_x[i] * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_x * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_010 : STEP: y
                        AUX_INT__p_s_s_s[m * 3 + 1] = P.PA_y[i] * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_y * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_001 : STEP: z
                        AUX_INT__p_s_s_s[m * 3 + 2] = P.PA_z[i] * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_z * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                    }

                    // Accumulating in contracted workspace
                    for(int i = 0; i < 3; i++)
                        PRIM_INT__p_s_s_s[i] += AUX_INT__p_s_s_s[i];


                    // Forming AUX_INT__d_s_s_s[1 * 6];
                    // Needed from this AM:
                    //    D_200
                    //    D_110
                    //    D_101
                    //    D_020
                    //    D_011
                    //    D_002
                    for(int m = 0; m < 1; m++)  // loop over orders of boys function
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
    // Steps: 9
    //////////////////////////////////////////////

    for(abcd = 0; abcd < nshell1234; ++abcd)
    {
        // form INT__p_p_s_s
        for(int iket = 0; iket < 1; ++iket)
        {
            // (P_100 P_100|_{i} = (D_200 S_000|_{t} + x_ab * (P_100 S_000|_{t}
            INT__p_p_s_s[abcd * 9 + 0 * 1 + iket] = INT__d_s_s_s[abcd * 6 + 0 * 1 + iket] + ( AB_x[abcd] * INT__p_s_s_s[abcd * 3 + 0 * 1 + iket] );

            // (P_100 P_010|_{i} = (D_110 S_000|_{t} + y_ab * (P_100 S_000|_{t}
            INT__p_p_s_s[abcd * 9 + 1 * 1 + iket] = INT__d_s_s_s[abcd * 6 + 1 * 1 + iket] + ( AB_y[abcd] * INT__p_s_s_s[abcd * 3 + 0 * 1 + iket] );

            // (P_100 P_001|_{i} = (D_101 S_000|_{t} + z_ab * (P_100 S_000|_{t}
            INT__p_p_s_s[abcd * 9 + 2 * 1 + iket] = INT__d_s_s_s[abcd * 6 + 2 * 1 + iket] + ( AB_z[abcd] * INT__p_s_s_s[abcd * 3 + 0 * 1 + iket] );

            // (P_010 P_100|_{i} = (D_110 S_000|_{t} + x_ab * (P_010 S_000|_{t}
            INT__p_p_s_s[abcd * 9 + 3 * 1 + iket] = INT__d_s_s_s[abcd * 6 + 1 * 1 + iket] + ( AB_x[abcd] * INT__p_s_s_s[abcd * 3 + 1 * 1 + iket] );

            // (P_010 P_010|_{i} = (D_020 S_000|_{t} + y_ab * (P_010 S_000|_{t}
            INT__p_p_s_s[abcd * 9 + 4 * 1 + iket] = INT__d_s_s_s[abcd * 6 + 3 * 1 + iket] + ( AB_y[abcd] * INT__p_s_s_s[abcd * 3 + 1 * 1 + iket] );

            // (P_010 P_001|_{i} = (D_011 S_000|_{t} + z_ab * (P_010 S_000|_{t}
            INT__p_p_s_s[abcd * 9 + 5 * 1 + iket] = INT__d_s_s_s[abcd * 6 + 4 * 1 + iket] + ( AB_z[abcd] * INT__p_s_s_s[abcd * 3 + 1 * 1 + iket] );

            // (P_001 P_100|_{i} = (D_101 S_000|_{t} + x_ab * (P_001 S_000|_{t}
            INT__p_p_s_s[abcd * 9 + 6 * 1 + iket] = INT__d_s_s_s[abcd * 6 + 2 * 1 + iket] + ( AB_x[abcd] * INT__p_s_s_s[abcd * 3 + 2 * 1 + iket] );

            // (P_001 P_010|_{i} = (D_011 S_000|_{t} + y_ab * (P_001 S_000|_{t}
            INT__p_p_s_s[abcd * 9 + 7 * 1 + iket] = INT__d_s_s_s[abcd * 6 + 4 * 1 + iket] + ( AB_y[abcd] * INT__p_s_s_s[abcd * 3 + 2 * 1 + iket] );

            // (P_001 P_001|_{i} = (D_002 S_000|_{t} + z_ab * (P_001 S_000|_{t}
            INT__p_p_s_s[abcd * 9 + 8 * 1 + iket] = INT__d_s_s_s[abcd * 6 + 5 * 1 + iket] + ( AB_z[abcd] * INT__p_s_s_s[abcd * 3 + 2 * 1 + iket] );

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

