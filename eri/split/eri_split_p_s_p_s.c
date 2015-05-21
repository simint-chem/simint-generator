#include <string.h>
#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"
#include "boys/boys_split.h"


int eri_split_p_s_p_s(struct multishell_pair const P,
                      struct multishell_pair const Q,
                      double * const restrict INT__p_s_p_s)
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

    ASSUME_ALIGN(INT__p_s_p_s);

    const int nshell1234 = P.nshell12 * Q.nshell12;

    memset(INT__p_s_p_s, 0, nshell1234*9*sizeof(double));

    int ab, cd, abcd;
    int i, j;



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
            double * const restrict PRIM_INT__p_s_p_s = INT__p_s_p_s + (abcd * 9);

            const int cdstart = Q.primstart[cd];
            const int cdend = Q.primend[cd];

            // this should have been set/aligned in fill_multishell_pair or something else
            ASSUME(cdstart%SIMD_ALIGN_DBL == 0);

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
                    double AUX_INT__s_s_s_s[3 * 1];

                    // AM = 1: Needed from this AM: 3
                    double AUX_INT__p_s_s_s[2 * 3];

                    // AM = 2: Needed from this AM: 6
                    double AUX_INT__d_s_s_s[1 * 6];



                    // Holds temporary integrals for electron transfer
                    double AUX_INT__p_s_p_s[9];


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
                        AUX_INT__p_s_s_s[m * 3 + 0] = P_PA_x * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_x * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_010 : STEP: y
                        AUX_INT__p_s_s_s[m * 3 + 1] = P_PA_y * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_y * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                        //P_001 : STEP: z
                        AUX_INT__p_s_s_s[m * 3 + 2] = P_PA_z * AUX_INT__s_s_s_s[m * 1 + 0] - a_over_p * PQ_z * AUX_INT__s_s_s_s[(m+1) * 1 + 0];

                    }


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




                    //////////////////////////////////////////////
                    // Primitive integrals: Electron transfer
                    //////////////////////////////////////////////

                    // ( P_100 S_000 | P_100 S_000 )^0_{t} = x * ( P_100 S_000 | S_000 S_000 )^0_{e} + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( D_200 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[0] = etfac[0] * AUX_INT__p_s_s_s[0] + 1 * one_over_2q * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__d_s_s_s[0];

                    // ( P_100 S_000 | P_010 S_000 )^0_{t} = y * ( P_100 S_000 | S_000 S_000 )^0_{e} - ( D_110 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[1] = etfac[1] * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__d_s_s_s[1];

                    // ( P_100 S_000 | P_001 S_000 )^0_{t} = z * ( P_100 S_000 | S_000 S_000 )^0_{e} - ( D_101 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[2] = etfac[2] * AUX_INT__p_s_s_s[0] - p_over_q * AUX_INT__d_s_s_s[2];

                    // ( P_010 S_000 | P_100 S_000 )^0_{t} = x * ( P_010 S_000 | S_000 S_000 )^0_{e} - ( D_110 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[3] = etfac[0] * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__d_s_s_s[1];

                    // ( P_010 S_000 | P_010 S_000 )^0_{t} = y * ( P_010 S_000 | S_000 S_000 )^0_{e} + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( D_020 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[4] = etfac[1] * AUX_INT__p_s_s_s[1] + 1 * one_over_2q * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__d_s_s_s[3];

                    // ( P_010 S_000 | P_001 S_000 )^0_{t} = z * ( P_010 S_000 | S_000 S_000 )^0_{e} - ( D_011 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[5] = etfac[2] * AUX_INT__p_s_s_s[1] - p_over_q * AUX_INT__d_s_s_s[4];

                    // ( P_001 S_000 | P_100 S_000 )^0_{t} = x * ( P_001 S_000 | S_000 S_000 )^0_{e} - ( D_101 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[6] = etfac[0] * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__d_s_s_s[2];

                    // ( P_001 S_000 | P_010 S_000 )^0_{t} = y * ( P_001 S_000 | S_000 S_000 )^0_{e} - ( D_011 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[7] = etfac[1] * AUX_INT__p_s_s_s[2] - p_over_q * AUX_INT__d_s_s_s[4];

                    // ( P_001 S_000 | P_001 S_000 )^0_{t} = z * ( P_001 S_000 | S_000 S_000 )^0_{e} + ( S_000 S_000 | S_000 S_000 )^0_{e} - ( D_002 S_000 | S_000 S_000 )^0_{e}
                    AUX_INT__p_s_p_s[8] = etfac[2] * AUX_INT__p_s_s_s[2] + 1 * one_over_2q * AUX_INT__s_s_s_s[0] - p_over_q * AUX_INT__d_s_s_s[5];


                    // Accumulating in contracted workspace
                    for(int n = 0; n < 9; n++)
                        PRIM_INT__p_s_p_s[n] += AUX_INT__p_s_p_s[n];

                 }
            }
        }
    }





    return nshell1234;
}

