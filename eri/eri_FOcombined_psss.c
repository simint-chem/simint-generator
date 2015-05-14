#include <string.h>
#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"


int eri_FOcombined_psss(struct multishell_pair const P,
                        struct multishell_pair const Q,
                        double * const restrict S_1_0_0_0)
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

    memset(S_1_0_0_0, 0, nshell1234*3*sizeof(double));

    // Holds AB_{xyz} and CD_{xyz} in a flattened fashion for later
    double AB_x[nshell1234];  double CD_x[nshell1234];
    double AB_y[nshell1234];  double CD_y[nshell1234];
    double AB_z[nshell1234];  double CD_z[nshell1234];

    int ab, cd, abcd;
    int i, j;

    // Workspace for contracted integrals


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
                    double AUX_S_0_0_0_0[2 * 1];

                    // AM = 1: Needed from this AM: 3
                    double AUX_S_1_0_0_0[1 * 3];



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
                    const double one_over_p = 1.0 / P.alpha[i];
                    const double one_over_q = 1.0 / Q.alpha[j];
                    const double a_over_p =  alpha * one_over_p;     // a/p from MEST
                    const double one_over_2p = 0.5 * one_over_p;  // gets multiplied by i in VRR
                    const double one_over_2q = 0.5 * one_over_q;
                    const double p_over_q = P.alpha[i] * one_over_q;

                    // for electron transfer
                    const double etfac[3] = {
                                             -(P.bAB_x[i] + Q.bAB_x[j]) * one_over_q,
                                             -(P.bAB_y[i] + Q.bAB_y[j]) * one_over_q,
                                             -(P.bAB_z[i] + Q.bAB_z[j]) * one_over_q,
                                            };


                    //////////////////////////////////////////////
                    // Boys function section
                    // Maximum v value: 1
                    //////////////////////////////////////////////
                    // The paremeter to the boys function
                    const double F_x = R2 * alpha;


                    AUX_S_0_0_0_0[0] = allprefac
                             * pow(
                                     (
                                       (
                                                   1.
                                         + F_x * ( 0.414016243000866299
                                         + F_x * ( 0.130448682735044324
                                         + F_x * ( 0.0281490811816026161
                                         + F_x * ( 0.00462868463720416281
                                         + F_x * ( 0.00062025147610678493
                                         + F_x * ( 0.0000686770885390617984
                                         + F_x * ( 6.28488230669978749e-6
                                         + F_x * ( 5.01986197619830788e-7
                                         + F_x * ( 3.96915046153987083e-8
                                         + F_x * ( 1.14619057438675389e-9
                                         + F_x * ( 2.21422857239286206e-10
                                         + F_x * ( -3.47087628137958658e-12
                                         + F_x * ( 5.26907054399378694e-13
                                                 )))))))))))))
                                       )
                                       /
                                       (
                                                   1.0
                                         + F_x * ( 1.08068290966581548
                                         + F_x * ( 0.539792844780494916
                                         + F_x * ( 0.166084230769033217
                                         + F_x * ( 0.035790298646556066
                                         + F_x * ( 0.00589703243578234382
                                         + F_x * ( 0.000789592116312117277
                                         + F_x * ( 0.0000874456178357265776
                                         + F_x * ( 8.00211107603171383e-6
                                         + F_x * ( 6.39149165582055646e-7
                                         + F_x * ( 5.05367903369666356e-8
                                         + F_x * ( 1.45937517486452486e-9
                                         + F_x * ( 2.81924337930412797e-10
                                         + F_x * ( -4.4192569363292127e-12
                                         + F_x * ( 6.70878898061210528e-13
                                                 ))))))))))))))
                                       )
                                     ), 0.0+0.5);

                    AUX_S_0_0_0_0[1] = allprefac
                             * pow(
                                     (
                                       (
                                                   0.480749856769136127
                                         + F_x * ( 0.0757107453935371611
                                         + F_x * ( 0.0207544733443622468
                                         + F_x * ( 0.00296686757159093428
                                         + F_x * ( 0.000385086850988198076
                                         + F_x * ( 0.0000396245291118678106
                                         + F_x * ( 3.44653527568129186e-6
                                         + F_x * ( 2.60728584781887378e-7
                                         + F_x * ( 1.64651276793705422e-8
                                         + F_x * ( 1.24924114646889903e-9
                                         + F_x * ( -4.06609033385506782e-13
                                         + F_x * ( 6.69244577908522819e-12
                                         + F_x * ( -1.74969643368118421e-13
                                         + F_x * ( 9.69901476816967128e-15
                                                 )))))))))))))
                                       )
                                       /
                                       (
                                                   1.0
                                         + F_x * ( 0.557484696706444012
                                         + F_x * ( 0.16330778039059299
                                         + F_x * ( 0.0332854355906960737
                                         + F_x * ( 0.00522222434942978777
                                         + F_x * ( 0.000658578219691024581
                                         + F_x * ( 0.0000682687148667475058
                                         + F_x * ( 5.92820704830126442e-6
                                         + F_x * ( 4.486040686836368e-7
                                         + F_x * ( 2.83282406028920981e-8
                                         + F_x * ( 2.14933012694134515e-9
                                         + F_x * ( -6.99576077110757546e-13
                                         + F_x * ( 1.15144066901510615e-11
                                         + F_x * ( -3.01036676011995688e-13
                                         + F_x * ( 1.66872327689919498e-14
                                                 ))))))))))))))
                                       )
                                     ), 1.0+0.5);


                    //////////////////////////////////////////////
                    // Primitive integrals: Vertical recurrance
                    //////////////////////////////////////////////

                    int idx = 0;

                    // Forming AUX_S_1_0_0_0[1 * 3];
                    // Needed from this AM:
                    //    P_100
                    //    P_010
                    //    P_001
                    idx = 0;
                    for(int m = 0; m < 1; m++)  // loop over orders of boys function
                    {
                        //P_100 : STEP: x
                        AUX_S_1_0_0_0[idx++] = P.PA_x[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_x * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                        //P_010 : STEP: y
                        AUX_S_1_0_0_0[idx++] = P.PA_y[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_y * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                        //P_001 : STEP: z
                        AUX_S_1_0_0_0[idx++] = P.PA_z[i] * AUX_S_0_0_0_0[m * 1 + 0] - a_over_p * PQ_z * AUX_S_0_0_0_0[(m+1) * 1 + 0];

                    }

                    // Accumulating in contracted workspace
                    idx = 0;
                    S_1_0_0_0[abcd * 3 + 0] += AUX_S_1_0_0_0[0];
                    S_1_0_0_0[abcd * 3 + 1] += AUX_S_1_0_0_0[1];
                    S_1_0_0_0[abcd * 3 + 2] += AUX_S_1_0_0_0[2];




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


    return nshell1234;
}

