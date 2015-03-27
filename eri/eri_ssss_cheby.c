#include <math.h>
#include <string.h> // for memset

#include "constants.h"
#include "boys/boys.h"
#include "eri/shell.h"


int eri_ssss_cheby(struct shell_pair const P,
                   struct shell_pair const Q,
                   double * const restrict integrals,
                   double * const restrict integralwork1,
                   double * const restrict integralwork2)
{
    ASSUME_ALIGN(P.x);
    ASSUME_ALIGN(P.y);
    ASSUME_ALIGN(P.z);
    ASSUME_ALIGN(P.alpha);
    ASSUME_ALIGN(P.prefac);
    ASSUME_ALIGN(Q.x);
    ASSUME_ALIGN(Q.y);
    ASSUME_ALIGN(Q.z);
    ASSUME_ALIGN(Q.alpha);
    ASSUME_ALIGN(Q.prefac);

    ASSUME_ALIGN(integrals);
    ASSUME_ALIGN(integralwork1);
    ASSUME_ALIGN(integralwork2);

    int nint = 0;
    int i, j;
    int ab, cd;

    const int nshell1234 = P.nshell12 * Q.nshell12;

    memset(integrals, 0, nshell1234*sizeof(double));
    for(ab = 0; ab < P.nshell12; ++ab)
    {
        const int abstart = P.primstart[ab];
        const int abend = P.primend[ab];

        // this should have been set in fill_shell_pair or something else
        ASSUME(abstart%SIMD_ALIGN_DBL == 0);

        for(cd = 0; cd < Q.nshell12; ++cd)
        {
            const int cdstart = Q.primstart[cd];
            const int cdend = Q.primend[cd];

            // this should have been set in fill_shell_pair or something else
            ASSUME(cdstart%SIMD_ALIGN_DBL == 0);

            for(i = abstart; i < abend; ++i)
            {
                #pragma simd
                for(j = cdstart; j < cdend; ++j)
                {
                    const double PQalpha_mul = P.alpha[i] * Q.alpha[j];
                    const double PQalpha_sum = P.alpha[i] + Q.alpha[j];

                    const double pfac = 1.0 / (PQalpha_mul * sqrt(PQalpha_sum));

                    /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
                    const double PQ_x = P.x[i] - Q.x[j];
                    const double PQ_y = P.y[i] - Q.y[j];
                    const double PQ_z = P.z[i] - Q.z[j];
                    const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;

                    // The paremeter to the boys function
                    const double x = R2 * PQalpha_mul/PQalpha_sum;
                    const unsigned int idx = cheby_idx(x);

                    const double bshift = (int)( ((1<<idx)-1) + ((1<<(idx+1))-1) );
                    const double x2 = x - (bshift * 0.5);
   
                    const double F =        boys_chebygrid_F0[idx][0]
                                   + x2 * ( boys_chebygrid_F0[idx][1]
                                   + x2 * ( boys_chebygrid_F0[idx][2]
                                   + x2 * ( boys_chebygrid_F0[idx][3]
                                   + x2 * ( boys_chebygrid_F0[idx][4]
                                   + x2 * ( boys_chebygrid_F0[idx][5]
                                   + x2 * ( boys_chebygrid_F0[idx][6]
                                   + x2 * ( boys_chebygrid_F0[idx][7]
                                   + x2 * ( boys_chebygrid_F0[idx][8]
                                   + x2 * ( boys_chebygrid_F0[idx][9]
                                   + x2 * ( boys_chebygrid_F0[idx][10]
                                   + x2 * ( boys_chebygrid_F0[idx][11]
                                   + x2 * ( boys_chebygrid_F0[idx][12]
                                   + x2 * ( boys_chebygrid_F0[idx][13]
                                   + x2 * ( boys_chebygrid_F0[idx][14]
                                   ))))))))))))));
                    integrals[nint] += pfac * P.prefac[i] * Q.prefac[j] * F;
                 }
            }

            ++nint;

        }
    }

    // apply constants to integrals
    // also heavily vectorized
    for(i = 0; i < nshell1234; ++i)
        integrals[i] *= ONESIX_OVER_SQRT_PI;

    return nshell1234;

}

