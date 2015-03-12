#include <math.h>

#include "constants.h"
#include "eri/shell.h"
#include "vectorization.h"

#include <string.h> // for memset

int eri_ssss_flat(struct shell_pair const P,
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

    int i, j;
    int ab, cd;
    int nint;
    int idx = 0;
    const int nshell1234 = P.nshell12 * Q.nshell12;

    // we don't want to align integralwork1,2 by shell quartet, etc, since
    // we just want to plow through then later in a 'flat' way

    for(ab = 0; ab < P.nshell12; ++ab)
    {
        const int abstart = P.primstart[ab];
        const int abend = P.primend[ab];

        // this should have been set in fill_shell_pair or something else
        ASSUME(abstart%SIMD_ALIGN_DBL == 0);

        for(cd = 0; cd < Q.nshell12; cd++)
        {
            for(i = abstart; i < abend; ++i)
            {
                const int cdstart = Q.primstart[cd];
                const int cdend = Q.primend[cd];

                // this should have been set in fill_shell_pair or something else
                ASSUME(cdstart%SIMD_ALIGN_DBL == 0);

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

                    // store the paremeter to the boys function in integralwork1
                    integralwork1[idx] = R2 * PQalpha_mul/PQalpha_sum;

                    // store the prefactors in integralwork2
                    integralwork2[idx] = pfac * P.prefac[i] * Q.prefac[j];

                    ++idx;
                 }
            }
        }
    }

    // rip through the integral work arrays and store result back in integralwork1
    // This is the loop that should be heavily vectorized
    for(i = 0; i < idx; ++i)
    {
        const double x2 = sqrt(integralwork1[i]);
        integralwork1[i] = integralwork2[i] * erf(x2) / x2;
    }

    // now sum them, forming the contracted integrals
    memset(integrals, 0, nshell1234);
    idx = 0;
    nint = 0;
    for(ab = 0; ab < P.nshell12; ++ab)
    for(cd = 0; cd < Q.nshell12; cd++)
    {
        const int nprim1234 = P.nprim12[ab] * Q.nprim12[cd];
        for(i = 0; i < nprim1234; ++i)
        {
            integrals[nint] += integralwork1[idx];
            ++idx;
        }

        nint++;
    }

    // apply constants to integrals
    // also heavily vectorized
    for(i = 0; i < nshell1234; ++i)
        integrals[i] *= F0_KFAC * ONESIX_OVER_SQRT_PI;

    return nshell1234;

}


