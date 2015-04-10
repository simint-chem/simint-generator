#include <math.h>
#include <string.h> // for memset

#include "constants.h"
#include "boys/boys.h"
#include "eri/shell.h"

#define F0_KFAC 0.88622692545275801364908374  // sqrt(pi)/2

int eri_ssss_combined(struct multishell_pair const P,
                      struct multishell_pair const Q,
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

    const int nshell1234 = P.nshell12 * Q.nshell12;

    memset(integrals, 0, nshell1234*sizeof(double));

    int ab, cd;
    int i, j;
    int nint = 0;

    for(ab = 0; ab < P.nshell12; ++ab)
    {
        const int abstart = P.primstart[ab];
        const int abend = P.primend[ab];

        // this should have been set/aligned in fill_multishell_pair or something else
        ASSUME(abstart%SIMD_ALIGN_DBL == 0);

        for(cd = 0; cd < Q.nshell12; ++cd)
        {
            const int cdstart = Q.primstart[cd];
            const int cdend = Q.primend[cd];

            // this should have been set/aligned in fill_multishell_pair or something else
            ASSUME(cdstart%SIMD_ALIGN_DBL == 0);

            for(i = abstart; i < abend; ++i)
            {
                for(j = cdstart; j < cdend; ++j)
                {
                    const double PQalpha_mul = P.alpha[i] * Q.alpha[j];
                    const double PQalpha_sum = P.alpha[i] + Q.alpha[j];

                    const double pfac = TWO_PI_52 / (PQalpha_mul * sqrt(PQalpha_sum));

                    /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
                    const double PQ_x = P.x[i] - Q.x[j];
                    const double PQ_y = P.y[i] - Q.y[j];
                    const double PQ_z = P.z[i] - Q.z[j];
                    const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;

                    // The paremeter to the boys function
                    const double x = R2 * PQalpha_mul/PQalpha_sum;
                    integrals[nint] += pfac * P.prefac[i] * Q.prefac[j] * Boys_F0_FO(x);
                 }
            }

            ++nint;

        }
    }

    return nshell1234;
}
