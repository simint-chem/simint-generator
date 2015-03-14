#include <math.h>

#include "constants.h"
#include "vectorization.h"

#include "eri/shell.h"
#include "boys/boys_grid.h"

#include <stdio.h>

extern double boys_grid[BOYS_GRID_NPOINT][BOYS_GRID_MAXN + 1];

#include <string.h> // for memset

int eri_ssss_combined(struct shell_pair const P,
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
    int nint = 0;
    const int nshell1234 = P.nshell12 * Q.nshell12;

    // we don't want to align integralwork1,2 by shell quartet, etc, since
    // we just want to plow through then later in a 'flat' way

    memset(integrals, 0, nshell1234*sizeof(double));
    for(ab = 0; ab < P.nshell12; ++ab)
    {
        const int abstart = P.primstart[ab];
        const int abend = P.primend[ab];

        // this should have been set in fill_shell_pair or something else
        ASSUME(abstart%SIMD_ALIGN_DBL == 0);

        for(cd = 0; cd < Q.nshell12; ++cd)
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
                    const double F0param = sqrt(R2 * PQalpha_mul/PQalpha_sum);
                    integrals[nint] += pfac * P.prefac[i] * Q.prefac[j] * erf(F0param) / F0param; 
                 }
            }
            ++nint;
        }
    }

    // apply constants to integrals
    // also heavily vectorized
    for(i = 0; i < nshell1234; ++i)
        integrals[i] *= F0_KFAC * ONESIX_OVER_SQRT_PI;

    return nshell1234;

}

