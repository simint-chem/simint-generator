#include <math.h>

#include "constants.h"
#include "vectorization.h"

#include "eri/shell.h"
#include "boys/boys_grid.h"

#include <stdio.h>

extern double boys_grid[BOYS_GRID_NPOINT][BOYS_GRID_MAXN + 1];

#include <string.h> // for memset

int eri_ssss_taylor(struct shell_pair const P,
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
        // lookup F7(xi) for xi nearest x
        // 7 = interpolation order, so use the define
        const double x = integralwork1[i];

        const int lookup_idx = (int)(BOYS_GRID_LOOKUPFAC*(x+BOYS_GRID_LOOKUPFAC2));
        const double xi = ((double)lookup_idx / BOYS_GRID_LOOKUPFAC);
        const double dx = xi-x;   // -delta x

        const double f0xi = boys_grid[lookup_idx][0];
        const double f1xi = boys_grid[lookup_idx][1];
        const double f2xi = boys_grid[lookup_idx][2];
        const double f3xi = boys_grid[lookup_idx][3];
        const double f4xi = boys_grid[lookup_idx][4];
        const double f5xi = boys_grid[lookup_idx][5];
        const double f6xi = boys_grid[lookup_idx][6];
        const double f7xi = boys_grid[lookup_idx][7];

        integralwork1[i] = f0xi
                         + dx * (                  f1xi
                         + dx * ( (1.0/2.0   )   * f2xi
                         + dx * ( (1.0/6.0   )   * f3xi
                         + dx * ( (1.0/24.0  )   * f4xi
                         + dx * ( (1.0/120.0 )   * f5xi
                         + dx * ( (1.0/720.0 )   * f6xi
                         + dx * ( (1.0/5040.0)   * f7xi 
                         )))))));

        // apply coefficients and prefactors
        integralwork1[i] *= integralwork2[i];
    }

    // now sum them, forming the contracted integrals
    memset(integrals, 0, nshell1234*sizeof(double));
    idx = 0;
    nint = 0;
    for(ab = 0; ab < P.nshell12; ++ab)
    for(cd = 0; cd < Q.nshell12; ++cd)
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

