#include <math.h>
#include <string.h> // for memset

#include "constants.h"
#include "boys/boys.h"
#include "eri/shell.h"


int eri_ssss_split(struct shell_pair const P,
                   struct shell_pair const Q,
                   double * const restrict integrals,
                   double * const integralwork1,
                   double * const integralwork2)
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

    // we need to split integralwork1 and integralwork2 into three separate spaces
    // One holds the whole thing, one holds the short regime, and one holds the long regime
    const int nprim1234 = P.nprim * Q.nprim;

    double * const restrict integralwork1_init = integralwork1;
    double * const restrict integralwork1_short = integralwork1 + SIMD_ROUND_DBL(nprim1234);
    double * const restrict integralwork1_long = integralwork1 + 2*SIMD_ROUND_DBL(nprim1234);

    double * const restrict integralwork2_init = integralwork2;
    double * const restrict integralwork2_short = integralwork2 + SIMD_ROUND_DBL(nprim1234);
    double * const restrict integralwork2_long = integralwork2 + 2*SIMD_ROUND_DBL(nprim1234);


    ASSUME_ALIGN(integralwork1_short);
    ASSUME_ALIGN(integralwork1_long);
    ASSUME_ALIGN(integralwork2_short);
    ASSUME_ALIGN(integralwork2_long);


    int i, j;
    int ab, cd;
    int idx = 0;
    int idx_short = 0;
    int idx_long = 0;
    const int nshell1234 = P.nshell12 * Q.nshell12;
    int blockends_short[nshell1234];
    int blockends_long[nshell1234];

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
                    integralwork1_init[idx] = R2 * PQalpha_mul/PQalpha_sum;
                    integralwork2_init[idx] = pfac * P.prefac[i] * Q.prefac[j];
                    ++idx;
                }
            }
        }
    }

    // partition into the long and short regimes 
    int iblock = 0;
    idx = 0;
    for(ab = 0; ab < P.nshell12; ++ab)
    {
        for(cd = 0; cd < Q.nshell12; ++cd)
        {
            const int abcd12 = P.nprim12[ab] * Q.nprim12[cd];

            for(i = 0; i < abcd12; ++i)
            {
                const double Fparam = integralwork1_init[idx];
                const double prefac = integralwork2_init[idx];
                if(Fparam < BOYS_SHORTGRID_MAXX)
                {
                    integralwork1_short[idx_short] = Fparam;
                    integralwork2_short[idx_short] = prefac;
                    ++idx_short;
                }
                else
                {
                    integralwork1_long[idx_long] = Fparam;
                    integralwork2_long[idx_long] = prefac;
                    ++idx_long;
                }
                ++idx;
            }

            blockends_short[iblock] = idx_short;
            blockends_long[iblock] = idx_long;
            ++iblock;
        }
    }

    // rip through the integral work arrays and store result back in integralwork1
    // large x
    for(i = 0; i < idx_long; ++i)
    {
        // will multiply by F0_KFAC later
        // but apply coefficients and prefactors here
        integralwork1_long[i] = integralwork2_long[i] / sqrt(integralwork1_long[i]); 
    }

    // short x
    for(i = 0; i < idx_short; ++i)
    {
        integralwork1_short[i] = integralwork2_short[i] * Boys_F0_taylor(integralwork1_short[i]);
    }

    // now sum them, forming the contracted integrals
    memset(integrals, 0, nshell1234*sizeof(double));

    // long first
    j = 0;
    for(i = 0; i < nshell1234; ++i)
    {
        const int iend = blockends_long[i];
        for(; j < iend; ++j)
            integrals[i] += integralwork1_long[j];
    }

    // apply factor, which isn't included above
    for(i = 0; i < nshell1234; ++i)
        integrals[i] *= F0_KFAC;

    // now short
    j = 0;
    for(i = 0; i < nshell1234; ++i)
    {
        const int iend = blockends_short[i];
        for(; j < iend; ++j)
        {
            integrals[i] += integralwork1_short[j];
        }
    }

    // apply constants to integrals
    // also heavily vectorized
    for(i = 0; i < nshell1234; ++i)
        integrals[i] *= ONESIX_OVER_SQRT_PI;

    return nshell1234;

}
