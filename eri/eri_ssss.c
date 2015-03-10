#include <math.h>

#include "constants.h"
#include "eri/shell.h"
#include "vectorization.h"

#include "string.h" // for memset

int eri_1pair_ssss_single(struct gaussian_shell const A,
                          struct gaussian_shell const B,
                          struct shell_pair const Q,
                          double * const restrict integrals)
{
    ASSUME_ALIGN(Q.x);
    ASSUME_ALIGN(Q.y);
    ASSUME_ALIGN(Q.z);
    ASSUME_ALIGN(Q.alpha);
    ASSUME_ALIGN(Q.prefac);

    int i, j, k;
    double sum = 0.0;

    // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
    const double Xab_x = A.x - B.x;
    const double Xab_y = A.y - B.y;
    const double Xab_z = A.z - B.z;
    const double Xab = Xab_x*Xab_x + Xab_y*Xab_y + Xab_z*Xab_z;


    for(i = 0; i < A.nprim; ++i)
    {
        const double AxAa = A.x * A.alpha[i];
        const double AyAa = A.y * A.alpha[i];
        const double AzAa = A.z * A.alpha[i];

        for(j = 0; j < B.nprim; ++j)
        {
            const double p_ab = A.alpha[i] + B.alpha[j];
            const double alpha_ab = A.alpha[i] * B.alpha[j];

            const double prefac_ab = A.coef[i] * B.coef[j]
                                     * pow(alpha_ab, 0.75)
                                     * exp(-Xab * alpha_ab / p_ab);

            const double Px = (AxAa + B.alpha[j]*B.x)/p_ab;
            const double Py = (AyAa + B.alpha[j]*B.y)/p_ab;
            const double Pz = (AzAa + B.alpha[j]*B.z)/p_ab;


            for(k = 0; k < Q.nprim; ++k)
            {
                const double PQalpha_mul = p_ab * Q.alpha[k];
                const double PQalpha_sum = p_ab + Q.alpha[k];

                const double pfac = prefac_ab * Q.prefac[k]
                                    / (PQalpha_mul * sqrt(PQalpha_sum));

                /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
                const double PQ_x = Px - Q.x[k];
                const double PQ_y = Py - Q.y[k];
                const double PQ_z = Pz - Q.z[k];
                const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;

                const double x2 = sqrt(R2 * PQalpha_mul/PQalpha_sum);

                sum += pfac * F0_KFAC * erf(x2) / x2;
            }
        }
    }

    integrals[0] = sum * ONESIX_OVER_SQRT_PI;

    return 1;
}


int eri_1pair_ssss_multi(int na, struct gaussian_shell const * const restrict A,
                         int nb, struct gaussian_shell const * const restrict B,
                         struct shell_pair const Q,
                         double * const restrict integrals,
                         double * const restrict integralwork)
{
    ASSUME_ALIGN(Q.x);
    ASSUME_ALIGN(Q.y);
    ASSUME_ALIGN(Q.z);
    ASSUME_ALIGN(Q.alpha);
    ASSUME_ALIGN(Q.prefac);

    ASSUME_ALIGN(integralwork);

    int i, j, k, l, sa, sb;
    int idx = 0;

    for(sa = 0; sa < na; ++sa)
    {
        for(sb = 0; sb < nb; ++sb)
        {
            // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
            const double Xab_x = A[sa].x - B[sb].x;
            const double Xab_y = A[sa].y - B[sb].y;
            const double Xab_z = A[sa].z - B[sb].z;
            const double Xab = Xab_x*Xab_x + Xab_y*Xab_y + Xab_z*Xab_z;

            for(i = 0; i < A[sa].nprim; ++i)
            {
                const double AxAa = A[sa].x * A[sa].alpha[i];
                const double AyAa = A[sa].y * A[sa].alpha[i];
                const double AzAa = A[sa].z * A[sa].alpha[i];

                for(j = 0; j < B[sb].nprim; ++j)
                {
                    const double p_ab = A[sa].alpha[i] + B[sb].alpha[j];
                    const double alpha_ab = A[sa].alpha[i] * B[sb].alpha[j];

                    const double prefac_ab = pow(alpha_ab, 0.75)
                                             * exp(-Xab * alpha_ab / p_ab);

                    const double Px = (AxAa + B[sb].alpha[j]*B[sb].x)/p_ab;
                    const double Py = (AyAa + B[sb].alpha[j]*B[sb].y)/p_ab;
                    const double Pz = (AzAa + B[sb].alpha[j]*B[sb].z)/p_ab;

                    ASSUME(idx%SIMD_ALIGN_DBL == 0);
                    for(k = 0; k < Q.nprim; ++k)
                    {
                        const double PQalpha_mul = p_ab * Q.alpha[k];
                        const double PQalpha_sum = p_ab + Q.alpha[k];

                        const double pfac = prefac_ab / (PQalpha_mul * sqrt(PQalpha_sum));

                        /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
                        const double PQ_x = Px - Q.x[k];
                        const double PQ_y = Py - Q.y[k];
                        const double PQ_z = Pz - Q.z[k];
                        const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;

                        const double x2 = sqrt(R2 * PQalpha_mul/PQalpha_sum);
                        integralwork[idx++] = pfac * F0_KFAC * erf(x2) / x2;
                    }
                    idx = SIMD_ROUND_DBL(idx);   

                }
            }
        }
    }

    // apply contraction coefficients and prefactor
    idx = 0;
    for(sa = 0; sa < na; ++sa)
    for(sb = 0; sb < nb; ++sb)
    for(i = 0; i < A[sa].nprim; ++i)
    for(j = 0; j < B[sb].nprim; ++j)
    {
        const double prefac_ab = ONESIX_OVER_SQRT_PI * A[sa].coef[i] * B[sb].coef[j];

        ASSUME(idx%SIMD_ALIGN_DBL == 0);
        for(k = 0; k < Q.nprim; ++k)
            integralwork[idx++] *= prefac_ab * Q.prefac[k];
        idx = SIMD_ROUND_DBL(idx);   
    }


    // now sum them
    idx = 0;
    int nintstart = 0;
    memset(integrals, 0, na*nb*Q.nshell12);
    for(sa = 0; sa < na; ++sa)
    for(sb = 0; sb < nb; ++sb)
    {
        const int sasb = A[sa].nprim * B[sb].nprim;
        for(i = 0; i < sasb; ++i)
        {
            int nint = nintstart;

            // this line isn't really needed, but why not
            ASSUME(idx%SIMD_ALIGN_DBL == 0);
            for(j = 0; j < Q.nshell12; ++j)
            {
                for(l = 0; l < Q.nprim12[j]; ++l)
                    integrals[nint] += integralwork[idx++];

                nint++;
            }
            idx = SIMD_ROUND_DBL(idx);

        }
        nintstart += Q.nshell12;
    }

    return na*nb*Q.nshell12;
}



// assuming both P and Q contain only pairs between one shell
int eri_2pair_ssss_single(struct shell_pair const P,
                          struct shell_pair const Q,
                          double * const restrict integrals)
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

    int i, j;
    double sum = 0.0;

    for(i = 0; i < P.nprim; ++i)
    {
        const double pfac_p = P.prefac[i];

        for(j = 0; j < Q.nprim; ++j)
        {
            const double PQalpha_mul = P.alpha[i] * Q.alpha[j];
            const double PQalpha_sum = P.alpha[i] + Q.alpha[j];

            const double pfac = pfac_p * Q.prefac[j]
                                / (PQalpha_mul * sqrt(PQalpha_sum));

            /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
            const double PQ_x = P.x[i] - Q.x[j];
            const double PQ_y = P.y[i] - Q.y[j];
            const double PQ_z = P.z[i] - Q.z[j];
            const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;

            const double x2 = sqrt(R2 * PQalpha_mul/PQalpha_sum);

            sum += pfac * F0_KFAC * erf(x2) / x2;
        }
    }

    integrals[0] = sum * ONESIX_OVER_SQRT_PI;

    return 1;
}


int eri_2pair_ssss_multi(struct shell_pair const P,
                         struct shell_pair const Q,
                         double * const restrict integrals,
                         double * const restrict integralwork)
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

    ASSUME_ALIGN(integralwork);

    int i, j, k, l;
    int idx = 0;

    for(i = 0; i < P.nprim; ++i)
    {
        ASSUME(idx%SIMD_ALIGN_DBL == 0);
        for(j = 0; j < Q.nprim; ++j)
        {
            const double PQalpha_mul = P.alpha[i] * Q.alpha[j];
            const double PQalpha_sum = P.alpha[i] + Q.alpha[j];

            const double pfac = 1.0 / (PQalpha_mul * sqrt(PQalpha_sum));

            /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
            const double PQ_x = P.x[i] - Q.x[j];
            const double PQ_y = P.y[i] - Q.y[j];
            const double PQ_z = P.z[i] - Q.z[j];
            const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;

            const double x2 = sqrt(R2 * PQalpha_mul/PQalpha_sum);

            integralwork[idx] = pfac * F0_KFAC * erf(x2) / x2;
            idx++;
        }
        idx = SIMD_ROUND_DBL(idx);
    }

    // apply contraction coefficients and prefactor
    idx = 0;
    for(i = 0; i < P.nprim; ++i)
    {
        ASSUME(idx%SIMD_ALIGN_DBL == 0);
        for(j = 0; j < Q.nprim; ++j)
            integralwork[idx++] *= ONESIX_OVER_SQRT_PI * P.prefac[i] * Q.prefac[j];
        idx = SIMD_ROUND_DBL(idx);
    }


    // now sum them
    idx = 0;
    int nintstart = 0;
    memset(integrals, 0, P.nshell12*Q.nshell12);
    for(l = 0; l < P.nshell12; l++)
    {
        for(k = 0; k < P.nprim12[l]; ++k)
        {
            int nint = nintstart;

            // this line isn't really needed, but why not
            ASSUME(idx%SIMD_ALIGN_DBL == 0);
            for(i = 0; i < Q.nshell12; ++i)
            {
                for(j = 0; j < Q.nprim12[i]; ++j)
                    integrals[nint] += integralwork[idx++];
                nint++;
            }
            idx = SIMD_ROUND_DBL(idx);

        }

        nintstart += Q.nshell12;
    }

    return idx;
}
