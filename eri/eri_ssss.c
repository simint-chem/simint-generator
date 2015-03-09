#include <math.h>

#include "constants.h"
#include "shell.h"


int eri__ssss(struct gaussian_shell const * restrict A,
              struct gaussian_shell const * restrict B,
              struct gaussian_shell const * restrict C,
              struct gaussian_shell const * restrict D,
              struct shell_pair * restrict P_tmp,
              struct shell_pair * restrict Q_tmp,
              double * restrict integrals);

int eri_0pair_ssss(int na, struct gaussian_shell const * const restrict A,
                   int nb, struct gaussian_shell const * const restrict B,
                   int nc, struct gaussian_shell const * const restrict C,
                   int nd, struct gaussian_shell const * const restrict D,
                   double * restrict integrals)
{
    int i, j, k, l;
    int sa, sb, sc, sd;
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

            for(sc = 0, idx = 0; sc < nc; ++sc)
            {
                for(sd = 0; sd < nd; ++sd)
                {
                    const double Xcd_x = C[sc].x - D[sb].x;
                    const double Xcd_y = C[sc].y - D[sb].y;
                    const double Xcd_z = C[sc].z - D[sb].z;
                    const double Xcd = Xcd_x*Xcd_x + Xcd_y*Xcd_y + Xcd_z*Xcd_z;

                    for(i = 0; i < A[sa].nprim; ++i)
                    {
                        const double AxAa = A[sa].x * A[sa].alpha[i];
                        const double AyAa = A[sa].y * A[sa].alpha[i];
                        const double AzAa = A[sa].z * A[sa].alpha[i];

                        for(j = 0; j < B[sb].nprim; ++j)
                        {
                            const double p_ab = A[sa].alpha[i] + B[sb].alpha[j];
                            const double alpha_ab = A[sa].alpha[i] * B[sb].alpha[j];

                            // throw the 16/sqrt(pi) in here and reduce
                            // the number of multiplies
                            // in the next loops
                            const double prefac_ab = A[sa].coef[i] * B[sb].coef[j]
                                                     * pow(alpha_ab, 0.75)
                                                     * exp(-Xab * alpha_ab / p_ab)
                                                     * ONESIX_OVER_SQRT_PI;

                            const double Px = (AxAa + B[sb].alpha[j]*B[sb].x)/p_ab;
                            const double Py = (AyAa + B[sb].alpha[j]*B[sb].y)/p_ab;
                            const double Pz = (AzAa + B[sb].alpha[j]*B[sb].z)/p_ab;


                            for(k = 0; k < C[sc].nprim; ++k)
                            {
                                const double CxCa = C[sc].x * C[sc].alpha[k];
                                const double CyCa = C[sc].y * C[sc].alpha[k];
                                const double CzCa = C[sc].z * C[sc].alpha[k];

                                for(l = 0; l < D[sd].nprim; ++l)
                                {
                                    const double p_cd = C[sc].alpha[k] + D[sd].alpha[l];
                                    const double alpha_cd = C[sc].alpha[k] * D[sd].alpha[l];

                                    const double prefac_cd = C[sc].coef[k] * D[sd].coef[l]
                                                             * pow(alpha_cd, 0.75)
                                                             * exp(-Xcd * alpha_cd / p_cd);

                                    const double Qx = (CxCa + D[sd].alpha[l]*D[sd].x)/p_cd;
                                    const double Qy = (CyCa + D[sd].alpha[l]*D[sd].y)/p_cd;
                                    const double Qz = (CzCa + D[sd].alpha[l]*D[sd].z)/p_cd;


                                    const double PQalpha_mul = p_ab * p_cd;
                                    const double PQalpha_sum = p_ab + p_cd;

                                    const double pfac = prefac_ab * prefac_cd
                                                        / (PQalpha_mul * sqrt(PQalpha_sum));

                                    /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
                                    const double PQ_x = Px - Qx;
                                    const double PQ_y = Py - Qy;
                                    const double PQ_z = Pz - Qz;
                                    const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;

                                    //double F0;
                                    //Boys_F(&F0, 0, R2 * (P->alpha[i] * Q.alpha[j])/(P->alpha[i] + Q.alpha[j]));
                                    // F0 = K * erf(sqrt(x)) / sqrt(x)
                                    //      with x = R2 * Palpha * Qalpha / (Palpha + Qalpha)
                                    const double x2 = sqrt(R2 * PQalpha_mul/PQalpha_sum);

                                    integrals[idx++] = pfac * F0_KFAC * erf(x2) / x2;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return idx;
}



int eri_1pair_ssss(int na, struct gaussian_shell const * const restrict A,
                   int nb, struct gaussian_shell const * const restrict B,
                   struct shell_pair const Q,
                   double * restrict integrals)
{
    int i, j, k, sa, sb;
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

                    const double prefac_ab = A[sa].coef[i] * B[sb].coef[j]
                                             * pow(alpha_ab, 0.75)
                                             * exp(-Xab * alpha_ab / p_ab)
                                             * ONESIX_OVER_SQRT_PI;  // throw this here and reduce
                    // the number of multiplies
                    // in the next loop

                    const double Px = (AxAa + B[sb].alpha[j]*B[sb].x)/p_ab;
                    const double Py = (AyAa + B[sb].alpha[j]*B[sb].y)/p_ab;
                    const double Pz = (AzAa + B[sb].alpha[j]*B[sb].z)/p_ab;


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

                        //double F0;
                        //Boys_F(&F0, 0, R2 * (P->alpha[i] * Q.alpha[j])/(P->alpha[i] + Q.alpha[j]));
                        // F0 = K * erf(sqrt(x)) / sqrt(x)
                        //      with x = R2 * Palpha * Qalpha / (Palpha + Qalpha)
                        const double x2 = sqrt(R2 * PQalpha_mul/PQalpha_sum);

                        integrals[idx++] = pfac * F0_KFAC * erf(x2) / x2;
                    }
                }
            }
        }
    }

    return idx;
}

int eri_2pair_ssss(struct shell_pair const P,
                   struct shell_pair const Q,
                   double * restrict integrals)
{
    int i, j;
    int idx = 0;

    for(i = 0; i < P.nprim; ++i)
    {
        const double pfac_p = ONESIX_OVER_SQRT_PI * P.prefac[i];

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

            //double F0;
            //Boys_F(&F0, 0, R2 * (P.alpha[i] * Q.alpha[j])/(P.alpha[i] + Q.alpha[j]));
            // F0 = K * erf(sqrt(x)) / sqrt(x)
            //      with x = R2 * Palpha * Qalpha / (Palpha + Qalpha)
            const double x2 = sqrt(R2 * PQalpha_mul/PQalpha_sum);

            integrals[idx++] = pfac * F0_KFAC * erf(x2) / x2;
        }
    }

    return idx;
}








