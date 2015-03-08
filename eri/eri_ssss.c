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


int eri_0pair_ssss(struct gaussian_shell const A,
                   struct gaussian_shell const B,
                   struct gaussian_shell const C,
                   struct gaussian_shell const D,
                   double * restrict integrals);


int eri_1pair_ssss(struct gaussian_shell const A,
                   struct gaussian_shell const B,
                   struct shell_pair const Q,
                   double * restrict integrals)
{
    int i, j, k, idx;

    // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
    double Xab_tmp;
    double Xab = 0;
    Xab_tmp = A.x - B.x;  Xab += Xab_tmp * Xab_tmp;
    Xab_tmp = A.y - B.y;  Xab += Xab_tmp * Xab_tmp;
    Xab_tmp = A.z - B.z;  Xab += Xab_tmp * Xab_tmp;

    for(i = 0, idx = 0; i < A.nprim; ++i)
    {
        const double AxAa = A.x * A.alpha[i];
        const double AyAa = A.y * A.alpha[i];
        const double AzAa = A.z * A.alpha[i];

        for(j = 0; j < B.nprim; ++j)
        {
            const double p_ab = A.alpha[i] + B.alpha[j];
            const double ABalpha = A.alpha[i] * B.alpha[j];
    
            const double prefac_ab = A.coef[i] * B.coef[j]
                                   * pow(ABalpha, 0.75)
                                   * exp(-Xab * ABalpha / p_ab)
                                   * ONESIX_OVER_SQRT_PI;  // throw this here and reduce
                                                           // the number of multiplies
                                                           // in the next loop

            const double Px = (AxAa + B.alpha[j]*B.x)/p_ab;
            const double Py = (AyAa + B.alpha[j]*B.y)/p_ab;
            const double Pz = (AzAa + B.alpha[j]*B.z)/p_ab;


            for(k = 0; k < Q.n; ++k, ++idx)
            { 
                double tmp;
       
                const double PQalpha_mul = p_ab * Q.alpha[k];
                const double PQalpha_sum = p_ab + Q.alpha[k];
     
                const double pfac = prefac_ab * Q.prefac[k]
                                    / (PQalpha_mul * sqrt(PQalpha_sum));
        
                /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
                double R2 = 0.0;
                tmp = Px - Q.x[k]; R2 += tmp * tmp;
                tmp = Py - Q.y[k]; R2 += tmp * tmp;
                tmp = Pz - Q.z[k]; R2 += tmp * tmp;
        
                //double F0;
                //Boys_F(&F0, 0, R2 * (P->alpha[i] * Q.alpha[j])/(P->alpha[i] + Q.alpha[j]));
                // F0 = K * erf(sqrt(x)) / sqrt(x)
                //      with x = R2 * Palpha * Qalpha / (Palpha + Qalpha)
                const double x2 = sqrt(R2 * PQalpha_mul/PQalpha_sum);
    
                integrals[idx] = pfac * F0_KFAC * erf(x2) / x2;
            }
        }
    }

    return idx;
}

int eri_2pair_ssss(struct shell_pair const P,
                   struct shell_pair const Q,
                   double * restrict integrals)
{
    int i, j, idx;

    for(i = 0, idx = 0; i < P.n; ++i)
    {
        const double pfac_p = ONESIX_OVER_SQRT_PI * P.prefac[i];

        for(j = 0; j < Q.n; ++j, ++idx)
        {
            double tmp;
   
            const double PQalpha_mul = P.alpha[i] * Q.alpha[j];
            const double PQalpha_sum = P.alpha[i] + Q.alpha[j];
 
            const double pfac = pfac_p * Q.prefac[j]
                                / (PQalpha_mul * sqrt(PQalpha_sum));
    
            /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
            double R2 = 0.0;
            tmp = P.x[i] - Q.x[j]; R2 += tmp * tmp;
            tmp = P.y[i] - Q.y[j]; R2 += tmp * tmp;
            tmp = P.z[i] - Q.z[j]; R2 += tmp * tmp;
    
            //double F0;
            //Boys_F(&F0, 0, R2 * (P.alpha[i] * Q.alpha[j])/(P.alpha[i] + Q.alpha[j]));
            // F0 = K * erf(sqrt(x)) / sqrt(x)
            //      with x = R2 * Palpha * Qalpha / (Palpha + Qalpha)
            const double x2 = sqrt(R2 * PQalpha_mul/PQalpha_sum);

            integrals[idx] = pfac * F0_KFAC * erf(x2) / x2;
        }
    }

    return idx;
}
