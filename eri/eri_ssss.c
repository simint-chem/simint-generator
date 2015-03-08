#include <math.h>

#include "constants.h"
#include "shell.h"

int eri_ssss(struct shell_pair const * restrict P,
             struct shell_pair const * restrict Q,
             double * restrict integrals)
{
    int i, j, idx;

    for(i = 0, idx = 0; i < P->n; ++i)
    {
        for(j = 0; j < Q->n; ++j, ++idx)
        {
            double tmp;
   
            const double PQalpha_mul = P->alpha[i] * Q->alpha[j];
            const double PQalpha_sum = P->alpha[i] + Q->alpha[j];
 
            const double pfac = P->prefac[i] * Q->prefac[j]
                                / (PQalpha_mul * sqrt(PQalpha_sum));
    
            /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
            double R2 = 0.0;
            tmp = P->x[i] - Q->x[j]; R2 += tmp * tmp;
            tmp = P->y[i] - Q->y[j]; R2 += tmp * tmp;
            tmp = P->z[i] - Q->z[j]; R2 += tmp * tmp;
    
            //double F0;
            //Boys_F(&F0, 0, R2 * (P->alpha[i] * Q->alpha[j])/(P->alpha[i] + Q->alpha[j]));
            // F0 = K * erf(sqrt(x)) / sqrt(x)
            //      with x = R2 * Palpha * Qalpha / (Palpha + Qalpha)
            const double x2 = sqrt(R2 * PQalpha_mul/PQalpha_sum);

            integrals[idx] = pfac * F0_KFAC * erf(x2) / x2;
        }
    }

    return idx;
}
