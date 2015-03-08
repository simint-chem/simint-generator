#include <math.h>
#include <stdio.h>

#include "eri/eri.h"
#include "boys/boys.h"
#include "valeev/valeev.h"


int eri_ssss(const struct shell_pair * P,
             const struct shell_pair * Q,
             double * integrals)
{
    int i, j, idx;

    for(i = 0, idx = 0; i < P->n; ++i)
    {
        for(j = 0; j < Q->n; ++j, ++idx)
        {
            double tmp;
            double R2;
            double F0;
            double pfac;
            double Rpq_x, Rpq_y, Rpq_z;
    
            pfac = 2 * P->prefac[i] * Q->prefac[j]
                  / (P->alpha[i] * Q->alpha[j] * sqrt(P->alpha[i] + Q->alpha[j]));
    
            Rpq_x = P->x[i] - Q->x[j];
            Rpq_y = P->y[i] - Q->y[j];
            Rpq_z = P->z[i] - Q->z[j];
            R2 = Rpq_x*Rpq_x + Rpq_y*Rpq_y + Rpq_z*Rpq_z;
    
            /* Same as above, but no more extra variables */
            R2 = 0;
            tmp = P->x[i] - Q->x[j]; R2 += tmp * tmp;
            tmp = P->y[i] - Q->y[j]; R2 += tmp * tmp;
            tmp = P->z[i] - Q->z[j]; R2 += tmp * tmp;
    
            Boys_F(&F0, 0, R2 * (P->alpha[i] * Q->alpha[j])/(P->alpha[i] + Q->alpha[j]));
            integrals[idx] = pfac * F0;
        }
    }

    return idx;
}
