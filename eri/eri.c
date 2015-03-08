#include <math.h>

#include "eri/eri.h"

int create_shell_pair(const struct gaussian_shell * A,
                      const struct gaussian_shell * B,
                      struct shell_pair * P)
{
    int i, j, idx;
    double Xab_tmp;

    double AxAa, AyAa, AzAa;

    // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
    double Xab = 0;
    Xab_tmp = A->x - B->x;  Xab += Xab_tmp * Xab_tmp;
    Xab_tmp = A->y - B->y;  Xab += Xab_tmp * Xab_tmp;
    Xab_tmp = A->z - B->z;  Xab += Xab_tmp * Xab_tmp;

    for(i = 0, idx = 0; i < A->nprim; ++i)
    {
        AxAa = A->x * A->alpha[i];
        AyAa = A->y * A->alpha[i];
        AzAa = A->z * A->alpha[i];

        for(j = 0; j < B->nprim; ++j, ++idx)
        {
            double p_ab;
            double pfac;
    
            p_ab = A->alpha[i] + B->alpha[j];
    
            // This accumulates the normalization and prefactors, coefficients, etc
            // pi^1.25 is part of the normalization
            pfac = A->coef[i]*B->coef[j] * pow(M_PI, 1.25);
            pfac *= pow(2.0 * A->alpha[i] / M_PI, 0.75);
            pfac *= pow(2.0 * B->alpha[j] / M_PI, 0.75);
    
            pfac *= exp(-((A->alpha[i] * B->alpha[j]) / (p_ab)) * Xab);
    
            P->x[idx] = (AxAa + B->alpha[j]*B->x)/p_ab;
            P->y[idx] = (AyAa + B->alpha[j]*B->y)/p_ab;
            P->z[idx] = (AzAa + B->alpha[j]*B->z)/p_ab;
            P->alpha[idx] = p_ab;
            P->prefac[idx] = pfac;
        }
    }
    return idx; // should be na*nb 
}
