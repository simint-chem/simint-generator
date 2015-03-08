#include <math.h>

#include "constants.h"
#include "shell.h"

int create_ss_shell_pair(struct gaussian_shell const * restrict A,
                         struct gaussian_shell const * restrict B,
                         struct shell_pair * restrict P)
{
    int i, j, idx;
    double Xab_tmp;

    // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
    double Xab = 0;
    Xab_tmp = A->x - B->x;  Xab += Xab_tmp * Xab_tmp;
    Xab_tmp = A->y - B->y;  Xab += Xab_tmp * Xab_tmp;
    Xab_tmp = A->z - B->z;  Xab += Xab_tmp * Xab_tmp;

    for(i = 0, idx = 0; i < A->nprim; ++i)
    {
        const double AxAa = A->x * A->alpha[i];
        const double AyAa = A->y * A->alpha[i];
        const double AzAa = A->z * A->alpha[i];

        for(j = 0; j < B->nprim; ++j, ++idx)
        {
            const double p_ab = A->alpha[i] + B->alpha[j];
            const double ABalpha = A->alpha[i] * B->alpha[j];
    
            P->prefac[idx] = A->coef[i] * B->coef[j]
                           * pow(ABalpha, 0.75)
                           * exp(-Xab * ABalpha / p_ab);

            P->x[idx] = (AxAa + B->alpha[j]*B->x)/p_ab;
            P->y[idx] = (AyAa + B->alpha[j]*B->y)/p_ab;
            P->z[idx] = (AzAa + B->alpha[j]*B->z)/p_ab;
            P->alpha[idx] = p_ab;
        }
    }

    return idx; // should be na*nb 
}
