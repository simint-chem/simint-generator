#include <math.h>

#include "constants.h"
#include "shell.h"

void create_ss_shell_pair(struct gaussian_shell const A,
                          struct gaussian_shell const B,
                          struct shell_pair * restrict P)
{
    int i, j, idx;

    // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
    const double Xab_x = A.x - B.x;
    const double Xab_y = A.y - B.y;
    const double Xab_z = A.z - B.z;
    const double Xab = Xab_x*Xab_x + Xab_y*Xab_y + Xab_z*Xab_z;

    P->am1 = A.am;
    P->am2 = B.am;
    P->nshell1 = P->nshell2 = 1;
    P->nprim1[0] = A.nprim;
    P->nprim2[0] = B.nprim;

    for(i = 0, idx = 0; i < A.nprim; ++i)
    {
        const double AxAa = A.x * A.alpha[i];
        const double AyAa = A.y * A.alpha[i];
        const double AzAa = A.z * A.alpha[i];

        for(j = 0; j < B.nprim; ++j, ++idx)
        {
            const double p_ab = A.alpha[i] + B.alpha[j];
            const double ABalpha = A.alpha[i] * B.alpha[j];

            P->prefac[idx] = A.coef[i] * B.coef[j]
                             * pow(ABalpha, 0.75)
                             * exp(-Xab * ABalpha / p_ab);

            P->x[idx] = (AxAa + B.alpha[j]*B.x)/p_ab;
            P->y[idx] = (AyAa + B.alpha[j]*B.y)/p_ab;
            P->z[idx] = (AzAa + B.alpha[j]*B.z)/p_ab;
            P->alpha[idx] = p_ab;
        }
    }

    P->nprim = idx;
}


void create_ss_shell_pair_multi(int na, struct gaussian_shell const * const A,
                                int nb, struct gaussian_shell const * const B,
                                struct shell_pair * restrict P)
{
    int i, j, sa, sb, idx;

    P->nshell1 = na;
    P->nshell2 = nb;
    P->nprim = 0;

    for(sa = 0, idx = 0; sa < na; ++sa)
    {
        P->am1 = A[sa].am;
        P->nprim1[sa] = A[sa].nprim;

        for(sb = 0; sb < nb; ++sb)
        {
            P->am2 = B[sb].am;
            P->nprim2[sb] = B[sb].nprim;

            // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
            const double Xab_x = A[sa].x - B[sb].x;
            const double Xab_y = A[sa].y - B[sb].y;
            const double Xab_z = A[sa].z - B[sb].z;
            const double Xab = Xab_x*Xab_x + Xab_y*Xab_y + Xab_z*Xab_z;

            for(i = 0, idx = 0; i < A[sa].nprim; ++i)
            {
                const double AxAa = A[sa].x * A[sa].alpha[i];
                const double AyAa = A[sa].y * A[sa].alpha[i];
                const double AzAa = A[sa].z * A[sa].alpha[i];

                for(j = 0; j < B[sb].nprim; ++j, ++idx)
                {
                    const double p_ab = A[sa].alpha[i] + B[sb].alpha[j];
                    const double ABalpha = A[sa].alpha[i] * B[sb].alpha[j];

                    P->prefac[idx] = A[sa].coef[i] * B[sb].coef[j]
                                     * pow(ABalpha, 0.75)
                                     * exp(-Xab * ABalpha / p_ab);

                    P->x[idx] = (AxAa + B[sb].alpha[j]*B[sb].x)/p_ab;
                    P->y[idx] = (AyAa + B[sb].alpha[j]*B[sb].y)/p_ab;
                    P->z[idx] = (AzAa + B[sb].alpha[j]*B[sb].z)/p_ab;
                    P->alpha[idx] = p_ab;
                }
            }
        }
    }

    P->nprim = idx;
}

