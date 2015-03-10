#include <math.h>

#include "constants.h"
#include "vectorization.h"
#include "eri/shell.h"

// Allocates a shell pair with correct alignment
struct shell_pair
allocate_shell_pair(int na, struct gaussian_shell const * const restrict A,
                    int nb, struct gaussian_shell const * const restrict B)
{
    struct shell_pair P;
    P.nprim = 0;
    for(int i = 0; i < na; ++i)
    for(int j = 0; j < nb; ++j)
        P.nprim += A[i].nprim * B[j].nprim;

    // round up to the nearest boundary
    int size = P.nprim * sizeof(double);
    int simd_size = SIMD_ROUND(size);
    int simd_per = simd_size / sizeof(double);

    // allocate one large space
    double * mem = ALLOC(simd_size * 5);   // 5 = x, y, z, alpha, prefac
    P.x = mem;
    P.y = mem + simd_per;
    P.z = mem + 2*simd_per;
    P.alpha = mem + 3*simd_per;
    P.prefac = mem + 4*simd_per;

    /* Should this be aligned? I don't think so */
    int * intmem = malloc((na+nb+na*nb)*sizeof(int));
    P.nprim1 = intmem;
    P.nprim2 = intmem + na;
    P.nprim12 = intmem + na + nb;
    
    return P; 
}


void free_shell_pair(struct shell_pair P)
{
   // Only need to free P.x since that points to the beginning of mem
   FREE(P.x);

   // similar with nprim1 and nprim2
   free(P.nprim1);
}




void fill_ss_shell_pair(int na, struct gaussian_shell const * const restrict A,
                        int nb, struct gaussian_shell const * const restrict B,
                        struct shell_pair * const restrict P)
{
    int i, j, sa, sb, sasb, idx;

    P->nshell1 = na;
    P->nshell2 = nb;
    P->nshell12 = na * nb;
    P->nprim = 0;

    sasb = 0;

    idx = 0;
    for(sa = 0; sa < na; ++sa)
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

            for(i = 0; i < A[sa].nprim; ++i)
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

            P->nprim12[sasb] = A[sa].nprim*B[sb].nprim;
            P->nprim += A[sa].nprim*B[sb].nprim;
            sasb++;
        }
    }

}


struct shell_pair
create_ss_shell_pair(int na, struct gaussian_shell const * const restrict A,
                     int nb, struct gaussian_shell const * const restrict B)
{
    struct shell_pair P = allocate_shell_pair(na, A, nb, B);
    fill_ss_shell_pair(na, A, nb, B, &P);
    return P; 
}

