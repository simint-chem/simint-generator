#include <math.h>
#include <string.h>

#include "vectorization/vectorization.h"
#include "constants.h"
#include "shell/shell.h"
#include "shell/shell_constants.h"

extern double const norm_fac[SHELL_PRIM_NORMFAC_MAXL+1];

// Allocate a gaussian shell with correct alignment
void allocate_gaussian_shell(int nprim, struct gaussian_shell * const restrict G)
{
    const int prim_size = SIMINT_SIMD_ROUND(nprim);
    const int size = prim_size * sizeof(double);   

    double * mem = ALLOC(2*size);
    G->alpha = mem;
    G->coef = mem + prim_size;
}


void free_gaussian_shell(struct gaussian_shell G)
{
    // only need to free G.alpha since only one memory space was used
    FREE(G.alpha);
}

struct gaussian_shell copy_gaussian_shell(const struct gaussian_shell G)
{
    struct gaussian_shell G_copy;
    allocate_gaussian_shell(G.nprim, &G_copy);

    G_copy.nprim = G.nprim;
    G_copy.am = G.am;
    G_copy.x = G.x;
    G_copy.y = G.y;
    G_copy.z = G.z;

    memcpy(G_copy.alpha, G.alpha, G.nprim * sizeof(double));
    memcpy(G_copy.coef, G.coef, G.nprim * sizeof(double));
    return G_copy;
}


void normalize_gaussian_shells(int n, struct gaussian_shell * const restrict G)
{
    for(int i = 0; i < n; ++i)
    {
        const int iam = G[i].am;
        const double am = (double)iam;
        const double m = am + 1.5;
        const double m2 = 0.5 * m;

        double sum = 0.0;

        for(int j = 0; j < G[i].nprim; j++)
        {
            const double a1 = G[i].alpha[j];
            const double c1 = G[i].coef[j];

            for(int k = 0; k < G[i].nprim; k++)
            {
                const double a2 = G[i].alpha[k];
                const double c2 = G[i].coef[k];
                sum += ( c1 * c2 *  pow(a1*a2, m2) ) / ( pow(a1+a2, m) );
            }
        }

        const double norm = 1.0 / sqrt(sum * norm_fac[iam]);

        // apply the rest of the normalization
        for (int j = 0; j < G[i].nprim; ++j)
            G[i].coef[j] *= norm * pow(G[i].alpha[j], m2);
    }

}




// Allocates a shell pair with correct alignment
// Only fills in {i,d}memsize and nprim_length members
void allocate_multishell_pair(int na, struct gaussian_shell const * const restrict A,
                              int nb, struct gaussian_shell const * const restrict B,
                              struct multishell_pair * const restrict P)
{
    ASSUME_ALIGN_DBL(A);
    ASSUME_ALIGN_DBL(B);

    int prim_size = 0;

    // with rounding up to the nearest boundary
    for(int i = 0; i < na; ++i)
    for(int j = 0; j < nb; ++j)
        prim_size += SIMINT_SIMD_ROUND(A[i].nprim * B[j].nprim);

    const int shell12_size = na*nb;


    const int dmemsize = (prim_size*11 + shell12_size*3)*sizeof(double);
    const int imemsize = (shell12_size)*sizeof(int);
    P->dmemsize = dmemsize;
    P->imemsize = imemsize;
    P->nprim_length = prim_size;

    // allocate one large space
    double * mem = ALLOC(dmemsize); 
    P->x      = mem;
    P->y      = mem +   prim_size;
    P->z      = mem + 2*prim_size;
    P->PA_x   = mem + 3*prim_size;
    P->PA_y   = mem + 4*prim_size;
    P->PA_z   = mem + 5*prim_size;
    P->bAB_x  = mem + 6*prim_size;
    P->bAB_y  = mem + 7*prim_size;
    P->bAB_z  = mem + 8*prim_size;
    P->alpha  = mem + 9*prim_size;
    P->prefac = mem + 10*prim_size;

    P->AB_x   = mem + 11*prim_size;
    P->AB_y   = mem + 11*prim_size +   shell12_size;
    P->AB_z   = mem + 11*prim_size + 2*shell12_size;

    /* Should this be aligned? I don't think so */
    int * intmem = malloc(imemsize);
    P->nprim12   = intmem;
}


void free_multishell_pair(struct multishell_pair P)
{
   // Only need to free P.x since that points to the beginning of mem
   FREE(P.x);

   // similar with nprim12
   free(P.nprim12);
}


void fill_multishell_pair(int na, struct gaussian_shell const * const restrict A,
                          int nb, struct gaussian_shell const * const restrict B,
                          struct multishell_pair * const restrict P)
{
    ASSUME_ALIGN_DBL(A);
    ASSUME_ALIGN_DBL(B);
    ASSUME_ALIGN_DBL(P->x);
    ASSUME_ALIGN_DBL(P->y);
    ASSUME_ALIGN_DBL(P->z);
    ASSUME_ALIGN_DBL(P->PA_x);
    ASSUME_ALIGN_DBL(P->PA_y);
    ASSUME_ALIGN_DBL(P->PA_z);
    ASSUME_ALIGN_DBL(P->bAB_x);
    ASSUME_ALIGN_DBL(P->bAB_y);
    ASSUME_ALIGN_DBL(P->bAB_z);
    ASSUME_ALIGN_DBL(P->alpha);
    ASSUME_ALIGN_DBL(P->prefac);

    int i, j, sa, sb, sasb, idx;

    P->nshell1 = na;
    P->nshell2 = nb;
    P->nshell12 = na * nb;
    P->nprim = 0;

    // zero out
    memset(P->x, 0, P->dmemsize);
    memset(P->nprim12, 0, P->imemsize);

    // set alpha to one for all
    for(i = 0; i < P->nprim_length; i++)
        P->alpha[i] = 1.0;

    sasb = 0;
    idx = 0;
    for(sa = 0; sa < na; ++sa)
    {
        P->am1 = A[sa].am;

        for(sb = 0; sb < nb; ++sb)
        {
            P->am2 = B[sb].am;

            // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
            const double Xab_x = A[sa].x - B[sb].x;
            const double Xab_y = A[sa].y - B[sb].y;
            const double Xab_z = A[sa].z - B[sb].z;
            const double Xab = Xab_x*Xab_x + Xab_y*Xab_y + Xab_z*Xab_z;

            ASSUME(idx%SIMINT_SIMD_LEN == 0);

            ASSUME_ALIGN_DBL(A[sa].alpha);
            ASSUME_ALIGN_DBL(A[sa].coef);
            ASSUME_ALIGN_DBL(B[sb].alpha);
            ASSUME_ALIGN_DBL(B[sb].coef);

            for(i = 0; i < A[sa].nprim; ++i)
            {
                const double AxAa = A[sa].x * A[sa].alpha[i];
                const double AyAa = A[sa].y * A[sa].alpha[i];
                const double AzAa = A[sa].z * A[sa].alpha[i];

                for(j = 0; j < B[sb].nprim; ++j)
                {
                    const double p_ab = A[sa].alpha[i] + B[sb].alpha[j];
                    const double ABalpha = A[sa].alpha[i] * B[sb].alpha[j];

                    P->prefac[idx] = A[sa].coef[i] * B[sb].coef[j]
                                     * exp(-Xab * ABalpha / p_ab)
                                     * SQRT_TWO_PI_52 / p_ab;

                    P->x[idx] = (AxAa + B[sb].alpha[j]*B[sb].x)/p_ab;
                    P->y[idx] = (AyAa + B[sb].alpha[j]*B[sb].y)/p_ab;
                    P->z[idx] = (AzAa + B[sb].alpha[j]*B[sb].z)/p_ab;
                    P->PA_x[idx] = P->x[idx] - A[sa].x;
                    P->PA_y[idx] = P->y[idx] - A[sa].y;
                    P->PA_z[idx] = P->z[idx] - A[sa].z;
                    P->bAB_x[idx] = B[sb].alpha[j]*Xab_x;
                    P->bAB_y[idx] = B[sb].alpha[j]*Xab_y;
                    P->bAB_z[idx] = B[sb].alpha[j]*Xab_z;
                    P->alpha[idx] = p_ab;
                    ++idx;
                }

            }

            // align to the next boundary
            idx = SIMINT_SIMD_ROUND(idx);

            P->nprim12[sasb] = A[sa].nprim*B[sb].nprim;
            P->nprim += A[sa].nprim*B[sb].nprim;

            P->AB_x[sasb] = Xab_x;
            P->AB_y[sasb] = Xab_y;
            P->AB_z[sasb] = Xab_z;

            sasb++;
        }
    }
}


struct multishell_pair
create_multishell_pair(int na, struct gaussian_shell const * const restrict A,
                       int nb, struct gaussian_shell const * const restrict B)
{
    struct multishell_pair P;
    allocate_multishell_pair(na, A, nb, B, &P);
    fill_multishell_pair(na, A, nb, B, &P);
    return P; 
}

