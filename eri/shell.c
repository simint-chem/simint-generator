#include <math.h>

#include "vectorization.h"
#include "constants.h"
#include "eri/shell.h"
#include "eri/shell_constants.h"

extern double const norm_fac[SHELL_PRIM_NORMFAC_MAXL+1];

// Allocate a gaussian shell with correct alignment
void allocate_gaussian_shell(int nprim, struct gaussian_shell * const restrict G)
{
    const int prim_size = SIMD_ROUND_DBL(nprim);
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


void normalize_gaussian_shells(int n, struct gaussian_shell * const restrict G)
{
    /* 
     * Normalizes both primitives and contractions. This is Eq. 2.25 in
     * "Fundamentals of Molecular Integrals Evaluation"
     * by Justin T. Fermann & Edward F. Valeev
     */
    for(int i = 0; i < n; ++i)
    {
        const double am = (double)(G[i].am);
        const double m = am + 1.5;

        double sum = 0.0;

        for(int j = 0; j < G[i].nprim; j++)
        {
            const double a1 = G[i].alpha[j];

            for(int k = 0; k < G[i].nprim; k++)
            {
                const double a2 = G[i].alpha[k];
                const double s = (G[i].coef[j] * G[i].coef[k]) / pow(a1+a2, m);
                sum += s;
            }
        }

        const double norm = 1.0 / sqrt(sum * norm_fac[G[i].am]);

        // apply the normalization
        for (int j = 0; j < G[i].nprim; ++j)
            G[i].coef[j] *= norm; 
    }

}


// Allocates a shell pair with correct alignment
// Only fills in nprim_length member
void allocate_shell_pair(struct gaussian_shell const A,
                         struct gaussian_shell const B,
                         struct shell_pair * const restrict P)
{
    const int prim_size = SIMD_ROUND_DBL(A.nprim * B.nprim);
    P->nprim_length = prim_size;

    const int size = prim_size * sizeof(double);

    // allocate one large space
    double * mem = ALLOC(size * 8);
    P->x      = mem;
    P->y      = mem +   prim_size;
    P->z      = mem + 2*prim_size;
    P->PA_x   = mem + 3*prim_size;
    P->PA_y   = mem + 4*prim_size;
    P->PA_z   = mem + 5*prim_size;
    P->alpha  = mem + 6*prim_size;
    P->prefac = mem + 7*prim_size;
}


// Allocates a shell pair with correct alignment
// Only fills in nprim_length member
void allocate_multishell_pair(int na, struct gaussian_shell const * const restrict A,
                              int nb, struct gaussian_shell const * const restrict B,
                              struct multishell_pair * const restrict P)
{
    ASSUME_ALIGN(A);
    ASSUME_ALIGN(B);

    int prim_size = 0;

    // with rounding up to the nearest boundary
    for(int i = 0; i < na; ++i)
    for(int j = 0; j < nb; ++j)
        prim_size += SIMD_ROUND_DBL(A[i].nprim * B[j].nprim);

    const int shell12_size = na*nb;

    P->nprim_length = prim_size;

    const int size = prim_size * sizeof(double);
    const int size2 = shell12_size*sizeof(double); // for holding Xab, etc

    // allocate one large space
    double * mem = ALLOC(size * 11 + size2 * 3 ); 
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
    P->AB_y   = mem + 11*prim_size + shell12_size;
    P->AB_z   = mem + 11*prim_size + 2*shell12_size;

    /* Should this be aligned? I don't think so */
    int * intmem = malloc((na+nb+3*na*nb)*sizeof(int));
    P->nprim1 = intmem;
    P->nprim2 = intmem + na;
    P->nprim12 = intmem + na + nb;
    P->primstart = intmem + na + nb + na*nb;
    P->primend = intmem + na + nb + 2*na*nb;

}


void free_shell_pair(struct shell_pair P)
{
   // Only need to free P.x since that points to the beginning of mem
   FREE(P.x);
}


void free_multishell_pair(struct multishell_pair P)
{
   // Only need to free P.x since that points to the beginning of mem
   FREE(P.x);

   // similar with nprim1 and nprim2
   free(P.nprim1);
}



void fill_shell_pair(struct gaussian_shell const A,
                        struct gaussian_shell const B,
                        struct shell_pair * const restrict P)
{
    ASSUME_ALIGN(A.alpha);
    ASSUME_ALIGN(A.coef);

    ASSUME_ALIGN(B.alpha);
    ASSUME_ALIGN(B.coef);

    ASSUME_ALIGN(P->x);
    ASSUME_ALIGN(P->y);
    ASSUME_ALIGN(P->z);
    ASSUME_ALIGN(P->PA_x);
    ASSUME_ALIGN(P->PA_y);
    ASSUME_ALIGN(P->PA_z);
    ASSUME_ALIGN(P->alpha);
    ASSUME_ALIGN(P->prefac);

    int i, j;

    P->nprim = A.nprim * B.nprim;

    // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
    const double Xab_x = A.x - B.x;
    const double Xab_y = A.y - B.y;
    const double Xab_z = A.z - B.z;
    const double Xab = Xab_x*Xab_x + Xab_y*Xab_y + Xab_z*Xab_z;

    int idx = 0;
    for(i = 0; i < A.nprim; ++i)
    {
        const double AxAa = A.x * A.alpha[i];
        const double AyAa = A.y * A.alpha[i];
        const double AzAa = A.z * A.alpha[i];

        for(j = 0; j < B.nprim; ++j)
        {
            const double p_ab = A.alpha[i] + B.alpha[j];
            const double ABalpha = A.alpha[i] * B.alpha[j];

            P->prefac[idx] = A.coef[i] * B.coef[j]
                           * exp(-Xab * ABalpha / p_ab);

            P->x[idx] = (AxAa + B.alpha[j]*B.x)/p_ab;
            P->y[idx] = (AyAa + B.alpha[j]*B.y)/p_ab;
            P->z[idx] = (AzAa + B.alpha[j]*B.z)/p_ab;
            P->PA_x[idx] = P->x[idx] - A.x;
            P->PA_y[idx] = P->y[idx] - A.y;
            P->PA_z[idx] = P->z[idx] - A.z;
            P->alpha[idx] = p_ab;
            ++idx;
        }
    }

    P->AB_x = Xab_x;
    P->AB_y = Xab_y;
    P->AB_z = Xab_z;
}


void fill_multishell_pair(int na, struct gaussian_shell const * const restrict A,
                             int nb, struct gaussian_shell const * const restrict B,
                             struct multishell_pair * const restrict P)
{
    ASSUME_ALIGN(A);
    ASSUME_ALIGN(B);
    ASSUME_ALIGN(P->x);
    ASSUME_ALIGN(P->y);
    ASSUME_ALIGN(P->z);
    ASSUME_ALIGN(P->PA_x);
    ASSUME_ALIGN(P->PA_y);
    ASSUME_ALIGN(P->PA_z);
    ASSUME_ALIGN(P->bAB_x);
    ASSUME_ALIGN(P->bAB_Y);
    ASSUME_ALIGN(P->bAB_Z);
    ASSUME_ALIGN(P->alpha);
    ASSUME_ALIGN(P->prefac);

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
            // align to the next boundary
            idx = SIMD_ROUND_DBL(idx);

            P->primstart[sasb] = idx;

            P->am2 = B[sb].am;
            P->nprim2[sb] = B[sb].nprim;

            // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
            const double Xab_x = A[sa].x - B[sb].x;
            const double Xab_y = A[sa].y - B[sb].y;
            const double Xab_z = A[sa].z - B[sb].z;
            const double Xab = Xab_x*Xab_x + Xab_y*Xab_y + Xab_z*Xab_z;

            ASSUME(idx%SIMD_ALIGN_DBL == 0);

            ASSUME_ALIGN(A[sa].alpha);
            ASSUME_ALIGN(A[sa].coef);
            ASSUME_ALIGN(B[sb].alpha);
            ASSUME_ALIGN(B[sb].coef);

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
                                     * exp(-Xab * ABalpha / p_ab);

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

            // store the end of this shell pair
            P->primend[sasb] = idx;

            P->nprim12[sasb] = A[sa].nprim*B[sb].nprim;
            P->nprim += A[sa].nprim*B[sb].nprim;

            P->AB_x[sasb] = Xab_x;
            P->AB_y[sasb] = Xab_y;
            P->AB_z[sasb] = Xab_z;

            sasb++;
        }
    }
}

struct shell_pair
create_shell_pair(struct gaussian_shell const A,
                     struct gaussian_shell const B)
{
    struct shell_pair P;
    allocate_shell_pair(A, B, &P);
    fill_shell_pair(A, B, &P);
    return P; 
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
