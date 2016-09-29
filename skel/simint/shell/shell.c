#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "simint/vectorization/vectorization.h"
#include "simint/constants.h"
#include "simint/shell/shell.h"
#include "simint/shell/shell_screen.h"
#include "simint/shell/shell_constants.h"


#if defined(__ICC) || defined(__INTEL_COMPILER)
    #pragma warning(disable:1338)                  // Pointer arithmetic on void
#elif defined(__GNUC__) || defined(__GNUG__)
    #pragma GCC diagnostic ignored -Wpointer-arith // Pointer arithmetic on void
#endif


#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

extern double const norm_fac[SHELL_PRIM_NORMFAC_MAXL+1];


#define SWAP_D(a) tmp_d = ((a)[idxi]); ((a)[idxi]) = ((a)[idxi+1]); ((a)[idxi+1]) = tmp_d;

static void
simint_sort_multi_shellpair(struct simint_multi_shellpair * P)
{
    if(P->nprim == 1)
        return;

    int idx = 0;
    double tmp_d;

    for(int ab = 0; ab < P->nshell12; ab++)
    {
        int swapped;

        do { 
            swapped = 0;
            for(int i = 0; i < (P->nprim12[ab]-1); i++)
            {
                const int idxi = idx + i;
                if(P->screen[idxi] > P->screen[idxi+1])
                {
                    SWAP_D(P->x)
                    SWAP_D(P->y)
                    SWAP_D(P->z)
                    SWAP_D(P->PA_x)
                    SWAP_D(P->PA_y)
                    SWAP_D(P->PA_z)
                    SWAP_D(P->PB_x)
                    SWAP_D(P->PB_y)
                    SWAP_D(P->PB_z)
                    SWAP_D(P->alpha)
                    SWAP_D(P->prefac)
                    SWAP_D(P->screen)
                    swapped = 1;
                }
            }
        } while(swapped);

        
        idx += P->nprim12[ab];
    }    
}



void simint_initialize_shell(struct simint_shell * G)
{
    G->ptr = NULL;
    G->memsize = 0;
}


// Allocate a gaussian shell with correct alignment
void simint_allocate_shell(int nprim, struct simint_shell * G)
{
    const size_t memsize = 2 * nprim * sizeof(double);   

    if(G->memsize < memsize)
    {
        simint_free_shell(G);
        G->ptr = malloc(memsize);
        memset(G->ptr, 0, memsize); // need to zero padding
        G->memsize = memsize;
    }

    G->alpha = G->ptr;
    G->coef = G->ptr + nprim*sizeof(double);
}


void simint_free_shell(struct simint_shell * G)
{
    if(G->ptr != NULL)
        free(G->ptr);

    simint_initialize_shell(G);
}

void simint_copy_shell(struct simint_shell const * src,
                       struct simint_shell * dest)
{
    simint_allocate_shell(src->nprim, dest);

    dest->nprim = src->nprim;
    dest->am = src->am;
    dest->x = src->x;
    dest->y = src->y;
    dest->z = src->z;

    memcpy(dest->ptr, src->ptr, dest->memsize);
}


void simint_normalize_shells(int n, struct simint_shell * G)
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


void simint_create_shell(int nprim, int am, double x, double y, double z,
                         double const * alpha,
                         double const * coef,
                         struct simint_shell * G)
{
    simint_allocate_shell(nprim, G);
    G->am = am;
    G->nprim = nprim;
    G->x = x;
    G->y = y;
    G->z = z;
    memcpy(G->alpha, alpha, nprim * sizeof(double));
    memcpy(G->coef, coef, nprim * sizeof(double));
}


void simint_initialize_multi_shellpair(struct simint_multi_shellpair * P)
{
    P->ptr = NULL;
    P->memsize = 0;
}


void simint_allocate_multi_shellpair(int na, struct simint_shell const * A,
                                     int nb, struct simint_shell const * B,
                                     struct simint_multi_shellpair * P)
{
    struct simint_shell AB[2*na*nb];

    int ij = 0;
    for(int i = 0; i < na; ++i)
    for(int j = 0; j < nb; ++j)
    {
        AB[ij] = A[i];
        AB[ij+1] = B[j];
        ij += 2;
    }

    simint_allocate_multi_shellpair2(na*nb, AB, P);
}


void simint_allocate_multi_shellpair2(int npair, struct simint_shell const * AB,
                                      struct simint_multi_shellpair * P)
{
    int nprim = 0;

    int batchprim = 0;
    int ij = 0;
    for(int i = 0; i < npair; i++)
    {
        if(compare_shell(&AB[ij], &AB[ij+1]))
            batchprim += ((AB[ij].nprim)*(AB[ij].nprim+1))/2;
        else
            batchprim += AB[ij].nprim * AB[ij+1].nprim;

        int ip1 = i+1;
        if((ip1 % SIMINT_NSHELL_SIMD) == 0 || ip1 >= npair)
        {        
            nprim += SIMINT_SIMD_ROUND(batchprim);
            batchprim = 0;
        }

        ij += 2;
    }


    int nshell12 = npair;

    const size_t dprim_size = nprim * sizeof(double);
    const size_t ishell12_size = nshell12 * sizeof(int);
    const size_t dshell12_size = nshell12 * sizeof(double);

    const size_t memsize = dprim_size*12 + dshell12_size*3 + ishell12_size;

    // Allocate one large space.
    // Only allocate if the currently allocated memory is too small
    if(P->memsize < memsize)
    {
        simint_free_multi_shellpair(P);
        P->ptr = ALLOC(memsize); 
        memset(P->ptr, 0, memsize);
        P->memsize = memsize;
    }

    int dcount = 0;
    P->x          = P->ptr + dprim_size*(dcount++);
    P->y          = P->ptr + dprim_size*(dcount++);
    P->z          = P->ptr + dprim_size*(dcount++);
    P->PA_x       = P->ptr + dprim_size*(dcount++);
    P->PA_y       = P->ptr + dprim_size*(dcount++);
    P->PA_z       = P->ptr + dprim_size*(dcount++);
    P->PB_x       = P->ptr + dprim_size*(dcount++);
    P->PB_y       = P->ptr + dprim_size*(dcount++);
    P->PB_z       = P->ptr + dprim_size*(dcount++);
    P->alpha      = P->ptr + dprim_size*(dcount++);
    P->prefac     = P->ptr + dprim_size*(dcount++);
    P->screen     = P->ptr + dprim_size*(dcount++);

    // below are unaligned
    P->AB_x       = P->ptr + 12*dprim_size;
    P->AB_y       = P->ptr + 12*dprim_size +   dshell12_size;
    P->AB_z       = P->ptr + 12*dprim_size + 2*dshell12_size;
    P->nprim12    = P->ptr + 12*dprim_size + 3*dshell12_size;
}


void simint_free_multi_shellpair(struct simint_multi_shellpair * P)
{
    if(P->ptr != NULL)
        FREE(P->ptr);

    simint_initialize_multi_shellpair(P);
}


void simint_fill_multi_shellpair(int na, struct simint_shell const * A,
                                 int nb, struct simint_shell const * B,
                                 struct simint_multi_shellpair * P,
                                 int screen)
{
    struct simint_shell AB[2*na*nb];

    int ij = 0;
    for(int i = 0; i < na; ++i)
    for(int j = 0; j < nb; ++j)
    {
        AB[ij] = A[i];
        AB[ij+1] = B[j];
        ij += 2;
    }

    simint_fill_multi_shellpair2(na*nb, AB, P, screen);
}


void simint_fill_multi_shellpair2(int npair, struct simint_shell const * AB,
                                  struct simint_multi_shellpair * P,
                                  int screen)
{
    int i, j, ij, sasb, idx;

    P->nshell12 = npair;
    P->nshell12_clip = npair; // by default, it's the same
    P->nprim = 0;

    // zero out
    memset(P->ptr, 0, P->memsize);

    sasb = 0;
    idx = 0;

    ij = 0;
    for(sasb = 0; sasb < npair; sasb++)
    {
        struct simint_shell const * A = &AB[ij];
        struct simint_shell const * B = &AB[ij+1];

        // compute the screening information
        if(screen)
            P->screen_max = simint_primscreen_schwarz_max(A, B, P->screen + idx);
        else
            P->screen_max = 0.0;

        // are these the same shells?
        const int same_shell = compare_shell(A, B);

        P->am1 = A->am;
        P->am2 = B->am;

        // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
        const double Xab_x = A->x - B->x;
        const double Xab_y = A->y - B->y;
        const double Xab_z = A->z - B->z;
        const double Xab = Xab_x*Xab_x + Xab_y*Xab_y + Xab_z*Xab_z;


        for(i = 0; i < A->nprim; ++i)
        {
            const double alpha_i = A->alpha[i];
            const double AxAa = A->x * alpha_i;
            const double AyAa = A->y * alpha_i;
            const double AzAa = A->z * alpha_i;

            int jend = B->nprim;
            if(same_shell)
                jend = (i+1);

            for(j = 0; j < jend; ++j)
            {
                const double alpha_j = B->alpha[j];
                const double ab_sum = alpha_i + alpha_j;
                const double oo_ab_sum = 1.0 / ab_sum;
                const double ab_mul = alpha_i * alpha_j;

                P->prefac[idx] = A->coef[i] * B->coef[j]
                                 * exp(-Xab * ab_mul * oo_ab_sum)
                                 * SQRT_TWO_PI_52 * oo_ab_sum;

                if(same_shell && (i != j))
                    P->prefac[idx] *= 2.0;

                P->x[idx] = (AxAa + alpha_j * B->x) * oo_ab_sum;
                P->y[idx] = (AyAa + alpha_j * B->y) * oo_ab_sum;
                P->z[idx] = (AzAa + alpha_j * B->z) * oo_ab_sum;
                P->PA_x[idx] = P->x[idx] - A->x;
                P->PA_y[idx] = P->y[idx] - A->y;
                P->PA_z[idx] = P->z[idx] - A->z;
                P->PB_x[idx] = P->x[idx] - B->x;
                P->PB_y[idx] = P->y[idx] - B->y;
                P->PB_z[idx] = P->z[idx] - B->z;

                P->alpha[idx] = ab_sum;
                ++idx;
            }
        }

        
        


        int nprim = (same_shell ? (A->nprim * (A->nprim+1))/2 : A->nprim * B->nprim);

        P->AB_x[sasb] = Xab_x;
        P->AB_y[sasb] = Xab_y;
        P->AB_z[sasb] = Xab_z;

        int sasbp1 = sasb + 1;

        if((sasbp1 % SIMINT_NSHELL_SIMD) == 0 || sasbp1 >= npair)
        {
            // fill in alpha = 1 until next boundary
            while(idx < SIMINT_SIMD_ROUND(idx))
            {
                P->alpha[idx++] = 1.0;
                nprim++;
            }
        }

        P->nprim12[sasb] = nprim;
        P->nprim += nprim;

        ij += 2;
    }


    // Sort, if we are screening
    if(screen)
        simint_sort_multi_shellpair(P);
}


void simint_create_multi_shellpair(int na, struct simint_shell const * A,
                                   int nb, struct simint_shell const * B,
                                   struct simint_multi_shellpair * P,
                                   int screen)
{
    simint_allocate_multi_shellpair(na, A, nb, B, P);
    simint_fill_multi_shellpair(na, A, nb, B, P, screen);
}


void simint_create_multi_shellpair2(int npair,
                                    struct simint_shell const * AB,
                                    struct simint_multi_shellpair * P,
                                    int screen)
{
    simint_allocate_multi_shellpair2(npair, AB, P);
    simint_fill_multi_shellpair2(npair, AB, P, screen);
}

