#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "simint/constants.h"
#include "simint/shell/shell.h"
#include "simint/shell/shell_screen.h"
#include "simint/shell/shell_constants.h"
#include "simint/vectorization/vectorization.h"


#if defined(__ICC) || defined(__INTEL_COMPILER)
    #pragma warning(disable:1338)                  // Pointer arithmetic on void
#elif defined(__GNUC__) || defined(__GNUG__)
    #pragma GCC diagnostic ignored "-Wpointer-arith"
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
                if(P->screen[idxi] < P->screen[idxi+1])
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

                    #if SIMINT_OSTEI_MAXDER > 0
                    SWAP_D(P->alpha2);
                    SWAP_D(P->beta2);
                    #endif

                    swapped = 1;
                }
            }
        } while(swapped);

        idx += P->nprim12[ab];

        int ab1 = ab+1;
        if( (ab1 % SIMINT_NSHELL_SIMD) == 0 || (ab1 >= P->nshell12) )
            idx = SIMINT_SIMD_ROUND(idx);

    }
}


static void simint_allocate_multi_shellpair_base(int npair, int nprim,
                                                 struct simint_multi_shellpair * P,
                                                 int screen_method)
{
    const size_t dprim_size = nprim * sizeof(double);
    const size_t ishell12_size = npair * sizeof(int);
    const size_t dshell12_size = npair * sizeof(double);

    int nprim_arr = 11;
    if(screen_method)
        nprim_arr++;
        
    #if SIMINT_OSTEI_MAXDER > 0
    nprim_arr += 2;
    #endif

    const size_t memsize = dprim_size*nprim_arr + dshell12_size*3 + ishell12_size;

    // Allocate one large space.
    // Only allocate if the currently allocated memory is too small
    if(P->memsize < memsize)
    {
        simint_free_multi_shellpair(P);
        P->ptr = SIMINT_ALLOC(memsize);
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

    if(screen_method)
        P->screen = P->ptr + dprim_size*(dcount++);
    else
        P->screen = NULL;

    #if SIMINT_OSTEI_MAXDER > 0
    P->alpha2     = P->ptr + dprim_size*(dcount++);
    P->beta2      = P->ptr + dprim_size*(dcount++);
    #endif

    P->AB_x       = P->ptr + dcount*dprim_size;
    P->AB_y       = P->ptr + dcount*dprim_size +   dshell12_size;
    P->AB_z       = P->ptr + dcount*dprim_size + 2*dshell12_size;
    P->nprim12    = P->ptr + dcount*dprim_size + 3*dshell12_size;
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
                                     struct simint_multi_shellpair * P,
                                     int screen_method)
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

    simint_allocate_multi_shellpair2(na*nb, AB, P, screen_method);
}


void simint_allocate_multi_shellpair2(int npair, struct simint_shell const * AB,
                                      struct simint_multi_shellpair * P,
                                      int screen_method)
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

    simint_allocate_multi_shellpair_base(npair, nprim, P, screen_method);
}


void simint_free_multi_shellpair(struct simint_multi_shellpair * P)
{
    if(P->ptr != NULL)
        SIMINT_FREE(P->ptr);

    simint_initialize_multi_shellpair(P);
}


void simint_fill_multi_shellpair(int na, struct simint_shell const * A,
                                 int nb, struct simint_shell const * B,
                                 struct simint_multi_shellpair * P,
                                 int screen_method)
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

    simint_fill_multi_shellpair2(na*nb, AB, P, screen_method);
}


void simint_fill_multi_shellpair2(int npair, struct simint_shell const * AB,
                                  struct simint_multi_shellpair * P,
                                  int screen_method)
{
    int i, j, ij, sasb, idx;

    P->nshell12 = npair;
    P->nshell12_clip = npair; // by default, it's the same
    P->nprim = 0;

    // zero out
    // This is not really needed, and can be expensive in
    // direct code
    // It's not needed since everything will be taken care of
    // with a prefactor of 0.0
    //memset(P->ptr, 0, P->memsize);

    sasb = 0;
    idx = 0;

    ij = 0;
    for(sasb = 0; sasb < npair; sasb++)
    {
        struct simint_shell const * A = &AB[ij];
        struct simint_shell const * B = &AB[ij+1];

        // compute the screening information
        if(screen_method)
        {
            double m = simint_primscreen(A, B, P->screen + idx, screen_method);
            P->screen_max = (m > P->screen_max ? m : P->screen_max);
        }
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
                const double ab_mul = alpha_i * alpha_j;

                // multiplying by reciprocal of ab_sum resulted
                // in small numerical differences
                //const double oo_ab_sum = 1.0 / ab_sum;

                P->prefac[idx] = A->coef[i] * B->coef[j]
                                 * exp(-Xab * ab_mul / ab_sum)
                                 * SQRT_TWO_PI_52 / ab_sum;

                if(same_shell && (i != j))
                    P->prefac[idx] *= 2.0;

                P->x[idx] = (AxAa + alpha_j * B->x) / ab_sum;
                P->y[idx] = (AyAa + alpha_j * B->y) / ab_sum;
                P->z[idx] = (AzAa + alpha_j * B->z) / ab_sum;
                P->PA_x[idx] = P->x[idx] - A->x;
                P->PA_y[idx] = P->y[idx] - A->y;
                P->PA_z[idx] = P->z[idx] - A->z;
                P->PB_x[idx] = P->x[idx] - B->x;
                P->PB_y[idx] = P->y[idx] - B->y;
                P->PB_z[idx] = P->z[idx] - B->z;

                P->alpha[idx] = ab_sum;

                #if SIMINT_OSTEI_MAXDER > 0
                if(same_shell && (i != j))
                {
                    // there is already a factor of 2.0 in the prefac,
                    // so we don't need it here (work it out and see
                    // for yourself why we don't need it)
                    P->alpha2[idx] = (alpha_i + alpha_j);
                    P->beta2[idx] = (alpha_i + alpha_j);
                }
                else
                {
                    P->alpha2[idx] = 2.0 * alpha_i;
                    P->beta2[idx] = 2.0 * alpha_j;
                }
                #endif

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
            // fill in some members until next boundary
            while(idx < SIMINT_SIMD_ROUND(idx))
            {
                P->alpha[idx] = 1.0;
                P->prefac[idx] = 0.0;
                P->x[idx] = 0.0;
                P->y[idx] = 0.0;
                P->z[idx] = 0.0;
                P->PA_x[idx] = 0.0;
                P->PA_y[idx] = 0.0;
                P->PA_z[idx] = 0.0;
                P->PB_x[idx] = 0.0;
                P->PB_y[idx] = 0.0;
                P->PB_z[idx] = 0.0;

                #if SIMINT_OSTEI_MAXDER > 0
                P->alpha2[idx] = 1.0;
                P->beta2[idx] = 1.0;
                #endif

                idx++;
            }
        }

        P->nprim12[sasb] = nprim;
        P->nprim += nprim;

        ij += 2;
    }


    // If we are screening, sort the primitives within each shell
    if(screen_method)
        simint_sort_multi_shellpair(P);
}


void simint_create_multi_shellpair(int na, struct simint_shell const * A,
                                   int nb, struct simint_shell const * B,
                                   struct simint_multi_shellpair * P,
                                   int screen_method)
{
    simint_allocate_multi_shellpair(na, A, nb, B, P, screen_method);
    simint_fill_multi_shellpair(na, A, nb, B, P, screen_method);
}


void simint_create_multi_shellpair2(int nshell12,
                                    struct simint_shell const * AB,
                                    struct simint_multi_shellpair * P,
                                    int screen_method)
{
    simint_allocate_multi_shellpair2(nshell12, AB, P, screen_method);
    simint_fill_multi_shellpair2(nshell12, AB, P, screen_method);
}


void simint_cat_multi_shellpair(int nmpair,
                                struct simint_multi_shellpair const ** Pin,
                                struct simint_multi_shellpair * Pout,
                                int screen_method)
{
    if(nmpair <= 0)
        return; //??

    // Total number of shell pair
    int nshell12 = 0;
    for(int i = 0; i < nmpair; i++)
        nshell12 += Pin[i]->nshell12;

    // How many primitives
    int nprim = 0;
    int batchprim = 0;

    int ipair = 0;
    for(int i = 0; i < nmpair; i++)
    for(int j = 0; j < Pin[i]->nshell12; j++)
    {
        batchprim += Pin[i]->nprim12[j];

        int ip1 = ipair+1;
        if((ip1 % SIMINT_NSHELL_SIMD) == 0 || ip1 >= nshell12)
        {
            nprim += SIMINT_SIMD_ROUND(batchprim);
            batchprim = 0;
        }

        ipair++;
    }

    // (re)allocate Pout
    simint_allocate_multi_shellpair_base(nshell12, nprim, Pout, screen_method);

    // these should all have the same AM
    Pout->am1 = Pin[0]->am1;
    Pout->am2 = Pin[0]->am2;
    Pout->screen_max = 0.0;

    // now copy data
    int idx = 0;
    int sasb = 0;
    for(int i = 0; i < nmpair; i++)
    {
        for(int j = 0; j < Pin[i]->nshell12; j++)
        {
            for(int p = 0; p < Pin[i]->nprim12[j]; p++)
            {
                Pout->x[idx] = Pin[i]->x[p];
                Pout->y[idx] = Pin[i]->y[p];
                Pout->z[idx] = Pin[i]->z[p];
                Pout->PA_x[idx] = Pin[i]->PA_x[p];
                Pout->PA_y[idx] = Pin[i]->PA_y[p];
                Pout->PA_z[idx] = Pin[i]->PA_z[p];
                Pout->PB_x[idx] = Pin[i]->PB_x[p];
                Pout->PB_y[idx] = Pin[i]->PB_y[p];
                Pout->PB_z[idx] = Pin[i]->PB_z[p];
                Pout->x[idx] = Pin[i]->x[p];
                Pout->y[idx] = Pin[i]->y[p];
                Pout->z[idx] = Pin[i]->z[p];

                Pout->alpha[idx] = Pin[i]->alpha[p];
                Pout->prefac[idx] = Pin[i]->prefac[p];

                #if SIMINT_OSTEI_MAXDER > 0
                Pout->alpha2[idx] = Pin[i]->alpha2[p];
                Pout->beta2[idx] = Pin[i]->beta2[p];
                #endif

                if(screen_method)
                    Pout->screen[idx] = Pin[i]->screen[p];

                if(Pin[i]->screen_max > Pout->screen_max)
                    Pout->screen_max = Pin[i]->screen_max;

                idx++;
            }

            
            Pout->AB_x[sasb] = Pin[i]->AB_x[j];
            Pout->AB_y[sasb] = Pin[i]->AB_y[j];
            Pout->AB_z[sasb] = Pin[i]->AB_z[j];
            Pout->nprim12[sasb] = Pin[i]->nprim12[j];

            int sasbp1 = sasb + 1;

            if((sasbp1 % SIMINT_NSHELL_SIMD) == 0 || sasbp1 >= nshell12)
            {
                // fill in some members until next boundary
                while(idx < SIMINT_SIMD_ROUND(idx))
                {
                    Pout->alpha[idx] = 1.0;
                    Pout->prefac[idx] = 0.0;
                    Pout->x[idx] = 0.0;
                    Pout->y[idx] = 0.0;
                    Pout->z[idx] = 0.0;
                    Pout->PA_x[idx] = 0.0;
                    Pout->PA_y[idx] = 0.0;
                    Pout->PA_z[idx] = 0.0;
                    Pout->PB_x[idx] = 0.0;
                    Pout->PB_y[idx] = 0.0;
                    Pout->PB_z[idx] = 0.0;

                    #if SIMINT_OSTEI_MAXDER > 0
                    Pout->alpha2[idx] = 1.0;
                    Pout->beta2[idx] = 1.0;
                    #endif

                    idx++;
                }
            }

            sasb++;
        }

        Pout->nprim += Pin[i]->nprim; 

    }

    Pout->nshell12 = nshell12;
    Pout->nshell12_clip = nshell12;

}



/*
void
simint_prune_multi_shellpair(struct simint_multi_shellpair const * P,
                             struct simint_multi_shellpair * out,
                             double screen_max, double screen_tol)
{
    const double screen_tol2 = screen_tol * screen_tol;
    const double screen_val = screen_tol2 / screen_max;

    simint_allocate_multi_shellpair_base(P->nshell12, P->nprim, out);
    memset(out->ptr, 0, out->memsize);
    out->am1 = P->am1;
    out->am2 = P->am2;
    out->nshell12 = P->nshell12;
    out->nshell12_clip = P->nshell12_clip;
    out->screen_max = P->screen_max;

    int read_idx = 0;  // index we are currently reading from
    int write_idx = 0; // index we are currently writing to

    for(int ab = 0; ab < P->nshell12; ab++)
    {
        out->AB_x[ab] = P->AB_x[ab];
        out->AB_y[ab] = P->AB_y[ab];
        out->AB_z[ab] = P->AB_z[ab];

        int shell_nprim = 0;
        for(int i = 0; i < P->nprim12[ab]; i++)
        {
            if(P->screen[read_idx] > screen_val)
            {
                out->x[write_idx] = P->x[read_idx];
                out->y[write_idx] = P->y[read_idx];
                out->z[write_idx] = P->z[read_idx];
                out->PA_x[write_idx] = P->PA_x[read_idx];
                out->PA_y[write_idx] = P->PA_y[read_idx];
                out->PA_z[write_idx] = P->PA_z[read_idx];
                out->PB_x[write_idx] = P->PB_x[read_idx];
                out->PB_y[write_idx] = P->PB_y[read_idx];
                out->PB_z[write_idx] = P->PB_z[read_idx];
                out->alpha[write_idx] = P->alpha[read_idx];
                out->prefac[write_idx] = P->prefac[read_idx];
                out->screen[write_idx] = P->screen[read_idx];
                write_idx++;
                shell_nprim++;
            }
            read_idx++;
        }

        int ab1 = ab+1;
        if((ab1 % SIMINT_NSHELL_SIMD) == 0 || ab1 >= P->nshell12)
        {
            // coefficients should be zero due to the memset above
            while(write_idx < SIMINT_SIMD_ROUND(write_idx))
            {
                out->alpha[write_idx++] = 1.0;
                shell_nprim++;
            }
        }

        out->nprim12[ab] = shell_nprim;
    }

    out->nprim = write_idx;
}
*/
