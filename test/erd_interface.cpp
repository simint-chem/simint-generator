#include <stdlib.h>
#include <string.h> // for memcpy
#include <math.h>

#include "eri/shell.h"

// 2 * pi**(2.5) * sqrt(pi) / 16
#define ERD_NORM_FIX 3.875784585037477521934539383387674400278161070735638461768067262975799364683

extern "C" {

void erd__gener_eri_batch_(const int *imax, const int *zmax, const int *nalpha, const int *ncoeff,
                           const int *ncsum, const int *ncgto1, const int *ncgto2,
                           const int *ncgto3, const int *ncgto4, const int *npgto1,
                           const int *npgto2, const int *npgto3, const int *npgto4,
                           const int *shell1, const int *shell2, const int *shell3,
                           const int *shell4, const double *x1, const double *y1, const double *z1,
                           const double *x2,const double *y2,const double *z2, const double *x3,
                           const double *y3,const double *z3,const double *x4, const double *y4, const double *z4,
                           const double *alpha, const double *cc, const int *ccbeg, const int *ccend,
                           const int *spheric,  const int *screen, int *icore,
                           int *nbatch, int * nfirst, double *zcore );

void erd__memory_eri_batch_(const int *nalpha, const int *ncoeff,
                            const int *ncgto1, const int *ncgto2, const int *ncgto3, const int *ncgto4,
                            const int *npgto1, const int *npgto2, const int *npgto3, const int *npgto4,
                            const int *shell1, const int *shell2, const int *shell3, const int *shell4,
                            const double *x1, const double *y1, const double *z1, const double *x2, const double *y2,
                            const double *z2, const double *x3, const double *y3, const double *z3, const double *x4,
                            const double *y4, const double *z4, const double *alpha, const double *cc, const int *spheric,
                            int *imin, int *iopt, int *zmin, int *zopt);
}

double * dscratch;
int * iscratch;

int i_buffer_size;
int d_buffer_size;

double * cc;
double * alpha;

void ERD_Init(int na, struct gaussian_shell const * const restrict A,
              int nb, struct gaussian_shell const * const restrict B,
              int nc, struct gaussian_shell const * const restrict C,
              int nd, struct gaussian_shell const * const restrict D)
{
    int maxam = 0;

    // get scratch requirements
    int npgto1 = 0;
    int sh1 = 0;
    for(int i = 0; i < na; ++i)
    {
        if(npgto1 < A[i].nprim)
        {
            npgto1 = A[i].nprim;
            sh1 = i;
        }
        if(maxam < A[i].am)
            maxam = A[i].am;
    }

    int npgto2 = 0;
    int sh2 = 0;
    for(int i = 0; i < nb; ++i)
    {
        if(npgto2 < B[i].nprim)
        {
            npgto2 = B[i].nprim;
            sh2 = i;
        }
        if(maxam < B[i].am)
            maxam = B[i].am;
    }

    int npgto3 = 0;
    int sh3 = 0;
    for(int i = 0; i < nc; ++i)
    {
        if(npgto3 < C[i].nprim)
        {
            npgto3 = C[i].nprim;
            sh3 = i;
        }
        if(maxam < C[i].am)
            maxam = C[i].am;
    }

    int npgto4 = 0;
    int sh4 = 0;
    for(int i = 0; i < nd; ++i)
    {
        if(npgto4 < D[i].nprim)
        {
            npgto4 = D[i].nprim;
            sh4 = i;
        }
        if(maxam < D[i].am)
            maxam = D[i].am;
    }
  
    cc    = (double *)malloc((npgto1 + npgto2 + npgto3 + npgto4) * sizeof(double));
    alpha = (double *)malloc((npgto1 + npgto2 + npgto3 + npgto4) * sizeof(double));
    
    memcpy(cc                           , A[sh1].coef, A[sh1].nprim * sizeof(double)); 
    memcpy(cc + npgto1                  , B[sh2].coef, B[sh2].nprim * sizeof(double)); 
    memcpy(cc + npgto1 + npgto2         , C[sh3].coef, C[sh3].nprim * sizeof(double)); 
    memcpy(cc + npgto1 + npgto2 + npgto3, D[sh4].coef, D[sh4].nprim * sizeof(double)); 
    memcpy(alpha                           , A[sh1].alpha, A[sh1].nprim * sizeof(double)); 
    memcpy(alpha + npgto1                  , B[sh2].alpha, B[sh2].nprim * sizeof(double)); 
    memcpy(alpha + npgto1 + npgto2         , C[sh3].alpha, C[sh3].nprim * sizeof(double)); 
    memcpy(alpha + npgto1 + npgto2 + npgto3, D[sh4].alpha, D[sh4].nprim * sizeof(double)); 
    
    // todo - general contraction?
    int ncgto1 = 1;
    int ncgto2 = 1;
    int ncgto3 = 1;
    int ncgto4 = 1;

    double X[4] = {1.0, 1.0, 1.0, 1.0};
    double Y[4] = {1.0, 1.0, 1.0, 1.0};
    double Z[4] = {1.0, 1.0, 1.0, 1.0};
    int zmin = 0;
    int zopt = 0;
    int imin = 0;
    int iopt = 0;
    int spheric = 0;
    int npgto = npgto1 + npgto2 + npgto3 + npgto4;

    erd__memory_eri_batch_(&npgto, &npgto,
                           &ncgto4, &ncgto3, &ncgto2, &ncgto1,
                           &npgto4, &npgto3, &npgto2, &npgto1,
                           &maxam, &maxam, &maxam, &maxam,
                           &X[3], &Y[3], &Z[3],
                           &X[2], &Y[2], &Z[2],
                           &X[1], &Y[1], &Z[1],
                           &X[0], &Y[0], &Z[0],
                           alpha, cc, &spheric, &imin, &iopt, &zmin, &zopt);
                             

    iscratch = (int *)malloc(iopt * sizeof(int)); 
    dscratch = (double *)malloc(zopt * sizeof(double)); 

    i_buffer_size = iopt;
    d_buffer_size = zopt;
   
}





void ERD_Compute_shell(struct gaussian_shell const A,
                       struct gaussian_shell const B,
                       struct gaussian_shell const C,
                       struct gaussian_shell const D,
                       double * integrals)
{
    // todo - general contraction?
    int ncgto1 = 1;
    int ncgto2 = 1;
    int ncgto3 = 1;
    int ncgto4 = 1;
    int ncgto = ncgto1 + ncgto2 + ncgto3 + ncgto4;
    int ccbeg[4] = { 1, 1, 1, 1 };
    int ccend[4] = { D.nprim, C.nprim, B.nprim, A.nprim };

    int offset_i = 0;
    int offset_j = offset_i + D.nprim; 
    int offset_k = offset_j + C.nprim; 
    int offset_l = offset_k + B.nprim; 

    int npgto = A.nprim + B.nprim + C.nprim + D.nprim;

    memcpy(alpha + offset_i, D.alpha, D.nprim * sizeof(double));
    memcpy(alpha + offset_j, C.alpha, C.nprim * sizeof(double));
    memcpy(alpha + offset_k, B.alpha, B.nprim * sizeof(double));
    memcpy(alpha + offset_l, A.alpha, A.nprim * sizeof(double));
    memcpy(cc + offset_i, D.coef, D.nprim * sizeof(double));
    memcpy(cc + offset_j, C.coef, C.nprim * sizeof(double));
    memcpy(cc + offset_k, B.coef, B.nprim * sizeof(double));
    memcpy(cc + offset_l, A.coef, A.nprim * sizeof(double));

    int nbatch = 0;
    int spheric = 0;
    int screen = 0;
    int buffer_offset = 0;

    erd__gener_eri_batch_(&i_buffer_size, &d_buffer_size,
                          &npgto, &npgto, &ncgto,
                          &ncgto4, &ncgto3, &ncgto2, &ncgto1,
                          &D.nprim, &C.nprim, &B.nprim, &A.nprim,
                          &D.am, &C.am, &B.am, &A.am,
                          &D.x, &D.y, &D.z,
                          &C.x, &C.y, &C.z,
                          &B.x, &B.y, &B.z,
                          &A.x, &A.y, &A.z,
                          alpha, cc, ccbeg, ccend, &spheric, &screen,
                          iscratch, &nbatch, &buffer_offset, dscratch);

    // remember, fortran has 1-based indexing
    //printf("ERD: %d integrals computed\n", nbatch); 
    //printf("Offset: %d\n", buffer_offset-1); 

        
    memcpy(integrals, dscratch + buffer_offset - 1, nbatch * sizeof(double));

    // fix normalization
    for(int i = 0; i < nbatch; i++)
        integrals[i] *= ERD_NORM_FIX;
}


void ERD_Finalize(void)
{
    free(dscratch);
    free(iscratch);
    free(cc);
    free(alpha);
}


