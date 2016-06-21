#include <stdio.h>
#include <stdlib.h>
#include <string.h> // for memcpy
#include <math.h>

#include "test/ERD.hpp"
#include "simint/vectorization/vectorization.h"

#define MAX(x, y) ((x) > (y) ? (x) : (y))


ERD_ERI::ERD_ERI(int am1, int nprim1, int ncgto1,
                 int am2, int nprim2, int ncgto2,
                 int am3, int nprim3, int ncgto3,
                 int am4, int nprim4, int ncgto4)
{
    Init_(am1, nprim1, ncgto1,
          am2, nprim2, ncgto2,
          am3, nprim3, ncgto3,
          am4, nprim4, ncgto4);
}



ERD_ERI::ERD_ERI(int am, int nprim, int ncgto)
{
    Init_(am, nprim, ncgto,
          am, nprim, ncgto,
          am, nprim, ncgto,
          am, nprim, ncgto);
}



ERD_ERI::ERD_ERI(int na, struct gaussian_shell const * const restrict A,
                 int nb, struct gaussian_shell const * const restrict B,
                 int nc, struct gaussian_shell const * const restrict C,
                 int nd, struct gaussian_shell const * const restrict D)
{
    int am1 = A[0].am;
    int am2 = B[0].am;
    int am3 = C[0].am;
    int am4 = D[0].am;

    int npgto1 = 0;
    for(int i = 0; i < na; ++i)
        npgto1 = MAX(npgto1, A[i].nprim);

    int npgto2 = 0;
    for(int i = 0; i < nb; ++i)
        npgto2 = MAX(npgto2, B[i].nprim);

    int npgto3 = 0;
    for(int i = 0; i < nc; ++i)
        npgto3 = MAX(npgto3, C[i].nprim);

    int npgto4 = 0;
    for(int i = 0; i < nd; ++i)
        npgto4 = MAX(npgto4, D[i].nprim);

    // todo - general contraction?
    int ncgto1 = 1;
    int ncgto2 = 1;
    int ncgto3 = 1;
    int ncgto4 = 1;
 
    Init_(am1, npgto1, ncgto1,
              am2, npgto2, ncgto2,
              am3, npgto3, ncgto3,
              am4, npgto4, ncgto4);
}



ERD_ERI::~ERD_ERI(void)
{
    FREE(dscratch);
    FREE(iscratch);
    FREE(cc);
    FREE(alpha);
}



void ERD_ERI::Init_(int am1, int nprim1, int ncgto1,
               int am2, int nprim2, int ncgto2,
               int am3, int nprim3, int ncgto3,
               int am4, int nprim4, int ncgto4)
{
    int nprim = nprim1 + nprim2 + nprim3 + nprim4;

    cc    = (double *)ALLOC(nprim * sizeof(double));
    alpha = (double *)ALLOC(nprim * sizeof(double));

    for(int i = 0; i < nprim; i++)
    {
        cc[i] = 2.0;
        alpha[i] = 10.0;
    }
 
    double X[4] = {1.0, 2.0, 3.0, 4.0};
    double Y[4] = {1.0, 2.0, 3.0, 4.0};
    double Z[4] = {1.0, 2.0, 3.0, 4.0};
    int zmin = 0;
    int zopt = 0;
    int imin = 0;
    int iopt = 0;
    int spheric = 0;
    int maxam = MAX(MAX(am1, am2), MAX(am3, am4));

    erd__memory_eri_batch_(&nprim, &nprim,
                           &ncgto1, &ncgto2, &ncgto3, &ncgto4,
                           &nprim1, &nprim2, &nprim3, &nprim4,
                           &maxam,  &maxam,  &maxam,  &maxam,
                           &X[0], &Y[0], &Z[0],
                           &X[1], &Y[1], &Z[1],
                           &X[2], &Y[2], &Z[2],
                           &X[3], &Y[3], &Z[3],
                           alpha, cc, &spheric, &imin, &iopt, &zmin, &zopt);
                             
    iscratch = (int *)ALLOC(iopt * sizeof(int)); 
    dscratch = (double *)ALLOC(zopt * sizeof(double)); 

    i_buffer_size = iopt;
    d_buffer_size = zopt;
}






TimeContrib ERD_ERI::Compute_shell_(struct gaussian_shell const A,
                                    struct gaussian_shell const B,
                                    struct gaussian_shell const C,
                                    struct gaussian_shell const D,
                                    double * integrals)
{
    TimeContrib time;
    TimerType ticks0, ticks1;

    CLOCK(ticks0);
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
    CLOCK(ticks1);
    time.copy_data = ticks1 - ticks0;


    CLOCK(ticks0);
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

    CLOCK(ticks1);
    time.integrals = ticks1 - ticks0;

    // remember, fortran has 1-based indexing
    //printf("ERD: %d integrals computed\n", nbatch); 
    //printf("Offset: %d\n", buffer_offset-1); 

    CLOCK(ticks0);
    if(nbatch > 0)    
        memcpy(integrals, dscratch + buffer_offset - 1, nbatch * sizeof(double));
    CLOCK(ticks1);
    time.permute = ticks1 - ticks0;

    return time;
}


TimeContrib ERD_ERI::Integrals(struct gaussian_shell const * const restrict A, int nshell1,
                               struct gaussian_shell const * const restrict B, int nshell2,
                               struct gaussian_shell const * const restrict C, int nshell3,
                               struct gaussian_shell const * const restrict D, int nshell4,
                               double * const integrals)
{
    const int nshell1234 = nshell1 * nshell2 * nshell3 * nshell4;

    const int am1 = A[0].am;
    const int am2 = B[0].am;
    const int am3 = C[0].am;
    const int am4 = D[0].am;
    const int ncart1234 = NCART(am1) * NCART(am2) * NCART(am3) * NCART(am4);

    const int ncart = nshell1234 * ncart1234;
    std::fill(integrals, integrals + ncart, 0.0);

    TimeContrib totaltime;

    int idx = 0;
    for(int i = 0; i < nshell1; i++)
    for(int j = 0; j < nshell2; j++)
    for(int k = 0; k < nshell3; k++)
    for(int l = 0; l < nshell4; l++)
    {
        totaltime += Compute_shell_(A[i], B[j], C[k], D[l], integrals + idx);
        idx += ncart1234;
    }

    return totaltime;
}


void normalize_gaussian_shells_erd(int n, struct gaussian_shell * const restrict G)
{
    //df[i] = (i-1)!!
    const double df[] = {
                          1.0,
                          1.0,
                          1.0,
                          2.0,
                          3.0,
                          8.0,
                          15.0,
                          48.0,
                          105.0,
                          384.0,
                          945.0,
                          3840.0,
                          10395.0,
                          46080.0,
                          135135.0,
                          645120.0,
                          2027025.0,
                          10321920.0,
                          34459425.0,
                          185794560.0,
                          654729075.0,
                          3715891200.0,
                          13749310575.0,
                          81749606400.0,
                          316234143225.0,
                          1961990553600.0,
                          7905853580625.0
                        };
 
    /* 
     * libERD version of normalization 
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

                double temp = G[i].coef[j] * G[i].coef[k];
                double temp2 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
                temp2 = pow(temp2, m);
                temp = temp * temp2;
                sum += temp;
            }
        }

        double prefac = 1.0;
        if(G[i].am > 1)
            prefac = pow(2.0, 2*G[i].am) / df[2*G[i].am];

        const double norm = sqrt(prefac / sum);

        // apply the normalization
        for (int j = 0; j < G[i].nprim; ++j)
            G[i].coef[j] *= norm * pow(G[i].alpha[j], 0.5*m);
    }

}


