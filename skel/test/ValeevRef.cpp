///////////////////////////////////////////
// Reference calculation of ERI
// Adapted from libint 2.0.5, by E. Valeev
// https://github.com/evaleev/libint
// http://www.valeevgroup.chem.vt.edu
//
// Originally released under GPL v2
///////////////////////////////////////////
//
// LIBINT (version 2) - a library for the evaluation of molecular integrals of many-body
// operators over Gaussian functions
// Copyright (C) 2004-2013 Edward F. Valeev
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program (see file LICENSE); if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
///////////////////////////////////////////


#include <math.h>
#include <stdlib.h>

#include "test/ValeevRef.hpp"
#include "test/Common.hpp"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define NUMDIFF_DELTA 1e-8

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static long double * fac;
static long double * df;
static long double ** bc;



static long double norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
                         long double alpha1)
{
    return powl(2 * alpha1 / M_PI, 0.75) * powl(4 * alpha1, 0.5 * (l1 + m1 + n1))
           / sqrtl(df[2 * l1] * df[2 * m1] * df[2 * n1]);
}

static long double* init_array(unsigned long int size)
{
    long double* result = (long double *)malloc(size * sizeof(long double));
    for (unsigned long int i = 0; i < size; i++)
        result[i] = 0.0;
    return result;
}

static void free_array(long double* array)
{
    free(array);
}


// copies a shell and moves it by delta in direction dir
// (used in numerical differentiation)
static void copy_gaussian_move(simint_shell const * src,
                               simint_shell * dest,
                               int dir, double delta)
{
    simint_copy_shell(src, dest);
    if(dir == 0)
        dest->x += delta;
    else if(dir == 1)
        dest->y += delta;
    else if(dir == 2)
        dest->z += delta;
}


static void FormDeriv(int ncart_abcd,
                      double const * buf_p, double const * buf_m,
                      double * dest, double delta)
{
    for(int i = 0; i < ncart_abcd; i++)
        dest[i] = (buf_p[i] - buf_m[i]) / (2.0*delta);
}



static void Valeev_F(long double *F, int n, long double x)
{
    int i, m;
    int m2;
    long double t2;
    long double num;
    long double sum;
    long double term1;
    const long double K = 0.8862269254527580136490837416705725913987747280611935641069038949264556422955160906874753283692723327l;
    long double et;


    if (x > 20.0)   /* For big t's do upward recursion */
    {
        t2 = 2 * x;
        et = expl(-x);
        x = sqrtl(x);
        F[0] = K * erfl(x) / x;
        for (m = 0; m <= n - 1; m++)
        {
            F[m + 1] = ((2 * m + 1) * F[m] - et) / (t2);
        }
    }
    else
    {
        /* For smaller t's compute F with highest n using
         asymptotic series (see I. Shavitt in
         Methods in Computational Physics, ed. B. Alder eta l,
         vol 2, 1963, page 8) */
        et = expl(-x);
        t2 = 2 * x;
        m2 = 2 * n;
        num = df[m2];
        i = 0;
        sum = 1.0 / (m2 + 1);
        do
        {
            i++;
            num = num * t2;
            term1 = num / df[m2 + 2 * i + 2];
            sum += term1;
        }
        while (fabsl(term1) > EPS && i < MAXFAC);
        F[n] = sum * et;
        for (m = n - 1; m >= 0; m--)   /* And then do downward recursion */
        {
            F[m] = (t2 * F[m + 1] + et) / (2 * m + 1);
        }
    }
}

static long double ValeevRef_eri(int l1, int m1, int n1, long double alpha1,
                                 const long double* A, int l2, int m2, int n2,
                                 long double alpha2, const long double* B, int l3, int m3,
                                 int n3, long double alpha3, const long double* C, int l4,
                                 int m4, int n4, long double alpha4, const long double* D,
                                 int norm_flag)
{

    const long double gammap = alpha1 + alpha2;
    const long double Px = (alpha1 * A[0] + alpha2 * B[0]) / gammap;
    const long double Py = (alpha1 * A[1] + alpha2 * B[1]) / gammap;
    const long double Pz = (alpha1 * A[2] + alpha2 * B[2]) / gammap;
    const long double PAx = Px - A[0];
    const long double PAy = Py - A[1];
    const long double PAz = Pz - A[2];
    const long double PBx = Px - B[0];
    const long double PBy = Py - B[1];
    const long double PBz = Pz - B[2];
    const long double AB2 = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1]
                       - B[1]) + (A[2] - B[2]) * (A[2] - B[2]);

    const long double gammaq = alpha3 + alpha4;
    const long double gammapq = gammap * gammaq / (gammap + gammaq);
    const long double Qx = (alpha3 * C[0] + alpha4 * D[0]) / gammaq;
    const long double Qy = (alpha3 * C[1] + alpha4 * D[1]) / gammaq;
    const long double Qz = (alpha3 * C[2] + alpha4 * D[2]) / gammaq;
    const long double QCx = Qx - C[0];
    const long double QCy = Qy - C[1];
    const long double QCz = Qz - C[2];
    const long double QDx = Qx - D[0];
    const long double QDy = Qy - D[1];
    const long double QDz = Qz - D[2];
    const long double CD2 = (C[0] - D[0]) * (C[0] - D[0]) + (C[1] - D[1]) * (C[1]
                       - D[1]) + (C[2] - D[2]) * (C[2] - D[2]);

    const long double PQx = Px - Qx;
    const long double PQy = Py - Qy;
    const long double PQz = Pz - Qz;
    const long double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

    int u1, u2, v1, v2, w1, w2, tx, ty, tz, txmax, tymax, tzmax;
    int i, j, k;
    int lp, lq, mp, mq, np, nq;
    int zeta;
    long double *flp, *flq, *fmp, *fmq, *fnp, *fnq;
    long double *F;
    long double K1, K2;
    long double Gx, Gy, Gz;
    long double pfac;
    long double result = 0.0;
    long double tmp;
    int u1max, u2max, v1max, v2max, w1max, w2max;

    K1 = expl(-alpha1 * alpha2 * AB2 / gammap);
    K2 = expl(-alpha3 * alpha4 * CD2 / gammaq);
    pfac = 2 * powl(M_PI, 2.5) * K1 * K2 / (gammap * gammaq
            * sqrtl(gammap + gammaq));

    if (norm_flag > 0)
    {
        pfac *= norm_const(l1, m1, n1, alpha1);
        pfac *= norm_const(l2, m2, n2, alpha2);
        pfac *= norm_const(l3, m3, n3, alpha3);
        pfac *= norm_const(l4, m4, n4, alpha4);
    }

    F = init_array(l1 + l2 + l3 + l4 + m1 + m2 + m3 + m4 + n1 + n2 + n3 + n4 + 1);
    Valeev_F(F, l1 + l2 + l3 + l4 + m1 + m2 + m3 + m4 + n1 + n2 + n3 + n4,
           PQ2 * gammapq);

    flp = init_array(l1 + l2 + 1);
    for (k = 0; k <= l1 + l2; k++)
        for (i = 0; i <= MIN(k,l1); i++)
        {
            j = k - i;
            if (j > l2)
                continue;
            tmp = bc[l1][i] * bc[l2][j];
            if (l1 - i > 0)
                tmp *= powl(PAx, l1 - i);
            if (l2 - j > 0)
                tmp *= powl(PBx, l2 - j);
            flp[k] += tmp;
        }
    fmp = init_array(m1 + m2 + 1);
    for (k = 0; k <= m1 + m2; k++)
        for (i = 0; i <= MIN(k,m1); i++)
        {
            j = k - i;
            if (j > m2)
                continue;
            tmp = bc[m1][i] * bc[m2][j];
            if (m1 - i > 0)
                tmp *= powl(PAy, m1 - i);
            if (m2 - j > 0)
                tmp *= powl(PBy, m2 - j);
            fmp[k] += tmp;
        }
    fnp = init_array(n1 + n2 + 1);
    for (k = 0; k <= n1 + n2; k++)
        for (i = 0; i <= MIN(k,n1); i++)
        {
            j = k - i;
            if (j > n2)
                continue;
            tmp = bc[n1][i] * bc[n2][j];
            if (n1 - i > 0)
                tmp *= powl(PAz, n1 - i);
            if (n2 - j > 0)
                tmp *= powl(PBz, n2 - j);
            fnp[k] += tmp;
        }
    flq = init_array(l3 + l4 + 1);
    for (k = 0; k <= l3 + l4; k++)
        for (i = 0; i <= MIN(k,l3); i++)
        {
            j = k - i;
            if (j > l4)
                continue;
            tmp = bc[l3][i] * bc[l4][j];
            if (l3 - i > 0)
                tmp *= powl(QCx, l3 - i);
            if (l4 - j > 0)
                tmp *= powl(QDx, l4 - j);
            flq[k] += tmp;
        }
    fmq = init_array(m3 + m4 + 1);
    for (k = 0; k <= m3 + m4; k++)
        for (i = 0; i <= MIN(k,m3); i++)
        {
            j = k - i;
            if (j > m4)
                continue;
            tmp = bc[m3][i] * bc[m4][j];
            if (m3 - i > 0)
                tmp *= powl(QCy, m3 - i);
            if (m4 - j > 0)
                tmp *= powl(QDy, m4 - j);
            fmq[k] += tmp;
        }
    fnq = init_array(n3 + n4 + 1);
    for (k = 0; k <= n3 + n4; k++)
        for (i = 0; i <= MIN(k,n3); i++)
        {
            j = k - i;
            if (j > n4)
                continue;
            tmp = bc[n3][i] * bc[n4][j];
            if (n3 - i > 0)
                tmp *= powl(QCz, n3 - i);
            if (n4 - j > 0)
                tmp *= powl(QDz, n4 - j);
            fnq[k] += tmp;
        }

    for (lp = 0; lp <= l1 + l2; lp++)
        for (lq = 0; lq <= l3 + l4; lq++)
        {
            u1max = lp / 2;
            u2max = lq / 2;
            for (u1 = 0; u1 <= u1max; u1++)
                for (u2 = 0; u2 <= u2max; u2++)
                {
                    Gx = powl(-1, lp) * flp[lp] * flq[lq] * fac[lp] * fac[lq]
                         * powl(gammap, u1 - lp) * powl(gammaq, u2 - lq) * fac[lp + lq - 2
                                 * u1 - 2 * u2] * powl(gammapq, lp + lq - 2 * u1 - 2 * u2)
                         / (fac[u1] * fac[u2] * fac[lp - 2 * u1] * fac[lq - 2 * u2]);
                    for (mp = 0; mp <= m1 + m2; mp++)
                        for (mq = 0; mq <= m3 + m4; mq++)
                        {
                            v1max = mp / 2;
                            v2max = mq / 2;
                            for (v1 = 0; v1 <= v1max; v1++)
                                for (v2 = 0; v2 <= v2max; v2++)
                                {
                                    Gy = powl(-1, mp) * fmp[mp] * fmq[mq] * fac[mp] * fac[mq]
                                         * powl(gammap, v1 - mp) * powl(gammaq, v2 - mq) * fac[mp
                                                 + mq - 2 * v1 - 2 * v2] * powl(gammapq,
                                                         mp + mq - 2 * v1 - 2 * v2)
                                         / (fac[v1] * fac[v2] * fac[mp - 2 * v1]
                                            * fac[mq - 2 * v2]);
                                    for (np = 0; np <= n1 + n2; np++)
                                        for (nq = 0; nq <= n3 + n4; nq++)
                                        {
                                            w1max = np / 2;
                                            w2max = nq / 2;
                                            for (w1 = 0; w1 <= w1max; w1++)
                                                for (w2 = 0; w2 <= w2max; w2++)
                                                {
                                                    Gz = powl(-1, np) * fnp[np] * fnq[nq] * fac[np]
                                                         * fac[nq] * powl(gammap, w1 - np) * powl(gammaq,
                                                                 w2 - nq)
                                                         * fac[np + nq - 2 * w1 - 2 * w2]
                                                         * powl(gammapq, np + nq - 2 * w1 - 2 * w2)
                                                         / (fac[w1] * fac[w2] * fac[np - 2 * w1] * fac[nq
                                                                 - 2 * w2]);
                                                    txmax = (lp + lq - 2 * u1 - 2 * u2) / 2;
                                                    tymax = (mp + mq - 2 * v1 - 2 * v2) / 2;
                                                    tzmax = (np + nq - 2 * w1 - 2 * w2) / 2;
                                                    for (tx = 0; tx <= txmax; tx++)
                                                        for (ty = 0; ty <= tymax; ty++)
                                                            for (tz = 0; tz <= tzmax; tz++)
                                                            {
                                                                zeta = lp + lq + mp + mq + np + nq - 2 * u1 - 2
                                                                       * u2 - 2 * v1 - 2 * v2 - 2 * w1 - 2 * w2
                                                                       - tx - ty - tz;
                                                                result += Gx * Gy * Gz * F[zeta]
                                                                          * powl(-1, tx + ty + tz) * powl(
                                                                              PQx,
                                                                              lp + lq - 2
                                                                              * u1 - 2
                                                                              * u2 - 2
                                                                              * tx)
                                                                          * powl(PQy,
                                                                                mp + mq - 2 * v1 - 2 * v2 - 2 * ty)
                                                                          * powl(PQz,
                                                                                np + nq - 2 * w1 - 2 * w2 - 2 * tz)
                                                                          / (powl(
                                                                                 4,
                                                                                 u1 + u2 + tx + v1 + v2 + ty + w1
                                                                                 + w2 + tz) * powl(gammapq, tx)
                                                                             * powl(gammapq, ty) * powl(gammapq, tz)
                                                                             * fac[lp + lq - 2 * u1 - 2 * u2 - 2
                                                                                   * tx] * fac[tx] * fac[mp + mq - 2
                                                                                           * v1 - 2 * v2 - 2 * ty] * fac[ty]
                                                                             * fac[np + nq - 2 * w1 - 2 * w2 - 2
                                                                                   * tz] * fac[tz]);
                                                            }
                                                }
                                        }
                                }
                        }
                }
        }

    free_array(F);
    free_array(flp);
    free_array(fmp);
    free_array(fnp);
    free_array(flq);
    free_array(fmq);
    free_array(fnq);

    return result * pfac;
}





void ValeevRef_Init(void)
{
    int i = 0;
    int j = 0;

    fac = (long double *)malloc(MAXFAC * sizeof(long double));
    fac[0] = 1.0;
    for(i = 1; i < MAXFAC; i++)
        fac[i] = fac[i-1] * i;

    df = (long double *)malloc(2 * MAXFAC * sizeof(long double));
    df[0] = 1.0;
    df[1] = 1.0;
    df[2] = 1.0;
    for(i = 3; i < MAXFAC*2; i++)
       df[i] = (i - 1) * df[i - 2];

    bc = (long double **)malloc(MAXFAC * sizeof(long double*));
    for(i = 0; i < MAXFAC; i++)
    {
        bc[i] = (long double *)malloc(MAXFAC * sizeof(long double));
        for(j = 0; j <= i; j++)
            bc[i][j] = fac[i] / (fac[i-j] * fac[j]);
    }
}

void ValeevRef_Finalize(void)
{
    int i = 0;
    free(fac);
    free(df);

    for(i = 0; i < MAXFAC; i++)
        free(bc[i]);
    free(bc);
}


// Helper for calculating via simint_shell structures
static void ValeevRef_Deriv0(simint_shell const * const A, int nshellA,
                             simint_shell const * const B, int nshellB,
                             simint_shell const * const C, int nshellC,
                             simint_shell const * const D, int nshellD,
                             double * const integrals, bool normalize)
{
    int inorm = (normalize ? 1 : 0);

    const int am1 = A[0].am;
    const int am2 = B[0].am;
    const int am3 = C[0].am;
    const int am4 = D[0].am;

    const int n234 = nshellB * nshellC * nshellD * NCART(am1) * NCART(am2) * NCART(am3) * NCART(am4);

    for(int i = 0; i < nshellA; i++)
    {
        const int idxstart = i * n234;
        int idx = 0;

        for(int j = 0; j < nshellB; j++)
        for(int k = 0; k < nshellC; k++)
        for(int l = 0; l < nshellD; l++)
        {
            long double vA[3] = { A[i].x, A[i].y, A[i].z };
            long double vB[3] = { B[j].x, B[j].y, B[j].z };
            long double vC[3] = { C[k].x, C[k].y, C[k].z };
            long double vD[3] = { D[l].x, D[l].y, D[l].z };

            std::array<int, 3> g1 = {{am1, 0, 0}};
            do
            {
                std::array<int, 3> g2 = {{am2, 0, 0}};
                do
                {
                    std::array<int, 3> g3 = {{am3, 0, 0}};
                    do
                    {
                        std::array<int, 3> g4 = {{am4, 0, 0}};
                        do
                        {
                            double myint = 0.0;

                            for(int m = 0; m < A[i].nprim; m++)
                            for(int n = 0; n < B[j].nprim; n++)
                            for(int o = 0; o < C[k].nprim; o++)
                            for(int p = 0; p < D[l].nprim; p++)
                            {

                                double val = (double)ValeevRef_eri(g1[0], g1[1], g1[2], A[i].alpha[m], vA,
                                                                   g2[0], g2[1], g2[2], B[j].alpha[n], vB,
                                                                   g3[0], g3[1], g3[2], C[k].alpha[o], vC,
                                                                   g4[0], g4[1], g4[2], D[l].alpha[p], vD, inorm);
                                myint += val * A[i].coef[m] * B[j].coef[n] * C[k].coef[o] * D[l].coef[p];

                            }

                            integrals[idxstart + idx] = myint;
                            idx++;

                        } while(IterateGaussian(g4));
                    } while(IterateGaussian(g3));
                } while(IterateGaussian(g2));
            } while(IterateGaussian(g1));
        }
    }
}


static void ValeevRef_Deriv1(simint_shell const * const A, int nshellA,
                             simint_shell const * const B, int nshellB,
                             simint_shell const * const C, int nshellC,
                             simint_shell const * const D, int nshellD,
                             double * const integrals, bool normalize)
{
    int inorm = (normalize ? 1 : 0);

    const int am1 = A[0].am;
    const int am2 = B[0].am;
    const int am3 = C[0].am;
    const int am4 = D[0].am;

    const int n234 = nshellB * nshellC * nshellD * NCART(am1) * NCART(am2) * NCART(am3) * NCART(am4);

    for(int i = 0; i < nshellA; i++)
    {
        const int idxstart = i * n234;
        int idx = 0;

        for(int j = 0; j < nshellB; j++)
        for(int k = 0; k < nshellC; k++)
        for(int l = 0; l < nshellD; l++)
        {
            long double vA[3] = { A[i].x, A[i].y, A[i].z };
            long double vB[3] = { B[j].x, B[j].y, B[j].z };
            long double vC[3] = { C[k].x, C[k].y, C[k].z };
            long double vD[3] = { D[l].x, D[l].y, D[l].z };

            std::array<int, 3> g1 = {{am1, 0, 0}};
            do
            {
                std::array<int, 3> g2 = {{am2, 0, 0}};
                do
                {
                    std::array<int, 3> g3 = {{am3, 0, 0}};
                    do
                    {
                        std::array<int, 3> g4 = {{am4, 0, 0}};
                        do
                        {
                            {
                                // the 12 total components
                                double myints[12] = {0.0};
                                
                                for(int m = 0; m < A[i].nprim; m++)
                                for(int n = 0; n < B[j].nprim; n++)
                                for(int o = 0; o < C[k].nprim; o++)
                                for(int p = 0; p < D[l].nprim; p++)
                                {
                                    // total prefactor
                                    const double prefac = A[i].coef[m] * B[j].coef[n] * C[k].coef[o] * D[l].coef[p];

                                    // temporary gaussians
                                    std::array<int, 3> gp;
                                    std::array<int, 3> gm;
                                    double valp, valm;

                                    // First center
                                    for(int d = 0; d < 3; d++)
                                    {
                                        gp = g1;  gp[d]++;
                                        gm = g1;  gm[d]--;
                                        valp = 0.0;
                                        valm = 0.0;

                                        valp = (double)ValeevRef_eri(gp[0], gp[1], gp[2], A[i].alpha[m], vA,
                                                                     g2[0], g2[1], g2[2], B[j].alpha[n], vB,
                                                                     g3[0], g3[1], g3[2], C[k].alpha[o], vC,
                                                                     g4[0], g4[1], g4[2], D[l].alpha[p], vD, inorm);

                                        if(ValidGaussian(gm))
                                            valm = (double)ValeevRef_eri(gm[0], gm[1], gm[2], A[i].alpha[m], vA,
                                                                         g2[0], g2[1], g2[2], B[j].alpha[n], vB,
                                                                         g3[0], g3[1], g3[2], C[k].alpha[o], vC,
                                                                         g4[0], g4[1], g4[2], D[l].alpha[p], vD, inorm);

                                        valp *= 2 * A[i].alpha[m];
                                        valm *= g1[d];  // ok, may be zero
                    
                                        myints[0*3+d] += (valp - valm) * prefac;
                                    }


                                    // Second Center
                                    for(int d = 0; d < 3; d++)
                                    {
                                        gp = g2;  gp[d]++;
                                        gm = g2;  gm[d]--;
                                        valp = 0.0;
                                        valm = 0.0;

                                        valp = (double)ValeevRef_eri(g1[0], g1[1], g1[2], A[i].alpha[m], vA,
                                                                     gp[0], gp[1], gp[2], B[j].alpha[n], vB,
                                                                     g3[0], g3[1], g3[2], C[k].alpha[o], vC,
                                                                     g4[0], g4[1], g4[2], D[l].alpha[p], vD, inorm);
                                        if(ValidGaussian(gm))
                                            valm = (double)ValeevRef_eri(g1[0], g1[1], g1[2], A[i].alpha[m], vA,
                                                                         gm[0], gm[1], gm[2], B[j].alpha[n], vB,
                                                                         g3[0], g3[1], g3[2], C[k].alpha[o], vC,
                                                                         g4[0], g4[1], g4[2], D[l].alpha[p], vD, inorm);

                                        valp *= 2 * B[j].alpha[n];
                                        valm *= g2[d];  // ok, may be zero
                    
                                        myints[1*3+d] += (valp - valm) * prefac;
                                    }


                                    // Third Center
                                    for(int d = 0; d < 3; d++)
                                    {
                                        gp = g3;  gp[d]++;
                                        gm = g3;  gm[d]--;
                                        valp = 0.0;
                                        valm = 0.0;

                                        valp = (double)ValeevRef_eri(g1[0], g1[1], g1[2], A[i].alpha[m], vA,
                                                                     g2[0], g2[1], g2[2], B[j].alpha[n], vB,
                                                                     gp[0], gp[1], gp[2], C[k].alpha[o], vC,
                                                                     g4[0], g4[1], g4[2], D[l].alpha[p], vD, inorm);
                                        if(ValidGaussian(gm))
                                            valm = (double)ValeevRef_eri(g1[0], g1[1], g1[2], A[i].alpha[m], vA,
                                                                         g2[0], g2[1], g2[2], B[j].alpha[n], vB,
                                                                         gm[0], gm[1], gm[2], C[k].alpha[o], vC,
                                                                         g4[0], g4[1], g4[2], D[l].alpha[p], vD, inorm);

                                        valp *= 2 * C[k].alpha[o];
                                        valm *= g3[d];  // ok, may be zero
                    
                                        myints[2*3+d] += (valp - valm) * prefac;
                                    }

                                    // The fourth center can be determined from the other three
                                    myints[3*3+0] = -(myints[0*3+0] + myints[1*3+0] + myints[2*3+0]);
                                    myints[3*3+1] = -(myints[0*3+1] + myints[1*3+1] + myints[2*3+1]);
                                    myints[3*3+2] = -(myints[0*3+2] + myints[1*3+2] + myints[2*3+2]);
                                }

                                // copy to the final integrals
                                std::copy(myints, myints+12, integrals + idxstart + idx);
                                idx += 12;

                            }
                        } while(IterateGaussian(g4));
                    } while(IterateGaussian(g3));
                } while(IterateGaussian(g2));
            } while(IterateGaussian(g1));
        }
    }
}



static
void ValeevRef_NumDeriv1(simint_shell const * const A, int nshellA,
                         simint_shell const * const B, int nshellB,
                         simint_shell const * const C, int nshellC,
                         simint_shell const * const D, int nshellD,
                         double * const integrals, bool normalize)
{
    const double delta = 1e-8;
    const int ncart_a = NCART(A->am);
    const int ncart_b = NCART(B->am);
    const int ncart_c = NCART(C->am);
    const int ncart_d = NCART(D->am);
    const int ncart_abcd = ncart_a * ncart_b * ncart_c * ncart_d;
    double buf_p[ncart_abcd];
    double buf_m[ncart_abcd];

    simint_shell tmp_p, tmp_m;
    simint_initialize_shell(&tmp_p);
    simint_initialize_shell(&tmp_m);


    double * integrals_ptr = integrals;

    for(int i = 0; i < nshellA; i++)
    for(int j = 0; j < nshellB; j++)
    for(int k = 0; k < nshellC; k++)
    for(int l = 0; l < nshellD; l++)
    {
        // Ax
        copy_gaussian_move(A+i, &tmp_p, 0,  delta);
        copy_gaussian_move(A+i, &tmp_m, 0, -delta);
        ValeevRef_Deriv0(&tmp_p, 1, B+j, 1, C+k, 1, D+l, 1, buf_p, normalize);
        ValeevRef_Deriv0(&tmp_m, 1, B+j, 1, C+k, 1, D+l, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 0*ncart_abcd, delta);

        // Ay
        copy_gaussian_move(A+i, &tmp_p, 1,  delta);
        copy_gaussian_move(A+i, &tmp_m, 1, -delta);
        ValeevRef_Deriv0(&tmp_p, 1, B+j, 1, C+k, 1, D+l, 1, buf_p, normalize);
        ValeevRef_Deriv0(&tmp_m, 1, B+j, 1, C+k, 1, D+l, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 1*ncart_abcd, delta);
                         
        // Az
        copy_gaussian_move(A+i, &tmp_p, 2,  delta);
        copy_gaussian_move(A+i, &tmp_m, 2, -delta);
        ValeevRef_Deriv0(&tmp_p, 1, B+j, 1, C+k, 1, D+l, 1, buf_p, normalize);
        ValeevRef_Deriv0(&tmp_m, 1, B+j, 1, C+k, 1, D+l, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 2*ncart_abcd, delta);

        // Bx
        copy_gaussian_move(B+j, &tmp_p, 0,  delta);
        copy_gaussian_move(B+j, &tmp_m, 0, -delta);
        ValeevRef_Deriv0(A+i, 1, &tmp_p, 1, C+k, 1, D+l, 1, buf_p, normalize);
        ValeevRef_Deriv0(A+i, 1, &tmp_m, 1, C+k, 1, D+l, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 3*ncart_abcd, delta);

        // By
        copy_gaussian_move(B+j, &tmp_p, 1,  delta);
        copy_gaussian_move(B+j, &tmp_m, 1, -delta);
        ValeevRef_Deriv0(A+i, 1, &tmp_p, 1, C+k, 1, D+l, 1, buf_p, normalize);
        ValeevRef_Deriv0(A+i, 1, &tmp_m, 1, C+k, 1, D+l, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 4*ncart_abcd, delta);

        // Bz
        copy_gaussian_move(B+j, &tmp_p, 2,  delta);
        copy_gaussian_move(B+j, &tmp_m, 2, -delta);
        ValeevRef_Deriv0(A+i, 1, &tmp_p, 1, C+k, 1, D+l, 1, buf_p, normalize);
        ValeevRef_Deriv0(A+i, 1, &tmp_m, 1, C+k, 1, D+l, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 5*ncart_abcd, delta);

        // Cx
        copy_gaussian_move(C+k, &tmp_p, 0,  delta);
        copy_gaussian_move(C+k, &tmp_m, 0, -delta);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, &tmp_p, 1, D+l, 1, buf_p, normalize);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, &tmp_m, 1, D+l, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 6*ncart_abcd, delta);

        // Cy
        copy_gaussian_move(C+k, &tmp_p, 1,  delta);
        copy_gaussian_move(C+k, &tmp_m, 1, -delta);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, &tmp_p, 1, D+l, 1, buf_p, normalize);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, &tmp_m, 1, D+l, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 7*ncart_abcd, delta);

        // Cz
        copy_gaussian_move(C+k, &tmp_p, 2,  delta);
        copy_gaussian_move(C+k, &tmp_m, 2, -delta);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, &tmp_p, 1, D+l, 1, buf_p, normalize);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, &tmp_m, 1, D+l, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 8*ncart_abcd, delta);

        // Dx
        copy_gaussian_move(D+l, &tmp_p, 0,  delta);
        copy_gaussian_move(D+l, &tmp_m, 0, -delta);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, C+k, 1, &tmp_p, 1, buf_p, normalize);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, C+k, 1, &tmp_m, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 9*ncart_abcd, delta);

        // Dy
        copy_gaussian_move(D+l, &tmp_p, 1,  delta);
        copy_gaussian_move(D+l, &tmp_m, 1, -delta);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, C+k, 1, &tmp_p, 1, buf_p, normalize);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, C+k, 1, &tmp_m, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 10*ncart_abcd, delta);

        // Dz
        copy_gaussian_move(D+l, &tmp_p, 2,  delta);
        copy_gaussian_move(D+l, &tmp_m, 2, -delta);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, C+k, 1, &tmp_p, 1, buf_p, normalize);
        ValeevRef_Deriv0(A+i, 1, B+j, 1, C+k, 1, &tmp_m, 1, buf_m, normalize);
        FormDeriv(ncart_abcd, buf_p, buf_m, integrals_ptr + 11*ncart_abcd, delta);
        
        integrals_ptr += ncart_abcd * 12;
    }

    simint_free_shell(&tmp_p);
    simint_free_shell(&tmp_m);

}

void ValeevRef_Integrals(simint_shell const * const A, int nshellA,
                         simint_shell const * const B, int nshellB,
                         simint_shell const * const C, int nshellC,
                         simint_shell const * const D, int nshellD,
                         double * const integrals, int deriv, bool normalize)
{
    if(deriv == 0)
        ValeevRef_Deriv0(A, nshellA, B, nshellB, C, nshellC, D, nshellD,
                         integrals, normalize);
    else if(deriv == 1)
        ValeevRef_Deriv1(A, nshellA, B, nshellB, C, nshellC, D, nshellD,
                         integrals, normalize);
    else if(deriv == -1)
        ValeevRef_NumDeriv1(A, nshellA, B, nshellB, C, nshellC, D, nshellD,
                            integrals, normalize);
    else
        throw std::runtime_error("Invlalid derivative for ValeevRef");
}

