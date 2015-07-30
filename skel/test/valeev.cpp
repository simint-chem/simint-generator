
#include <math.h>
#include <stdlib.h>

#include "test/valeev.hpp"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static long double * fac;
static long double * df;
static long double ** bc;


void Valeev_Init(void)
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

void Valeev_Finalize(void)
{
    int i = 0;
    free(fac);
    free(df);

    for(i = 0; i < MAXFAC; i++)
        free(bc[i]);
    free(bc);
}

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


void Valeev_F(long double *F, int n, long double x)
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

long double Valeev_eri(int l1, int m1, int n1, long double alpha1,
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
