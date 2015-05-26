#include <math.h>
#include <stdlib.h>

#define EPS 1e-17
#define MAXFAC 200

static double * df;

void Boys_F_VRef_Init(void)
{
    int i = 0;

    df = malloc(2 * MAXFAC * sizeof(double));
    df[0] = 1.0;
    df[1] = 1.0;
    df[2] = 1.0;
    for(i = 3; i < MAXFAC*2; i++)
       df[i] = (i - 1) * df[i - 2];
}

void Boys_F_VRef_Finalize(void)
{
    free(df);
}

void Boys_F_VRef(double *F, int n, double x)
{
    int i, m;
    int m2;
    double t2;
    double num;
    double sum;
    double term1;
    const double K = 1.0 / 1.12837916709551257390;
    double et;


    if (x > 20.0)   /* For big t's do upward recursion */
    {
        t2 = 2 * x;
        et = exp(-x);
        x = sqrt(x);
        F[0] = K * erf(x) / x;
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
        et = exp(-x);
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
        while (fabs(term1) > EPS && i < MAXFAC);
        F[n] = sum * et;
        for (m = n - 1; m >= 0; m--)   /* And then do downward recursion */
        {
            F[m] = (t2 * F[m + 1] + et) / (2 * m + 1);
        }
    }
}
