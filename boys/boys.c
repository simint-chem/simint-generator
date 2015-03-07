#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "boys/boys.h"
#include "boys/boys_grid.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern double boys_grid[BOYS_GRID_NPOINT][BOYS_GRID_MAXN + 1];

double dfac[BOYS_MAX_DFAC]; // (2n-1)!!

void Boys_Init(void)
{
    int i;

    dfac[0] = 1.0;   // (-1)!!
    dfac[1] = 1.0;   // ( 1)!!
    dfac[2] = 3.0;   // ( 3)!!
    for(i = 3; i < BOYS_MAX_DFAC; ++i)
        dfac[i] = (2*i-1)*dfac[i-1];
    
}

void Boys_Finalize(void)
{
}



void Boys_F(double * F, int n, double x)
{
    if(x < BOYS_GRID_MAXX)
      Boys_F_short(F, n, x);
    else
      Boys_F_long(F, n, x);
}



#ifdef SIMINT_SIMD
  #pragma omp declare simd simdlen(SIMD_LEN)
#endif
void Boys_F_short(double * F, int n, double x)
{
    int i;
    // rounding, assuming x is positive
    double fac = 1.0;
    const double ex = exp(-x);
    const double x2 = 2 * x;

    // The boys_grid is stored with x as the first index
    // therefore, the n index is contiguous
    // This lookup should be the only part that can't be vectorized
    const int idx = (int)(BOYS_GRID_LOOKUPFAC*(x+1.0/(BOYS_GRID_LOOKUPFAC*2.0))); 
    const double dx = (((double)idx / BOYS_GRID_LOOKUPFAC)-x);
    double const * fval = boys_grid[idx] + n;

    F[n] = 0;
    for(i = 0; i <= BOYS_INTERP_ORDER; i++)
    {
        F[n] += (*fval)*pow(dx, i)/fac;
        fac *= (i+1);
        ++fval;
    }

    // recurse down
    for(i = n-1; i >= 0; --i)
        F[i] = (x2* F[i+1] + ex)/(2*i + 1);
}



#ifdef SIMINT_SIMD
  #pragma omp declare simd simdlen(SIMD_LEN)
#endif
void Boys_F_long(double * F, int n, double x)
{
    const double ex = exp(-x);
    const double x2 = 2 * x;
    int i;

    F[n] = dfac[n] / (1<<(n+1)) * sqrt(M_PI / pow(x, 2*n+1));

    // recurse down
    for(i = n-1; i >= 0; --i)
        F[i] = (x2 * F[i+1] + ex)/(2*i + 1);
}


// Maximum value needed for Boys function
double Boys_Max(const int ncenter,
                const double * X, const double * Y, const double * Z,
                const int * n_prim_per_center, const double * alpha)
{
    int i, j, a, b;

    /*
    int maxpa = 0;
    int maxpb = 0;
    */

    double maxd = 0;
    double d, d2;

    int pa = 0; // counts over alpha
    int pb = 0; // counts over alpha

    // these hold the start of this center in the alpha array
    int paorig = 0;
    int pborig = 0;

    for(i = 0; i < ncenter; i++)
    {
        pborig = 0;
        for(j = 0; j < ncenter; j++)
        {
            d = (X[i]-X[j])*(X[i]-X[j])
              + (Y[i]-Y[j])*(Y[i]-Y[j])
              + (Z[i]-Z[j])*(Z[i]-Z[j]);

            pa = paorig;
            for(a = 0; a < n_prim_per_center[i]; a++)
            {
                pb = pborig;
                for(b = 0; b < n_prim_per_center[j]; b++)
                {

                    d2 = d * (alpha[pa]*alpha[pb]) / (alpha[pa] + alpha[pb]);
                    if(d2 > maxd)
                    {
                        /*
                        maxpa = pa;
                        maxpb = pb;
                        */
                        maxd = d2;
                    }
                    pb++;
                }
                pa++;
            }
            pborig += n_prim_per_center[j];
        }
        paorig += n_prim_per_center[i];
    }

    //printf("Max i, j: %d %d\n", maxpa, maxpb);

    // missing a factor of 2
    return 2.0 * maxd;;
}
