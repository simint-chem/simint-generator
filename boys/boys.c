#include <math.h>
#include <string.h> // for memcpy

#include "vectorization.h"
#include "boys/boys.h"
#include "boys/boys_longfac.h"

extern const double boys_shortgrid[BOYS_SHORTGRID_NPOINT][BOYS_SHORTGRID_MAXN + 1];
extern const double boys_longfac[BOYS_LONGFAC_MAXN];

double * boys_grid_flat;
double ** boys_grid;

double boys_grid_max_x = 0;
int boys_grid_max_n = 0;

#define M_PI 3.14159265358979323846

void Boys_Init(double max_x, int max_n)
{
    // max_x should fit in the 15 digits of precision. It should
    // usually be less than 1e6 or so. If not, you've got some problems.

    // + 1 for fencepost
    // + 2 for rounding error + safety
    int nelements_x = 3 + (int)((max_x / BOYS_GRID_SPACE));
    int nelements_n = SIMD_ROUND_DBL(max_n+1);
    boys_grid_flat = ALLOC(nelements_x * nelements_n * sizeof(double));

    // set up the pointers
    boys_grid = ALLOC(nelements_x * sizeof(double*));

    int ii = 0;
    for(int i = 0; i < nelements_x; ++i, ii += nelements_n)
        boys_grid[i] = boys_grid_flat + ii;
    
    // copy over the shortgrid elements for all n
    // note that we may be requesting a shorter range on x than is provided
    int shortend = (max_x > BOYS_SHORTGRID_MAXX ? BOYS_SHORTGRID_NPOINT : nelements_x);
    for(int i = 0; i < shortend; ++i)
        memcpy(boys_grid[i], boys_shortgrid[i], max_n * sizeof(double));

    // calculate the rest with the asymptotic formula
    // (won't run if nelements_x < BOYS_SHORTGRID_NPOINT)
    double x = BOYS_SHORTGRID_MAXX + BOYS_GRID_SPACE;

    // vectorize the outer loop
    #ifdef SIMINT_SIMD
    #pragma simd
    #endif
    for(int i = BOYS_SHORTGRID_NPOINT; i < nelements_x; ++i)
    {
        double powx = 1.0/sqrt(x);
        for(int n = 0; n < max_n; ++n)
        {
            boys_grid[i][n] = boys_longfac[n] * powx;
            powx /= x;
        }
        x += BOYS_GRID_SPACE;
    }

    boys_grid_max_x = max_x;
    boys_grid_max_n = max_n;
}

void Boys_Finalize(void)
{
    FREE(boys_grid_flat);
    FREE(boys_grid);
}


void Boys_F(double * const restrict F, int n, double x)
{
    ASSUME_ALIGN(boys_grid);

    const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(x+BOYS_SHORTGRID_LOOKUPFAC2));
    const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
    const double dx = xi-x;   // -delta x

    double * gridpts = &(boys_grid[lookup_idx][0]);

    for(int i = 0; i <= n; ++i)
    {
        const double f0xi = gridpts[i];
        const double f1xi = gridpts[i+1];
        const double f2xi = gridpts[i+2];
        const double f3xi = gridpts[i+3];
        const double f4xi = gridpts[i+4];
        const double f5xi = gridpts[i+5];
        const double f6xi = gridpts[i+6];
        const double f7xi = gridpts[i+7];

        F[i] = f0xi
               + dx * (                  f1xi
               + dx * ( (1.0/2.0   )   * f2xi
               + dx * ( (1.0/6.0   )   * f3xi
               + dx * ( (1.0/24.0  )   * f4xi
               + dx * ( (1.0/120.0 )   * f5xi
               + dx * ( (1.0/720.0 )   * f6xi
               + dx * ( (1.0/5040.0)   * f7xi
               )))))));
    }
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
