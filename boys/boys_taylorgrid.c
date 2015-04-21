#include <string.h> // for memcpy

#include "vectorization.h"
#include "boys/boys_shortgrid.h"
#include "boys/boys_long.h"

extern const double boys_shortgrid[BOYS_SHORTGRID_NPOINT][BOYS_SHORTGRID_MAXN + 1];

double * boys_grid_flat;
double ** boys_grid;

double boys_grid_max_x = 0;
int boys_grid_max_n = 0;


void Boys_taylorgrid_Init(double max_x, int max_n)
{
    // max_x should fit in the 15 digits of precision. It should
    // usually be less than 1e6 or so. If not, you've got some problems.

    // + 1 for fencepost
    // + 2 for rounding error + safety
    int nelements_x = 3 + (int)(max_x / BOYS_SHORTGRID_SPACE);
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
    double x = BOYS_SHORTGRID_MAXX + BOYS_SHORTGRID_SPACE;
    for(int i = BOYS_SHORTGRID_NPOINT; i < nelements_x; ++i)
    {
        Boys_F_long(boys_grid[i], max_n, x);
        x += BOYS_SHORTGRID_SPACE;
    }

    boys_grid_max_x = max_x;
    boys_grid_max_n = max_n;
}

void Boys_taylorgrid_Finalize(void)
{
    // deallocates data
    FREE(boys_grid_flat);

    // deallocates pointers
    FREE(boys_grid);
}

