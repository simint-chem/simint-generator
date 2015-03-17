#ifndef BOYS_H
#define BOYS_H

#include "vectorization.h"
#include "boys_shortgrid.h"

// These don't only apply to the short-range grid
// so sopy the definitions
#define BOYS_GRID_MAXN       BOYS_SHORTGRID_MAXN
#define BOYS_GRID_SPACE      BOYS_SHORTGRID_SPACE
#define BOYS_GRID_LOOKUPFAC  BOYS_SHORTGRID_LOOKUPFAC
#define BOYS_GRID_LOOKUPFAC2 BOYS_SHORTGRID_LOOKUPFAC2


#define BOYS_INTERP_ORDER 7


void Boys_Init(double max_x, int max_n);
void Boys_Finalize(void);

void Boys_F(double* const restrict F, int n, double x);

double Boys_Max(const int ncenter,
                const double * X, const double * Y, const double * Z,
                const int * n_prim_per_center, const double * alpha);


#endif
