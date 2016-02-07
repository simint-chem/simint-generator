#ifndef BOYS_CHEBY_H
#define BOYS_CHEBY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "simint/boys/boys_chebygrid.h"

//extern double const * const restrict * const restrict boys_chebygrid_F0; // size [BOYS_CHEBY_NBIN][BOYS_CHEBY_ORDER+1]
extern double const boys_chebygrid_F0[BOYS_CHEBY_NBIN][BOYS_CHEBY_ORDER+1];

#ifdef __cplusplus
}
#endif

#endif
