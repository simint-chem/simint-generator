#pragma once

#include "simint/boys/boys_taylor.h"
#include "simint/boys/boys_long.h"

#ifdef __cplusplus
extern "C" {
#endif

void boys_F_split(double * restrict Farr, int n, double const * restrict xvec);


#ifdef __cplusplus
}
#endif

