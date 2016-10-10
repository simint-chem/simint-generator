#pragma once

#ifdef __cplusplus
extern "C" {
#endif

static inline
void contract(int ncart,
              double const * restrict PRIM_INT,
              double * restrict PRIM_PTR)
{
    for(int np = 0; np < ncart; ++np)
        PRIM_PTR[np] += PRIM_INT[np]; 
}

#ifdef __cplusplus
}
#endif


