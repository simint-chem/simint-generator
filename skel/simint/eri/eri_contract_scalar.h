#pragma once

inline void contract(int ncart,
                     double const * restrict PRIM_INT,
                     double * const restrict PRIM_PTR)
{
    for(int np = 0; np < ncart; ++np)
        PRIM_PTR[np] += PRIM_INT[np]; 
}

