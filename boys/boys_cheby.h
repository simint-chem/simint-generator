#ifndef BOYS_CHEBY_H
#define BOYS_CHEBY_H

#include "boys/boys_chebygrid.h"

//extern double const * const restrict * const restrict boys_chebygrid_F0; // size [BOYS_CHEBY_NBIN][BOYS_CHEBY_ORDER+1]
extern double const boys_chebygrid_F0[BOYS_CHEBY_NBIN][BOYS_CHEBY_ORDER+1];

inline double Boys_F0_cheby(double x)
{
    const unsigned int idx = cheby_idx(x);

    double bshift = (int)( ((1<<idx)-1) + ((1<<(idx+1))-1) );
    x -= (bshift * 0.5);
   
    return         boys_chebygrid_F0[idx][0]
           + x * ( boys_chebygrid_F0[idx][1]
           + x * ( boys_chebygrid_F0[idx][2]
           + x * ( boys_chebygrid_F0[idx][3]
           + x * ( boys_chebygrid_F0[idx][4]
           + x * ( boys_chebygrid_F0[idx][5]
           + x * ( boys_chebygrid_F0[idx][6]
           + x * ( boys_chebygrid_F0[idx][7]
           + x * ( boys_chebygrid_F0[idx][8]
           + x * ( boys_chebygrid_F0[idx][9]
           + x * ( boys_chebygrid_F0[idx][10]
           + x * ( boys_chebygrid_F0[idx][11]
           + x * ( boys_chebygrid_F0[idx][12]
           + x * ( boys_chebygrid_F0[idx][13]
           + x * ( boys_chebygrid_F0[idx][14]
           ))))))))))))));
}

#endif
