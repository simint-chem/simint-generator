#ifndef BOYS_H
#define BOYS_H

#include <string.h> // for memcpy
#include "vectorization.h"
#include "boys/boys_param.h"

#define ASM __asm

// For various Boys_F functions
extern double const * const restrict * const restrict boys_grid;
extern const double boys_grid_max_x;
extern const int boys_grid_max_n;

extern double const boys_longfac[BOYS_LONGFAC_MAXN];

//extern double const * const restrict * const restrict boys_chebygrid_F0; // size [BOYS_CHEBY_NBIN][BOYS_CHEBY_ORDER+1]
extern double const boys_chebygrid_F0[BOYS_CHEBY_NBIN][BOYS_CHEBY_ORDER+1];

// Prototypes
void Boys_Init(double max_x, int max_n);
void Boys_Finalize(void);

double Boys_Max(const int ncenter,
                const double * X, const double * Y, const double * Z,
                const int * n_prim_per_center, const double * alpha);


/////////////////////////
// Inline functions
/////////////////////////

// Calculates the index for cheby
inline unsigned int cheby_idx(double d)
{
#ifdef ASM_CHEBYIDX
    unsigned int i = (d+1.0);
    unsigned int l = 0;
    ASM("bsrl %1,%0" : "=r"(l) : "r"(i));
    return l;
#else
    unsigned int r = 0;
    unsigned int i = (d+1);

    while (i >>= 1)
        r++;
    return r;
#endif
}



inline double Boys_F0_taylor(double x)
{
    ASSUME_ALIGN(boys_grid);

    const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(x+BOYS_SHORTGRID_LOOKUPFAC2));
    const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
    const double dx = xi-x;   // -delta x

    double const * const restrict gridpts = &(boys_grid[lookup_idx][0]);

/*
    const double f0xi = gridpts[0];
    const double f1xi = gridpts[1];
    const double f2xi = gridpts[2];
    const double f3xi = gridpts[3];
    const double f4xi = gridpts[4];
    const double f5xi = gridpts[5];
    const double f6xi = gridpts[6];
    const double f7xi = gridpts[7];

    return f0xi
           + dx * (                  f1xi
           + dx * ( (1.0/2.0   )   * f2xi
           + dx * ( (1.0/6.0   )   * f3xi
           + dx * ( (1.0/24.0  )   * f4xi
           + dx * ( (1.0/120.0 )   * f5xi
           + dx * ( (1.0/720.0 )   * f6xi
           + dx * ( (1.0/5040.0)   * f7xi
           )))))));
*/
    const double f0xi = gridpts[0];
    const double f1xi = gridpts[1];
    const double f2xi = gridpts[2];
    const double f3xi = gridpts[3];
    const double f4xi = gridpts[4];
    const double f5xi = gridpts[5];
    const double f6xi = gridpts[6];

    return f0xi
           + dx * (                  f1xi
           + dx * ( (1.0/2.0   )   * f2xi
           + dx * ( (1.0/6.0   )   * f3xi
           + dx * ( (1.0/24.0  )   * f4xi
           + dx * ( (1.0/120.0 )   * f5xi
           + dx * ( (1.0/720.0 )   * f6xi
           ))))));
}

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


// Values from this are missing F0_KFAC
inline double Boys_F0_erf(double x)
{
    const double x2 = sqrt(x);
    return erf(x2) / x2;   
}


// fast sqrt on the interval [0,1]
inline double Boys_fast_sqrt(double x)
{
    double s = 0.5;

    // may be able to get away with 11?
    for(int i = 0; i < 12; i++)
        s = 0.5 * (s + (x/s));
    return s;
}

inline double Boys_F0_FO(double x)
{
    const double num =        1.000000000000000e+00
                      + x * ( 0.310455397905288e+00
                      + x * ( 0.928325083325019e-01
                      + x * ( 0.164499706024780e-01 
                      + x * ( 0.227657096317876e-02
                      + x * ( 0.249999914708882e-03
                      + x * ( 0.219655434719760e-04
                      + x * ( 0.117515516961103e-05 
                      + x * ( 0.114088446739842e-06 )))))))); 

    const double den =        1.000000000000000e+00
                      + x * ( 0.977122076203642e+00
                      + x * ( 0.433135963320387e+00
                      + x * ( 0.115500167344045e+00
                      + x * ( 0.210603755673049e-01
                      + x * ( 0.289562352101294e-02
                      + x * ( 0.318354593651793e-03
                      + x * ( 0.279670555233004e-04
                      + x * ( 0.149625501987285e-05
                      + x * ( 0.145261921986579e-06 ))))))))); 

    const double frac = num/den;
    return sqrt(frac);
}

inline double Boys_F0_FO2(double x)
{
    const double num =        1.000000000000000e+00
                      + x * ( 3.09893858976444114e-01
                      + x * ( 9.26949023938073990e-02
                      + x * ( 1.64090380400552030e-02
                      + x * ( 2.27022457430278929e-03
                      + x * ( 2.49480075487682406e-04
                      + x * ( 2.17903871761419160e-05
                      + x * ( 1.19111485941561993e-06
                      + x * ( 1.11405368915781790e-07 )))))))); 

    const double den =        1.000000000000000e+00
                      + x * ( 9.76560532709625506e-01 
                      + x * ( 4.32624051808767562e-01
                      + x * ( 1.15292412795871560e-01
                      + x * ( 2.10111743064058006e-02
                      + x * ( 2.88739013575455309e-03
                      + x * ( 3.17697268656422444e-04
                      + x * ( 2.77439702781294230e-05
                      + x * ( 1.51657593378322051e-06
                      + x * ( 1.41845721199382853e-07 ))))))))); 

    const double frac = num/den;
    return sqrt(frac);
}



inline double Boys_F0_long(double x)
{
    return boys_longfac[0] / sqrt(x);
}

inline void Boys_F_taylor(double * const restrict F, int n, double x)
{
    ASSUME_ALIGN(boys_grid);

    const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(x+BOYS_SHORTGRID_LOOKUPFAC2));
    const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
    const double dx = xi-x;   // -delta x

    double const * const restrict gridpts = &(boys_grid[lookup_idx][0]);

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


#endif
