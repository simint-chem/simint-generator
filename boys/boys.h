#ifndef BOYS_H
#define BOYS_H

#include <string.h> // for memcpy
#include "vectorization.h"
#include "boys/boys_param.h"

#define ASM __asm

// For various Boys_F functions
//extern double const * const restrict * const restrict boys_grid;
extern double const * const restrict * const restrict boys_grid;
extern const double boys_grid_max_x;
extern const int boys_grid_max_n;

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
inline cheby_idx(double d)
{
#ifdef SAFE_CHEBYIDX
    unsigned int r = 0;
    unsigned int i = (d+1);

    while (i >>= 1)
        r++;
    return r;
#else
    unsigned int i = (d+1.0);
    unsigned int l = 0;
    ASM("bsrl %1,%0" : "=r"(l) : "r"(i));
    return l;
#endif
}




// This includes F0_KFAC
inline double Boys_F0_taylor(double x)
{
    ASSUME_ALIGN(boys_grid);

    const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(x+BOYS_SHORTGRID_LOOKUPFAC2));
    const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
    const double dx = xi-x;   // -delta x

    double const * const restrict gridpts = &(boys_grid[lookup_idx][0]);

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
}

// Values include F0_KFAC
inline double Boys_F0_cheby(double x)
{
    double polynomial[BOYS_CHEBY_ORDER+1];

    const unsigned int idx = cheby_idx(x);

    double bshift = (int)( ((1<<idx)-1) + ((1<<(idx+1))-1) );
    x -= (bshift * 0.5);
    
    memcpy(polynomial, boys_chebygrid_F0[idx], (BOYS_CHEBY_ORDER+1)*sizeof(double));

    double F = polynomial[0];
    double xpow = x;
    for(int i = 1; i < BOYS_CHEBY_ORDER+1; i++)
    {
        F += polynomial[i]*xpow;
        xpow *= x;
    }

    return F;
 
}


// Values from this are missing F0_KFAC
inline double Boys_F0_erf(double x)
{
    const double x2 = sqrt(x);
    return erf(x2) / x2;   
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
