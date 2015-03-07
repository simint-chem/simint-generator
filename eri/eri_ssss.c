#include <math.h>

#include "boys/boys.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Computes (A B | C D) when A,B,C,D are all s-type orbitals */
int eri_ssss(int n,
             double * Ax, double * Ay, double * Az,
             double * Bx, double * By, double * Bz,
             double * Cx, double * Cy, double * Cz,
             double * Dx, double * Dy, double * Dz,
             double * Aexp, double * Bexp, double * Cexp, double * Dexp,
             double * res)
{
    int i;
    double tmp, tmp2;
    double p_ab, p_cd;
    double R2;
    double F0;

    double pfac;

    #ifdef SIMINT_SIMD
      #pragma simd
    #endif
    for(i = 0; i < n; i++)
    {
        p_ab = Aexp[i] + Bexp[i];
        p_cd = Cexp[i] + Dexp[i];

        // This accumulates the normalization and prefactors
        pfac = 2 * pow(M_PI, 2.5) / (p_ab * p_cd * sqrt(p_ab + p_cd));

        // Rewritten to use only one 2 tmp variable
        // do Xab = (Xab_x **2 + Xab_y ** 2 + Xab_z **2)
        // and apply to the prefactor
        tmp2 = 0;
        tmp =  Ax[i] - Bx[i];  tmp2 += tmp * tmp;
        tmp =  Ay[i] - By[i];  tmp2 += tmp * tmp;
        tmp =  Az[i] - Bz[i];  tmp2 += tmp * tmp;
        pfac *= exp(-((Aexp[i] * Bexp[i]) / (p_ab)) * tmp2);
    
        // and similarly for Xcd
        tmp2 = 0;
        tmp =  Cx[i] - Dx[i];  tmp2 += tmp * tmp;
        tmp =  Cy[i] - Dy[i];  tmp2 += tmp * tmp;
        tmp =  Cz[i] - Dz[i];  tmp2 += tmp * tmp;
        pfac *= exp(-((Cexp[i] * Dexp[i]) / (p_cd)) * tmp2);

        /*
        Pab_x = (Aexp[i]*Ax[i] + Bexp[i]*Bx[i])/p_ab;
        Pab_y = (Aexp[i]*Ay[i] + Bexp[i]*By[i])/p_ab;
        Pab_z = (Aexp[i]*Az[i] + Bexp[i]*Bz[i])/p_ab;
        Pcd_x = (Cexp[i]*Cx[i] + Dexp[i]*Dx[i])/p_cd;
        Pcd_y = (Cexp[i]*Cy[i] + Dexp[i]*Dy[i])/p_cd;
        Pcd_z = (Cexp[i]*Cz[i] + Dexp[i]*Dz[i])/p_cd;

        // now the two electron integral
        Rpq_x = Pab_x - Pcd_x;
        Rpq_y = Pab_y - Pcd_y;
        Rpq_z = Pab_z - Pcd_z;
        R2 = Rpq_x*Rpq_x + Rpq_y*Rpq_y + Rpq_z*Rpq_z;
        */
        /* Same as above, but no more extra variables */
        R2 = 0;
        tmp = (Aexp[i]*Ax[i] + Bexp[i]*Bx[i])/p_ab
              - (Cexp[i]*Cx[i] + Dexp[i]*Dx[i])/p_cd;
        R2 += tmp * tmp;

        tmp = (Aexp[i]*Ay[i] + Bexp[i]*By[i])/p_ab
              - (Cexp[i]*Cy[i] + Dexp[i]*Dy[i])/p_cd;
        R2 += tmp * tmp;

        tmp = (Aexp[i]*Az[i] + Bexp[i]*Bz[i])/p_ab
              - (Cexp[i]*Cz[i] + Dexp[i]*Dz[i])/p_cd;
        R2 += tmp * tmp;

        //alpha = (p_ab * p_cd)/(p_ab + p_cd);
        //Valeev_F(&F0, 0, alpha * R2);
        Boys_F_long(&F0, 0, R2 * (p_ab * p_cd)/(p_ab + p_cd));
 
        pfac *= pow(2.0 * Aexp[i] / M_PI, 0.75);
        pfac *= pow(2.0 * Bexp[i] / M_PI, 0.75);
        pfac *= pow(2.0 * Cexp[i] / M_PI, 0.75);
        pfac *= pow(2.0 * Dexp[i] / M_PI, 0.75);

        res[i] = pfac * F0;
    }

    return n;
}
