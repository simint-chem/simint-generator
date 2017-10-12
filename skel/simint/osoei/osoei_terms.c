#include <math.h>
#include "simint/constants.h"
#include "simint/osoei/osoei.h"

// Get a value of S_IJ
#define S_IJ(i,j) (s_ij[((i)*(nam2) + j)])

void simint_osoei_overlap_terms(const double alpha1, const double * xyz1,
                                const double alpha2, const double * xyz2,
                                int nam1, int nam2,
                                double * restrict terms)
{
    const double AB[3] = { xyz1[0] - xyz2[0], xyz1[1] - xyz2[1], xyz1[2] - xyz2[2] };
    const double AB2[3] = { AB[0]*AB[0], AB[1]*AB[1], AB[2]*AB[2] };


    /////////////////////////////////////////////////////////
    // General notes about the following
    //
    // This is the OS algorithm taken pretty much verbatim
    // from Helgaker, Jorgensen, & Olsen. For each x,y,z
    // direction, we construct an array of Sij of length nam1*nam2.
    // These are placed in the pre-allocated workspace.
    //
    // The indexing of these arrays is pretty straightforward.
    // Sij = ptr[i*nam2+j]. This is in the S_IJ macro, hopefully
    // to make the code clearer.
    /////////////////////////////////////////////////////////
    const double a1xyz[3] = { alpha1*xyz1[0], alpha1*xyz1[1], alpha1*xyz1[2] };
    const double oop = 1.0/(alpha1 + alpha2); // = 1/p = 1/(a1 + a2)
    const double mu = alpha1*alpha2*oop;      // (a1*a2)/(a1+a2)

    const double oo2p = 0.5*oop;
    const double a2xyz[3] = { alpha2*xyz2[0], alpha2*xyz2[1], alpha2*xyz2[2] };

    const double P[3] = { (a1xyz[0]+a2xyz[0])*oop,
                          (a1xyz[1]+a2xyz[1])*oop,
                          (a1xyz[2]+a2xyz[2])*oop };

    const double PA[3] = { P[0] - xyz1[0], P[1] - xyz1[1], P[2] - xyz1[2] };
    const double PB[3] = { P[0] - xyz2[0], P[1] - xyz2[1], P[2] - xyz2[2] };


    // three cartesian directions
    for(int d = 0; d < 3; d++)
    {
        // the workspace for this direction
        // THIS IS THEN ACCESSED THROUGH THE S_IJ MACRO
        double * const s_ij = terms + d*nam1*nam2;

        S_IJ(0,0) = sqrt(PI * oop) * exp(-mu * AB2[d]);

        // do j = 0 for all remaining i
        for(int i = 1; i < nam1; i++)
        {
            S_IJ(i,0) = PA[d]*S_IJ(i-1,0);
            if(i > 1)
                S_IJ(i,0) += (i-1)*oo2p*S_IJ(i-2,0);
        }


        // now do i = 0 for all remaining j
        for(int j = 1; j < nam2; j++)
        {
            S_IJ(0,j) = PB[d]*S_IJ(0,j-1);
            if(j > 1)
                S_IJ(0,j) += (j-1)*oo2p*S_IJ(0,j-2);
        }

        // now all the rest
        for(int i = 1; i < nam1; i++)
        for(int j = 1; j < nam2; j++)
        {
            S_IJ(i,j) = PB[d]*S_IJ(i,j-1) + oo2p*i*S_IJ(i-1,j-1);
            if(j > 1)
                S_IJ(i,j) += oo2p*(j-1)*S_IJ(i,j-2);
        }
    }
}
