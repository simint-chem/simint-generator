#include <math.h>
#include <string.h>

#include "simint/boys/boys.h"
#include "simint/constants.h"
#include "simint/recur_lookup.h"
#include "simint/osoei/osoei.h"


int simint_compute_osoei_potential(int ncenter,
                                   double * Z, double * x, double * y, double * z,
                                   struct simint_shell const * sh1,
                                   struct simint_shell const * sh2,
                                   double * restrict integrals)
{
    const int am1 = sh1->am;
    const int am2 = sh2->am;
    const int am12 = am1 + am2;

    const int nam1 = am1 + 1;
    const int nam2 = am2 + 1;
    const int nam12 = nam1*nam2;

    // Workspace for calculating terms
    double work[6*nam12];

    const double xyz1[3] = { sh1->x, sh1->y, sh1->z };
    const double xyz2[3] = { sh2->x, sh2->y, sh2->z };

    const double AB[3] = { xyz1[0] - xyz2[0], xyz1[1] - xyz2[1], xyz1[2] - xyz2[2] };
    const double AB2 = AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2];

    const int ncart1 = NCART(am1);
    const int ncart2 = NCART(am2);
    const int ncart12 = ncart1*ncart2;

    //TODO - too big
    double amwork[nam1][nam2][NCART(am1+1)*NCART(am2+1)*(nam1+nam2+1)];

    for(int n = 0; n < ncenter; n++)
    {
        for(int a = 0; a < sh1->nprim; a++)
        {
            const double a1 = sh1->alpha[a];
            const double a1xyz[3] = { a1*xyz1[0], a1*xyz1[1], a1*xyz1[2] };

            for(int b = 0; b < sh2->nprim; b++)
            {
                const double a2 = sh2->alpha[b];
                const double p = a1 + a2;
                const double oop = 1.0/(a1 + a2); // = 1/p = 1/(a1 + a2)
                const double mu = a1*a2*oop; // (a1+a2)/(a1*a2)

                const double oo2p = 0.5*oop;
                const double a2xyz[3] = { a2*xyz2[0], a2*xyz2[1], a2*xyz2[2] };

                const double P[3] = { (a1xyz[0]+a2xyz[0])*oop,
                                      (a1xyz[1]+a2xyz[1])*oop,
                                      (a1xyz[2]+a2xyz[2])*oop };

                const double PA[3] = { P[0] - xyz1[0], P[1] - xyz1[1], P[2] - xyz1[2] };
                const double PB[3] = { P[0] - xyz2[0], P[1] - xyz2[1], P[2] - xyz2[2] };
                const double PC[3] = { P[0] - x[n], P[1] - y[n], P[2] - z[n] };
                const double PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];

                // boys function
                const double T = PC2 * p;

                // TODO - HACK
                SIMINT_DBLTYPE Ftmp[am12+1];
                SIMINT_DBLTYPE Ttmp = SIMINT_DBLSET1(T);
                boys_F_split(Ftmp, Ttmp, am12);
                
                for(int i = 0; i <= am12; i++)
                {
                    //TODO - hacky type pun
                    double * Fval = (double*)(&Ftmp[i]);
                    amwork[0][0][i] = Fval[0] * 2*PI*oop*exp(-mu*AB2);
                }

                // nested recurrence
                // we skip (0,0) since that is the boys function
                for(int i = 0; i <= am1; i++)
                {
                    const int arrstart1 = am_recur_map[sh1->am];
                    struct RecurInfo const * aminfo1 =  &recurinfo_array[arrstart1];

                    // number of cartesians in the previous two shells
                    const int incart   = NCART(i);
                    const int incart_1 = (i > 0) ? NCART(i-1) : 0;
                    const int incart_2 = (i > 1) ? NCART(i-2) : 0;

                    // only form if i != 0 (ie, don't do (0,0)
                    if(i > 0)
                    {
                        // form (i,0)
                        double * iwork = amwork[i][0];
                        double * iwork14 = amwork[i-1][0];                      // location of the 1st and 4th terms
                        double * iwork25 = (i > 1) ? amwork[i-2][0] : NULL;  // location of the 2nd and 5th terms

     
                        // maximum value of (m) to calculate
                        // we need [0, am12-i] inclusive
                        const int max_m = am12 - i; 
                        int idx = 0;
                        for(int m = 0; m <= max_m; m++)
                        {
                            // commonly used in dimensioning
                            const int offset_1 = m*incart_1;
                            const int offset_4 = offset_1 + incart_1;  // (m+1)*incart_1
                            const int offset_2 = m*incart_2;
                            const int offset_5 = offset_2 + incart_2;  // (m+1)*incart_1

                            for(int n = 0; n < incart; n++)
                            {
                                // get the recurrence information for this cartesian
                                const int8_t d = aminfo1[n].dir;
                                const int8_t i_ijk = aminfo1[n].ijk[d];
                                const int idx1 = offset_1 + aminfo1[n].idx[d][0]; // index for 1st term
                                const int idx4 = offset_4 + aminfo1[n].idx[d][0]; // index for 4th term
                                const int idx2 = offset_2 + aminfo1[n].idx[d][1]; // index for 2nd term
                                const int idx5 = offset_5 + aminfo1[n].idx[d][1]; // index for 5th term

                                iwork[idx] = PA[d]*iwork14[idx1] - PC[d]*iwork14[idx4];  // 1st and 4th terms

                                if(i_ijk > 1)
                                    iwork[idx] += oo2p*(i_ijk-1)*(iwork25[idx2] - iwork25[idx5]); // 2nd and 5th terms

                                idx++;
                            }
                        }
                    }

                    // now (i,j) via second vertical recurrence
                    for(int j = 1; j <= am2; j++)
                    {
                        const int arrstart2 = am_recur_map[sh2->am];
                        struct RecurInfo const * aminfo2 =  &recurinfo_array[arrstart2];

                        // number of cartesians in the previous two shells
                        const int jncart   = NCART(j);
                        const int jncart_1 = NCART(j-1);   // j can't be zero (loop starts at 1)
                        const int jncart_2 = (j > 1) ? NCART(j-2) : 0;

                        double * jwork = amwork[i][j];

                        double * jwork14 = amwork[i][j-1];                        // location of the 1st and 4th terms
                        double * jwork36 = (j > 1) ? amwork[i][j-2] : NULL;    // location of the 3rd and 6th terms
                        double * jwork25 = (i > 0) ? amwork[i-1][j-1] : NULL;  // location of the 2nd and 5th terms

                        const int max_m2 = am2 - j + 1; // need [0, am2-1] inclusive

                        int cartidx = 0; // index of the pair of cartesians
                        for(int m = 0; m <= max_m2; m++)
                        {
                            int cartidx_1 = 0;  // index of just the first cartesian
                            for(int n = 0; n < incart; n++)
                            {
                                // precompute some of the offsets
                                // storage is  m, cart1, cart2
                                // so total index would be (m*ncart1*ncart2 + cart1*ncart2 + cart2)
                                const int offset1 = jncart_1*(m*incart + cartidx_1);                     // m*incart*jncart_1 + cartidx_1*jncart_1 
                                const int offset4 = jncart_1*((m+1)*incart + cartidx_1);                 // (m+1)*incart*jncart_1 + cartidx_1*jncart_1
                                const int offset3 = jncart_2*(m*incart + cartidx_1);                     // m*incart*jncart_2 + cartidx_1*jncart_2
                                const int offset6 = jncart_2*((m+1)*incart + cartidx_1);                 // (m+1)*incart*jncart_2 + cartidx_1*jncart_2
                                const int offset2[3] = { jncart_1*(m*incart_1 + aminfo1[n].idx[0][0]),        // m*incart_1*jncart_1 + aminfo1[n].idx[0][0]*jncart_1 
                                                            jncart_1*(m*incart_1 + aminfo1[n].idx[1][0]),        // m*incart_1*jncart_1 + aminfo1[n].idx[1][0]*jncart_1 
                                                            jncart_1*(m*incart_1 + aminfo1[n].idx[2][0]) };      // m*incart_1*jncart_1 + aminfo1[n].idx[2][0]*jncart_1 
                                const int offset5[3] = { jncart_1*((m+1)*incart_1 + aminfo1[n].idx[0][0]),    // (m+1)*incart_1*jncart_1 + aminfo1[n].idx[0][0]*jncart_1 
                                                            jncart_1*((m+1)*incart_1 + aminfo1[n].idx[1][0]),    // (m+1)*incart_1*jncart_1 + aminfo1[n].idx[1][0]*jncart_1 
                                                            jncart_1*((m+1)*incart_1 + aminfo1[n].idx[2][0]) };  // (m+1)*incart_1*jncart_1 + aminfo1[n].idx[2][0]*jncart_1 

                                for(int o = 0; o < jncart; o++)
                                {
                                    const int8_t d = aminfo2[o].dir; // direction we should recurse
                                    const int8_t i_ijk = aminfo1[n].ijk[d];  // values of i and j in that direction
                                    const int8_t j_ijk = aminfo2[o].ijk[d];
                                    const int idx1 = offset1 + aminfo2[o].idx[d][0];     // 1st term
                                    const int idx4 = offset4 + aminfo2[o].idx[d][0];     // 4th term
                                    const int idx3 = offset3 + aminfo2[o].idx[d][1];     // 3rd term
                                    const int idx6 = offset6 + aminfo2[o].idx[d][1];     // 6th term
                                    const int idx2 = offset2[d] + aminfo2[o].idx[d][0];  // 2nd term
                                    const int idx5 = offset5[d] + aminfo2[o].idx[d][0];  // 5th term


                                    jwork[cartidx] = PB[d]*jwork14[idx1] - PC[d]*jwork14[idx4]; // terms 1 & 4

                                    if(i_ijk > 0)
                                        jwork[cartidx] += oo2p*(i_ijk)*jwork25[idx2] - oo2p*(i_ijk)*jwork25[idx5]; // terms 2 & 5

                                    if(j_ijk > 1)
                                        jwork[cartidx] += oo2p*(j_ijk-1)*(jwork36[idx3] - jwork36[idx6]); // terms 3 & 6

                                    cartidx++;
                                }

                                cartidx_1++;
                            }
                        } // end loop over m
                    } // end loop over j
                } // end loop over i

                // apply the coefficients
                double * amintegrals = amwork[sh1->am][sh2->am];
                for(int i = 0; i < ncart12; i++)
                    integrals[i] -= amintegrals[i] * Z[n] * sh1->coef[a] * sh2->coef[b];

            } // end loop over primitive a
        } // end loop over primitive b
    } // close loop over atoms

    return 1;
}

