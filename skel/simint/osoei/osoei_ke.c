#include <math.h>
#include <string.h>
#include "simint/recur_lookup.h"
#include "simint/osoei/osoei.h"

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)
#define S_IJ(i,j) (s_ij[((i)*(nam2) + j)])
#define T_IJ(i,j) (t_ij[((i)*(nam2) + j)])

int simint_compute_osoei_ke(struct simint_shell const * sh1,
                            struct simint_shell const * sh2,
                            double * restrict integrals)
{
    const int am1 = sh1->am;
    const int am2 = sh2->am;

    const int nam1 = am1 + 1;
    const int nam2 = am2 + 1;
    const int nam12 = nam1*nam2;

    // Workspace for calculating terms
    double work[6*nam12];

    const double xyz1[3] = { sh1->x, sh1->y, sh1->z };
    const double xyz2[3] = { sh2->x, sh2->y, sh2->z };

    const int ncart1 = NCART(am1);
    const int ncart2 = NCART(am2);

    // Mostly needed just for the ordering
    const int arrstart1 = am_recur_map[sh1->am];
    const int arrstart2 = am_recur_map[sh2->am];
    struct RecurInfo const * aminfo1 =  &recurinfo_array[arrstart1];
    struct RecurInfo const * aminfo2 =  &recurinfo_array[arrstart2];

    memset(integrals, 0, ncart1*ncart2*sizeof(double));

    for(int i = 0; i < sh1->nprim; i++)
    {
        for(int j = 0; j < sh2->nprim; j++)
        {
            const double prefac = sh1->coef[i] * sh2->coef[j];
            simint_osoei_ke_terms(sh1->alpha[i], xyz1, sh2->alpha[j], xyz2,
                                  sh1->am+1, sh2->am+1, work);

            size_t outidx = 0;

            for(int n = 0; n < ncart1; n++)
            for(int m = 0; m < ncart2; m++)
            {
                const int8_t * ijk1 = aminfo1[n].ijk;
                const int8_t * ijk2 = aminfo2[m].ijk;
                const int xidx = ijk1[0]*nam2 + ijk2[0];
                const int yidx = ijk1[1]*nam2 + ijk2[1];
                const int zidx = ijk1[2]*nam2 + ijk2[2];

                const double val = work[3*nam12 + xidx]*work[1*nam12 + yidx]*work[2*nam12 + zidx]  // Tij*Skl*Smn
                                 + work[0*nam12 + xidx]*work[4*nam12 + yidx]*work[2*nam12 + zidx]  // Sij*Tkl*Smn
                                 + work[0*nam12 + xidx]*work[1*nam12 + yidx]*work[5*nam12 + zidx]; // Sij*Skl*Tmn;

                // remember: a and b are indices of primitives
                integrals[outidx++] += prefac * val;
            }
        }
    }

    return 1;
}

