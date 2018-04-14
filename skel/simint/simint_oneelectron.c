#include "simint/simint_oneelectron.h"
#include "simint/osoei/osoei.h"

int simint_compute_overlap(struct simint_shell * a,
                           struct simint_shell * b, 
                           double * restrict integrals)
{
    return simint_compute_osoei_overlap(a, b, integrals);
}


int simint_compute_ke(struct simint_shell * a,
                           struct simint_shell * b, 
                           double * restrict integrals)
{
    return simint_compute_osoei_ke(a, b, integrals);
}


int simint_compute_potential(int ncenter,
                             double * Z, double * x, double * y, double * z,
                             struct simint_shell const * sh1,
                             struct simint_shell const * sh2,
                             double * restrict integrals)
{
    return simint_compute_osoei_potential(ncenter, Z, x, y, z, sh1, sh2, integrals);
}
