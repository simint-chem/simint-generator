#include "simint/simint_oneelectron.h"
#include "simint/osoei/osoei.h"

int simint_compute_overlap(struct simint_shell * a,
                           struct simint_shell * b, 
                           double * restrict integrals)
{
    return simint_compute_osoei_overlap(a, b, integrals);
}

