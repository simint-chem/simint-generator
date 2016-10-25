#include "simint/simint_eri.h"
#include "simint/ostei/ostei.h"
#include "simint/boys/boys.h"



int simint_compute_eri(struct simint_multi_shellpair const * P,
                       struct simint_multi_shellpair const * Q,
                       double screen_tol,
                       double * restrict integrals)
{
    return simint_compute_ostei(P, Q, boys_F_split, screen_tol, integrals);
}



int simint_compute_eri_sharedwork(struct simint_multi_shellpair const * P,
                                  struct simint_multi_shellpair const * Q,
                                  double screen_tol,
                                  double * restrict work,
                                  double * restrict integrals)
{
    return simint_compute_ostei_sharedwork(P, Q, boys_F_split, screen_tol, work, integrals);
}

