#include "simint/ostei/ostei.h"
#include "simint/ostei/ostei_config.h"

// This is the actual storage for this array
#define AMSIZE   SIMINT_OSTEI_MAXAM+1
#define DERSIZE  SIMINT_OSTEI_MAXDER+1

simint_osteifunc simint_osteifunc_array[DERSIZE][AMSIZE][AMSIZE][AMSIZE][AMSIZE];


int simint_compute_ostei(struct simint_multi_shellpair const * P,
                         struct simint_multi_shellpair const * Q,
                         double screen_tol,
                         double * restrict work,
                         double * restrict integrals)
{
    // don't forget that we don't include the square root in the screen values
    // stored in the shell pair
    double screen_tol2 = screen_tol * screen_tol;
    if(screen_tol > 0.0 && (P->screen_max * Q->screen_max) < screen_tol2 )
        return -1;

    return simint_osteifunc_array[0][P->am1][P->am2][Q->am1][Q->am2](*P, *Q,
                                                screen_tol2, work, integrals);
}


int simint_compute_ostei_deriv(int deriv,
                               struct simint_multi_shellpair const * P,
                               struct simint_multi_shellpair const * Q,
                               double screen_tol,
                               double * restrict work,
                               double * restrict integrals)
{
    // don't forget that we don't include the square root in the screen values
    // stored in the shell pair
    double screen_tol2 = screen_tol * screen_tol;
    if(screen_tol > 0.0 && (P->screen_max * Q->screen_max) < screen_tol2 )
        return -1;

    return simint_osteifunc_array[deriv][P->am1][P->am2][Q->am1][Q->am2](*P, *Q,
                                                  screen_tol2, work, integrals);
}

