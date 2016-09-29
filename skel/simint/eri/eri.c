#include "simint/eri/eri.h"

#include <string.h> // for memset

// Stores pointers to the eri functions
extern simint_erifunc_sharedwork  simint_erifunc_sharedwork_array[SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1];
extern simint_erifunc             simint_erifunc_array[SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1][SIMINT_ERI_MAXAM+1];


int simint_compute_eri(struct simint_multi_shellpair const * P,
                       struct simint_multi_shellpair const * Q,
                       double screen_tol,
                       double * restrict integrals)
{
    // don't forget that we don't include the square root in the screen values
    // stored in the shell pair
    double screen_tol2 = screen_tol * screen_tol;
    if(screen_tol > 0.0 && (P->screen_max * Q->screen_max) < screen_tol2 )
        return -1;

    return simint_erifunc_array[P->am1][P->am2][Q->am1][Q->am2](*P, *Q, screen_tol2, integrals);
}



/*! \brief Compute an eri given shell pair information
 *
 * \param [in] P The shell pairs for the bra side of the integral 
 * \param [in] Q The shell pairs for the ket side of the integral
 * \param [in] work Workspace to use in calculating the integrals
 * \param [inout] integrals Storage for the final integrals. Since size information
 *                          is not passed, you are expected to ensure that this buffer
 *                          is large enough
 */
int simint_compute_eri_sharedwork(struct simint_multi_shellpair const * P,
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

    return simint_erifunc_sharedwork_array[P->am1][P->am2][Q->am1][Q->am2](*P, *Q, screen_tol2, work, integrals);
}

