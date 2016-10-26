#pragma once

#include "simint/shell/shell.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Calculate the screening value for a shell pair */
double
simint_shellscreen(struct simint_shell const * A,
                   struct simint_shell const * B,
                   int screen_method);


/*! \brief Calculate primitive screening information for a shell pair
 *
 * The screening algorithm depends on the \p method paramter.
 */
double
simint_primscreen(struct simint_shell const * A,
                  struct simint_shell const * B,
                  double * out,
                  int screen_method);


/*! \brief Calculate the screening value for a shell pair */
double
simint_shellscreen_schwarz(struct simint_shell const * A,
                           struct simint_shell const * B);


/*! \brief Calculate the screening value for a shell pair */
double
simint_shellscreen_fastschwarz(struct simint_shell const * A,
                               struct simint_shell const * B);


/*! \brief Calculate primitive screening information for a shell pair
 *
 * Calculates screening info via Schwarz screening
 */
double
simint_primscreen_schwarz(struct simint_shell const * A,
                          struct simint_shell const * B,
                          double * out);


/*! \brief Calculate primitive screening information for a shell pair
 *
 * Similar to simint_primscreen_schwarz, however forces the angular
 * momentum to zero.  (00|00) >= (ab|cd) if all other parameters
 * are the same.
 */
double
simint_primscreen_fastschwarz(struct simint_shell const * A,
                              struct simint_shell const * B,
                              double * out);



#ifdef __cplusplus
}
#endif

