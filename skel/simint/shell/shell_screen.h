#pragma once

#include "simint/shell/shell.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Calculate the screening value for a shell pair */
double
simint_shellscreen_schwarz_max(struct simint_shell const * A,
                               struct simint_shell const * B);


/*! \brief Calculate primitive screening information for a shell pair */
double
simint_primscreen_schwarz_max(struct simint_shell const * A,
                              struct simint_shell const * B,
                              double * out);


#ifdef __cplusplus
}
#endif

