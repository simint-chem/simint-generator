#pragma once

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Calculate screening information for a shell pair */
void
simint_primscreen_schwarz_max(struct simint_shell const * A,
                              struct simint_shell const * B,
                              double * out);


#ifdef __cplusplus
}
#endif

