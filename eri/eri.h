#ifndef ERI_H
#define ERI_H

#ifdef __cplusplus
extern "C" {
#endif

#include "shell.h"

int eri_split_ssss(struct multishell_pair const P,
                   struct multishell_pair const Q,
                   double * const restrict integrals,
                   double * const integralwork1,
                   double * const integralwork2);

int eri_splitcombined_ssss(struct multishell_pair const P,
                           struct multishell_pair const Q,
                           double * const restrict integrals);

int eri_FOcombined_ssss(struct multishell_pair const P,
                        struct multishell_pair const Q,
                        double * const restrict integrals);

#ifdef __cplusplus
}
#endif

#endif
