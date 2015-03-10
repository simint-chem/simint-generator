#ifndef ERI_H
#define ERI_H

#include "shell.h"

int eri_1pair_ssss_single(struct gaussian_shell const A,
                          struct gaussian_shell const B,
                          struct shell_pair const Q,
                          double * const restrict integrals);

int eri_1pair_ssss_multi(int na, struct gaussian_shell const * const restrict A,
                         int nb, struct gaussian_shell const * const restrict B,
                         struct shell_pair const Q,
                         double * const restrict integrals,
                         double * const restrict integralwork);

int eri_2pair_ssss_single(struct shell_pair const P,
                          struct shell_pair const Q,
                          double * const restrict integrals);

int eri_2pair_ssss_multi(struct shell_pair const P,
                         struct shell_pair const Q,
                         double * const restrict integrals,
                         double * const restrict integralwork);

#endif
