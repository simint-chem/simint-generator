#ifndef ERI_H
#define ERI_H

#include "shell.h"

int eri_ssss_flat(struct shell_pair const P,
                  struct shell_pair const Q,
                  double * const restrict integrals,
                  double * const restrict integralwork1,
                  double * const restrict integralwork2);

#endif
