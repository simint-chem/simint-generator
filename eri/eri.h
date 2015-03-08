#ifndef ERI_H
#define ERI_H

#include "shell.h"

int create_shell_pair(struct gaussian_shell const * restrict A,
                      struct gaussian_shell const * restrict B,
                      struct shell_pair * restrict P);


int eri_ssss(struct shell_pair const * restrict P,
             struct shell_pair const * restrict Q,
             double * restrict integrals);

#endif
