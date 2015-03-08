#ifndef ERI_H
#define ERI_H

#include "shell.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


int create_shell_pair(const struct gaussian_shell * A,
                      const struct gaussian_shell * B,
                      struct shell_pair * P);


int eri_ssss(const struct shell_pair * P,
             const struct shell_pair * Q,
             double * integrals);

#endif
