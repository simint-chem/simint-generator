#ifndef ERI_H
#define ERI_H

#include "shell.h"

int create_ss_shell_pair(struct gaussian_shell const A,
                         struct gaussian_shell const B,
                         struct shell_pair * restrict P);

int eri__ssss(struct gaussian_shell const A,
              struct gaussian_shell const B,
              struct gaussian_shell const C,
              struct gaussian_shell const D,
              struct shell_pair * restrict P_tmp,
              struct shell_pair * restrict Q_tmp,
              double * restrict integrals);

int eri_0pair_ssss(struct gaussian_shell const A,
                   struct gaussian_shell const B,
                   struct gaussian_shell const C,
                   struct gaussian_shell const D,
                   double * restrict integrals);

int eri_1pair_ssss(struct gaussian_shell const A,
                   struct gaussian_shell const B,
                   struct shell_pair const Q,
                   double * restrict integrals);

int eri_2pair_ssss(struct shell_pair const P,
                   struct shell_pair const Q,
                   double * restrict integrals);


#endif
