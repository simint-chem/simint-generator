#ifndef SIMINT_ERD_INTERFACE_H
#define SIMINT_ERD_INTERFACE_H

#include "eri/shell.h"

void ERD_Init(int am1, int nprim1, int ncgto1,
              int am2, int nprim2, int ncgto2,
              int am3, int nprim3, int ncgto3,
              int am4, int nprim4, int ncgto4);

void ERD_Init(int na, struct gaussian_shell const * const restrict A,
              int nb, struct gaussian_shell const * const restrict B,
              int nc, struct gaussian_shell const * const restrict C,
              int nd, struct gaussian_shell const * const restrict D);

void ERD_Compute_shell(struct gaussian_shell const A,
                       struct gaussian_shell const B,
                       struct gaussian_shell const C,
                       struct gaussian_shell const D,
                       double * integrals);

void ERD_Finalize(void);

#endif
