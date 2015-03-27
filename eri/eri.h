#ifndef ERI_H
#define ERI_H

#include "shell.h"

int eri_ssss(struct shell_pair const P,
             struct shell_pair const Q,
             double * const restrict integrals,
             double * const restrict integralwork1,
             double * const restrict integralwork2);

int eri_ssss_cheby(struct shell_pair const P,
                   struct shell_pair const Q,
                   double * const restrict integrals,
                   double * const restrict integralwork1,
                   double * const restrict integralwork2);

int eri_ssss_combined(struct shell_pair const P,
                      struct shell_pair const Q,
                      double * const restrict integrals,
                      double * const restrict integralwork1,
                      double * const restrict integralwork2);

int eri_ssss_taylorcombined(struct shell_pair const P,
                            struct shell_pair const Q,
                            double * const restrict integrals,
                            double * const restrict integralwork1,
                            double * const restrict integralwork2);

int eri_ssss_taylor(struct shell_pair const P,
                    struct shell_pair const Q,
                    double * const restrict integrals,
                    double * const restrict integralwork1,
                    double * const restrict integralwork2);

int eri_ssss_split(struct shell_pair const P,
                   struct shell_pair const Q,
                   double * const restrict integrals,
                   double * const integralwork1,
                   double * const integralwork2);
#endif
