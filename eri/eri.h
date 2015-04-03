#ifndef ERI_H
#define ERI_H

#include "shell.h"

int eri_ssss(struct multishell_pair const P,
             struct multishell_pair const Q,
             double * const restrict integrals,
             double * const restrict integralwork1,
             double * const restrict integralwork2);

int eri_ssss_cheby(struct multishell_pair const P,
                   struct multishell_pair const Q,
                   double * const restrict integrals,
                   double * const restrict integralwork1,
                   double * const restrict integralwork2);

int eri_ssss_combined(struct multishell_pair const P,
                      struct multishell_pair const Q,
                      double * const restrict integrals,
                      double * const restrict integralwork1,
                      double * const restrict integralwork2);

int eri_ssss_taylorcombined(struct multishell_pair const P,
                            struct multishell_pair const Q,
                            double * const restrict integrals,
                            double * const restrict integralwork1,
                            double * const restrict integralwork2);

int eri_ssss_taylor(struct multishell_pair const P,
                    struct multishell_pair const Q,
                    double * const restrict integrals,
                    double * const restrict integralwork1,
                    double * const restrict integralwork2);

int eri_ssss_split(struct multishell_pair const P,
                   struct multishell_pair const Q,
                   double * const restrict integrals,
                   double * const integralwork1,
                   double * const integralwork2);
#endif
