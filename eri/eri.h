#ifndef ERI_H
#define ERI_H

#include "shell.h"

int eri_erf_ssss(struct shell_pair const P,
                 struct shell_pair const Q,
                 double * const restrict integrals,
                 double * const restrict integralwork1,
                 double * const restrict integralwork2);

int eri_erf_multi_ssss(struct multishell_pair const P,
                       struct multishell_pair const Q,
                       double * const restrict integrals,
                       double * const restrict integralwork1,
                       double * const restrict integralwork2);

int eri_erf_split_ssss(struct multishell_pair const P,
                       struct multishell_pair const Q,
                       double * const restrict integrals,
                       double * const integralwork1,
                       double * const integralwork2);


int eri_cheby_ssss(struct multishell_pair const P,
                   struct multishell_pair const Q,
                   double * const restrict integrals,
                   double * const restrict integralwork1,
                   double * const restrict integralwork2);


int eri_erf_combined_ssss(struct multishell_pair const P,
                          struct multishell_pair const Q,
                          double * const restrict integrals,
                          double * const restrict integralwork1,
                          double * const restrict integralwork2);


int eri_taylor_ssss(struct multishell_pair const P,
                    struct multishell_pair const Q,
                    double * const restrict integrals,
                    double * const restrict integralwork1,
                    double * const restrict integralwork2);


int eri_taylorcombined_ssss(struct multishell_pair const P,
                            struct multishell_pair const Q,
                            double * const restrict integrals,
                            double * const restrict integralwork1,
                            double * const restrict integralwork2);


int eri_taylorcombined_psss(struct multishell_pair const P,
                            struct multishell_pair const Q,
                            double * const restrict integrals,
                            double * const restrict integralwork1,
                            double * const restrict integralwork2);


#endif
