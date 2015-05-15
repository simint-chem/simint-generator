#ifndef ERI_FOCOMBINED_H
#define ERI_FOCOMBINED_H

#ifdef __cplusplus
extern "C" {
#endif

#include "eri/shell.h"

int eri_FOcombined_ssss(struct multishell_pair const P,
                        struct multishell_pair const Q,
                        double * const restrict S_0_0_0_0);

int eri_FOcombined_psss(struct multishell_pair const P,
                        struct multishell_pair const Q,
                        double * const restrict S_1_0_0_0);

int eri_FOcombined_spss(struct multishell_pair const P,
                        struct multishell_pair const Q,
                        double * const restrict S_0_1_0_0);

int eri_FOcombined_ssps(struct multishell_pair const P,
                        struct multishell_pair const Q,
                        double * const restrict S_0_0_1_0);

int eri_FOcombined_sssp(struct multishell_pair const P,
                        struct multishell_pair const Q,
                        double * const restrict S_0_0_0_1);


#ifdef __cplusplus
}
#endif

#endif
