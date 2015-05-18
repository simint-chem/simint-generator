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

int eri_FOcombined_ppss(struct multishell_pair const P,
                        struct multishell_pair const Q,
                        double * const restrict S_1_1_0_0);

int eri_FOcombined_psps(struct multishell_pair const P,
                        struct multishell_pair const Q,
                        double * const restrict S_1_0_1_0);

int eri_FOcombined_ppps(struct multishell_pair const P,
                        struct multishell_pair const Q,
                        double * const restrict S_1_1_1_0);



#ifdef __cplusplus
}
#endif

#endif
