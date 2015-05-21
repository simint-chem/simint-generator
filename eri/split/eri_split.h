#ifndef ERI_SPLIT_H
#define ERI_SPLIT_H

#ifdef __cplusplus
extern "C" {
#endif


#include "eri/shell.h"


#define ERI_SPLIT_MAXAM 2


int eri_split_s_s_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__s_s_s_s);
int eri_split_p_s_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__p_s_s_s);
int eri_split_p_s_p_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__p_s_p_s);
int eri_split_p_p_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__p_p_s_s);
int eri_split_p_p_p_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__p_p_p_s);
int eri_split_p_p_p_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__p_p_p_p);
int eri_split_d_s_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_s_s_s);
int eri_split_d_s_p_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_s_p_s);
int eri_split_d_s_p_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_s_p_p);
int eri_split_d_s_d_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_s_d_s);
int eri_split_d_p_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_p_s_s);
int eri_split_d_p_p_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_p_p_s);
int eri_split_d_p_p_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_p_p_p);
int eri_split_d_p_d_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_p_d_s);
int eri_split_d_p_d_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_p_d_p);
int eri_split_d_d_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_d_s_s);
int eri_split_d_d_p_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_d_p_s);
int eri_split_d_d_p_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_d_p_p);
int eri_split_d_d_d_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_d_d_s);
int eri_split_d_d_d_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_d_d_p);
int eri_split_d_d_d_d(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict INT__d_d_d_d);


#ifdef __cplusplus
}
#endif

#endif
