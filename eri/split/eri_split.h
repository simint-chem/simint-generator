#ifndef ERI_SPLIT_H
#define ERI_SPLIT_H

#ifdef __cplusplus
extern "C" {
#endif


#include "eri/shell.h"


#define ERI_SPLIT_MAXAM 2


int eri_split_s_s_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_0_0_0_0);
int eri_split_p_s_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_1_0_0_0);
int eri_split_p_s_p_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_1_0_1_0);
int eri_split_p_p_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_1_1_0_0);
int eri_split_p_p_p_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_1_1_1_0);
int eri_split_p_p_p_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_1_1_1_1);
int eri_split_d_s_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_0_0_0);
int eri_split_d_s_p_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_0_1_0);
int eri_split_d_s_p_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_0_1_1);
int eri_split_d_s_d_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_0_2_0);
int eri_split_d_p_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_1_0_0);
int eri_split_d_p_p_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_1_1_0);
int eri_split_d_p_p_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_1_1_1);
int eri_split_d_p_d_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_1_2_0);
int eri_split_d_p_d_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_1_2_1);
int eri_split_d_d_s_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_2_0_0);
int eri_split_d_d_p_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_2_1_0);
int eri_split_d_d_p_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_2_1_1);
int eri_split_d_d_d_s(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_2_2_0);
int eri_split_d_d_d_p(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_2_2_1);
int eri_split_d_d_d_d(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict S_2_2_2_2);


#ifdef __cplusplus
}
#endif

#endif
