#include "simint/ostei/ostei_general.h"
#include "simint/ostei/recur_lookup.h"


void ostei_general_hrr_I(int i, int j, int k, int l,
                         const double AB[3],
                         const double * theta1,
                         const double * theta2,
                         double * output)
{
    const int ncart_i = ((i+1)*(i+2))/2;
    const int ncart_j = ((j+1)*(j+2))/2;
    const int ncart_k = ((k+1)*(k+2))/2;
    const int ncart_l = ((l+1)*(l+2))/2;

    const int ncart_kl = ncart_k * ncart_l;

    const int ncart_j_p1 = ((j+2)*(j+3))/2;

    const int arrstart_i = am_recur_map[i];
    const int arrstart_j = am_recur_map[j];
    struct RecurInfo const * aminfo_i = &recurinfo_array[arrstart_i];
    struct RecurInfo const * aminfo_j = &recurinfo_array[arrstart_j];

    int outidx = 0;
    for(int idx_i = 0; idx_i < ncart_i; idx_i++)
    {
        struct RecurInfo const * ri = aminfo_i + idx_i;

        for(int idx_j = 0; idx_j < ncart_j; idx_j++)
        {
            struct RecurInfo const * rj = aminfo_j + idx_j;
            const int dir = ri->dir;

            const int idx_i_m1 = ri->idx[dir][0];
            const int idx_j_p1 = rj->idx[dir][2];

            const int idx1 = idx_i_m1 * ncart_j_p1 * ncart_kl +
                             idx_j_p1 * ncart_kl;

            const int idx2 = idx_i_m1 * ncart_j    * ncart_kl +
                             idx_j    * ncart_kl;

            for(int idx_kl = 0; idx_kl < ncart_kl; idx_kl++)
            {
                const int idx1_full = idx1 + idx_kl;
                const int idx2_full = idx2 + idx_kl;
                output[outidx] = theta1[idx1_full] - AB[dir]*theta2[idx2_full];
                outidx++;
            }
        }
    }
}



void ostei_general_hrr_J(int i, int j, int k, int l,
                         const double AB[3],
                         const double * theta1,
                         const double * theta2,
                         double * output)
{
    const int ncart_i = ((i+1)*(i+2))/2;
    const int ncart_j = ((j+1)*(j+2))/2;
    const int ncart_k = ((k+1)*(k+2))/2;
    const int ncart_l = ((l+1)*(l+2))/2;

    const int ncart_kl = ncart_k * ncart_l;

    const int ncart_j_m1 = ((j)*(j+1))/2;

    const int arrstart_i = am_recur_map[i];
    const int arrstart_j = am_recur_map[j];
    struct RecurInfo const * aminfo_i = &recurinfo_array[arrstart_i];
    struct RecurInfo const * aminfo_j = &recurinfo_array[arrstart_j];

    int outidx = 0;
    for(int idx_i = 0; idx_i < ncart_i; idx_i++)
    {
        struct RecurInfo const * ri = aminfo_i + idx_i;

        for(int idx_j = 0; idx_j < ncart_j; idx_j++)
        {
            struct RecurInfo const * rj = aminfo_j + idx_j;
            const int dir = rj->dir;

            const int idx_i_p1 = ri->idx[dir][2];
            const int idx_j_m1 = rj->idx[dir][0];

            const int idx1 = idx_i_p1 * ncart_j_m1 * ncart_kl +
                             idx_j_m1 * ncart_kl;

            const int idx2 = idx_i    * ncart_j_m1 * ncart_kl +
                             idx_j_m1 * ncart_kl;

            for(int idx_kl = 0; idx_kl < ncart_kl; idx_kl++)
            {
                const int idx1_full = idx1 + idx_kl;
                const int idx2_full = idx2 + idx_kl;
                output[outidx] = theta1[idx1_full] + AB[dir]*theta2[idx2_full];
                outidx++;
            }
        }
    }
}


void ostei_general_hrr_K(int i, int j, int k, int l,
                         const double CD[3],
                         const double * theta1,
                         const double * theta2,
                         double * output)
{
    const int ncart_i = ((i+1)*(i+2))/2;
    const int ncart_j = ((j+1)*(j+2))/2;
    const int ncart_k = ((k+1)*(k+2))/2;
    const int ncart_l = ((l+1)*(l+2))/2;

    const int ncart_ij = ncart_i * ncart_j;

    const int ncart_k_m1 = ((k)*(k+1))/2;
    const int ncart_l_p1 = ((l+2)*(l+3))/2;

    const int arrstart_k = am_recur_map[k];
    const int arrstart_l = am_recur_map[l];
    struct RecurInfo const * aminfo_k = &recurinfo_array[arrstart_k];
    struct RecurInfo const * aminfo_l = &recurinfo_array[arrstart_l];

    int outidx = 0;
    for(int idx_ij = 0; idx_ij < ncart_ij; idx_ij++)
    {
        for(int idx_k = 0; idx_k < ncart_k; idx_k++)
        {
            struct RecurInfo const * rk = aminfo_k + idx_k;

            for(int idx_l = 0; idx_l < ncart_l; idx_l++)
            {
                struct RecurInfo const * rl = aminfo_l + idx_l;
                const int dir = rk->dir;
                const int idx_k_m1 = rk->idx[dir][0];
                const int idx_l_p1 = rl->idx[dir][2];

                const int idx1_full = idx_ij   * ncart_k_m1 * ncart_l_p1 +
                                      idx_k_m1 * ncart_l_p1 +
                                      idx_l_p1;

                const int idx2_full = idx_ij   * ncart_k_m1 * ncart_l +
                                      idx_k_m1 * ncart_l +
                                      idx_l;
 
                output[outidx] = theta1[idx1_full] - CD[dir]*theta2[idx2_full];
                outidx++;
            }
        }
    }
}

void ostei_general_hrr_L(int i, int j, int k, int l,
                         const double CD[3],
                         const double * theta1,
                         const double * theta2,
                         double * output)
{
    const int ncart_i = ((i+1)*(i+2))/2;
    const int ncart_j = ((j+1)*(j+2))/2;
    const int ncart_k = ((k+1)*(k+2))/2;
    const int ncart_l = ((l+1)*(l+2))/2;

    const int ncart_ij = ncart_i * ncart_j;

    const int ncart_k_p1 = ((k+2)*(k+3))/2;
    const int ncart_l_m1 = ((l)*(l+1))/2;

    const int arrstart_k = am_recur_map[k];
    const int arrstart_l = am_recur_map[l];
    struct RecurInfo const * aminfo_k = &recurinfo_array[arrstart_k];
    struct RecurInfo const * aminfo_l = &recurinfo_array[arrstart_l];

    int outidx = 0;
    for(int idx_ij = 0; idx_ij < ncart_ij; idx_ij++)
    {
        for(int idx_k = 0; idx_k < ncart_k; idx_k++)
        {
            struct RecurInfo const * rk = aminfo_k + idx_k;

            for(int idx_l = 0; idx_l < ncart_l; idx_l++)
            {
                struct RecurInfo const * rl = aminfo_l + idx_l;
                const int dir = rl->dir;
                const int idx_k_p1 = rk->idx[dir][2];
                const int idx_l_m1 = rl->idx[dir][0];

                const int idx1_full = idx_ij   * ncart_k_p1 * ncart_l_m1 +
                                      idx_k_p1 * ncart_l_m1 +
                                      idx_l_m1;

                const int idx2_full = idx_ij   * ncart_k  * ncart_l_m1 +
                                      idx_k    * ncart_l_m1 +
                                      idx_l_m1;
 
                output[outidx] = theta1[idx1_full] + CD[dir]*theta2[idx2_full];
                outidx++;
            }
        }
    }
}
