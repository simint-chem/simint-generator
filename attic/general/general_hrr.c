#include "simint/eri/general_eri.h"
#include "simint/vectorization/vectorization.h"
#include "simint/shell/shell.h"

#include "simint/eri/recur_lookup.h"


void general_hrr_j(int i, int j, int k, int l,
                   __m256d AB[3],
                   __m256d * theta1, __m256d * theta2,
                   __m256d * output)
{
    const int ncart_i = ((i+1)*(i+2))/2;
    const int ncart_j = ((j+1)*(j+2))/2;
    const int ncart_k = ((k+1)*(k+2))/2;
    const int ncart_l = ((l+1)*(l+2))/2;

    const int ncart_kl = ncart_k * ncart_l;

    const int ncart_j_m1 = ((j)*(j+1))/2;

    const int arrstart = am_recur_map[j];
    struct RecurInfo const * aminfo = &recurinfo_array[arrstart]; 

    int outidx = 0;
    for(int16_t idx_i = 0; idx_i < ncart_i; idx_i++)
    for(int16_t idx_j = 0; idx_j < ncart_j; idx_j++)
    {
        struct RecurInfo const * rj = aminfo + idx_j;
        struct RecurInfo const * ri = aminfo + idx_i;
        const int dir = rj->dir;

        for(int16_t idx_kl = 0; idx_kl < ncart_kl; idx_kl++)
        {
            const int idx1 = (ri->idx[dir][2]) * ncart_j_m1 * ncart_kl
                           + (rj->idx[dir][0]) * ncart_kl
                           + idx_kl;

            const int idx2 = idx_i * ncart_j_m1 * ncart_kl
                           + (rj->idx[dir][0]) * ncart_kl
                           + idx_kl;
                           

            output[outidx] = theta1[idx1] + AB[dir]*theta2[idx2];
            outidx++;
        }
    }
}

