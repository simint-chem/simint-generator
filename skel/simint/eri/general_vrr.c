#include "simint/eri/general_eri.h"
#include "simint/vectorization/vectorization.h"
#include "simint/shell/shell.h"

#include "simint/eri/recur_lookup.h"


void general_vrr1(int i, int N,
                  __m256d one_over_2p, __m256d a_over_p,
                  __m256d PA[3], __m256d PQ[3],
                  __m256d * theta1, __m256d * theta2,
                  __m256d * output)
{
    // in this case, we want (i, 0, 0, 0), so i is equal
    // to the angular momentum.

    const int ncart = ((i+1)*(i+2))/2;
    const int ncart1 = ((i)*(i+1))/2;
    const int ncart2 = ((i-1)*(i))/2;

    const int arrstart = am_recur_map[i];
    struct RecurInfo const * aminfo = &recurinfo_array[arrstart]; 

    int curidx = 0;
    for(int16_t idx_i = 0; idx_i < ncart; idx_i++)
    {
        struct RecurInfo const * ri = aminfo + idx_i;
        int dir = ri->dir;

        for(int n = 0; n <= N; n++)
        {
            const int outidx = n * ncart + curidx;
            const int idx1 = n*(ncart1) + ri->idx[dir][0];
            const int idx2 = (n+1)*(ncart1) + ri->idx[dir][0];

            output[outidx] = PA[dir]*theta1[idx1] + a_over_p*PQ[dir]*theta1[idx2];

            if(ri->ijk[dir] > 1)
            {
                const int idx3 = n*(ncart2) + ri->idx[dir][1];
                const int idx4 = (n+1)*(ncart2) + ri->idx[dir][1];
                const __m256d i1 = _mm256_set1_pd(ri->ijk[dir]-1);
                output[outidx] += one_over_2p * i1 * (theta2[idx3] + a_over_p * theta2[idx4]);
            }
        }
        curidx++;
    }
}


void general_vrr2(int i, int k, int N,
                  __m256d one_over_2q, __m256d one_over_2pq,
                  __m256d a_over_q,
                  __m256d QC[3], __m256d PQ[3],
                  __m256d * theta1, __m256d * theta2, __m256d * theta3,
                  __m256d * output)
{
    const int ncart_i = ((i+1)*(i+2))/2;
    const int ncart_k = ((k+1)*(k+2))/2;
    const int ncart_i_1 = ((i)*(i+1))/2;
    const int ncart_k_1 = ((k)*(k+1))/2;
    const int ncart_k_2 = ((k-1)*(k))/2;

    const int ncart = ncart_i * ncart_k;
    const int ncart1 = ncart_i * ncart_k_1;
    const int ncart2 = ncart_i * ncart_k_2;
    const int ncart3 = ncart_i_1 * ncart_k_1;

    // we are mainly interested in k-1, etc. But we also need to know
    // some info about i-1 
    const int arrstart_i = am_recur_map[i];
    const int arrstart_k = am_recur_map[k];
    struct RecurInfo const * aminfo_i = &recurinfo_array[arrstart_i]; 
    struct RecurInfo const * aminfo_k = &recurinfo_array[arrstart_k]; 

    int curidx = 0;
    for(int16_t idx_i = 0; idx_i < ncart_i; idx_i++)
    for(int16_t idx_k = 0; idx_k < ncart_k; idx_k++)
    {
        struct RecurInfo const * rk = aminfo_k + idx_k;
        struct RecurInfo const * ri = aminfo_i + idx_i;
        int dir = rk->dir;

        for(int n = 0; n <= N; n++)
        {
            const int outidx = n * ncart + curidx;
            const int idx1 =     n*(ncart1) + idx_i * ncart_k_1 + rk->idx[dir][0];
            const int idx2 = (n+1)*(ncart1) + idx_i * ncart_k_1 + rk->idx[dir][0];

            output[outidx] = QC[dir]*theta1[idx1] - a_over_q*PQ[dir]*theta1[idx2];

            if(ri->idx[dir][0] > -1) // -1 signifies that it doesn't exist
            {
                const int idx5 = (n+1)*(ncart3) + ri->idx[dir][0] * ncart_k_1 + rk->idx[dir][0];
                const __m256d i1 = _mm256_set1_pd(ri->ijk[dir]);
                output[outidx] += i1 * one_over_2pq * theta3[idx5];
            }

            if(rk->ijk[dir] > 1)
            {
                const int idx3 = n*(ncart2) + idx_i * ncart_k_2 + rk->idx[dir][1];
                const int idx4 = (n+1)*(ncart2) + idx_i * ncart_k_2 + rk->idx[dir][1];
                const __m256d k1 = _mm256_set1_pd(rk->ijk[dir]-1);
                output[outidx] += one_over_2q * k1 * (theta2[idx3] + a_over_q * theta2[idx4]);
            }
        }
        curidx++;
    }
}
