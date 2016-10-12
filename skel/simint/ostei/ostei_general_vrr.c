#include "simint/ostei/ostei_general.h"
#include "simint/ostei/recur_lookup.h"


/**********************************************************/
/* All vrr1 are the same, with some swapping of variables */
/**********************************************************/
void ostei_general_vrr1(int i, int num_n,
                        __m256d one_over_2p, __m256d a_over_p,
                        __m256d aop_PQ[3], __m256d PA[3],
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
        const int dir = ri->dir;

        for(int n = 0; n < num_n; n++)
        {
            const int outidx = n * ncart + curidx;
            const int idx1 = n*(ncart1) + ri->idx[dir][0];
            const int idx2 = idx1 + ncart1;

            output[outidx] = PA[dir]*theta1[idx1] + aop_PQ[dir]*theta1[idx2];

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



void ostei_general_vrr2_i(int i, int j, int k, int l, int num_n,
                          __m256d one_over_2p, __m256d a_over_p,
                          __m256d one_over_2pq,
                          __m256d aop_PQ[3], __m256d PA[3],
                          __m256d * theta1, __m256d * theta2,
                          __m256d * theta3, __m256d * theta4,
                          __m256d * theta5, __m256d * output)
{
    const int ncart_i = ((i+1)*(i+2))/2;
    const int ncart_j = ((j+1)*(j+2))/2;
    const int ncart_k = ((k+1)*(k+2))/2;
    const int ncart_l = ((l+1)*(l+2))/2;

    const int ncart_i_1 = ((i)*(i+1))/2;
    const int ncart_j_1 = ((j)*(j+1))/2;
    const int ncart_k_1 = ((k)*(k+1))/2;
    const int ncart_l_1 = ((l)*(l+1))/2;

    const int ncart_i_2 = ((i-1)*(i))/2;

    const int ncart = ncart_i * ncart_j * ncart_k * ncart_l;
    const int ncart1 = ncart_i_1 * ncart_j   * ncart_k   * ncart_l;
    const int ncart2 = ncart_i_2 * ncart_j   * ncart_k   * ncart_l;
    const int ncart3 = ncart_i_1 * ncart_j_1 * ncart_k   * ncart_l;
    const int ncart4 = ncart_i_1 * ncart_j   * ncart_k_1 * ncart_l;
    const int ncart5 = ncart_i_1 * ncart_j   * ncart_k   * ncart_l_1;

    const int arrstart_i = am_recur_map[i];
    const int arrstart_j = am_recur_map[j];
    const int arrstart_k = am_recur_map[k];
    const int arrstart_l = am_recur_map[l];

    struct RecurInfo const * const aminfo_i = &recurinfo_array[arrstart_i]; 
    struct RecurInfo const * const aminfo_j = &recurinfo_array[arrstart_j]; 
    struct RecurInfo const * const aminfo_k = &recurinfo_array[arrstart_k]; 
    struct RecurInfo const * const aminfo_l = &recurinfo_array[arrstart_l]; 

    int curidx = 0;
    for(int idx_i = 0; idx_i < ncart_i; idx_i++)
    {
        struct RecurInfo const * const ri = aminfo_i + idx_i;
        const int dir = ri->dir;
        const int idx_i_1 = ri->idx[dir][0];
        const int idx_i_2 = ri->idx[dir][1];
        const int base_idx1_i = idx_i_1 * ncart_j   * ncart_k   * ncart_l;
        const int base_idx3_i = idx_i_2 * ncart_j   * ncart_k   * ncart_l;
        const int base_idx5_i = idx_i_1 * ncart_j_1 * ncart_k   * ncart_l;
        const int base_idx6_i = idx_i_1 * ncart_j   * ncart_k_1 * ncart_l;
        const int base_idx7_i = idx_i_1 * ncart_j   * ncart_k_1 * ncart_l_1;

        for(int16_t idx_j = 0; idx_j < ncart_j; idx_j++)
        {
            struct RecurInfo const * const rj = aminfo_j + idx_j;
            const int idx_j_1 = rj->idx[dir][0];
            const int base_idx1_ij = base_idx1_i + idx_j   * ncart_k   * ncart_l;
            const int base_idx3_ij = base_idx3_i + idx_j   * ncart_k   * ncart_l;
            const int base_idx5_ij = base_idx5_i + idx_j_1 * ncart_k   * ncart_l;
            const int base_idx6_ij = base_idx6_i + idx_j   * ncart_k_1 * ncart_l;
            const int base_idx7_ij = base_idx7_i + idx_j   * ncart_k   * ncart_l_1;

            for(int16_t idx_k = 0; idx_k < ncart_k; idx_k++)
            {
                struct RecurInfo const * const rk = aminfo_k + idx_k;
                const int idx_k_1 = rk->idx[dir][0];
                const int base_idx1_ijk = base_idx1_ij + idx_k   * ncart_l;
                const int base_idx3_ijk = base_idx3_ij + idx_k   * ncart_l;
                const int base_idx5_ijk = base_idx5_ij + idx_k   * ncart_l;
                const int base_idx6_ijk = base_idx6_ij + idx_k_1 * ncart_l;
                const int base_idx7_ijk = base_idx7_ij + idx_k   * ncart_l_1;

                for(int16_t idx_l = 0; idx_l < ncart_l; idx_l++)
                {
                    struct RecurInfo const * const rl = aminfo_l + idx_l;
                    const int idx_l_1 = rl->idx[dir][0];

                    const int base_idx1 = base_idx1_ijk + idx_l;
                    const int base_idx3 = base_idx3_ijk + idx_l;
                    const int base_idx5 = base_idx5_ijk + idx_l;
                    const int base_idx6 = base_idx6_ijk + idx_l;
                    const int base_idx7 = base_idx7_ijk + idx_l_1;

                    for(int n = 0; n < num_n; n++)
                    {
                        const int outidx = n * ncart + curidx;

                        const int idx1 = base_idx1 + n*ncart1;
                        const int idx2 = idx1 + ncart1;

                        output[outidx] = PA[dir]*theta1[idx1] + aop_PQ[dir]*theta1[idx2];

                        if(ri->ijk[dir] > 1)
                        {
                            const int idx3 = base_idx3 + n*ncart2;
                            const int idx4 = idx3 + ncart2;
                            const __m256d ival = _mm256_set1_pd(ri->ijk[dir]-1);
                            output[outidx] += one_over_2p * ival * (theta2[idx3] + a_over_p * theta2[idx4]);
                        }

                        if(rj->ijk[dir])
                        {
                            const int idx5 = base_idx5 + ncart3;
                            const int idx6 = idx5 + ncart3;
                            const __m256d jval = _mm256_set1_pd(rj->ijk[dir]);
                            output[outidx] += one_over_2p * jval * (theta3[idx5] + a_over_p * theta3[idx6]);
                        }

                        if(rk->ijk[dir])
                        {
                            const int idx6 = base_idx6 + (n+1)*ncart4;
                            const __m256d kval = _mm256_set1_pd(rk->ijk[dir]);
                            output[outidx] += one_over_2pq * kval * theta4[idx6];
                        }

                        if(rl->ijk[dir])
                        {
                            const int idx7 = base_idx7 + (n+1)*ncart5;
                            const __m256d lval = _mm256_set1_pd(rl->ijk[dir]);
                            output[outidx] += one_over_2pq * lval * theta5[idx7];
                        }
                    }

                    curidx++;
                }
            }
        }
    }
}


void ostei_general_vrr2_j(int i, int j, int k, int l, int num_n,
                          __m256d one_over_2p, __m256d a_over_p,
                          __m256d one_over_2pq,
                          __m256d aop_PQ[3], __m256d PB[3],
                          __m256d * theta1, __m256d * theta2,
                          __m256d * theta3, __m256d * theta4,
                          __m256d * theta5, __m256d * output)
{
    const int ncart_i = ((i+1)*(i+2))/2;
    const int ncart_j = ((j+1)*(j+2))/2;
    const int ncart_k = ((k+1)*(k+2))/2;
    const int ncart_l = ((l+1)*(l+2))/2;

    const int ncart_i_1 = ((i)*(i+1))/2;
    const int ncart_j_1 = ((j)*(j+1))/2;
    const int ncart_k_1 = ((k)*(k+1))/2;
    const int ncart_l_1 = ((l)*(l+1))/2;

    const int ncart_j_2 = ((i-1)*(i))/2;

    const int ncart = ncart_i * ncart_j * ncart_k * ncart_l;
    const int ncart1 = ncart_i   * ncart_j_1 * ncart_k   * ncart_l;
    const int ncart2 = ncart_i_1 * ncart_j_1 * ncart_k   * ncart_l;
    const int ncart3 = ncart_i   * ncart_j_2 * ncart_k   * ncart_l;
    const int ncart4 = ncart_i   * ncart_j_1 * ncart_k_1 * ncart_l;
    const int ncart5 = ncart_i   * ncart_j_1 * ncart_k   * ncart_l_1;

    const int arrstart_i = am_recur_map[i];
    const int arrstart_j = am_recur_map[j];
    const int arrstart_k = am_recur_map[k];
    const int arrstart_l = am_recur_map[l];

    struct RecurInfo const * const aminfo_i = &recurinfo_array[arrstart_i]; 
    struct RecurInfo const * const aminfo_j = &recurinfo_array[arrstart_j]; 
    struct RecurInfo const * const aminfo_k = &recurinfo_array[arrstart_k]; 
    struct RecurInfo const * const aminfo_l = &recurinfo_array[arrstart_l]; 

    int curidx = 0;
    for(int idx_i = 0; idx_i < ncart_i; idx_i++)
    {
        struct RecurInfo const * const ri = aminfo_i + idx_i;
        const int dir = ri->dir;
        const int idx_i_1 = ri->idx[dir][0];
        const int base_idx1_i = idx_i   * ncart_j_1 * ncart_k   * ncart_l;
        const int base_idx3_i = idx_i_1 * ncart_j_1 * ncart_k   * ncart_l;
        const int base_idx5_i = idx_i   * ncart_j_2 * ncart_k   * ncart_l;
        const int base_idx6_i = idx_i   * ncart_j_1 * ncart_k_1 * ncart_l;
        const int base_idx7_i = idx_i   * ncart_j_1 * ncart_k_1 * ncart_l_1;

        for(int16_t idx_j = 0; idx_j < ncart_j; idx_j++)
        {
            struct RecurInfo const * const rj = aminfo_j + idx_j;
            const int idx_j_1 = rj->idx[dir][0];
            const int idx_j_2 = rj->idx[dir][1];
            const int base_idx1_ij = base_idx1_i + idx_j_1 * ncart_k   * ncart_l;
            const int base_idx3_ij = base_idx3_i + idx_j_1 * ncart_k   * ncart_l;
            const int base_idx5_ij = base_idx5_i + idx_j_2 * ncart_k   * ncart_l;
            const int base_idx6_ij = base_idx6_i + idx_j_1 * ncart_k_1 * ncart_l;
            const int base_idx7_ij = base_idx7_i + idx_j_1 * ncart_k   * ncart_l_1;

            for(int16_t idx_k = 0; idx_k < ncart_k; idx_k++)
            {
                struct RecurInfo const * const rk = aminfo_k + idx_k;
                const int idx_k_1 = rk->idx[dir][0];
                const int base_idx1_ijk = base_idx1_ij + idx_k   * ncart_l;
                const int base_idx3_ijk = base_idx3_ij + idx_k   * ncart_l;
                const int base_idx5_ijk = base_idx5_ij + idx_k   * ncart_l;
                const int base_idx6_ijk = base_idx6_ij + idx_k_1 * ncart_l;
                const int base_idx7_ijk = base_idx7_ij + idx_k   * ncart_l_1;

                for(int16_t idx_l = 0; idx_l < ncart_l; idx_l++)
                {
                    struct RecurInfo const * const rl = aminfo_l + idx_l;
                    const int idx_l_1 = rl->idx[dir][0];

                    const int base_idx1 = base_idx1_ijk + idx_l;
                    const int base_idx3 = base_idx3_ijk + idx_l;
                    const int base_idx5 = base_idx5_ijk + idx_l;
                    const int base_idx6 = base_idx6_ijk + idx_l;
                    const int base_idx7 = base_idx7_ijk + idx_l_1;

                    for(int n = 0; n < num_n; n++)
                    {
                        const int outidx = n * ncart + curidx;

                        const int idx1 = base_idx1 + n*ncart1;
                        const int idx2 = idx1 + ncart1;

                        output[outidx] = PB[dir]*theta1[idx1] + aop_PQ[dir]*theta1[idx2];

                        if(ri->ijk[dir])
                        {
                            const int idx3 = base_idx3 + n*ncart2;
                            const int idx4 = idx3 + ncart2;
                            const __m256d ival = _mm256_set1_pd(ri->ijk[dir]);
                            output[outidx] += one_over_2p * ival * (theta2[idx3] + a_over_p * theta2[idx4]);
                        }

                        if(rj->ijk[dir] > 1)
                        {
                            const int idx5 = base_idx5 + ncart3;
                            const int idx6 = idx5 + ncart3;
                            const __m256d jval = _mm256_set1_pd(rj->ijk[dir]-1);
                            output[outidx] += one_over_2p * jval * (theta3[idx5] + a_over_p * theta3[idx6]);
                        }

                        if(rk->ijk[dir])
                        {
                            const int idx6 = base_idx6 + (n+1)*ncart4;
                            const __m256d kval = _mm256_set1_pd(rk->ijk[dir]);
                            output[outidx] += one_over_2pq * kval * theta4[idx6];
                        }

                        if(rl->ijk[dir])
                        {
                            const int idx7 = base_idx7 + (n+1)*ncart5;
                            const __m256d lval = _mm256_set1_pd(rl->ijk[dir]);
                            output[outidx] += one_over_2pq * lval * theta5[idx7];
                        }
                    }

                    curidx++;
                }
            }
        }
    }
}


/*
void ostei_general_vrr2_ket(int i, int j, int k, int l, int num_n,
                            __m256d one_over_2q, __m256d a_over_q,
                            __m256d one_over_2pq,
                            __m256d aoq_PQ[3], __m256d QC[3],
                            __m256d * theta1, __m256d * theta2,
                            __m256d * theta3, __m256d * theta4,
                            __m256d * theta5, __m256d * output)
{
    const int ncart_i = ((i+1)*(i+2))/2;
    const int ncart_j = ((j+1)*(j+2))/2;
    const int ncart_k = ((k+1)*(k+2))/2;
    const int ncart_l = ((l+1)*(l+2))/2;

    const int ncart_i_1 = ((i)*(i+1))/2;
    const int ncart_j_1 = ((j)*(j+1))/2;
    const int ncart_k_1 = ((k)*(k+1))/2;
    const int ncart_l_1 = ((l)*(l+1))/2;

    const int ncart_k_2 = ((k-1)*(k))/2;

    const int ncart = ncart_i * ncart_j * ncart_k * ncart_l;
    const int ncart1 = ncart_i   * ncart_j   * ncart_k_1 * ncart_l;
    const int ncart2 = ncart_i   * ncart_j   * ncart_k_2 * ncart_l;
    const int ncart3 = ncart_i   * ncart_j   * ncart_k_1 * ncart_l_1;
    const int ncart4 = ncart_i_1 * ncart_j   * ncart_k_1 * ncart_l;
    const int ncart5 = ncart_i   * ncart_j_1 * ncart_k_1 * ncart_l;

    const int arrstart_i = am_recur_map[i];
    const int arrstart_j = am_recur_map[j];
    const int arrstart_k = am_recur_map[k];
    const int arrstart_l = am_recur_map[l];

    struct RecurInfo const * const aminfo_i = &recurinfo_array[arrstart_i]; 
    struct RecurInfo const * const aminfo_j = &recurinfo_array[arrstart_j]; 
    struct RecurInfo const * const aminfo_k = &recurinfo_array[arrstart_k]; 
    struct RecurInfo const * const aminfo_l = &recurinfo_array[arrstart_l]; 

    int curidx = 0;
    for(int16_t idx_i = 0; idx_i < ncart_i; idx_i++)
    {
        struct RecurInfo const * const ri = aminfo_i + idx_i;

        for(int16_t idx_j = 0; idx_j < ncart_j; idx_j++)
        {
            struct RecurInfo const * const rj = aminfo_j + idx_j;

            for(int16_t idx_k = 0; idx_k < ncart_k; idx_k++)
            {
                struct RecurInfo const * const rk = aminfo_k + idx_k;
                const int dir = rk->dir;

                for(int16_t idx_l = 0; idx_l < ncart_l; idx_l++)
                {
                    struct RecurInfo const * const rl = aminfo_l + idx_l;

                    for(int n = 0; n < num_n; n++)
                    {
                        const int outidx = n * ncart + curidx;
                        const int idx1 =  n*ncart1 + 
                                          idx_i           * ncart_j   * ncart_k_1 * ncart_l +
                                          idx_j           * ncart_k_1 * ncart_l   +
                                          rk->idx[dir][0] * ncart_l   +
                                          idx_l;

                        const int idx2 = idx1 + ncart1;

                        output[outidx] = QC[dir]*theta1[idx1] + aoq_PQ[dir]*theta1[idx2];

                        if(rk->ijk[dir] > 1)
                        {
                            const int idx3 =  n*ncart2 +
                                              idx_i           * ncart_j   * ncart_k_2 * ncart_l +
                                              idx_j           * ncart_k_2 * ncart_l   +
                                              rk->idx[dir][1] * ncart_l   +
                                              idx_l;

                            const int idx4 = idx3 + ncart2;
                            const __m256d kval = _mm256_set1_pd(rk->ijk[dir]-1);
                            output[outidx] += one_over_2q * kval * (theta2[idx3] + a_over_q * theta2[idx4]);
                        }

                        if(rl->ijk[dir])
                        {
                            const int idx5 = n*ncart3 +
                                             idx_i           * ncart_j   * ncart_k_1 * ncart_l_1 +
                                             idx_j           * ncart_k_1 * ncart_l_1 +
                                             rk->idx[dir][0] * ncart_l_1 +
                                             rl->idx[dir][0];

                            const int idx6 = idx5 + ncart3;
                            const __m256d lval = _mm256_set1_pd(rl->ijk[dir]);
                            output[outidx] += one_over_2q * lval * (theta3[idx5] + a_over_q * theta3[idx6]);
                        }

                        if(ri->ijk[dir])
                        {
                            const int idx6 = (n+1)*ncart4 +
                                             ri->idx[dir][0] * ncart_j   * ncart_k_1 * ncart_l +
                                             idx_j           * ncart_k_1 * ncart_l   +
                                             rk->idx[dir][0] * ncart_l   +
                                             idx_l;

                            const __m256d ival = _mm256_set1_pd(ri->ijk[dir]);
                            output[outidx] += one_over_2pq * ival * theta4[idx6];
                        }

                        if(rj->ijk[dir])
                        {
                            const int idx7 = (n+1)*ncart5 +
                                             idx_i           * ncart_j_1 * ncart_k_1 * ncart_l +
                                             rj->idx[dir][0] * ncart_k_1 * ncart_l   +
                                             rk->idx[dir][0] * ncart_l   +
                                             idx_l;

                            const __m256d jval = _mm256_set1_pd(rj->ijk[dir]);
                            output[outidx] += one_over_2pq * jval * theta5[idx7];
                        }
                    }

                    curidx++;
                }
            }
        }
    }
}
*/
