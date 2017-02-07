#include "simint/ostei/ostei_general.h"
#include "simint/ostei/recur_lookup.h"


/**********************************************************/
/* All vrr1 are the same, with some swapping of variables */
/**********************************************************/
void ostei_general_vrr1(int i, int num_n,
                        SIMINT_DBLTYPE one_over_2p, SIMINT_DBLTYPE a_over_p,
                        const SIMINT_DBLTYPE aop_PQ[3],
                        const SIMINT_DBLTYPE PA[3],
                        SIMINT_DBLTYPE const * restrict theta1,
                        SIMINT_DBLTYPE const * restrict theta2,
                        SIMINT_DBLTYPE * restrict output)
{
    // in this case, we want (i, 0, 0, 0), so i is equal
    // to the angular momentum.

    const int ncart = ((i+1)*(i+2))/2;
    const int ncart1 = ((i)*(i+1))/2;
    const int ncart2 = ((i-1)*(i))/2;

    const int arrstart = am_recur_map[i];
    struct RecurInfo const * aminfo = &recurinfo_array[arrstart];

    int curidx = 0;
    for(int idx_i = 0; idx_i < ncart; idx_i++)
    {
        struct RecurInfo const * ri = aminfo + idx_i;
        const int dir = ri->dir;

        for(int n = 0; n < num_n; n++)
        {
            const int outidx = n * ncart + curidx;
            const int idx1 = n*(ncart1) + ri->idx[dir][0];
            const int idx2 = idx1 + ncart1;

            output[outidx] = SIMINT_FMADD(PA[dir], theta1[idx1], SIMINT_MUL(aop_PQ[dir], theta1[idx2]));

            if(ri->ijk[dir] > 1)
            {
                const int idx3 = n*(ncart2) + ri->idx[dir][1];
                const int idx4 = (n+1)*(ncart2) + ri->idx[dir][1];
                const SIMINT_DBLTYPE i1 = SIMINT_DBLSET1(ri->ijk[dir]-1);

                const SIMINT_DBLTYPE tmp1 = SIMINT_FMADD(a_over_p, theta2[idx4], theta2[idx3]);
                const SIMINT_DBLTYPE tmp2 = SIMINT_MUL(i1, one_over_2p);
                output[outidx] = SIMINT_FMADD(tmp1, tmp2, output[outidx]);
            }
        }
        curidx++;
    }
}



void ostei_general_vrr_I(int i, int j, int k, int l, int num_n,
                         SIMINT_DBLTYPE one_over_2p, SIMINT_DBLTYPE a_over_p,
                         SIMINT_DBLTYPE one_over_2pq,
                         const SIMINT_DBLTYPE aop_PQ[3],
                         const SIMINT_DBLTYPE PA[3],
                         SIMINT_DBLTYPE const * restrict theta1,
                         SIMINT_DBLTYPE const * restrict theta2,
                         SIMINT_DBLTYPE const * restrict theta3,
                         SIMINT_DBLTYPE const * restrict theta4,
                         SIMINT_DBLTYPE const * restrict theta5,
                         SIMINT_DBLTYPE * restrict output)
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

    const int ncart  = ncart_i   * ncart_j   * ncart_k   * ncart_l;
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

        for(int idx_j = 0; idx_j < ncart_j; idx_j++)
        {
            struct RecurInfo const * const rj = aminfo_j + idx_j;

            for(int idx_k = 0; idx_k < ncart_k; idx_k++)
            {
                struct RecurInfo const * const rk = aminfo_k + idx_k;

                for(int idx_l = 0; idx_l < ncart_l; idx_l++)
                {
                    struct RecurInfo const * const rl = aminfo_l + idx_l;
                    const int dir = ri->dir;
                    const int idx_i_1 = ri->idx[dir][0];

                    for(int n = 0; n < num_n; n++)
                    {
                        const int outidx = n * ncart + curidx;

                        const int idx1 = n * ncart1 +
                                         idx_i_1 * ncart_j   * ncart_k   * ncart_l   +
                                         idx_j   * ncart_k   * ncart_l   +
                                         idx_k   * ncart_l   +
                                         idx_l;

                        const int idx2 = idx1 + ncart1;

                        output[outidx] = SIMINT_FMADD(PA[dir], theta1[idx1], SIMINT_MUL(aop_PQ[dir], theta1[idx2]));

                        if(ri->ijk[dir] > 1)
                        {
                            const int idx_i_2 = ri->idx[dir][1];
                            const int idx3 = n * ncart2 +
                                             idx_i_2 * ncart_j   * ncart_k   * ncart_l   +
                                             idx_j   * ncart_k   * ncart_l   +
                                             idx_k   * ncart_l   +
                                             idx_l;

                            const int idx4 = idx3 + ncart2;

                            const SIMINT_DBLTYPE ival = SIMINT_DBLSET1(ri->ijk[dir]-1);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_FMADD(a_over_p, theta2[idx4], theta2[idx3]);
                            const SIMINT_DBLTYPE tmp2 = SIMINT_MUL(ival, one_over_2p);
                            output[outidx] = SIMINT_FMADD(tmp1, tmp2, output[outidx]);
                        }

                        if(rj->ijk[dir])
                        {
                            const int idx_j_1 = rj->idx[dir][0];
                            const int idx5 = n * ncart3 +
                                             idx_i_1 * ncart_j_1 * ncart_k   * ncart_l   +
                                             idx_j_1 * ncart_k   * ncart_l   +
                                             idx_k   * ncart_l   +
                                             idx_l;

                            const int idx6 = idx5 + ncart3;

                            const SIMINT_DBLTYPE jval = SIMINT_DBLSET1(rj->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_FMADD(a_over_p, theta3[idx6], theta3[idx5]);
                            const SIMINT_DBLTYPE tmp2 = SIMINT_MUL(jval, one_over_2p);
                            output[outidx] = SIMINT_FMADD(tmp1, tmp2, output[outidx]);
                        }

                        if(rk->ijk[dir])
                        {
                            const int idx_k_1 = rk->idx[dir][0];
                            const int idx7 = (n+1) * ncart4 +
                                             idx_i_1 * ncart_j   * ncart_k_1 * ncart_l   +
                                             idx_j   * ncart_k_1 * ncart_l   +
                                             idx_k_1 * ncart_l   +
                                             idx_l;

                            const SIMINT_DBLTYPE kval = SIMINT_DBLSET1(rk->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_MUL(one_over_2pq, kval);
                            output[outidx] = SIMINT_FMADD(tmp1, theta4[idx7], output[outidx]);
                        }

                        if(rl->ijk[dir])
                        {
                            const int idx_l_1 = rl->idx[dir][0];
                            const int idx8 = (n+1) * ncart5 +
                                             idx_i_1 * ncart_j   * ncart_k   * ncart_l_1 +
                                             idx_j   * ncart_k   * ncart_l_1 +
                                             idx_k   * ncart_l_1 +
                                             idx_l_1;

                            const SIMINT_DBLTYPE lval = SIMINT_DBLSET1(rl->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_MUL(one_over_2pq, lval);
                            output[outidx] = SIMINT_FMADD(tmp1, theta5[idx8], output[outidx]);
                        }
                    }

                    curidx++;
                }
            }
        }
    }
}


void ostei_general_vrr_J(int i, int j, int k, int l, int num_n,
                         SIMINT_DBLTYPE one_over_2p, SIMINT_DBLTYPE a_over_p,
                         SIMINT_DBLTYPE one_over_2pq,
                         const SIMINT_DBLTYPE aop_PQ[3],
                         const SIMINT_DBLTYPE PB[3],
                         SIMINT_DBLTYPE const * restrict theta1,
                         SIMINT_DBLTYPE const * restrict theta2,
                         SIMINT_DBLTYPE const * restrict theta3,
                         SIMINT_DBLTYPE const * restrict theta4,
                         SIMINT_DBLTYPE const * restrict theta5,
                         SIMINT_DBLTYPE * restrict output)
{
    const int ncart_i = ((i+1)*(i+2))/2;
    const int ncart_j = ((j+1)*(j+2))/2;
    const int ncart_k = ((k+1)*(k+2))/2;
    const int ncart_l = ((l+1)*(l+2))/2;

    const int ncart_i_1 = ((i)*(i+1))/2;
    const int ncart_j_1 = ((j)*(j+1))/2;
    const int ncart_k_1 = ((k)*(k+1))/2;
    const int ncart_l_1 = ((l)*(l+1))/2;

    const int ncart_j_2 = ((j-1)*(j))/2;

    const int ncart  = ncart_i   * ncart_j   * ncart_k   * ncart_l;
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

        for(int idx_j = 0; idx_j < ncart_j; idx_j++)
        {
            struct RecurInfo const * const rj = aminfo_j + idx_j;

            for(int idx_k = 0; idx_k < ncart_k; idx_k++)
            {
                struct RecurInfo const * const rk = aminfo_k + idx_k;

                for(int idx_l = 0; idx_l < ncart_l; idx_l++)
                {
                    struct RecurInfo const * const rl = aminfo_l + idx_l;
                    const int dir = rj->dir;
                    const int idx_j_1 = rj->idx[dir][0];

                    for(int n = 0; n < num_n; n++)
                    {
                        const int outidx = n * ncart + curidx;

                        const int idx1 = n * ncart1 +
                                         idx_i   * ncart_j_1 * ncart_k   * ncart_l   +
                                         idx_j_1 * ncart_k   * ncart_l   +
                                         idx_k   * ncart_l   +
                                         idx_l;

                        const int idx2 = idx1 + ncart1;

                        output[outidx] = SIMINT_FMADD(PB[dir], theta1[idx1], SIMINT_MUL(aop_PQ[dir], theta1[idx2]));

                        if(ri->ijk[dir])
                        {
                            const int idx_i_1 = ri->idx[dir][0];
                            const int idx3 = n * ncart2 +
                                             idx_i_1 * ncart_j_1 * ncart_k   * ncart_l   +
                                             idx_j_1 * ncart_k   * ncart_l   +
                                             idx_k   * ncart_l   +
                                             idx_l;

                            const int idx4 = idx3 + ncart2;

                            const SIMINT_DBLTYPE ival = SIMINT_DBLSET1(ri->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_FMADD(a_over_p, theta2[idx4], theta2[idx3]);
                            const SIMINT_DBLTYPE tmp2 = SIMINT_MUL(ival, one_over_2p);
                            output[outidx] = SIMINT_FMADD(tmp1, tmp2, output[outidx]);
                        }

                        if(rj->ijk[dir] > 1)
                        {
                            const int idx_j_2 = rj->idx[dir][1];
                            const int idx5 = n * ncart3 +
                                             idx_i   * ncart_j_2 * ncart_k   * ncart_l   +
                                             idx_j_2 * ncart_k   * ncart_l   +
                                             idx_k   * ncart_l   +
                                             idx_l;

                            const int idx6 = idx5 + ncart3;

                            const SIMINT_DBLTYPE jval = SIMINT_DBLSET1(rj->ijk[dir]-1);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_FMADD(a_over_p, theta3[idx6], theta3[idx5]);
                            const SIMINT_DBLTYPE tmp2 = SIMINT_MUL(jval, one_over_2p);
                            output[outidx] = SIMINT_FMADD(tmp1, tmp2, output[outidx]);
                        }

                        if(rk->ijk[dir])
                        {
                            const int idx_k_1 = rk->idx[dir][0];
                            const int idx7 = (n+1) * ncart4 +
                                             idx_i   * ncart_j_1 * ncart_k_1 * ncart_l   +
                                             idx_j_1 * ncart_k_1 * ncart_l   +
                                             idx_k_1 * ncart_l   +
                                             idx_l;

                            const SIMINT_DBLTYPE kval = SIMINT_DBLSET1(rk->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_MUL(one_over_2pq, kval);
                            output[outidx] = SIMINT_FMADD(tmp1, theta4[idx7], output[outidx]);
                        }

                        if(rl->ijk[dir])
                        {
                            const int idx_l_1 = rl->idx[dir][0];
                            const int idx8 = (n+1) * ncart5 +
                                             idx_i   * ncart_j_1 * ncart_k   * ncart_l_1 +
                                             idx_j_1 * ncart_k   * ncart_l_1 +
                                             idx_k   * ncart_l_1 +
                                             idx_l_1;

                            const SIMINT_DBLTYPE lval = SIMINT_DBLSET1(rl->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_MUL(one_over_2pq, lval);
                            output[outidx] = SIMINT_FMADD(tmp1, theta5[idx8], output[outidx]);
                        }
                    }

                    curidx++;
                }
            }
        }
    }
}


void ostei_general_vrr_K(int i, int j, int k, int l, int num_n,
                         SIMINT_DBLTYPE one_over_2q, SIMINT_DBLTYPE a_over_q,
                         SIMINT_DBLTYPE one_over_2pq,
                         const SIMINT_DBLTYPE aoq_PQ[3],
                         const SIMINT_DBLTYPE QC[3],
                         SIMINT_DBLTYPE const * restrict theta1,
                         SIMINT_DBLTYPE const * restrict theta2,
                         SIMINT_DBLTYPE const * restrict theta3,
                         SIMINT_DBLTYPE const * restrict theta4,
                         SIMINT_DBLTYPE const * restrict theta5,
                         SIMINT_DBLTYPE * restrict output)
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

    const int ncart  = ncart_i   * ncart_j   * ncart_k   * ncart_l;
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
    for(int idx_i = 0; idx_i < ncart_i; idx_i++)
    {
        struct RecurInfo const * const ri = aminfo_i + idx_i;

        for(int idx_j = 0; idx_j < ncart_j; idx_j++)
        {
            struct RecurInfo const * const rj = aminfo_j + idx_j;

            for(int idx_k = 0; idx_k < ncart_k; idx_k++)
            {
                struct RecurInfo const * const rk = aminfo_k + idx_k;

                for(int idx_l = 0; idx_l < ncart_l; idx_l++)
                {
                    struct RecurInfo const * const rl = aminfo_l + idx_l;
                    const int dir = rk->dir;
                    const int idx_k_1 = rk->idx[dir][0];

                    for(int n = 0; n < num_n; n++)
                    {
                        const int outidx = n * ncart + curidx;

                        const int idx1 = n * ncart1 +
                                         idx_i   * ncart_j   * ncart_k_1 * ncart_l   +
                                         idx_j   * ncart_k_1 * ncart_l   +
                                         idx_k_1 * ncart_l   +
                                         idx_l;

                        const int idx2 = idx1 + ncart1;

                        output[outidx] = SIMINT_FMADD(QC[dir], theta1[idx1], SIMINT_MUL(aoq_PQ[dir], theta1[idx2]));

                        if(rk->ijk[dir] > 1)
                        {
                            const int idx_k_2 = rk->idx[dir][1];
                            const int idx3 = n * ncart2 +
                                             idx_i   * ncart_j   * ncart_k_2 * ncart_l   +
                                             idx_j   * ncart_k_2 * ncart_l   +
                                             idx_k_2 * ncart_l   +
                                             idx_l;

                            const int idx4 = idx3 + ncart2;

                            const SIMINT_DBLTYPE kval = SIMINT_DBLSET1(rk->ijk[dir]-1);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_FMADD(a_over_q, theta2[idx4], theta2[idx3]);
                            const SIMINT_DBLTYPE tmp2 = SIMINT_MUL(one_over_2q, kval);
                            output[outidx] = SIMINT_FMADD(tmp1, tmp2, output[outidx]);
                        }

                        if(rl->ijk[dir])
                        {
                            const int idx_l_1 = rl->idx[dir][0];
                            const int idx5 = n * ncart3 +
                                             idx_i   * ncart_j   * ncart_k_1 * ncart_l_1 +
                                             idx_j   * ncart_k_1 * ncart_l_1 +
                                             idx_k_1 * ncart_l_1 +
                                             idx_l_1;

                            const int idx6 = idx5 + ncart3;

                            const SIMINT_DBLTYPE lval = SIMINT_DBLSET1(rl->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_FMADD(a_over_q, theta3[idx6], theta3[idx5]);
                            const SIMINT_DBLTYPE tmp2 = SIMINT_MUL(one_over_2q, lval);
                            output[outidx] = SIMINT_FMADD(tmp1, tmp2, output[outidx]);
                        }

                        if(ri->ijk[dir])
                        {
                            const int idx_i_1 = ri->idx[dir][0];
                            const int idx7 = (n+1) * ncart4 +
                                             idx_i_1 * ncart_j   * ncart_k_1 * ncart_l   +
                                             idx_j   * ncart_k_1 * ncart_l   +
                                             idx_k_1 * ncart_l   +
                                             idx_l;

                            const SIMINT_DBLTYPE ival = SIMINT_DBLSET1(ri->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_MUL(one_over_2pq, ival);
                            output[outidx] = SIMINT_FMADD(tmp1, theta4[idx7], output[outidx]);
                        }

                        if(rj->ijk[dir])
                        {
                            const int idx_j_1 = rj->idx[dir][0];
                            const int idx8 = (n+1) * ncart5 +
                                             idx_i   * ncart_j_1 * ncart_k_1 * ncart_l   +
                                             idx_j_1 * ncart_k_1 * ncart_l   +
                                             idx_k_1 * ncart_l +
                                             idx_l;

                            const SIMINT_DBLTYPE jval = SIMINT_DBLSET1(rj->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_MUL(one_over_2pq, jval);
                            output[outidx] = SIMINT_FMADD(tmp1, theta5[idx8], output[outidx]);
                        }
                    }

                    curidx++;
                }
            }
        }
    }
}

void ostei_general_vrr_L(int i, int j, int k, int l, int num_n,
                         SIMINT_DBLTYPE one_over_2q, SIMINT_DBLTYPE a_over_q,
                         SIMINT_DBLTYPE one_over_2pq,
                         const SIMINT_DBLTYPE aoq_PQ[3],
                         const SIMINT_DBLTYPE QD[3],
                         SIMINT_DBLTYPE const * restrict theta1,
                         SIMINT_DBLTYPE const * restrict theta2,
                         SIMINT_DBLTYPE const * restrict theta3,
                         SIMINT_DBLTYPE const * restrict theta4,
                         SIMINT_DBLTYPE const * restrict theta5,
                         SIMINT_DBLTYPE * restrict output)
{
    const int ncart_i = ((i+1)*(i+2))/2;
    const int ncart_j = ((j+1)*(j+2))/2;
    const int ncart_k = ((k+1)*(k+2))/2;
    const int ncart_l = ((l+1)*(l+2))/2;

    const int ncart_i_1 = ((i)*(i+1))/2;
    const int ncart_j_1 = ((j)*(j+1))/2;
    const int ncart_k_1 = ((k)*(k+1))/2;
    const int ncart_l_1 = ((l)*(l+1))/2;

    const int ncart_l_2 = ((l-1)*(l))/2;

    const int ncart  = ncart_i   * ncart_j   * ncart_k   * ncart_l;
    const int ncart1 = ncart_i   * ncart_j   * ncart_k   * ncart_l_1;
    const int ncart2 = ncart_i   * ncart_j   * ncart_k_1 * ncart_l_1;
    const int ncart3 = ncart_i   * ncart_j   * ncart_k   * ncart_l_2;
    const int ncart4 = ncart_i_1 * ncart_j   * ncart_k   * ncart_l_1;
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

        for(int idx_j = 0; idx_j < ncart_j; idx_j++)
        {
            struct RecurInfo const * const rj = aminfo_j + idx_j;

            for(int idx_k = 0; idx_k < ncart_k; idx_k++)
            {
                struct RecurInfo const * const rk = aminfo_k + idx_k;

                for(int idx_l = 0; idx_l < ncart_l; idx_l++)
                {
                    struct RecurInfo const * const rl = aminfo_l + idx_l;
                    const int dir = rl->dir;
                    const int idx_l_1 = rl->idx[dir][0];

                    for(int n = 0; n < num_n; n++)
                    {
                        const int outidx = n * ncart + curidx;

                        const int idx1 = n * ncart1 +
                                         idx_i   * ncart_j   * ncart_k   * ncart_l_1 +
                                         idx_j   * ncart_k   * ncart_l_1 +
                                         idx_k   * ncart_l_1 +
                                         idx_l_1;

                        const int idx2 = idx1 + ncart1;

                        output[outidx] = SIMINT_FMADD(QD[dir], theta1[idx1], SIMINT_MUL(aoq_PQ[dir], theta1[idx2]));

                        if(rk->ijk[dir])
                        {
                            const int idx_k_1 = rk->idx[dir][0];
                            const int idx3 = n * ncart2 +
                                             idx_i   * ncart_j   * ncart_k_1 * ncart_l_1 +
                                             idx_j   * ncart_k_1 * ncart_l_1 +
                                             idx_k_1 * ncart_l_1 +
                                             idx_l_1;

                            const int idx4 = idx3 + ncart2;

                            const SIMINT_DBLTYPE kval = SIMINT_DBLSET1(rk->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_FMADD(a_over_q, theta2[idx4], theta2[idx3]);
                            const SIMINT_DBLTYPE tmp2 = SIMINT_MUL(one_over_2q, kval);
                            output[outidx] = SIMINT_FMADD(tmp1, tmp2, output[outidx]);
                        }

                        if(rl->ijk[dir] > 1)
                        {
                            const int idx_l_2 = rl->idx[dir][1];
                            const int idx5 = n * ncart3 +
                                             idx_i   * ncart_j   * ncart_k   * ncart_l_2 +
                                             idx_j   * ncart_k   * ncart_l_2 +
                                             idx_k   * ncart_l_2 +
                                             idx_l_2;

                            const int idx6 = idx5 + ncart3;

                            const SIMINT_DBLTYPE lval = SIMINT_DBLSET1(rl->ijk[dir]-1);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_FMADD(a_over_q, theta3[idx6], theta3[idx5]);
                            const SIMINT_DBLTYPE tmp2 = SIMINT_MUL(one_over_2q, lval);
                            output[outidx] = SIMINT_FMADD(tmp1, tmp2, output[outidx]);
                        }

                        if(ri->ijk[dir])
                        {
                            const int idx_i_1 = ri->idx[dir][0];
                            const int idx7 = (n+1) * ncart4 +
                                             idx_i_1 * ncart_j   * ncart_k   * ncart_l_1 +
                                             idx_j   * ncart_k   * ncart_l_1 +
                                             idx_k   * ncart_l_1 +
                                             idx_l_1;

                            const SIMINT_DBLTYPE ival = SIMINT_DBLSET1(ri->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_MUL(one_over_2pq, ival);
                            output[outidx] = SIMINT_FMADD(tmp1, theta4[idx7], output[outidx]);
                        }

                        if(rj->ijk[dir])
                        {
                            const int idx_j_1 = rj->idx[dir][0];
                            const int idx8 = (n+1) * ncart5 +
                                             idx_i   * ncart_j_1 * ncart_k   * ncart_l_1   +
                                             idx_j_1 * ncart_k   * ncart_l_1 +
                                             idx_k   * ncart_l_1 +
                                             idx_l_1;

                            const SIMINT_DBLTYPE jval = SIMINT_DBLSET1(rj->ijk[dir]);
                            const SIMINT_DBLTYPE tmp1 = SIMINT_MUL(one_over_2pq, jval);
                            output[outidx] = SIMINT_FMADD(tmp1, theta5[idx8], output[outidx]);
                        }
                    }

                    curidx++;
                }
            }
        }
    }
}
