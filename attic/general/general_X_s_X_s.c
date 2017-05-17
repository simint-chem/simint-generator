#include "simint/boys/boys_split.h"
#include "simint/eri/eri.h"
#include "simint/eri/general_eri.h"
#include <math.h>
#include <alloca.h>
#include <string.h>


int eri_sharedwork_X_s_X_s(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__X_s_X_s)
{
    // Since this function can be called for X_s_X_s, s_X_s_X, X_s_s_X, and s_X_X_s,
    // we combined the am of the bra into a single L1 value, and the
    // am of the ket into L2. This allows a bit of generality

    const int ncart1 = ((P.am1+1)*(P.am1+2))/2;
    const int ncart2 = ((P.am2+1)*(P.am2+2))/2;
    const int ncart3 = ((Q.am1+1)*(Q.am1+2))/2;
    const int ncart4 = ((Q.am2+1)*(Q.am2+2))/2;
    const int ncart12 = ncart1 * ncart2;
    const int ncart34 = ncart3 * ncart4;
    const int ncart1234 = ncart12 * ncart34;
    const int L = P.am1 + P.am2 + Q.am1 + Q.am2;
    const int L1 = P.am1 + P.am2;
    const int L2 = Q.am1 + Q.am2;

    memset(INT__X_s_X_s, 0, P.nshell12 * Q.nshell12 * ncart1234 * sizeof(double));

    int ab, cd, cdbatch, abcd;
    int istart, jstart;
    int iprimcd, nprim_icd, icd;
    int i, j, n, n1, n2;
    int np;

    const __m256d const_1 = _mm256_set1_pd(1);
    const __m256d one_half = _mm256_set1_pd(0.5);

    // Size of arrays needed in the stack
    size_t arr_stack_size = 0;



    if(L1 >= L2)
    {
        // for ( i 0 | 0 0 )
        for(n1 = 0; n1 <= L1; n1++)
           arr_stack_size += ( ((n1+1)*(n1+2) )/2 ) * (L1+L2-n1+1);

        // for ( i 0 | k 0 )
        for(n1 = (L1-L2+1); n1 <= L1; n1++)
        for(n2 = 1; n2 <= (n1-(L1-L2)); n2++)
            arr_stack_size += ( ((n1+1)*(n1+2) )/2 ) * ( ((n2+1)*(n2+2))/2 ) * (L2 - n2 + 1);
    }
    else
    {
        // for ( 0 0 | (i+k) 0 )
        for(n2 = 0; n2 <= L2; n2++)
           arr_stack_size += ( ((n2+1)*(n2+2) )/2 ) * (L1+L2-n2+1);

        // for ( i 0 | k 0 )
        for(n2 = (L2-L1+1); n2 <= L2; n2++)
        for(n1 = 1; n1 <= (n2-(L2-L1)); n1++)
            arr_stack_size += ( ((n1+1)*(n1+2) )/2 ) * ( ((n2+1)*(n2+2))/2 ) * (L1 - n1 + 1);
    }


    __m256d * arr_stack = (__m256d *)(ALLOC(arr_stack_size * sizeof(__m256d)));
    __m256d * prim_ptr[L1+1][L2+1];
    __m256d * pos = arr_stack;


    if(L1 >= L2)
    {
        for(n1 = 0; n1 <= L1; n1++)
        {
            prim_ptr[n1][0] = pos;
            pos += ( ((n1+1)*(n1+2))/2 ) * (L1+L2-n1+1);
        }

        for(n1 = (L1-L2+1); n1 <= L1; n1++)
        for(n2 = 1; n2 <= (n1-(L1-L2)); n2++)
        {
            prim_ptr[n1][n2] = pos;
            pos += ( ((n1+1)*(n1+2) )/2 ) * ( ((n2+1)*(n2+2))/2 ) * (L2 - n2 + 1);
        }
    }
    else
    {
        for(n2 = 0; n2 <= L2; n2++)
        {
            prim_ptr[0][n2] = pos;
            pos += ( ((n2+1)*(n2+2))/2 ) * (L1+L2-n2+1);
        }

        for(n2 = (L2-L1+1); n2 <= L2; n2++)
        for(n1 = 1; n1 <= (n2-(L2-L1)); n1++)
        {
            prim_ptr[n1][n2] = pos;
            pos += ( ((n1+1)*(n1+2) )/2 ) * ( ((n2+1)*(n2+2))/2 ) * (L1 - n1 + 1);
        }
    }



    ////////////////////////////////////////
    // Loop over shells and primitives
    ////////////////////////////////////////

    abcd = 0;
    istart = 0;
    for(ab = 0; ab < P.nshell12; ++ab)
    {
        const int iend = istart + P.nprim12[ab];

        cd = 0;
        jstart = 0;

        for(cdbatch = 0; cdbatch < Q.nbatch; ++cdbatch)
        {
            const int nshellbatch = ((cd + SIMINT_NSHELL_SIMD) > Q.nshell12) ? Q.nshell12 - cd : SIMINT_NSHELL_SIMD;
            const int jend = jstart + Q.nbatchprim[cdbatch];
            for(i = istart; i < iend; ++i)
            {

                icd = 0;
                iprimcd = 0;
                nprim_icd = Q.nprim12[cd];
                double * restrict PRIM_PTR_INT__X_s_X_s = INT__X_s_X_s + abcd * ncart1234;

                // Load these one per loop over i
                const __m256d P_alpha = _mm256_set1_pd(P.alpha[i]);
                const __m256d P_prefac = _mm256_set1_pd(P.prefac[i]);
                const __m256d P_x = _mm256_set1_pd(P.x[i]);
                const __m256d P_y = _mm256_set1_pd(P.y[i]);
                const __m256d P_z = _mm256_set1_pd(P.z[i]);
                const __m256d P_PA_x = _mm256_set1_pd(P.PA_x[i]);
                const __m256d P_PA_y = _mm256_set1_pd(P.PA_y[i]);
                const __m256d P_PA_z = _mm256_set1_pd(P.PA_z[i]);


                for(j = jstart; j < jend; j += SIMINT_SIMD_LEN)
                {
                    // calculate the shell offsets
                    // these are the offset from the shell pointed to by cd
                    // for each element
                    int shelloffsets[SIMINT_SIMD_LEN] = {0};
                    if((iprimcd + SIMINT_SIMD_LEN) >= nprim_icd)
                    {
                        // Handle if the first element of the vector is a new shell
                        if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))
                        {
                            nprim_icd += Q.nprim12[cd + (++icd)];
                            PRIM_PTR_INT__X_s_X_s += ncart1234;
                        }
                        iprimcd++;
                        for(n = 1; n < SIMINT_SIMD_LEN; ++n)
                        {
                            if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))
                            {
                                shelloffsets[n] = shelloffsets[n-1] + 1;
                                nprim_icd += Q.nprim12[cd + (++icd)];
                            }
                            else
                                shelloffsets[n] = shelloffsets[n-1];
                            iprimcd++;
                        }
                    }
                    else
                        iprimcd += SIMINT_SIMD_LEN;


                    const __m256d Q_alpha = _mm256_load_pd(Q.alpha + j);
                    const __m256d PQalpha_mul = P_alpha * Q_alpha;
                    const __m256d PQalpha_sum = P_alpha + Q_alpha;
                    const __m256d one_over_PQalpha_sum = const_1 / PQalpha_sum;


                    /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
                    const __m256d PQ_x = P_x - _mm256_load_pd(Q.x + j);
                    const __m256d PQ_y = P_y - _mm256_load_pd(Q.y + j);
                    const __m256d PQ_z = P_z - _mm256_load_pd(Q.z + j);
                    const __m256d R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;

                    const __m256d alpha = PQalpha_mul * one_over_PQalpha_sum;   // alpha from MEST
                    const __m256d one_over_p = const_1 / P_alpha;
                    const __m256d one_over_q = const_1 / Q_alpha;
                    const __m256d one_over_2p = one_half * one_over_p;  // gets multiplied in VRR
                    const __m256d one_over_2q = one_half * one_over_q;  // gets multiplied in VRR
                    const __m256d one_over_2pq = one_half * one_over_PQalpha_sum;
                    const __m256d Q_PA_x = _mm256_load_pd(Q.PA_x + j);
                    const __m256d Q_PA_y = _mm256_load_pd(Q.PA_y + j);
                    const __m256d Q_PA_z = _mm256_load_pd(Q.PA_z + j);

                    // NOTE: Minus sign!
                    const __m256d a_over_p =  -alpha * one_over_p;     // a/p from MEST
                    const __m256d a_over_q =  -alpha * one_over_q;     // a/q from MEST


                    //////////////////////////////////////////////
                    // Boys function section
                    //////////////////////////////////////////////
                    // The parameter to the boys function
                    const __m256d F_x = R2 * alpha;


                    Boys_F_split_simd((double *)prim_ptr[0][0], L, (const double *)(&F_x));
                    const __m256d prefac = _mm256_sqrt_pd(one_over_PQalpha_sum) * P_prefac * _mm256_load_pd(Q.prefac + j);

                    for(n = 0; n <= L; n++)
                        prim_ptr[0][0][n] *= prefac;

                    //////////////////////////////////////////////
                    // Primitive integrals: Vertical recurrance
                    //////////////////////////////////////////////

                    __m256d PA[3] = {P_PA_x, P_PA_y, P_PA_z};
                    __m256d QC[3] = {Q_PA_x, Q_PA_y, Q_PA_z};

                    __m256d aop_PQ[3] = { a_over_p * PQ_x,
                                          a_over_p * PQ_y,
                                          a_over_p * PQ_z };

                    __m256d aoq_PQ[3] = { -a_over_q * PQ_x,
                                          -a_over_q * PQ_y,
                                          -a_over_q * PQ_z };

                    if( L1 >= L2 )
                    {
                        // do the first vrr, forming p_s_s_s
                        general_vrr1(1, L-1, one_over_2p, a_over_p, aop_PQ, PA,
                                     prim_ptr[0][0], NULL, prim_ptr[1][0]);

                        // now do all the rest, forming X_s_s_s
                        for(n1 = 2; n1 <= L1; n1++)
                            general_vrr1(n1, L-n1, one_over_2p, a_over_p, aop_PQ, PA,
                                         prim_ptr[n1-1][0], prim_ptr[n1-2][0], prim_ptr[n1][0]);


                        // up the ket
                        for(n1 = (L1-L2+1); n1 <= L1; n1++)
                        for(n2 = 1; n2 <= (n1-(L1-L2)); n2++)
                        {
                            if(n2 == 1)
                            {
                                general_vrr2_k(n1, 1, L2-1,
                                               one_over_2q, one_over_2pq, a_over_q, aoq_PQ, QC,
                                               prim_ptr[n1][0], NULL, prim_ptr[n1-1][0],
                                               prim_ptr[n1][1]);
                            }
                            else
                            {
                                general_vrr2_k(n1, n2, L2 - n2,
                                               one_over_2q, one_over_2pq, a_over_q, aoq_PQ, QC,
                                               prim_ptr[n1][n2-1], prim_ptr[n1][n2-2], prim_ptr[n1-1][n2-1],
                                               prim_ptr[n1][n2]);
                            }
                        }
                    }
                    else
                    {
                        // do the first vrr, forming s_s_p_s
                        general_vrr1(1, L-1, one_over_2q, a_over_q, aoq_PQ, QC,
                                     prim_ptr[0][0], NULL, prim_ptr[0][1]);

                        // now do all the rest, forming s_s_X_s
                        for(n2 = 2; n2 <= L2; n2++)
                            general_vrr1(n2, L-n2, one_over_2q, a_over_q, aoq_PQ, QC,
                                         prim_ptr[0][n2-1], prim_ptr[0][n2-2], prim_ptr[0][n2]);


                        // up the bra
                        for(n2 = (L2-L1+1); n2 <= L2; n2++)
                        for(n1 = 1; n1 <= (n2-(L2-L1)); n1++)
                        {
                            if(n1 == 1)
                            {
                                general_vrr2_i(1, n2, L1-1,
                                               one_over_2p, one_over_2pq, a_over_p, aop_PQ, PA,
                                               prim_ptr[0][n2], NULL, prim_ptr[0][n2-1],
                                               prim_ptr[1][n2]);
                            }
                            else
                            {
                                general_vrr2_i(n1, n2, L1 - n1,
                                               one_over_2p, one_over_2pq, a_over_p, aop_PQ, PA,
                                               prim_ptr[n1-1][n2], prim_ptr[n1-2][n2], prim_ptr[n1-1][n2-1],
                                               prim_ptr[n1][n2]);
                            }
                        }
                    }


                    ////////////////////////////////////
                    // Accumulate contracted integrals
                    ////////////////////////////////////
                    __m256d * result_ptr = prim_ptr[L1][L2];
                    for(np = 0; np < ncart1234; ++np)
                    {
                        const union simint_double4 tmp = (union simint_double4)result_ptr[np];
                        PRIM_PTR_INT__X_s_X_s[np] += tmp.d[0];   // first offset is always zero
                        for(n = 1; n < SIMINT_SIMD_LEN; ++n)
                            PRIM_PTR_INT__X_s_X_s[shelloffsets[n]*ncart1234+np] += tmp.d[n];
                    }

                    const int lastoffset = shelloffsets[SIMINT_SIMD_LEN-1];
                    PRIM_PTR_INT__X_s_X_s += lastoffset*ncart1234;

                }  // close loop over j
            }  // close loop over i

            //Advance to the next batch
            jstart = SIMINT_SIMD_ROUND(jend);
            abcd += nshellbatch;


            cd += nshellbatch;
        }   // close loop cdbatch

        istart = iend;
        // if this is the end of a batch in the bra part, skip the padding
        if( ((ab+1) % SIMINT_NSHELL_SIMD) == 0)
            istart = SIMINT_SIMD_ROUND(istart);

    }  // close loop over ab


    FREE(arr_stack);

    return P.nshell12 * Q.nshell12;
}


int eri_X_s_X_s(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__X_s_X_s)
{
    return eri_sharedwork_X_s_X_s(P, Q, NULL, INT__X_s_X_s);
}




int eri_sharedwork_s_X_s_X(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__s_X_s_X)
{
	// Can be accomplished by swapping some variables
	// and calling another function

	struct simint_multi_shellpair P_tmp = P;
	struct simint_multi_shellpair Q_tmp = Q;

	P_tmp.PA_x = P.PB_x;  P_tmp.PA_y = P.PB_y;  P_tmp.PA_z = P.PB_z;	
	P_tmp.PB_x = P.PA_x;  P_tmp.PB_y = P.PA_y;  P_tmp.PB_z = P.PA_z;	
	Q_tmp.PA_x = Q.PB_x;  Q_tmp.PA_y = Q.PB_y;  Q_tmp.PA_z = Q.PB_z;	
	Q_tmp.PB_x = Q.PA_x;  Q_tmp.PB_y = Q.PA_y;  Q_tmp.PB_z = Q.PA_z;	

	return eri_sharedwork_X_s_X_s(P_tmp, Q_tmp, contwork, INT__s_X_s_X);
}


int eri_s_X_s_X(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__s_X_s_X)
{
    return eri_sharedwork_s_X_s_X(P, Q, NULL, INT__s_X_s_X);
}


int eri_sharedwork_s_X_X_s(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__s_X_X_s)
{
	// Can be accomplished by swapping some variables
	// and calling another function

	struct simint_multi_shellpair P_tmp = P;

	P_tmp.PA_x = P.PB_x;  P_tmp.PA_y = P.PB_y;  P_tmp.PA_z = P.PB_z;	
	P_tmp.PB_x = P.PA_x;  P_tmp.PB_y = P.PA_y;  P_tmp.PB_z = P.PA_z;	

	return eri_sharedwork_X_s_X_s(P_tmp, Q, contwork, INT__s_X_X_s);
}


int eri_s_X_X_s(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__s_X_X_s)
{
    return eri_sharedwork_s_X_X_s(P, Q, NULL, INT__s_X_X_s);
}


int eri_sharedwork_X_s_s_X(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__X_s_s_X)
{
	// Can be accomplished by swapping some variables
	// and calling another function

	struct simint_multi_shellpair Q_tmp = Q;

	Q_tmp.PA_x = Q.PB_x;  Q_tmp.PA_y = Q.PB_y;  Q_tmp.PA_z = Q.PB_z;	
	Q_tmp.PB_x = Q.PA_x;  Q_tmp.PB_y = Q.PA_y;  Q_tmp.PB_z = Q.PA_z;	

	return eri_sharedwork_X_s_X_s(P, Q_tmp, contwork, INT__X_s_s_X);
}


int eri_X_s_s_X(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__X_s_s_X)
{
    return eri_sharedwork_X_s_s_X(P, Q, NULL, INT__X_s_s_X);
}
