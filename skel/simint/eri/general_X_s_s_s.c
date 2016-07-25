#include "simint/boys/boys_split.h"
#include "simint/eri/eri.h"
#include "simint/eri/general_eri.h"
#include <math.h>
#include <alloca.h>
#include <string.h>


int eri_sharedwork_X_s_s_s(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__X_s_s_s)
{
    // Note - In this case, L would generally equal P.am1. However,
    //        this function may be called for s_X_s_s, s_s_X_s, and s_s_X_s,
    //        So I use L in places where you would expect P.am1

    const int ncart1 = ((P.am1+1)*(P.am1+2))/2;
    const int ncart2 = ((P.am2+1)*(P.am2+2))/2;
    const int ncart3 = ((Q.am1+1)*(Q.am1+2))/2;
    const int ncart4 = ((Q.am2+1)*(Q.am2+2))/2;
    const int ncart12 = ncart1 * ncart2;
    const int ncart34 = ncart3 * ncart4;
    const int ncart1234 = ncart12 * ncart34;
    const int L = P.am1 + P.am2 + Q.am1 + Q.am2;

    memset(INT__X_s_s_s, 0, P.nshell12 * Q.nshell12 * ncart1234 * sizeof(double));

    int ab, cd, cdbatch, abcd;
    int istart, jstart;
    int iprimcd, nprim_icd, icd;
    int i, j, n;
    int np;

    const __m256d const_1 = _mm256_set1_pd(1);
    const __m256d one_half = _mm256_set1_pd(0.5);

    // Size of arrays needed in the stack
    size_t arr_stack_size = 0;
    for(n = 0; n <= L; n++)
        arr_stack_size += (( (n+1)*(n+2) )/2) * (L-n+1);

    __m256d * arr_stack = (__m256d *)(ALLOC(arr_stack_size * sizeof(__m256d)));
    __m256d * prim_ptr[L+1];

    __m256d * pos = arr_stack;
    for(n = 0; n <= L; n++)
    {
        prim_ptr[n] = pos;
        pos += ( ((n+1)*(n+2))/2 ) * (L-n+1);
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
                double * restrict PRIM_PTR_INT__X_s_s_s = INT__X_s_s_s + abcd * ncart1234;

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
                            PRIM_PTR_INT__X_s_s_s += ncart1234;
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
                    const __m256d one_over_2p = one_half * one_over_p;  // gets multiplied in VRR

                    // NOTE: Minus sign!
                    const __m256d a_over_p =  -alpha * one_over_p;     // a/p from MEST


                    //////////////////////////////////////////////
                    // Boys function section
                    //////////////////////////////////////////////
                    // The parameter to the boys function
                    const __m256d F_x = R2 * alpha;


                    Boys_F_split_simd((double *)prim_ptr[0], L, (const double *)(&F_x));
                    const __m256d prefac = _mm256_sqrt_pd(one_over_PQalpha_sum) * P_prefac * _mm256_load_pd(Q.prefac + j);

                    for(n = 0; n <= L; n++)
                        prim_ptr[0][n] *= prefac;

                    //////////////////////////////////////////////
                    // Primitive integrals: Vertical recurrance
                    //////////////////////////////////////////////

                    __m256d PA[3] = {P_PA_x, P_PA_y, P_PA_z};
                    __m256d PQ[3] = {PQ_x, PQ_y, PQ_z};

                    // do the first vrr, forming p_s_s_s
                    general_vrr1(1, L-1, one_over_2p, a_over_p, PA, PQ,
                                 prim_ptr[0], NULL, prim_ptr[1]);

                    // now do all the rest
                    for(n = 2; n <= L; n++)
                        general_vrr1(n, L-n, one_over_2p, a_over_p, PA, PQ,
                                     prim_ptr[n-1], prim_ptr[n-2], prim_ptr[n]);
                    
                                 

                    ////////////////////////////////////
                    // Accumulate contracted integrals
                    ////////////////////////////////////
                    __m256d * result_ptr = prim_ptr[L];
                    for(np = 0; np < ncart1234; ++np)
                    {
                        const union double4 tmp = (union double4)result_ptr[np];
                        PRIM_PTR_INT__X_s_s_s[np] += tmp.d[0];   // first offset is always zero
                        for(n = 1; n < SIMINT_SIMD_LEN; ++n)
                            PRIM_PTR_INT__X_s_s_s[shelloffsets[n]*ncart1234+np] += tmp.d[n];
                    }

                    const int lastoffset = shelloffsets[SIMINT_SIMD_LEN-1];
                    PRIM_PTR_INT__X_s_s_s += lastoffset*ncart1234; 

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



int eri_X_s_s_s(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__X_s_s_s)
{
    return eri_sharedwork_X_s_s_s(P, Q, NULL, INT__X_s_s_s);
}



int eri_sharedwork_s_X_s_s(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__s_X_s_s)
{
	// Can be accomplished by swapping some variables
	// and calling another function

	struct simint_multi_shellpair P_tmp = P;

	P_tmp.PA_x = P.PB_x;  P_tmp.PA_y = P.PB_y;  P_tmp.PA_z = P.PB_z;	
	P_tmp.PB_x = P.PA_x;  P_tmp.PB_y = P.PA_y;  P_tmp.PB_z = P.PA_z;	

	return eri_sharedwork_X_s_s_s(P_tmp, Q, contwork, INT__s_X_s_s);
}



int eri_s_X_s_s(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__s_X_s_s)
{
    return eri_sharedwork_s_X_s_s(P, Q, NULL, INT__s_X_s_s);
}



int eri_sharedwork_s_s_X_s(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__s_s_X_s)
{
	return eri_sharedwork_X_s_s_s(Q, P, contwork, INT__s_s_X_s);
}



int eri_s_s_X_s(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__s_s_X_s)
{
    return eri_sharedwork_s_s_X_s(P, Q, NULL, INT__s_s_X_s);
}



int eri_sharedwork_s_s_s_X(struct simint_multi_shellpair const P,
                           struct simint_multi_shellpair const Q,
                           double * const restrict contwork,
                           double * const restrict INT__s_s_s_X)
{
	// Can be accomplished by swapping some variables
	// and calling another function

	return eri_sharedwork_s_X_s_s(Q, P, contwork, INT__s_s_s_X);
}



int eri_s_s_s_X(struct simint_multi_shellpair const P,
                struct simint_multi_shellpair const Q,
                double * const restrict INT__s_s_s_X)
{
    return eri_sharedwork_s_s_s_X(P, Q, NULL, INT__s_s_s_X);
}

