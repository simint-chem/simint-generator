#include "simint/ostei/ostei_contract.h"
#include "simint/vectorization/vectorization.h"
#include "simint/ostei/ostei_general.h"
#include "simint/ostei/recur_lookup.h"

#include <math.h>
#include <string.h>

int ostei_general_sharedwork(struct simint_multi_shellpair const P,
                             struct simint_multi_shellpair const Q,
                             Fjtfunc fjt,
                             double screen_tol,
                             double * const restrict contwork,
                             double * const restrict output)
{
    int ab, cd, abcd;
    int istart, jstart;
    int iprimcd, nprim_icd, icd;
    const int check_screen = (screen_tol > 0.0);
    int i, j;
    int li, lk;
    int n;
    int not_screened;

    const int ncart1 = ((P.am1+1)*(P.am1+2))/2;
    const int ncart2 = ((P.am2+1)*(P.am2+2))/2;
    const int ncart3 = ((Q.am1+1)*(Q.am1+2))/2;
    const int ncart4 = ((Q.am2+1)*(Q.am2+2))/2;
    const int ncart12 = ncart1 * ncart2;
    const int ncart34 = ncart3 * ncart4;
    const int ncart1234 = ncart12 * ncart34;
    const int Lbra = P.am1 + P.am2;
    const int Lket = Q.am1 + Q.am2;
    const int ncart_lbra = ((Lbra+1)*(Lbra+2))/2;
    const int ncart_lket = ((Lket+1)*(Lket+2))/2;
    const int L = Lbra + Lket;

    // pointer arrays
    double * vrr_ptrs[P.am1+P.am2+1][Q.am1+Q.am2+1];

    // how much workspace do we need for VRR?
    size_t worksize = 0;
    for(li = 0; li <= Lbra; li++)
    {
        const int ncart_li = (((li+1)*(li+2))/2);
        worksize += ncart_li * (L-li+1); 

        for(lj = 1; lj <= Lket; lj++)
            worksize += (((lj+1)*(lj+2))/2) * ncart_li * (Lket-lj+1);
    }
        

    __m256d * work = SIMINT_ALLOC(worksize * SIMINT_SIMD_LEN * sizeof(double));
    __m256d * work_ptr = work;

    // now set up the VRR pointers
    for(li = 0; li <= Lbra; li++)
    {
        const int ncart_li = (((li+1)*(li+2))/2);
        vrr_ptrs[li][0] = work_ptr;

        work_ptr += ncart_li * (L-li+1);

        for(lj = 1; lj <= Lket; lj++)
        {
            vrr_ptrs[li][lj] = work_ptr;
            work_ptr += (((lj+1)*(lj+2))/2) * ncart_li * (Lket-lj+1);
        }
    }

    // zero the output and the work array
    memset(output, 0, P.nshell12_clip * Q.nshell12_clip * ncart1234 * sizeof(double));

    // Some constants
    const __m256d const_1 = _mm256_set1_pd(1);
    const __m256d one_half = _mm256_set1_pd(0.5);


    ////////////////////////////////////////
    // Loop over shells and primitives
    ////////////////////////////////////////

    abcd = 0;
    istart = 0;
    for(ab = 0; ab < P.nshell12_clip; ++ab)
    {
        const int iend = istart + P.nprim12[ab];

        cd = 0;
        jstart = 0;

        for(cd = 0; cd < Q.nshell12_clip; cd += SIMINT_NSHELL_SIMD)
        {
            const int nshellbatch = ((cd + SIMINT_NSHELL_SIMD) > Q.nshell12_clip) ? Q.nshell12_clip - cd : SIMINT_NSHELL_SIMD;
            int jend = jstart;
            for(i = 0; i < nshellbatch; i++)
                jend += Q.nprim12[cd+i];


            for(i = istart; i < iend; ++i)
            {

                // Skip this whole thing if always insignificant
                if(check_screen && (P.screen[i] * Q.screen_max) < screen_tol)
                    continue;

                icd = 0;
                iprimcd = 0;
                nprim_icd = Q.nprim12[cd];
                double * restrict PTR_output = INT__output + abcd * 1;



                // Load these one per loop over i
                const __m256d P_alpha = _mm256_set1_pd(P.alpha[i]);
                const __m256d P_prefac = _mm256_set1_pd(P.prefac[i]);
                const __m256d Pxyz[3] = { _mm256_set1_pd(P.x[i]), _mm256_set1_pd(P.y[i]), _mm256_set1_pd(P.z[i]) };
                const __m256d bra_screen_max = _mm256_set1_pd(P.screen[i]);


                for(j = jstart; j < jend; j += SIMINT_SIMD_LEN)
                {
                    // calculate the shell offsets
                    // these are the offset from the shell pointed to by cd
                    // for each element
                    int shelloffsets[SIMINT_SIMD_LEN] = {0};
                    int lastoffset = 0;
                    const int nlane = ( ((j + SIMINT_SIMD_LEN) < jend) ? SIMINT_SIMD_LEN : (jend - j));

                    if((iprimcd + SIMINT_SIMD_LEN) >= nprim_icd)
                    {
                        // Handle if the first element of the vector is a new shell
                        if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))
                        {
                            nprim_icd += Q.nprim12[cd + (++icd)];
                            PTR__output += ncart1234;
                        }
                        iprimcd++;
                        for(n = 1; n < SIMINT_SIMD_LEN; ++n)
                        {
                            if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))
                            {
                                shelloffsets[n] = shelloffsets[n-1] + 1;
                                lastoffset++;
                                nprim_icd += Q.nprim12[cd + (++icd)];
                            }
                            else
                                shelloffsets[n] = shelloffsets[n-1];
                            iprimcd++;
                        }
                    }
                    else
                        iprimcd += SIMINT_SIMD_LEN;

                    // Do we have to compute this vector (or has it been screened out)?
                    // (not_screened != 0 means we have to do this vector)
                    if(check_screen)
                    {
                        union simint_double4 screen_max = (union simint_double4)(bra_screen_max * _mm256_load_pd(Q.screen + j));
                        not_screened = 0;
                        for(n = 0; n < SIMINT_SIMD_LEN; n++)
                            not_screened = ( screen_max.d[n] >= screen_tol ? 1 : not_screened );
                        if(not_screened == 0)
                        {
                            PTR_output += lastoffset*ncart1234;
                            continue;
                        }
                    }

                    const __m256d Q_alpha = _mm256_load_pd(Q.alpha + j);
                    const __m256d PQalpha_mul = P_alpha * Q_alpha;
                    const __m256d PQalpha_sum = P_alpha + Q_alpha;
                    const __m256d one_over_PQalpha_sum = const_1 / PQalpha_sum;


                    /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
                    const __m256d PQ[3] = { Pxyz[0] - _mm256_load_pd(Q.x + j), Pxyz[1] - _mm256_load_pd(Q.y + j), Pxyz[2] - _mm256_load_pd(Q.z + j) };
                    const __m256d R2 = PQ[0]*PQ[0] + PQ[1]*PQ[1] + PQ[2]*PQ[2];

                    const __m256d alpha = PQalpha_mul * one_over_PQalpha_sum;   // alpha from MEST
                    const __m256d one_over_p = const_1 / P_alpha;
                    const __m256d one_over_q = const_1 / Q_alpha;
                    const __m256d one_over_2p = one_half * one_over_p;  // gets multiplied in VRR
                    const __m256d one_over_2q = one_half * one_over_q;  // gets multiplied in VRR
                    const __m256d one_over_2pq = one_half * one_over_PQalpha_sum;


                    //////////////////////////////////////////////
                    // Fjt function section
                    // Maximum v value: 0
                    //////////////////////////////////////////////
                    // The parameter to the Fjt function
                    const __m256d F_x = R2 * alpha;


                    union simint_double4 Q_prefac_u = (union simint_double4)_mm256_load_pd(Q.prefac + j);
                    for(n = nlane; n < SIMINT_SIMD_LEN; n++)
                        Q_prefac_u.d[n] = 0.0;
                    const __m256d Q_prefac = Q_prefac_u.d_256;


                    fjt((double *)vrr_work[0][0], 0, (const double *)(&F_x));
                    const __m256d prefac = _mm256_sqrt_pd(one_over_PQalpha_sum) * P_prefac * Q_prefac;
                    for(n = 0; n <= 0; n++)
                        vrr_work[0][0][n] *= prefac;

                    // Do VRR
                    for(li = 0; li <= Lbra; li++)
                    {
                        for(lj = 1; lj <= Lket; lj++)
                        {
                            // do something
                        }
                    }



                    ////////////////////////////////////
                    // Accumulate contracted integrals
                    ////////////////////////////////////
                    if(lastoffset == 0)
                    {
                        contract_all(ncart1234, vrr_work[Lbra][Lket], PTR_output);
                    }
                    else
                    {
                        contract(ncart1234, shelloffsets, vrr_work[Lbra][Lket], PTR_output);
                        PTR_output += lastoffset*ncart1234;
                    }

                }  // close loop over j
            }  // close loop over i
            
            //Advance to the next batch
            jstart = jend;
            abcd += nshellbatch;
            
        }   // close loop cdbatch

        istart = iend;
    }  // close loop over ab

#endif
    return P.nshell12_clip * Q.nshell12_clip;
}

int ostei_general(struct simint_multi_shellpair const P,
                  struct simint_multi_shellpair const Q,
                  Fjtfunc fjt,
                  double screen_tol,
                  double * const restrict INT__s_d_s_p)
{
    return 0;
}

