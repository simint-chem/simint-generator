#include "simint/boys/boys.h"
#include "simint/ostei/gen/ostei_generated.h"
#include "simint/vectorization/vectorization.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>


int ostei_s_s_s_s(struct simint_multi_shellpair const P,
                  struct simint_multi_shellpair const Q,
                  double screen_tol,
                  double * const restrict work,
                  double * const restrict INT__s_s_s_s)
{

    SIMINT_ASSUME_ALIGN_DBL(work);
    SIMINT_ASSUME_ALIGN_DBL(INT__s_s_s_s);
    memset(INT__s_s_s_s, 0, P.nshell12_clip * Q.nshell12_clip * 1 * sizeof(double));

    int ab, cd, abcd;
    int istart, jstart;
    int iprimcd, nprim_icd, icd;
    const int check_screen = (screen_tol > 0.0);
    int i, j;
    int n;
    int not_screened;

    // partition workspace
    SIMINT_DBLTYPE * const primwork = (SIMINT_DBLTYPE *)(work + SIMINT_NSHELL_SIMD*0);
    SIMINT_DBLTYPE * const restrict PRIM_INT__s_s_s_s = primwork + 0;
    double * const hrrwork = (double *)(primwork + 1);


    // Create constants
    const SIMINT_DBLTYPE const_1 = SIMINT_DBLSET1(1);
    const SIMINT_DBLTYPE one_half = SIMINT_DBLSET1(0.5);
    
    #if defined SIMINT_AVX512 || defined SIMINT_MICAVX512
    // Buffer for contract_all()
    SIMINT_DBLTYPE ca_buf[8];
    double ca_res[8];
    int ca_res_idx[8];
    int ca_cnt = 0;
    #endif
    
    // Create offset buffer
    int ivec;
    const int TopAM_size = 1;
    int n_info_vector = SIMINT_NSHELL_SIMD * 4;
    int *offset_info  = (int*) malloc(sizeof(int) * (SIMINT_SIMD_LEN + 1 + TopAM_size) * n_info_vector);
    assert(offset_info != NULL);

    #ifdef SIMINT_PRIM_SCREEN_STAT
    int calc_nprim = 0, skip_nprim = 0, calc_nvec = 0, skip_nvec = 0;
    #endif
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
            
            // Check if the offset buffer is large enough
            int j_vec = (jend - jstart + SIMINT_SIMD_LEN - 1) / SIMINT_SIMD_LEN;
            if (j_vec > n_info_vector)
            {
                n_info_vector = j_vec;
                if (offset_info != NULL) free(offset_info);
                offset_info  = (int*) malloc(sizeof(int) * (SIMINT_SIMD_LEN + 1 + TopAM_size) * n_info_vector);
            }
            // Calculate all the shell offsets in the j-loop
            // These are the offset from the shell pointed to by cd for each element
            ivec = 0;
            iprimcd = 0;
            nprim_icd = Q.nprim12[cd];
            icd = 0;
            for (j = jstart; j < jend; j += SIMINT_SIMD_LEN)
            {
                int *shelloffsets = offset_info + (SIMINT_SIMD_LEN + 1 + TopAM_size) * ivec;
                shelloffsets[0] = 0;
                shelloffsets[SIMINT_SIMD_LEN] = 0;  // for lastoffset 
                shelloffsets[SIMINT_SIMD_LEN + 1 + 0] = 0;
                ivec++;

                if((iprimcd + SIMINT_SIMD_LEN) >= nprim_icd)
                {
                    // Handle if the first element of the vector is a new shell
                    if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))
                    {
                        nprim_icd += Q.nprim12[cd + (++icd)];
                        shelloffsets[SIMINT_SIMD_LEN + 1 + 0] += 1;    // for PRIM_PTR_INT__s_s_s_s
                    }
                    iprimcd++;
                    for(n = 1; n < SIMINT_SIMD_LEN; ++n)
                    {
                        if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))
                        {
                            shelloffsets[n] = shelloffsets[n-1] + 1;
                            shelloffsets[SIMINT_SIMD_LEN]++;
                            nprim_icd += Q.nprim12[cd + (++icd)];
                        }
                        else
                            shelloffsets[n] = shelloffsets[n-1];
                        iprimcd++;
                    }
                }
                else iprimcd += SIMINT_SIMD_LEN; 
            }


            for(i = istart; i < iend; ++i)
            {
                SIMINT_DBLTYPE bra_screen_max;  // only used if check_screen
		bra_screen_max = SIMINT_DBLSET1(0.);

                if(check_screen)
                {
                    // Skip this whole thing if always insignificant
                    if((P.screen[i] * Q.screen_max) < screen_tol)
                    {
                        #ifdef SIMINT_PRIM_SCREEN_STAT
                        int j_len = jend - jstart;
                        skip_nprim += j_len;
                        skip_nvec  += (j_len + SIMINT_SIMD_LEN - 1) / SIMINT_SIMD_LEN;
                        #endif
                        continue;
                    }
                    bra_screen_max = SIMINT_DBLSET1(P.screen[i]);
                }

                icd = 0;
                iprimcd = 0;
                nprim_icd = Q.nprim12[cd];
                double * restrict PRIM_PTR_INT__s_s_s_s = INT__s_s_s_s + abcd * 1;



                // Load these one per loop over i
                const SIMINT_DBLTYPE P_alpha = SIMINT_DBLSET1(P.alpha[i]);
                const SIMINT_DBLTYPE P_prefac = SIMINT_DBLSET1(P.prefac[i]);
                const SIMINT_DBLTYPE Pxyz[3] = { SIMINT_DBLSET1(P.x[i]), SIMINT_DBLSET1(P.y[i]), SIMINT_DBLSET1(P.z[i]) };

                ivec = 0;
                for(j = jstart; j < jend; j += SIMINT_SIMD_LEN)
                {
                    // calculate the shell offsets
                    // these are the offset from the shell pointed to by cd
                    // for each element
                    //int shelloffsets[SIMINT_SIMD_LEN] = {0};
                    //int lastoffset = 0;
                    const int nlane = ( ((j + SIMINT_SIMD_LEN) < jend) ? SIMINT_SIMD_LEN : (jend - j));

                    /*
                    if((iprimcd + SIMINT_SIMD_LEN) >= nprim_icd)
                    {
                        // Handle if the first element of the vector is a new shell
                        if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))
                        {
                            nprim_icd += Q.nprim12[cd + (++icd)];
                            PRIM_PTR_INT__s_s_s_s += 1;
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
                    else iprimcd += SIMINT_SIMD_LEN;
                    */
                    
                    int *shelloffsets = offset_info + (SIMINT_SIMD_LEN + 1 + TopAM_size) * ivec;
                    int lastoffset = shelloffsets[SIMINT_SIMD_LEN];
                    if (shelloffsets[SIMINT_SIMD_LEN + 1]) PRIM_PTR_INT__s_s_s_s += 1;
                    ivec++;

                    // Do we have to compute this vector (or has it been screened out)?
                    // (not_screened != 0 means we have to do this vector)
                    SIMINT_DBLTYPE prim_screen_res = SIMINT_DBLSET1(0.);
                    if(check_screen)
                    {
		      if(jend -j < SIMINT_SIMD_LEN ){
			//initialize remainders when vlen < SIMD_LEN			
			for(n = jend; n < j+SIMINT_SIMD_LEN; n++)
			  Q.screen[n] = 0.;
		      }
		      prim_screen_res = SIMINT_MUL(bra_screen_max, SIMINT_DBLLOAD(Q.screen, j));
                        const double vmax = vector_max(prim_screen_res);
                        if(vmax < screen_tol)
                        {
                            #ifdef SIMINT_PRIM_SCREEN_STAT
                            skip_nvec++;
                            skip_nprim += SIMINT_SIMD_LEN;
                            #endif
                            PRIM_PTR_INT__s_s_s_s += lastoffset*1;    
                            continue;
                        }
                    }
                    #ifdef SIMINT_PRIM_SCREEN_STAT
                    calc_nvec++;
                    int calc_nprim_in_vec = count_prim_screen_survival(prim_screen_res, screen_tol);
                    calc_nprim += calc_nprim_in_vec;
                    skip_nprim += (SIMINT_SIMD_LEN - calc_nprim_in_vec);
                    #endif

                    const SIMINT_DBLTYPE Q_alpha = SIMINT_DBLLOAD(Q.alpha, j);
                    const SIMINT_DBLTYPE PQalpha_mul = SIMINT_MUL(P_alpha, Q_alpha);
                    const SIMINT_DBLTYPE PQalpha_sum = SIMINT_ADD(P_alpha, Q_alpha);
                    const SIMINT_DBLTYPE one_over_PQalpha_sum = SIMINT_DIV(const_1, PQalpha_sum);


                    /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */
                    SIMINT_DBLTYPE PQ[3];
                    PQ[0] = SIMINT_SUB(Pxyz[0], SIMINT_DBLLOAD(Q.x, j));
                    PQ[1] = SIMINT_SUB(Pxyz[1], SIMINT_DBLLOAD(Q.y, j));
                    PQ[2] = SIMINT_SUB(Pxyz[2], SIMINT_DBLLOAD(Q.z, j));
                    SIMINT_DBLTYPE R2 = SIMINT_MUL(PQ[0], PQ[0]);
                    R2 = SIMINT_FMADD(PQ[1], PQ[1], R2);
                    R2 = SIMINT_FMADD(PQ[2], PQ[2], R2);

                    const SIMINT_DBLTYPE alpha = SIMINT_MUL(PQalpha_mul, one_over_PQalpha_sum); // alpha from MEST
                    const SIMINT_DBLTYPE one_over_p = SIMINT_DIV(const_1, P_alpha);
                    const SIMINT_DBLTYPE one_over_q = SIMINT_DIV(const_1, Q_alpha);
                    const SIMINT_DBLTYPE one_over_2p = SIMINT_MUL(one_half, one_over_p);
                    const SIMINT_DBLTYPE one_over_2q = SIMINT_MUL(one_half, one_over_q);
                    const SIMINT_DBLTYPE one_over_2pq = SIMINT_MUL(one_half, one_over_PQalpha_sum);


                    //////////////////////////////////////////////
                    // Fjt function section
                    // Maximum v value: 0
                    //////////////////////////////////////////////
                    // The parameter to the Fjt function
                    const SIMINT_DBLTYPE F_x = SIMINT_MUL(R2, alpha);


                    const SIMINT_DBLTYPE Q_prefac = mask_load(nlane, Q.prefac + j);


                    boys_F_split(PRIM_INT__s_s_s_s, F_x, 0);
                    SIMINT_DBLTYPE prefac = SIMINT_SQRT(one_over_PQalpha_sum);
                    prefac = SIMINT_MUL(SIMINT_MUL(P_prefac, Q_prefac), prefac);
                    for(n = 0; n <= 0; n++)
                        PRIM_INT__s_s_s_s[n] = SIMINT_MUL(PRIM_INT__s_s_s_s[n], prefac);


                    ////////////////////////////////////
                    // Accumulate contracted integrals
                    ////////////////////////////////////
                    if(lastoffset == 0)
                    {
                        #if defined SIMINT_AVX512 || defined SIMINT_MICAVX512
                        int new_idx = PRIM_PTR_INT__s_s_s_s - INT__s_s_s_s;
                        if ((ca_cnt > 0) && (new_idx == ca_res_idx[ca_cnt - 1]))
                        {
                            ca_buf[ca_cnt - 1] = _mm512_add_pd(ca_buf[ca_cnt - 1], PRIM_INT__s_s_s_s[0]);
                        } else {
                            ca_buf[ca_cnt]     = PRIM_INT__s_s_s_s[0];
                            ca_res_idx[ca_cnt] = PRIM_PTR_INT__s_s_s_s - INT__s_s_s_s;
                            ca_cnt++;
                            if (ca_cnt == 8)
                            {
                                memset(ca_res, 0, sizeof(double) * 8);
                                contract_all(8, ca_buf, ca_res);
                                ca_cnt = 0;
                                
                                //#pragma simd  <--- should not use this, since for same i 
                                // different j PRIM_PTR_INT__s_s_s_s may repeat
                                for (int jj = 0; jj < 8; jj++)
                                    INT__s_s_s_s[ca_res_idx[jj]] += ca_res[jj];
                            }
                        }
                        #else
                        contract_all(1, PRIM_INT__s_s_s_s, PRIM_PTR_INT__s_s_s_s);
                        #endif
                    }
                    else
                    {
                        contract(1, shelloffsets, PRIM_INT__s_s_s_s, PRIM_PTR_INT__s_s_s_s);
                        PRIM_PTR_INT__s_s_s_s += lastoffset*1;
                    }

                }  // close loop over j
            }  // close loop over i
            
            //Advance to the next batch
            jstart = SIMINT_SIMD_ROUND(jend);
            abcd += nshellbatch;
            
        }   // close loop cdbatch

        istart = iend;
    }  // close loop over ab
    
    #if defined SIMINT_AVX512 || defined SIMINT_MICAVX512
    for (int jj = 0; jj < ca_cnt; jj++)
        INT__s_s_s_s[ca_res_idx[jj]] += _mm512_reduce_add_pd(ca_buf[jj]);
    #endif 
    
    if (offset_info != NULL) free(offset_info);
    
    #ifdef SIMINT_PRIM_SCREEN_STAT
    double *eri_res_end_pos = INT__s_s_s_s + 1 * P.nshell12_clip * Q.nshell12_clip;
    eri_res_end_pos[0] = (double) calc_nprim;
    eri_res_end_pos[1] = (double) skip_nprim;
    eri_res_end_pos[2] = (double) calc_nvec;
    eri_res_end_pos[3] = (double) skip_nvec;
    #endif

    return P.nshell12_clip * Q.nshell12_clip;
}

