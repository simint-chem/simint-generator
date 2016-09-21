#include <cstring> // memset
#include "simint/constants.h"
#include "test/Libint2.hpp"
#include "simint/boys/boys_split.h"
#include "test/Common.hpp"
#include "simint/vectorization/vectorization.h"


// Disable intel warnings
// 193 : zero used for undefined preprocessing identifier "LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18
//       Not much I can do about that 
#ifdef __INTEL_COMPILER
    #pragma warning(disable:193)
#endif


#define LIBINT2_SIMD_LEN SIMINT_SIMD_LEN
#define LIBINT2_LOAD(x) LIBINT2_REALTYPE(*(x+3), *(x+2), *(x+1), *(x+0))


static_assert(false, "THIS NEEDS TO BE UPDATED TO TAKE INTO ACCOUNT P,Q PERMUTATIONS. BENCHMARK_LIBINT2.CPP AS WELL (IE, NSHELL FOR P AND Q)")

Libint2_ERI::Libint2_ERI(int maxam, size_t maxnprim)
{
    size_t size = maxnprim * maxnprim * maxnprim * maxnprim;
    erival_.resize(size);
    LIBINT2_PREFIXED_NAME(libint2_init_eri)(erival_.data(), maxam, 0); 

    // temporary workspace
    // note the LIBINT2_MAX_AM_ERI + 1
    //      the +1 is needed since for a given AM, we need to store
    //      AM+1 values (since we store F[0])
    worksize_ = LIBINT2_SIMD_LEN * size * (23 + LIBINT2_MAX_AM_ERI + 1);
    work_ = new double[worksize_];

    tmp_AB_x_  = work_ +  0; 
    tmp_AB_y_  = work_ +  1*LIBINT2_SIMD_LEN*size; 
    tmp_AB_z_  = work_ +  2*LIBINT2_SIMD_LEN*size; 
               
    tmp_CD_x_  = work_ +  3*LIBINT2_SIMD_LEN*size; 
    tmp_CD_y_  = work_ +  4*LIBINT2_SIMD_LEN*size; 
    tmp_CD_z_  = work_ +  5*LIBINT2_SIMD_LEN*size; 
               
    tmp_PA_x_  = work_ +  6*LIBINT2_SIMD_LEN*size; 
    tmp_PA_y_  = work_ +  7*LIBINT2_SIMD_LEN*size; 
    tmp_PA_z_  = work_ +  8*LIBINT2_SIMD_LEN*size; 
               
    tmp_QC_x_  = work_ +  9*LIBINT2_SIMD_LEN*size; 
    tmp_QC_y_  = work_ + 10*LIBINT2_SIMD_LEN*size; 
    tmp_QC_z_  = work_ + 11*LIBINT2_SIMD_LEN*size; 
               
    tmp_WP_x_  = work_ + 12*LIBINT2_SIMD_LEN*size; 
    tmp_WP_y_  = work_ + 13*LIBINT2_SIMD_LEN*size; 
    tmp_WP_z_  = work_ + 14*LIBINT2_SIMD_LEN*size; 
               
    tmp_WQ_x_  = work_ + 15*LIBINT2_SIMD_LEN*size; 
    tmp_WQ_y_  = work_ + 16*LIBINT2_SIMD_LEN*size; 
    tmp_WQ_z_  = work_ + 17*LIBINT2_SIMD_LEN*size; 

    tmp_oo2z_  = work_ + 18*LIBINT2_SIMD_LEN*size;
    tmp_oo2e_  = work_ + 19*LIBINT2_SIMD_LEN*size;
    tmp_oo2ze_ = work_ + 20*LIBINT2_SIMD_LEN*size;
    tmp_roz_   = work_ + 21*LIBINT2_SIMD_LEN*size;
    tmp_roe_   = work_ + 22*LIBINT2_SIMD_LEN*size;

    tmp_vecF_  = work_ + 23*LIBINT2_SIMD_LEN*size;
}


Libint2_ERI::~Libint2_ERI()
{
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(erival_.data());

    delete [] work_;
}




TimerType Libint2_ERI::Integrals(struct simint_multi_shellpair P,
                                 struct simint_multi_shellpair Q,
                                 double * integrals)
{
    Libint_eri_t * erival = erival_.data();

    // according to libint, l(c)+l(d) >= l(a) + l(b)
    // so we switch P and Q, then permute
    const int M = P.am1 + P.am2 + Q.am1 + Q.am2;
    double F[M+1];


    size_t ncart1234 = NCART(P.am1) * NCART(P.am2) * NCART(Q.am1) * NCART(Q.am2);
    memset(integrals, 0, P.nshell12 * Q.nshell12 * ncart1234);

    // timing
    TimerType totaltime = 0;
    TimerType ticks0, ticks1;

    int istart = 0;
    for(int a = 0, ab = 0; a < Q.nshell1; a++)
    for(int b = 0;         b < Q.nshell2; b++, ab++)
    {
        int iend = istart + Q.nprim12[ab];

        int jstart = 0;

        for(int cd = 0; cd < P.nshell12; cd += LIBINT2_SIMD_LEN)
        {
            const int simdlen = std::min(LIBINT2_SIMD_LEN, P.nshell12 - cd);
            int maxnprim = 0;

            std::fill(work_, work_ + worksize_, 0.0);

            for(int scd = 0; scd < simdlen; scd++)
            {
                int jend = jstart + P.nprim12[cd+scd];
                int nprim = 0;

                for(int i = istart; i < iend; i++)
                for(int j = jstart; j < jend; j++)
                {
                    const int idx = nprim*LIBINT2_SIMD_LEN + scd;

                    const double PQalpha_sum = P.alpha[j] + Q.alpha[i];
                    const double PQalpha_mul = P.alpha[j] * Q.alpha[i];
                    const double one_over_PQalpha_sum = 1.0 / PQalpha_sum;
                    const double one_over_Palpha = 1.0 / P.alpha[j];
                    const double one_over_Qalpha = 1.0 / Q.alpha[i];

                    tmp_AB_x_[idx] = Q.AB_x[ab];
                    tmp_AB_y_[idx] = Q.AB_y[ab];
                    tmp_AB_z_[idx] = Q.AB_z[ab];

                    tmp_CD_x_[idx] = P.AB_x[cd+scd];
                    tmp_CD_y_[idx] = P.AB_y[cd+scd];
                    tmp_CD_z_[idx] = P.AB_z[cd+scd];

                    tmp_PA_x_[idx] = Q.PA_x[i];
                    tmp_PA_y_[idx] = Q.PA_y[i];
                    tmp_PA_z_[idx] = Q.PA_z[i];

                    tmp_QC_x_[idx] = P.PA_x[j];
                    tmp_QC_y_[idx] = P.PA_y[j];
                    tmp_QC_z_[idx] = P.PA_z[j];

                    double W[3];
                    W[0] = (P.alpha[j] * P.x[j] + Q.alpha[i] * Q.x[i]) * one_over_PQalpha_sum;
                    W[1] = (P.alpha[j] * P.y[j] + Q.alpha[i] * Q.y[i]) * one_over_PQalpha_sum;
                    W[2] = (P.alpha[j] * P.z[j] + Q.alpha[i] * Q.z[i]) * one_over_PQalpha_sum;
                    
                    tmp_WP_x_[idx] = W[0] - Q.x[i]; 
                    tmp_WP_y_[idx] = W[1] - Q.y[i]; 
                    tmp_WP_z_[idx] = W[2] - Q.z[i]; 

                    tmp_WQ_x_[idx] = W[0] - P.x[j]; 
                    tmp_WQ_y_[idx] = W[1] - P.y[j]; 
                    tmp_WQ_z_[idx] = W[2] - P.z[j]; 

                    double rho = PQalpha_mul * one_over_PQalpha_sum;
                    tmp_oo2z_[idx] = 0.5 * one_over_Qalpha; 
                    tmp_oo2e_[idx] = 0.5 * one_over_Palpha; 
                    tmp_oo2ze_[idx] = 0.5 * one_over_PQalpha_sum;
                    tmp_roz_[idx] = rho * one_over_Qalpha; 
                    tmp_roe_[idx] = rho * one_over_Palpha; 

                    double PQ2 = 0.0;
                    PQ2 += (P.x[j] - Q.x[i])*(P.x[j] - Q.x[i]);
                    PQ2 += (P.y[j] - Q.y[i])*(P.y[j] - Q.y[i]);
                    PQ2 += (P.z[j] - Q.z[i])*(P.z[j] - Q.z[i]);
                    double T = rho * PQ2; 


                    // calculate the boys function
                    Boys_F_split(F, M, T); 
                    double scale = sqrt(one_over_PQalpha_sum) * P.prefac[j] * Q.prefac[i];


                    // expand and scale
                    //      the +1 is needed since for a given AM, we need to store
                    //      AM+1 values (since we store F[0])
                    for(int fi = 0; fi <= M; fi++)
                        tmp_vecF_[nprim*LIBINT2_SIMD_LEN*(M+1) + fi*LIBINT2_SIMD_LEN + scd] = F[fi] * scale;

                    nprim++;
                }

                // handle padding in my shell pair
                jstart = jend;
                if(((cd+scd+1) % SIMINT_NSHELL_SIMD) == 0)
                    jstart = SIMINT_SIMD_ROUND(jstart);

                maxnprim = std::max(nprim, maxnprim);
            }

            erival[0].contrdepth = maxnprim;

            // now fill libint
            for(int nprim = 0; nprim < maxnprim; nprim++)
            {
                const int idx = nprim * LIBINT2_SIMD_LEN;

                erival[nprim].AB_x[0] = LIBINT2_LOAD(tmp_AB_x_ + idx);
                erival[nprim].AB_y[0] = LIBINT2_LOAD(tmp_AB_y_ + idx);
                erival[nprim].AB_z[0] = LIBINT2_LOAD(tmp_AB_z_ + idx);

                erival[nprim].CD_x[0] = LIBINT2_LOAD(tmp_CD_x_ + idx);
                erival[nprim].CD_y[0] = LIBINT2_LOAD(tmp_CD_y_ + idx);
                erival[nprim].CD_z[0] = LIBINT2_LOAD(tmp_CD_z_ + idx);

                erival[nprim].PA_x[0] = LIBINT2_LOAD(tmp_PA_x_ + idx);
                erival[nprim].PA_y[0] = LIBINT2_LOAD(tmp_PA_y_ + idx);
                erival[nprim].PA_z[0] = LIBINT2_LOAD(tmp_PA_z_ + idx);

                erival[nprim].QC_x[0] = LIBINT2_LOAD(tmp_QC_x_ + idx);
                erival[nprim].QC_y[0] = LIBINT2_LOAD(tmp_QC_y_ + idx);
                erival[nprim].QC_z[0] = LIBINT2_LOAD(tmp_QC_z_ + idx);

                erival[nprim].WP_x[0] = LIBINT2_LOAD(tmp_WP_x_ + idx);
                erival[nprim].WP_y[0] = LIBINT2_LOAD(tmp_WP_y_ + idx);
                erival[nprim].WP_z[0] = LIBINT2_LOAD(tmp_WP_z_ + idx);

                erival[nprim].WQ_x[0] = LIBINT2_LOAD(tmp_WQ_x_ + idx);
                erival[nprim].WQ_y[0] = LIBINT2_LOAD(tmp_WQ_y_ + idx);
                erival[nprim].WQ_z[0] = LIBINT2_LOAD(tmp_WQ_z_ + idx);

                erival[nprim].oo2z[0]  = LIBINT2_LOAD(tmp_oo2z_ + idx);
                erival[nprim].oo2e[0]  = LIBINT2_LOAD(tmp_oo2e_ + idx);
                erival[nprim].oo2ze[0] = LIBINT2_LOAD(tmp_oo2ze_ + idx);
                erival[nprim].roz[0]   = LIBINT2_LOAD(tmp_roz_ + idx);
                erival[nprim].roe[0]   = LIBINT2_LOAD(tmp_roe_ + idx);


                const int idx2 = nprim*LIBINT2_SIMD_LEN*(M+1);

                switch(M)
                {
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(28))
                    case 28:
                        erival[nprim].LIBINT_T_SS_EREP_SS(28)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 28*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(27))
                    case 27:
                        erival[nprim].LIBINT_T_SS_EREP_SS(27)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 27*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(26))
                    case 26:
                        erival[nprim].LIBINT_T_SS_EREP_SS(26)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 26*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(25))
                    case 25:
                        erival[nprim].LIBINT_T_SS_EREP_SS(25)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 25*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(24))
                    case 24:
                        erival[nprim].LIBINT_T_SS_EREP_SS(24)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 24*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(23))
                    case 23:
                        erival[nprim].LIBINT_T_SS_EREP_SS(23)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 23*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(22))
                    case 22:
                        erival[nprim].LIBINT_T_SS_EREP_SS(22)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 22*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(21))
                    case 21:
                        erival[nprim].LIBINT_T_SS_EREP_SS(21)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 21*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
                    case 20:
                        erival[nprim].LIBINT_T_SS_EREP_SS(20)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 20*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
                    case 19:
                        erival[nprim].LIBINT_T_SS_EREP_SS(19)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 19*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
                    case 18:
                        erival[nprim].LIBINT_T_SS_EREP_SS(18)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 18*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
                    case 17:
                        erival[nprim].LIBINT_T_SS_EREP_SS(17)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 17*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
                    case 16:
                        erival[nprim].LIBINT_T_SS_EREP_SS(16)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 16*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
                    case 15:
                        erival[nprim].LIBINT_T_SS_EREP_SS(15)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 15*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
                    case 14:
                        erival[nprim].LIBINT_T_SS_EREP_SS(14)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 14*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
                    case 13:
                        erival[nprim].LIBINT_T_SS_EREP_SS(13)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 13*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
                    case 12:
                        erival[nprim].LIBINT_T_SS_EREP_SS(12)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 12*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
                    case 11:
                        erival[nprim].LIBINT_T_SS_EREP_SS(11)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 11*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
                    case 10:
                        erival[nprim].LIBINT_T_SS_EREP_SS(10)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 10*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
                    case 9:
                        erival[nprim].LIBINT_T_SS_EREP_SS(9)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 9*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
                    case 8:
                        erival[nprim].LIBINT_T_SS_EREP_SS(8)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 8*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
                    case 7:
                        erival[nprim].LIBINT_T_SS_EREP_SS(7)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 7*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
                    case 6:
                        erival[nprim].LIBINT_T_SS_EREP_SS(6)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 6*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
                    case 5:
                        erival[nprim].LIBINT_T_SS_EREP_SS(5)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 5*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
                    case 4:
                        erival[nprim].LIBINT_T_SS_EREP_SS(4)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 4*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
                    case 3:
                        erival[nprim].LIBINT_T_SS_EREP_SS(3)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 3*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
                    case 2:
                        erival[nprim].LIBINT_T_SS_EREP_SS(2)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 2*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
                    case 1:
                        erival[nprim].LIBINT_T_SS_EREP_SS(1)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 1*LIBINT2_SIMD_LEN);
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
                    case 0:
                        erival[nprim].LIBINT_T_SS_EREP_SS(0)[0] = LIBINT2_LOAD(tmp_vecF_ + idx2 + 0*LIBINT2_SIMD_LEN);
                    #endif
                }   
            }

            // loop over vector elements
            for(int scd = 0; scd < simdlen; scd++)
            {
                double * intptr = integrals + ( (cd+scd) * Q.nshell12 + ab ) * ncart1234;
        
                if(M)
                {
                    CLOCK(ticks0);
                    LIBINT2_PREFIXED_NAME(libint2_build_eri)[Q.am1][Q.am2][P.am1][P.am2](erival);
                    CLOCK(ticks1);
                    totaltime += (ticks1 - ticks0);

                    // permute
                    int nn = 0;
                    for(int n1 = 0; n1 < NCART(Q.am1); n1++)
                    for(int n2 = 0; n2 < NCART(Q.am2); n2++)
                    for(int n3 = 0; n3 < NCART(P.am1); n3++)
                    for(int n4 = 0; n4 < NCART(P.am2); n4++)
                    {
                        #ifdef TESTS_LIBINT2_SIMD
                        double const * const splt = reinterpret_cast<double *>(&erival[0].targets[0][nn++].d);
                        //std::cout << "nn: " << (nn-1) << "  ->  " << splt[0] << "  " << splt[1] << "  " << splt[2] << "  " << splt[3] << "\n";

                        intptr[ n3*NCART(P.am2)*NCART(Q.am1)*NCART(Q.am2)
                               +n4*NCART(Q.am1)*NCART(Q.am2)
                               +n1*NCART(Q.am2)
                               +n2] = splt[scd];
                        #else
                        intptr[ n3*NCART(P.am2)*NCART(Q.am1)*NCART(Q.am2)
                               +n4*NCART(Q.am1)*NCART(Q.am2)
                               +n1*NCART(Q.am2)
                               +n2] = erival[0].targets[0][nn++]; 
                        #endif
                    }
                }
                else
                {
                    intptr[0] = 0.0;
                    #pragma novector
                    for(int i = 0; i < P.nprim12[cd+scd]; ++i)
                    {
                        #pragma novector
                        for(int j = 0; j < Q.nprim12[ab]; ++j)
                        {
                            int idx = j * P.nprim12[cd+scd] + i;

                            #ifdef TESTS_LIBINT2_SIMD
                            double const * const splt = reinterpret_cast<double *>(&erival[idx].LIBINT_T_SS_EREP_SS(0)[0].d);
                            intptr[0] += splt[scd];
                            #else
                            intptr[0] += erival[idx].LIBINT_T_SS_EREP_SS(0)[0]; 
                            #endif
                        }
                    }
                }
            }

        }

        istart = iend;
        if(((ab+1) % SIMINT_NSHELL_SIMD) == 0)
            istart = SIMINT_SIMD_ROUND(istart);
        
    }

    return totaltime;
}

