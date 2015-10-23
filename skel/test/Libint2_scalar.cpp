#include "constants.h"
#include "test/Libint2.hpp"
#include "boys/boys_split.h"
#include "test/common.hpp"
#include "vectorization/vectorization.h"


// Disable intel warnings
// 193 : zero used for undefined preprocessing identifier "LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18
//       Not much I can do about that 
#ifdef __INTEL_COMPILER
    #pragma warning(disable:193)
#endif



Libint2_ERI::Libint2_ERI(int maxam, size_t maxnprim)
{
    size_t size = maxnprim * maxnprim * maxnprim * maxnprim;
    erival_.resize(size);
    LIBINT2_PREFIXED_NAME(libint2_init_eri)(erival_.data(), maxam, 0); 
}


Libint2_ERI::~Libint2_ERI()
{
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(erival_.data());
}



TimerType Libint2_ERI::Integrals(struct multishell_pair P,
                                 struct multishell_pair Q,
                                 double * integrals)
{
    Libint_eri_t * erival = erival_.data();

    // according to libint, l(c)+l(d) >= l(a) + l(b)
    // so we switch P and Q, then permute
    const int M = P.am1 + P.am2 + Q.am1 + Q.am2;
    double F[M+1];

    size_t ncart1234 = NCART(P.am1) * NCART(P.am2) * NCART(Q.am1) * NCART(Q.am2);
    std::fill(integrals, integrals + P.nshell12 * Q.nshell12 * ncart1234, 0.0);

    // timing
    TimerType totaltime = 0;
    TimerType ticks0, ticks1;

    int istart = 0;
    for(int a = 0, ab = 0; a < Q.nshell1; a++)
    for(int b = 0;         b < Q.nshell2; b++, ab++)
    {
        int iend = istart + Q.nprim12[ab];

        int jstart = 0;

        for(int c = 0, cd = 0; c < P.nshell1; c++)
        for(int d = 0;         d < P.nshell2; d++, cd++)
        {
            int jend = jstart + P.nprim12[cd];

            size_t nprim1234 = Q.nprim12[ab] * P.nprim12[cd];
            erival[0].contrdepth = nprim1234;
            int nprim = 0;

            CLOCK(ticks0);
            for(int i = istart; i < iend; i++)
            for(int j = jstart; j < jend; j++)
            {
                const double PQalpha_sum = P.alpha[j] + Q.alpha[i];
                const double PQalpha_mul = P.alpha[j] * Q.alpha[i];
                const double one_over_PQalpha_sum = 1.0 / PQalpha_sum;
                const double one_over_Palpha = 1.0 / P.alpha[j];
                const double one_over_Qalpha = 1.0 / Q.alpha[i];

                erival[nprim].AB_x[0] = Q.AB_x[ab];
                erival[nprim].AB_y[0] = Q.AB_y[ab];
                erival[nprim].AB_z[0] = Q.AB_z[ab];

                erival[nprim].CD_x[0] = P.AB_x[cd];
                erival[nprim].CD_y[0] = P.AB_y[cd];
                erival[nprim].CD_z[0] = P.AB_z[cd];

                erival[nprim].PA_x[0] = Q.PA_x[i];
                erival[nprim].PA_y[0] = Q.PA_y[i];
                erival[nprim].PA_z[0] = Q.PA_z[i];

                erival[nprim].QC_x[0] = P.PA_x[j];
                erival[nprim].QC_y[0] = P.PA_y[j];
                erival[nprim].QC_z[0] = P.PA_z[j];

                double W[3];
                W[0] = (P.alpha[j] * P.x[j] + Q.alpha[i] * Q.x[i]) * one_over_PQalpha_sum;
                W[1] = (P.alpha[j] * P.y[j] + Q.alpha[i] * Q.y[i]) * one_over_PQalpha_sum;
                W[2] = (P.alpha[j] * P.z[j] + Q.alpha[i] * Q.z[i]) * one_over_PQalpha_sum;
                
                erival[nprim].WP_x[0] = W[0] - Q.x[i]; 
                erival[nprim].WP_y[0] = W[1] - Q.y[i]; 
                erival[nprim].WP_z[0] = W[2] - Q.z[i]; 

                erival[nprim].WQ_x[0] = W[0] - P.x[j]; 
                erival[nprim].WQ_y[0] = W[1] - P.y[j]; 
                erival[nprim].WQ_z[0] = W[2] - P.z[j]; 

                double rho = PQalpha_mul * one_over_PQalpha_sum;
                erival[nprim].oo2z[0] = 0.5 * one_over_Qalpha; 
                erival[nprim].oo2e[0] = 0.5 * one_over_Palpha; 
                erival[nprim].oo2ze[0] = 0.5 * one_over_PQalpha_sum;
                erival[nprim].roz[0] = rho * one_over_Qalpha; 
                erival[nprim].roe[0] = rho * one_over_Palpha; 
               

                double PQ2 = 0.0;
                PQ2 += (P.x[j] - Q.x[i])*(P.x[j] - Q.x[i]);
                PQ2 += (P.y[j] - Q.y[i])*(P.y[j] - Q.y[i]);
                PQ2 += (P.z[j] - Q.z[i])*(P.z[j] - Q.z[i]);
                double T = rho * PQ2; 


                // calculate the boys function
                Boys_F_split(F, M, T); 

                // temporarily disable so that it matches bit-for-bit my code
                double scale = sqrt(one_over_PQalpha_sum) * P.prefac[j] * Q.prefac[i];

                switch(M)
                {
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(28))
                    case 28:
                        erival[nprim].LIBINT_T_SS_EREP_SS(28)[0] = F[28] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(27))
                    case 27:
                        erival[nprim].LIBINT_T_SS_EREP_SS(27)[0] = F[27] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(26))
                    case 26:
                        erival[nprim].LIBINT_T_SS_EREP_SS(26)[0] = F[26] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(25))
                    case 25:
                        erival[nprim].LIBINT_T_SS_EREP_SS(25)[0] = F[25] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(24))
                    case 24:
                        erival[nprim].LIBINT_T_SS_EREP_SS(24)[0] = F[24] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(23))
                    case 23:
                        erival[nprim].LIBINT_T_SS_EREP_SS(23)[0] = F[23] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(22))
                    case 22:
                        erival[nprim].LIBINT_T_SS_EREP_SS(22)[0] = F[22] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(21))
                    case 21:
                        erival[nprim].LIBINT_T_SS_EREP_SS(21)[0] = F[21] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
                    case 20:
                        erival[nprim].LIBINT_T_SS_EREP_SS(20)[0] = F[20] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
                    case 19:
                        erival[nprim].LIBINT_T_SS_EREP_SS(19)[0] = F[19] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
                    case 18:
                        erival[nprim].LIBINT_T_SS_EREP_SS(18)[0] = F[18] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
                    case 17:
                        erival[nprim].LIBINT_T_SS_EREP_SS(17)[0] = F[17] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
                    case 16:
                        erival[nprim].LIBINT_T_SS_EREP_SS(16)[0] = F[16] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
                    case 15:
                        erival[nprim].LIBINT_T_SS_EREP_SS(15)[0] = F[15] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
                    case 14:
                        erival[nprim].LIBINT_T_SS_EREP_SS(14)[0] = F[14] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
                    case 13:
                        erival[nprim].LIBINT_T_SS_EREP_SS(13)[0] = F[13] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
                    case 12:
                        erival[nprim].LIBINT_T_SS_EREP_SS(12)[0] = F[12] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
                    case 11:
                        erival[nprim].LIBINT_T_SS_EREP_SS(11)[0] = F[11] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
                    case 10:
                        erival[nprim].LIBINT_T_SS_EREP_SS(10)[0] = F[10] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
                    case 9:
                        erival[nprim].LIBINT_T_SS_EREP_SS(9)[0] = F[9] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
                    case 8:
                        erival[nprim].LIBINT_T_SS_EREP_SS(8)[0] = F[8] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
                    case 7:
                        erival[nprim].LIBINT_T_SS_EREP_SS(7)[0] = F[7] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
                    case 6:
                        erival[nprim].LIBINT_T_SS_EREP_SS(6)[0] = F[6] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
                    case 5:
                        erival[nprim].LIBINT_T_SS_EREP_SS(5)[0] = F[5] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
                    case 4:
                        erival[nprim].LIBINT_T_SS_EREP_SS(4)[0] = F[4] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
                    case 3:
                        erival[nprim].LIBINT_T_SS_EREP_SS(3)[0] = F[3] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
                    case 2:
                        erival[nprim].LIBINT_T_SS_EREP_SS(2)[0] = F[2] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
                    case 1:
                        erival[nprim].LIBINT_T_SS_EREP_SS(1)[0] = F[1] * scale;
                    #endif
                    #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
                    case 0:
                        erival[nprim].LIBINT_T_SS_EREP_SS(0)[0] = F[0] * scale;
                    #endif
                }   

                nprim++;
            }
            CLOCK(ticks1);
            totaltime += (ticks1 - ticks0);

            double * intptr = integrals + ( cd * Q.nshell12 + ab ) * ncart1234;
    
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
                    intptr[ n3*NCART(P.am2)*NCART(Q.am1)*NCART(Q.am2)
                           +n4*NCART(Q.am1)*NCART(Q.am2)
                           +n1*NCART(Q.am2)
                           +n2] = erival[0].targets[0][nn++];
            }
            else
            {
                intptr[0] = 0.0;
                // sum in the same order as would happen in my code
                #pragma novector
                for(int i = 0; i < P.nprim12[cd]; ++i)
                {
                    #pragma novector
                    for(int j = 0; j < Q.nprim12[ab]; ++j)
                    {
                        int idx = j * P.nprim12[cd] + i;
                        intptr[0] += erival[idx].LIBINT_T_SS_EREP_SS(0)[0];
                    }
                }
            }

            jstart = jend;
            if(((cd+1) % SIMINT_NSHELL_SIMD) == 0)
                jstart = SIMINT_SIMD_ROUND(jstart);
        }

        istart = iend;
        if(((ab+1) % SIMINT_NSHELL_SIMD) == 0)
            istart = SIMINT_SIMD_ROUND(istart);
        
    }

    return totaltime;
}


