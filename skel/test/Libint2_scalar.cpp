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



TimeContrib Libint2_ERI::Integrals(struct multishell_pair P,
                                   struct multishell_pair Q,
                                   double * integrals)
{
    Libint_eri_t * erival = erival_.data();

    // permute?
    // Libint requires Q.am >= P.am
    struct multishell_pair * libint_P = &P;
    struct multishell_pair * libint_Q = &Q;
    
    bool permutePQ = false;
    if( (P.am1 + P.am2) > (Q.am1 + Q.am2) )
    {
        libint_P = &Q;
        libint_Q = &P;
        permutePQ = true;
    }
    

    const int M = P.am1 + P.am2 + Q.am1 + Q.am2;
    double F[M+1];

    size_t ncart1234 = NCART(P.am1) * NCART(P.am2) * NCART(Q.am1) * NCART(Q.am2);
    std::fill(integrals, integrals + P.nshell12 * Q.nshell12 * ncart1234, 0.0);

    // timing
    TimeContrib times;
    TimerType ticks_copy_0, ticks_copy_1;
    TimerType ticks_boys_0, ticks_boys_1;

    int istart = 0;
    for(int a = 0, ab = 0; a < libint_P->nshell1; a++)
    for(int b = 0;         b < libint_P->nshell2; b++, ab++)
    {
        int iend = istart + libint_P->nprim12[ab];

        int jstart = 0;

        for(int c = 0, cd = 0; c < libint_Q->nshell1; c++)
        for(int d = 0;         d < libint_Q->nshell2; d++, cd++)
        {
            int jend = jstart + libint_Q->nprim12[cd];

            size_t nprim1234 = libint_P->nprim12[ab] * libint_Q->nprim12[cd];
            erival[0].contrdepth = nprim1234;
            int nprim = 0;

            CLOCK(ticks_copy_0);
            for(int i = istart; i < iend; i++)
            for(int j = jstart; j < jend; j++)
            {
                const double PQalpha_sum = libint_Q->alpha[j] + libint_P->alpha[i];
                const double PQalpha_mul = libint_Q->alpha[j] * libint_P->alpha[i];
                const double one_over_PQalpha_sum = 1.0 / PQalpha_sum;
                const double one_over_Palpha = 1.0 / libint_Q->alpha[j];
                const double one_over_Qalpha = 1.0 / libint_P->alpha[i];
                const double rho = PQalpha_mul * one_over_PQalpha_sum;

                if(M)
                {
                    erival[nprim].AB_x[0] = libint_P->AB_x[ab];
                    erival[nprim].AB_y[0] = libint_P->AB_y[ab];
                    erival[nprim].AB_z[0] = libint_P->AB_z[ab];

                    erival[nprim].CD_x[0] = libint_Q->AB_x[cd];
                    erival[nprim].CD_y[0] = libint_Q->AB_y[cd];
                    erival[nprim].CD_z[0] = libint_Q->AB_z[cd];

                    erival[nprim].PA_x[0] = libint_P->PA_x[i];
                    erival[nprim].PA_y[0] = libint_P->PA_y[i];
                    erival[nprim].PA_z[0] = libint_P->PA_z[i];

                    erival[nprim].QC_x[0] = libint_Q->PA_x[j];
                    erival[nprim].QC_y[0] = libint_Q->PA_y[j];
                    erival[nprim].QC_z[0] = libint_Q->PA_z[j];

                    double W[3];
                    W[0] = (libint_Q->alpha[j] * libint_Q->x[j] + libint_P->alpha[i] * libint_P->x[i]) * one_over_PQalpha_sum;
                    W[1] = (libint_Q->alpha[j] * libint_Q->y[j] + libint_P->alpha[i] * libint_P->y[i]) * one_over_PQalpha_sum;
                    W[2] = (libint_Q->alpha[j] * libint_Q->z[j] + libint_P->alpha[i] * libint_P->z[i]) * one_over_PQalpha_sum;
                    
                    erival[nprim].WP_x[0] = W[0] - libint_P->x[i]; 
                    erival[nprim].WP_y[0] = W[1] - libint_P->y[i]; 
                    erival[nprim].WP_z[0] = W[2] - libint_P->z[i]; 

                    erival[nprim].WQ_x[0] = W[0] - libint_Q->x[j]; 
                    erival[nprim].WQ_y[0] = W[1] - libint_Q->y[j]; 
                    erival[nprim].WQ_z[0] = W[2] - libint_Q->z[j]; 

                    erival[nprim].oo2z[0] = 0.5 * one_over_Qalpha; 
                    erival[nprim].oo2e[0] = 0.5 * one_over_Palpha; 
                    erival[nprim].oo2ze[0] = 0.5 * one_over_PQalpha_sum;
                    erival[nprim].roz[0] = rho * one_over_Qalpha; 
                    erival[nprim].roe[0] = rho * one_over_Palpha; 
                }

               

                double PQ2 = 0.0;
                PQ2 += (libint_Q->x[j] - libint_P->x[i])*(libint_Q->x[j] - libint_P->x[i]);
                PQ2 += (libint_Q->y[j] - libint_P->y[i])*(libint_Q->y[j] - libint_P->y[i]);
                PQ2 += (libint_Q->z[j] - libint_P->z[i])*(libint_Q->z[j] - libint_P->z[i]);
                double T = rho * PQ2; 


                // calculate the boys function
                // NOTE - the time is also included in ticks_copy_*, so
                // we have to remember to take care of that later
                CLOCK(ticks_boys_0);
                Boys_F_split(F, M, T); 
                CLOCK(ticks_boys_1);

                const double scale = sqrt(one_over_PQalpha_sum) * libint_Q->prefac[j] * libint_P->prefac[i];

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
            CLOCK(ticks_copy_1);
            times.boys += ticks_boys_1 - ticks_boys_0;
            times.copy_data += (ticks_copy_1 - ticks_copy_0) - (ticks_boys_1 - ticks_boys_0); // time for calculating boys is included in copy

            double * intptr;
            if(permutePQ)  // ab loops through Q, cd loops through P
                intptr = integrals + ( cd * Q.nshell12 + ab ) * ncart1234;
            else
                intptr = integrals + ( ab * Q.nshell12 + cd ) * ncart1234;
   
            TimerType ticks_integrals_0, ticks_integrals_1;
 
            if(M)
            {
                CLOCK(ticks_integrals_0);
                LIBINT2_PREFIXED_NAME(libint2_build_eri)[libint_P->am1][libint_P->am2][libint_Q->am1][libint_Q->am2](erival);
                CLOCK(ticks_integrals_1);
                times.integrals += (ticks_integrals_1 - ticks_integrals_0);

                // permute
                TimerType ticks_permute_0, ticks_permute_1;
                CLOCK(ticks_permute_0);
                if(permutePQ)
                {
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
                    std::copy(erival[0].targets[0], erival[0].targets[0] + ncart1234, intptr);
                CLOCK(ticks_permute_1); // also include "not permuting" (ie copying)

                times.permute += ticks_permute_1 - ticks_permute_0;
            }
            else
            {
                CLOCK(ticks_integrals_0);
                intptr[0] = 0.0;
                for(int i = 0; i < P.nprim12[ab]*Q.nprim12[cd]; ++i)
                    intptr[0] += erival[i].LIBINT_T_SS_EREP_SS(0)[0];
                CLOCK(ticks_integrals_1);
                times.integrals += (ticks_integrals_1 - ticks_integrals_0);
            }

            jstart = jend;
            if(((cd+1) % SIMINT_NSHELL_SIMD) == 0)
                jstart = SIMINT_SIMD_ROUND(jstart);
        }

        istart = iend;
        if(((ab+1) % SIMINT_NSHELL_SIMD) == 0)
            istart = SIMINT_SIMD_ROUND(istart);
        
    }

    return times;
}


