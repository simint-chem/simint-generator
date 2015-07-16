#include <libint2.h>

#include "test/Libint2.hpp"
#include "boys/boys_split.h"
#include "test/common.hpp"

#define MAX(a,b) (((a)>(b))?(a):(b))


TimerInfo Libint2_ERI::Integrals(struct multishell_pair P,
                                 struct multishell_pair Q,
                                 double * integrals)
{

    int maxam = MAX(MAX(P.am1, P.am2), MAX(Q.am1, Q.am2));
    size_t ncart1234 = NCART(P.am1) * NCART(P.am2) * NCART(Q.am1) * NCART(Q.am2);

    // compute the parameters
    int nprim = 0;
    int ab = 0;
    int cd = 0;

    // timing
    TimerInfo totaltime = {0, 0.0};

    // we will stride with this
    double * intptr = integrals;

    for(int a = 0; a < P.nshell1; a++)
    for(int b = 0; b < P.nshell2; b++, ab++)
    for(int c = 0; c < Q.nshell1; c++)
    for(int d = 0; d < Q.nshell2; d++, cd++)
    {
        // assumes proper input ordering
        std::vector<Libint_eri_t> erival(P.nprim12[ab]*Q.nprim12[cd]);
        LIBINT2_PREFIXED_NAME(libint2_init_eri)(erival.data(), maxam, 0); 
        erival[0].contrdepth = erival.size();

        for(int i = P.primstart[ab]; i < P.primend[ab]; i++)
        for(int j = Q.primstart[cd]; j < Q.primend[cd]; j++)
        {
            erival[nprim].AB_x[0] = P.AB_x[i];
            erival[nprim].AB_y[0] = P.AB_y[i];
            erival[nprim].AB_z[0] = P.AB_z[i];

            erival[nprim].CD_x[0] = Q.AB_x[j];
            erival[nprim].CD_y[0] = Q.AB_y[j];
            erival[nprim].CD_z[0] = Q.AB_z[j];

            erival[nprim].PA_x[0] = P.PA_x[i];
            erival[nprim].PA_y[0] = P.PA_y[i];
            erival[nprim].PA_z[0] = P.PA_z[i];

            erival[nprim].QC_x[0] = Q.PA_x[j];
            erival[nprim].QC_y[0] = Q.PA_y[j];
            erival[nprim].QC_z[0] = Q.PA_z[j];

            double W[3], WP[3], WQ[3];
            W[0] = (P.alpha[i] * P.x[i] + Q.alpha[j] * Q.x[j]) / (P.alpha[i] + Q.alpha[j]);
            W[1] = (P.alpha[i] * P.y[i] + Q.alpha[j] * Q.y[j]) / (P.alpha[i] + Q.alpha[j]);
            W[2] = (P.alpha[i] * P.z[i] + Q.alpha[j] * Q.z[j]) / (P.alpha[i] + Q.alpha[j]);
            WP[0] = W[0] - P.x[i];
            WP[1] = W[1] - P.y[i];
            WP[2] = W[2] - P.z[i];
            WQ[0] = W[0] - Q.x[j];
            WQ[1] = W[1] - Q.y[j];
            WQ[2] = W[2] - Q.z[j];
            
            erival[nprim].WP_x[0] = WP[0]; 
            erival[nprim].WP_y[0] = WP[1]; 
            erival[nprim].WP_z[0] = WP[2]; 
            erival[nprim].WQ_x[0] = WP[0]; 
            erival[nprim].WQ_y[0] = WP[1]; 
            erival[nprim].WQ_z[0] = WP[2]; 

            double rho = (P.alpha[i] * Q.alpha[j])/(P.alpha[i] + Q.alpha[j]);
            erival[nprim].oo2z[0] = 1.0 / (2.0 * P.alpha[i]); 
            erival[nprim].oo2e[0] = 1.0 / (2.0 * Q.alpha[j]); 
            erival[nprim].oo2ze[0] = 1.0 / (2.0 * (P.alpha[i] + Q.alpha[j])); 
            erival[nprim].roz[0] = rho / P.alpha[i]; 
            erival[nprim].roe[0] = rho / Q.alpha[j]; 
           

            double PQ2 = 0;
            PQ2 += (P.x[i] - Q.x[j])*(P.x[i] - Q.x[j]);
            PQ2 += (P.y[i] - Q.y[j])*(P.y[i] - Q.y[j]);
            PQ2 += (P.z[i] - Q.z[j])*(P.z[i] - Q.z[j]);
            double T = rho * PQ2; 

            // calculate the boys function
            int M = P.am1 + P.am2 + Q.am1 + Q.am2;
            std::vector<double> F(M+1);
            Boys_F_split(F.data(), M, T); 

            double scale = P.prefac[i] * Q.prefac[j];

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

        // timing
        unsigned long long ticks0, ticks1;
        double walltime0, walltime1;

        CLOCK(ticks0, walltime0);
        LIBINT2_PREFIXED_NAME(libint2_build_eri)[P.am1][P.am2][Q.am1][Q.am2](erival.data());
        CLOCK(ticks1, walltime1);
        std::copy(erival[0].targets[0], erival[0].targets[0] + ncart1234, intptr);
        intptr += ncart1234;


        totaltime += {ticks1 - ticks0, walltime1 - walltime0};
        LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(erival.data());
    }


    return totaltime;
}


