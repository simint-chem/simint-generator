#include "constants.h"
#include "test/Libint2.hpp"
#include "boys/boys_FO.h"
#include "test/common.hpp"


#define MAX(a,b) (((a)>(b))?(a):(b))


Libint2_ERI::Libint2_ERI(int maxam, size_t maxnprim, size_t maxsize)
{
    size_t size = maxnprim * maxnprim * maxnprim * maxnprim;

    erival_.resize(size);
    LIBINT2_PREFIXED_NAME(libint2_init_eri)(erival_.data(), maxam, 0); 
}


Libint2_ERI::~Libint2_ERI()
{
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(erival_.data());
}



TimerInfo Libint2_ERI::Integrals(struct multishell_pair P,
                                 struct multishell_pair Q,
                                 double * integrals)
{
    // according to libint, l(c)+l(d) >= l(a) + l(b)
    // so we switch P and Q, then permute

    int M = P.am1 + P.am2 + Q.am1 + Q.am2;

    size_t ncart1234 = NCART(P.am1) * NCART(P.am2) * NCART(Q.am1) * NCART(Q.am2);

    std::fill(integrals, integrals + (P.nshell12 * Q.nshell12 * ncart1234), 999.0);

    // timing
    TimerInfo totaltime = {0, 0.0};

    for(int a = 0, ab = 0; a < Q.nshell1; a++)
    for(int b = 0;         b < Q.nshell2; b++, ab++)
    for(int c = 0, cd = 0; c < P.nshell1; c++)
    for(int d = 0;         d < P.nshell2; d++, cd++)
    {
        size_t nprim1234 = Q.nprim12[ab] * P.nprim12[cd];
        erival_[0].contrdepth = nprim1234;
        int nprim = 0;

        // timing
        unsigned long long ticks0, ticks1;

        for(int i = Q.primstart[ab]; i < Q.primstart[ab]+Q.nprim12[ab]; i++)
        for(int j = P.primstart[cd]; j < P.primstart[cd]+P.nprim12[cd]; j++)
        {
            erival_[nprim].AB_x[0] = Q.AB_x[ab];
            erival_[nprim].AB_y[0] = Q.AB_y[ab];
            erival_[nprim].AB_z[0] = Q.AB_z[ab];

            erival_[nprim].CD_x[0] = P.AB_x[cd];
            erival_[nprim].CD_y[0] = P.AB_y[cd];
            erival_[nprim].CD_z[0] = P.AB_z[cd];

            erival_[nprim].PA_x[0] = Q.PA_x[i];
            erival_[nprim].PA_y[0] = Q.PA_y[i];
            erival_[nprim].PA_z[0] = Q.PA_z[i];

            erival_[nprim].QC_x[0] = P.PA_x[j];
            erival_[nprim].QC_y[0] = P.PA_y[j];
            erival_[nprim].QC_z[0] = P.PA_z[j];

            double W[3];
            W[0] = (Q.alpha[i] * Q.x[i] + P.alpha[j] * P.x[j]) / (Q.alpha[i] + P.alpha[j]);
            W[1] = (Q.alpha[i] * Q.y[i] + P.alpha[j] * P.y[j]) / (Q.alpha[i] + P.alpha[j]);
            W[2] = (Q.alpha[i] * Q.z[i] + P.alpha[j] * P.z[j]) / (Q.alpha[i] + P.alpha[j]);
            
            erival_[nprim].WP_x[0] = W[0] - Q.x[i]; 
            erival_[nprim].WP_y[0] = W[1] - Q.y[i]; 
            erival_[nprim].WP_z[0] = W[2] - Q.z[i]; 

            erival_[nprim].WQ_x[0] = W[0] - P.x[j]; 
            erival_[nprim].WQ_y[0] = W[1] - P.y[j]; 
            erival_[nprim].WQ_z[0] = W[2] - P.z[j]; 

            double rho = (Q.alpha[i] * P.alpha[j])/(Q.alpha[i] + P.alpha[j]);
            erival_[nprim].oo2z[0] = 1.0 / (2.0 * Q.alpha[i]); 
            erival_[nprim].oo2e[0] = 1.0 / (2.0 * P.alpha[j]); 
            erival_[nprim].oo2ze[0] = 1.0 / (2.0 * (Q.alpha[i] + P.alpha[j])); 
            erival_[nprim].roz[0] = rho / Q.alpha[i]; 
            erival_[nprim].roe[0] = rho / P.alpha[j]; 
           

            double PQ2 = 0.0;
            PQ2 += (Q.x[i] - P.x[j])*(Q.x[i] - P.x[j]);
            PQ2 += (Q.y[i] - P.y[j])*(Q.y[i] - P.y[j]);
            PQ2 += (Q.z[i] - P.z[j])*(Q.z[i] - P.z[j]);
            double T = rho * PQ2; 


            // calculate the boys function
            std::vector<double> F(M+1);
            //CLOCK(ticks0);
            Boys_F_FO(F.data(), M, T); 
            //CLOCK(ticks1);
            //totaltime += {ticks1 - ticks0, (ticks1 - ticks0)/(1.0e9*PROC_CYCLES_PER_SECOND)};

            double scale = Q.prefac[i] * P.prefac[j] * sqrt(rho) * sqrt( (1.0 / Q.alpha[i] ) * (1.0 / P.alpha[j]));

            switch(M)
            {
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(28))
                case 28:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(28)[0] = F[28] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(27))
                case 27:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(27)[0] = F[27] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(26))
                case 26:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(26)[0] = F[26] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(25))
                case 25:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(25)[0] = F[25] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(24))
                case 24:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(24)[0] = F[24] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(23))
                case 23:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(23)[0] = F[23] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(22))
                case 22:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(22)[0] = F[22] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(21))
                case 21:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(21)[0] = F[21] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
                case 20:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(20)[0] = F[20] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
                case 19:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(19)[0] = F[19] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
                case 18:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(18)[0] = F[18] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
                case 17:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(17)[0] = F[17] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
                case 16:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(16)[0] = F[16] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
                case 15:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(15)[0] = F[15] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
                case 14:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(14)[0] = F[14] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
                case 13:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(13)[0] = F[13] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
                case 12:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(12)[0] = F[12] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
                case 11:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(11)[0] = F[11] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
                case 10:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(10)[0] = F[10] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
                case 9:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(9)[0] = F[9] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
                case 8:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(8)[0] = F[8] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
                case 7:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(7)[0] = F[7] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
                case 6:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(6)[0] = F[6] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
                case 5:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(5)[0] = F[5] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
                case 4:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(4)[0] = F[4] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
                case 3:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(3)[0] = F[3] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
                case 2:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(2)[0] = F[2] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
                case 1:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(1)[0] = F[1] * scale;
                #endif
                #if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
                case 0:
                    erival_[nprim].LIBINT_T_SS_EREP_SS(0)[0] = F[0] * scale;
                #endif
            } 

            nprim++;
        }

        double * intptr = integrals + (  c*P.nshell2*Q.nshell1*Q.nshell2
                                       + d*Q.nshell1*Q.nshell2
                                       + a*Q.nshell2
                                       + b )*ncart1234;

        if(M)
        {
            CLOCK(ticks0);
            LIBINT2_PREFIXED_NAME(libint2_build_eri)[Q.am1][Q.am2][P.am1][P.am2](erival_.data());
            CLOCK(ticks1);
            totaltime += {ticks1 - ticks0, (ticks1 - ticks0)/(1.0e9*PROC_CYCLES_PER_SECOND)};

            // permute
            int nn = 0;
            for(int n1 = 0; n1 < NCART(Q.am1); n1++)
            for(int n2 = 0; n2 < NCART(Q.am2); n2++)
            for(int n3 = 0; n3 < NCART(P.am1); n3++)
            for(int n4 = 0; n4 < NCART(P.am2); n4++)
                intptr[ n3*NCART(P.am2)*NCART(Q.am1)*NCART(Q.am2)
                       +n4*NCART(Q.am1)*NCART(Q.am2)
                       +n1*NCART(Q.am2)
                       +n2] = erival_[0].targets[0][nn++];
        }
        else
        {
            intptr[0] = 0.0;
            for(size_t n = 0; n < nprim1234; n++)
                intptr[0] += erival_[n].LIBINT_T_SS_EREP_SS(0)[0];
        }
        
    }


    return totaltime;
}


