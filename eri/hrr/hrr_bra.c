
#include "vectorization.h"


    //////////////////////////////////////////////
    // BRA: ( p p |
    // Steps: 9
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_BRA_p_p(
                   double * const restrict BRA_p_s,
                   double * const restrict BRA_p_p,
                   double * const restrict BRA_d_s,
                   const double bAB_x, const double bAB_y, const double bAB_z, const int ncart_ket
                 )
{
    int iket;


        for(iket = 0; iket < ncart_ket; ++iket)
        {
            // (p_1_0_0  p_1_0_0|_{i} = (d_2_0_0  s_0_0_0|_{t} + x_ab * (p_1_0_0  s_0_0_0|_{t}
            BRA_p_p[0 * ncart_ket + iket] = BRA_d_s[0 * ncart_ket + iket] + ( bAB_x * BRA_p_s[0 * ncart_ket + iket] );

            // (p_1_0_0  p_0_1_0|_{i} = (d_1_1_0  s_0_0_0|_{t} + y_ab * (p_1_0_0  s_0_0_0|_{t}
            BRA_p_p[1 * ncart_ket + iket] = BRA_d_s[1 * ncart_ket + iket] + ( bAB_y * BRA_p_s[0 * ncart_ket + iket] );

            // (p_1_0_0  p_0_0_1|_{i} = (d_1_0_1  s_0_0_0|_{t} + z_ab * (p_1_0_0  s_0_0_0|_{t}
            BRA_p_p[2 * ncart_ket + iket] = BRA_d_s[2 * ncart_ket + iket] + ( bAB_z * BRA_p_s[0 * ncart_ket + iket] );

            // (p_0_1_0  p_1_0_0|_{i} = (d_1_1_0  s_0_0_0|_{t} + x_ab * (p_0_1_0  s_0_0_0|_{t}
            BRA_p_p[3 * ncart_ket + iket] = BRA_d_s[1 * ncart_ket + iket] + ( bAB_x * BRA_p_s[1 * ncart_ket + iket] );

            // (p_0_1_0  p_0_1_0|_{i} = (d_0_2_0  s_0_0_0|_{t} + y_ab * (p_0_1_0  s_0_0_0|_{t}
            BRA_p_p[4 * ncart_ket + iket] = BRA_d_s[3 * ncart_ket + iket] + ( bAB_y * BRA_p_s[1 * ncart_ket + iket] );

            // (p_0_1_0  p_0_0_1|_{i} = (d_0_1_1  s_0_0_0|_{t} + z_ab * (p_0_1_0  s_0_0_0|_{t}
            BRA_p_p[5 * ncart_ket + iket] = BRA_d_s[4 * ncart_ket + iket] + ( bAB_z * BRA_p_s[1 * ncart_ket + iket] );

            // (p_0_0_1  p_1_0_0|_{i} = (d_1_0_1  s_0_0_0|_{t} + x_ab * (p_0_0_1  s_0_0_0|_{t}
            BRA_p_p[6 * ncart_ket + iket] = BRA_d_s[2 * ncart_ket + iket] + ( bAB_x * BRA_p_s[2 * ncart_ket + iket] );

            // (p_0_0_1  p_0_1_0|_{i} = (d_0_1_1  s_0_0_0|_{t} + y_ab * (p_0_0_1  s_0_0_0|_{t}
            BRA_p_p[7 * ncart_ket + iket] = BRA_d_s[4 * ncart_ket + iket] + ( bAB_y * BRA_p_s[2 * ncart_ket + iket] );

            // (p_0_0_1  p_0_0_1|_{i} = (d_0_0_2  s_0_0_0|_{t} + z_ab * (p_0_0_1  s_0_0_0|_{t}
            BRA_p_p[8 * ncart_ket + iket] = BRA_d_s[5 * ncart_ket + iket] + ( bAB_z * BRA_p_s[2 * ncart_ket + iket] );

        }


}


    //////////////////////////////////////////////
    // BRA: ( d p |
    // Steps: 18
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_BRA_d_p(
                   double * const restrict BRA_d_s,
                   double * const restrict BRA_d_p,
                   double * const restrict BRA_f_s,
                   const double bAB_x, const double bAB_y, const double bAB_z, const int ncart_ket
                 )
{
    int iket;


        for(iket = 0; iket < ncart_ket; ++iket)
        {
            // (d_2_0_0  p_1_0_0|_{i} = (f_3_0_0  s_0_0_0|_{t} + x_ab * (d_2_0_0  s_0_0_0|_{t}
            BRA_d_p[0 * ncart_ket + iket] = BRA_f_s[0 * ncart_ket + iket] + ( bAB_x * BRA_d_s[0 * ncart_ket + iket] );

            // (d_2_0_0  p_0_1_0|_{i} = (f_2_1_0  s_0_0_0|_{t} + y_ab * (d_2_0_0  s_0_0_0|_{t}
            BRA_d_p[1 * ncart_ket + iket] = BRA_f_s[1 * ncart_ket + iket] + ( bAB_y * BRA_d_s[0 * ncart_ket + iket] );

            // (d_2_0_0  p_0_0_1|_{i} = (f_2_0_1  s_0_0_0|_{t} + z_ab * (d_2_0_0  s_0_0_0|_{t}
            BRA_d_p[2 * ncart_ket + iket] = BRA_f_s[2 * ncart_ket + iket] + ( bAB_z * BRA_d_s[0 * ncart_ket + iket] );

            // (d_1_1_0  p_1_0_0|_{i} = (f_2_1_0  s_0_0_0|_{t} + x_ab * (d_1_1_0  s_0_0_0|_{t}
            BRA_d_p[3 * ncart_ket + iket] = BRA_f_s[1 * ncart_ket + iket] + ( bAB_x * BRA_d_s[1 * ncart_ket + iket] );

            // (d_1_1_0  p_0_1_0|_{i} = (f_1_2_0  s_0_0_0|_{t} + y_ab * (d_1_1_0  s_0_0_0|_{t}
            BRA_d_p[4 * ncart_ket + iket] = BRA_f_s[3 * ncart_ket + iket] + ( bAB_y * BRA_d_s[1 * ncart_ket + iket] );

            // (d_1_1_0  p_0_0_1|_{i} = (f_1_1_1  s_0_0_0|_{t} + z_ab * (d_1_1_0  s_0_0_0|_{t}
            BRA_d_p[5 * ncart_ket + iket] = BRA_f_s[4 * ncart_ket + iket] + ( bAB_z * BRA_d_s[1 * ncart_ket + iket] );

            // (d_1_0_1  p_1_0_0|_{i} = (f_2_0_1  s_0_0_0|_{t} + x_ab * (d_1_0_1  s_0_0_0|_{t}
            BRA_d_p[6 * ncart_ket + iket] = BRA_f_s[2 * ncart_ket + iket] + ( bAB_x * BRA_d_s[2 * ncart_ket + iket] );

            // (d_1_0_1  p_0_1_0|_{i} = (f_1_1_1  s_0_0_0|_{t} + y_ab * (d_1_0_1  s_0_0_0|_{t}
            BRA_d_p[7 * ncart_ket + iket] = BRA_f_s[4 * ncart_ket + iket] + ( bAB_y * BRA_d_s[2 * ncart_ket + iket] );

            // (d_1_0_1  p_0_0_1|_{i} = (f_1_0_2  s_0_0_0|_{t} + z_ab * (d_1_0_1  s_0_0_0|_{t}
            BRA_d_p[8 * ncart_ket + iket] = BRA_f_s[5 * ncart_ket + iket] + ( bAB_z * BRA_d_s[2 * ncart_ket + iket] );

            // (d_0_2_0  p_1_0_0|_{i} = (f_1_2_0  s_0_0_0|_{t} + x_ab * (d_0_2_0  s_0_0_0|_{t}
            BRA_d_p[9 * ncart_ket + iket] = BRA_f_s[3 * ncart_ket + iket] + ( bAB_x * BRA_d_s[3 * ncart_ket + iket] );

            // (d_0_2_0  p_0_1_0|_{i} = (f_0_3_0  s_0_0_0|_{t} + y_ab * (d_0_2_0  s_0_0_0|_{t}
            BRA_d_p[10 * ncart_ket + iket] = BRA_f_s[6 * ncart_ket + iket] + ( bAB_y * BRA_d_s[3 * ncart_ket + iket] );

            // (d_0_2_0  p_0_0_1|_{i} = (f_0_2_1  s_0_0_0|_{t} + z_ab * (d_0_2_0  s_0_0_0|_{t}
            BRA_d_p[11 * ncart_ket + iket] = BRA_f_s[7 * ncart_ket + iket] + ( bAB_z * BRA_d_s[3 * ncart_ket + iket] );

            // (d_0_1_1  p_1_0_0|_{i} = (f_1_1_1  s_0_0_0|_{t} + x_ab * (d_0_1_1  s_0_0_0|_{t}
            BRA_d_p[12 * ncart_ket + iket] = BRA_f_s[4 * ncart_ket + iket] + ( bAB_x * BRA_d_s[4 * ncart_ket + iket] );

            // (d_0_1_1  p_0_1_0|_{i} = (f_0_2_1  s_0_0_0|_{t} + y_ab * (d_0_1_1  s_0_0_0|_{t}
            BRA_d_p[13 * ncart_ket + iket] = BRA_f_s[7 * ncart_ket + iket] + ( bAB_y * BRA_d_s[4 * ncart_ket + iket] );

            // (d_0_1_1  p_0_0_1|_{i} = (f_0_1_2  s_0_0_0|_{t} + z_ab * (d_0_1_1  s_0_0_0|_{t}
            BRA_d_p[14 * ncart_ket + iket] = BRA_f_s[8 * ncart_ket + iket] + ( bAB_z * BRA_d_s[4 * ncart_ket + iket] );

            // (d_0_0_2  p_1_0_0|_{i} = (f_1_0_2  s_0_0_0|_{t} + x_ab * (d_0_0_2  s_0_0_0|_{t}
            BRA_d_p[15 * ncart_ket + iket] = BRA_f_s[5 * ncart_ket + iket] + ( bAB_x * BRA_d_s[5 * ncart_ket + iket] );

            // (d_0_0_2  p_0_1_0|_{i} = (f_0_1_2  s_0_0_0|_{t} + y_ab * (d_0_0_2  s_0_0_0|_{t}
            BRA_d_p[16 * ncart_ket + iket] = BRA_f_s[8 * ncart_ket + iket] + ( bAB_y * BRA_d_s[5 * ncart_ket + iket] );

            // (d_0_0_2  p_0_0_1|_{i} = (f_0_0_3  s_0_0_0|_{t} + z_ab * (d_0_0_2  s_0_0_0|_{t}
            BRA_d_p[17 * ncart_ket + iket] = BRA_f_s[9 * ncart_ket + iket] + ( bAB_z * BRA_d_s[5 * ncart_ket + iket] );

        }


}


    //////////////////////////////////////////////
    // BRA: ( d d |
    // Steps: 79
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_BRA_d_d(
                   double * const restrict BRA_d_s,
                   double * const restrict BRA_d_d,
                   double * const restrict BRA_f_s,
                   double * const restrict BRA_g_s,
                   const double bAB_x, const double bAB_y, const double bAB_z, const int ncart_ket
                 )
{
    int iket;


        for(iket = 0; iket < ncart_ket; ++iket)
        {
            // (d_2_0_0  p_1_0_0| = (f_3_0_0  s_0_0_0|_{t} + x_ab * (d_2_0_0  s_0_0_0|_{t}
            const double B_2_0_0_1_0_0 = BRA_f_s[0 * ncart_ket + iket] + ( bAB_x * BRA_d_s[0 * ncart_ket + iket] );

            // (d_2_0_0  p_0_1_0| = (f_2_1_0  s_0_0_0|_{t} + y_ab * (d_2_0_0  s_0_0_0|_{t}
            const double B_2_0_0_0_1_0 = BRA_f_s[1 * ncart_ket + iket] + ( bAB_y * BRA_d_s[0 * ncart_ket + iket] );

            // (d_2_0_0  p_0_0_1| = (f_2_0_1  s_0_0_0|_{t} + z_ab * (d_2_0_0  s_0_0_0|_{t}
            const double B_2_0_0_0_0_1 = BRA_f_s[2 * ncart_ket + iket] + ( bAB_z * BRA_d_s[0 * ncart_ket + iket] );

            // (d_1_1_0  p_1_0_0| = (f_2_1_0  s_0_0_0|_{t} + x_ab * (d_1_1_0  s_0_0_0|_{t}
            const double B_1_1_0_1_0_0 = BRA_f_s[1 * ncart_ket + iket] + ( bAB_x * BRA_d_s[1 * ncart_ket + iket] );

            // (d_1_1_0  p_0_1_0| = (f_1_2_0  s_0_0_0|_{t} + y_ab * (d_1_1_0  s_0_0_0|_{t}
            const double B_1_1_0_0_1_0 = BRA_f_s[3 * ncart_ket + iket] + ( bAB_y * BRA_d_s[1 * ncart_ket + iket] );

            // (d_1_1_0  p_0_0_1| = (f_1_1_1  s_0_0_0|_{t} + z_ab * (d_1_1_0  s_0_0_0|_{t}
            const double B_1_1_0_0_0_1 = BRA_f_s[4 * ncart_ket + iket] + ( bAB_z * BRA_d_s[1 * ncart_ket + iket] );

            // (d_1_0_1  p_1_0_0| = (f_2_0_1  s_0_0_0|_{t} + x_ab * (d_1_0_1  s_0_0_0|_{t}
            const double B_1_0_1_1_0_0 = BRA_f_s[2 * ncart_ket + iket] + ( bAB_x * BRA_d_s[2 * ncart_ket + iket] );

            // (d_1_0_1  p_0_1_0| = (f_1_1_1  s_0_0_0|_{t} + y_ab * (d_1_0_1  s_0_0_0|_{t}
            const double B_1_0_1_0_1_0 = BRA_f_s[4 * ncart_ket + iket] + ( bAB_y * BRA_d_s[2 * ncart_ket + iket] );

            // (d_1_0_1  p_0_0_1| = (f_1_0_2  s_0_0_0|_{t} + z_ab * (d_1_0_1  s_0_0_0|_{t}
            const double B_1_0_1_0_0_1 = BRA_f_s[5 * ncart_ket + iket] + ( bAB_z * BRA_d_s[2 * ncart_ket + iket] );

            // (d_0_2_0  p_1_0_0| = (f_1_2_0  s_0_0_0|_{t} + x_ab * (d_0_2_0  s_0_0_0|_{t}
            const double B_0_2_0_1_0_0 = BRA_f_s[3 * ncart_ket + iket] + ( bAB_x * BRA_d_s[3 * ncart_ket + iket] );

            // (d_0_2_0  p_0_1_0| = (f_0_3_0  s_0_0_0|_{t} + y_ab * (d_0_2_0  s_0_0_0|_{t}
            const double B_0_2_0_0_1_0 = BRA_f_s[6 * ncart_ket + iket] + ( bAB_y * BRA_d_s[3 * ncart_ket + iket] );

            // (d_0_2_0  p_0_0_1| = (f_0_2_1  s_0_0_0|_{t} + z_ab * (d_0_2_0  s_0_0_0|_{t}
            const double B_0_2_0_0_0_1 = BRA_f_s[7 * ncart_ket + iket] + ( bAB_z * BRA_d_s[3 * ncart_ket + iket] );

            // (d_0_1_1  p_1_0_0| = (f_1_1_1  s_0_0_0|_{t} + x_ab * (d_0_1_1  s_0_0_0|_{t}
            const double B_0_1_1_1_0_0 = BRA_f_s[4 * ncart_ket + iket] + ( bAB_x * BRA_d_s[4 * ncart_ket + iket] );

            // (d_0_1_1  p_0_1_0| = (f_0_2_1  s_0_0_0|_{t} + y_ab * (d_0_1_1  s_0_0_0|_{t}
            const double B_0_1_1_0_1_0 = BRA_f_s[7 * ncart_ket + iket] + ( bAB_y * BRA_d_s[4 * ncart_ket + iket] );

            // (d_0_1_1  p_0_0_1| = (f_0_1_2  s_0_0_0|_{t} + z_ab * (d_0_1_1  s_0_0_0|_{t}
            const double B_0_1_1_0_0_1 = BRA_f_s[8 * ncart_ket + iket] + ( bAB_z * BRA_d_s[4 * ncart_ket + iket] );

            // (d_0_0_2  p_1_0_0| = (f_1_0_2  s_0_0_0|_{t} + x_ab * (d_0_0_2  s_0_0_0|_{t}
            const double B_0_0_2_1_0_0 = BRA_f_s[5 * ncart_ket + iket] + ( bAB_x * BRA_d_s[5 * ncart_ket + iket] );

            // (d_0_0_2  p_0_1_0| = (f_0_1_2  s_0_0_0|_{t} + y_ab * (d_0_0_2  s_0_0_0|_{t}
            const double B_0_0_2_0_1_0 = BRA_f_s[8 * ncart_ket + iket] + ( bAB_y * BRA_d_s[5 * ncart_ket + iket] );

            // (d_0_0_2  p_0_0_1| = (f_0_0_3  s_0_0_0|_{t} + z_ab * (d_0_0_2  s_0_0_0|_{t}
            const double B_0_0_2_0_0_1 = BRA_f_s[9 * ncart_ket + iket] + ( bAB_z * BRA_d_s[5 * ncart_ket + iket] );

            // (f_3_0_0  p_1_0_0| = (g_4_0_0  s_0_0_0|_{t} + x_ab * (f_3_0_0  s_0_0_0|_{t}
            const double B_3_0_0_1_0_0 = BRA_g_s[0 * ncart_ket + iket] + ( bAB_x * BRA_f_s[0 * ncart_ket + iket] );

            // (f_2_1_0  p_1_0_0| = (g_3_1_0  s_0_0_0|_{t} + x_ab * (f_2_1_0  s_0_0_0|_{t}
            const double B_2_1_0_1_0_0 = BRA_g_s[1 * ncart_ket + iket] + ( bAB_x * BRA_f_s[1 * ncart_ket + iket] );

            // (f_2_1_0  p_0_1_0| = (g_2_2_0  s_0_0_0|_{t} + y_ab * (f_2_1_0  s_0_0_0|_{t}
            const double B_2_1_0_0_1_0 = BRA_g_s[3 * ncart_ket + iket] + ( bAB_y * BRA_f_s[1 * ncart_ket + iket] );

            // (f_2_0_1  p_1_0_0| = (g_3_0_1  s_0_0_0|_{t} + x_ab * (f_2_0_1  s_0_0_0|_{t}
            const double B_2_0_1_1_0_0 = BRA_g_s[2 * ncart_ket + iket] + ( bAB_x * BRA_f_s[2 * ncart_ket + iket] );

            // (f_2_0_1  p_0_1_0| = (g_2_1_1  s_0_0_0|_{t} + y_ab * (f_2_0_1  s_0_0_0|_{t}
            const double B_2_0_1_0_1_0 = BRA_g_s[4 * ncart_ket + iket] + ( bAB_y * BRA_f_s[2 * ncart_ket + iket] );

            // (f_2_0_1  p_0_0_1| = (g_2_0_2  s_0_0_0|_{t} + z_ab * (f_2_0_1  s_0_0_0|_{t}
            const double B_2_0_1_0_0_1 = BRA_g_s[5 * ncart_ket + iket] + ( bAB_z * BRA_f_s[2 * ncart_ket + iket] );

            // (f_1_2_0  p_1_0_0| = (g_2_2_0  s_0_0_0|_{t} + x_ab * (f_1_2_0  s_0_0_0|_{t}
            const double B_1_2_0_1_0_0 = BRA_g_s[3 * ncart_ket + iket] + ( bAB_x * BRA_f_s[3 * ncart_ket + iket] );

            // (f_1_2_0  p_0_1_0| = (g_1_3_0  s_0_0_0|_{t} + y_ab * (f_1_2_0  s_0_0_0|_{t}
            const double B_1_2_0_0_1_0 = BRA_g_s[6 * ncart_ket + iket] + ( bAB_y * BRA_f_s[3 * ncart_ket + iket] );

            // (f_1_1_1  p_1_0_0| = (g_2_1_1  s_0_0_0|_{t} + x_ab * (f_1_1_1  s_0_0_0|_{t}
            const double B_1_1_1_1_0_0 = BRA_g_s[4 * ncart_ket + iket] + ( bAB_x * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_1_1  p_0_1_0| = (g_1_2_1  s_0_0_0|_{t} + y_ab * (f_1_1_1  s_0_0_0|_{t}
            const double B_1_1_1_0_1_0 = BRA_g_s[7 * ncart_ket + iket] + ( bAB_y * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_1_1  p_0_0_1| = (g_1_1_2  s_0_0_0|_{t} + z_ab * (f_1_1_1  s_0_0_0|_{t}
            const double B_1_1_1_0_0_1 = BRA_g_s[8 * ncart_ket + iket] + ( bAB_z * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_0_2  p_1_0_0| = (g_2_0_2  s_0_0_0|_{t} + x_ab * (f_1_0_2  s_0_0_0|_{t}
            const double B_1_0_2_1_0_0 = BRA_g_s[5 * ncart_ket + iket] + ( bAB_x * BRA_f_s[5 * ncart_ket + iket] );

            // (f_1_0_2  p_0_1_0| = (g_1_1_2  s_0_0_0|_{t} + y_ab * (f_1_0_2  s_0_0_0|_{t}
            const double B_1_0_2_0_1_0 = BRA_g_s[8 * ncart_ket + iket] + ( bAB_y * BRA_f_s[5 * ncart_ket + iket] );

            // (f_1_0_2  p_0_0_1| = (g_1_0_3  s_0_0_0|_{t} + z_ab * (f_1_0_2  s_0_0_0|_{t}
            const double B_1_0_2_0_0_1 = BRA_g_s[9 * ncart_ket + iket] + ( bAB_z * BRA_f_s[5 * ncart_ket + iket] );

            // (f_0_3_0  p_1_0_0| = (g_1_3_0  s_0_0_0|_{t} + x_ab * (f_0_3_0  s_0_0_0|_{t}
            const double B_0_3_0_1_0_0 = BRA_g_s[6 * ncart_ket + iket] + ( bAB_x * BRA_f_s[6 * ncart_ket + iket] );

            // (f_0_3_0  p_0_1_0| = (g_0_4_0  s_0_0_0|_{t} + y_ab * (f_0_3_0  s_0_0_0|_{t}
            const double B_0_3_0_0_1_0 = BRA_g_s[10 * ncart_ket + iket] + ( bAB_y * BRA_f_s[6 * ncart_ket + iket] );

            // (f_0_2_1  p_1_0_0| = (g_1_2_1  s_0_0_0|_{t} + x_ab * (f_0_2_1  s_0_0_0|_{t}
            const double B_0_2_1_1_0_0 = BRA_g_s[7 * ncart_ket + iket] + ( bAB_x * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_2_1  p_0_1_0| = (g_0_3_1  s_0_0_0|_{t} + y_ab * (f_0_2_1  s_0_0_0|_{t}
            const double B_0_2_1_0_1_0 = BRA_g_s[11 * ncart_ket + iket] + ( bAB_y * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_2_1  p_0_0_1| = (g_0_2_2  s_0_0_0|_{t} + z_ab * (f_0_2_1  s_0_0_0|_{t}
            const double B_0_2_1_0_0_1 = BRA_g_s[12 * ncart_ket + iket] + ( bAB_z * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_1_2  p_1_0_0| = (g_1_1_2  s_0_0_0|_{t} + x_ab * (f_0_1_2  s_0_0_0|_{t}
            const double B_0_1_2_1_0_0 = BRA_g_s[8 * ncart_ket + iket] + ( bAB_x * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_1_2  p_0_1_0| = (g_0_2_2  s_0_0_0|_{t} + y_ab * (f_0_1_2  s_0_0_0|_{t}
            const double B_0_1_2_0_1_0 = BRA_g_s[12 * ncart_ket + iket] + ( bAB_y * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_1_2  p_0_0_1| = (g_0_1_3  s_0_0_0|_{t} + z_ab * (f_0_1_2  s_0_0_0|_{t}
            const double B_0_1_2_0_0_1 = BRA_g_s[13 * ncart_ket + iket] + ( bAB_z * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_0_3  p_1_0_0| = (g_1_0_3  s_0_0_0|_{t} + x_ab * (f_0_0_3  s_0_0_0|_{t}
            const double B_0_0_3_1_0_0 = BRA_g_s[9 * ncart_ket + iket] + ( bAB_x * BRA_f_s[9 * ncart_ket + iket] );

            // (f_0_0_3  p_0_1_0| = (g_0_1_3  s_0_0_0|_{t} + y_ab * (f_0_0_3  s_0_0_0|_{t}
            const double B_0_0_3_0_1_0 = BRA_g_s[13 * ncart_ket + iket] + ( bAB_y * BRA_f_s[9 * ncart_ket + iket] );

            // (f_0_0_3  p_0_0_1| = (g_0_0_4  s_0_0_0|_{t} + z_ab * (f_0_0_3  s_0_0_0|_{t}
            const double B_0_0_3_0_0_1 = BRA_g_s[14 * ncart_ket + iket] + ( bAB_z * BRA_f_s[9 * ncart_ket + iket] );

            // (d_2_0_0  d_2_0_0|_{i} = (f_3_0_0  p_1_0_0| + x_ab * (d_2_0_0  p_1_0_0|
            BRA_d_d[0 * ncart_ket + iket] = B_3_0_0_1_0_0 + ( bAB_x * B_2_0_0_1_0_0 );

            // (d_2_0_0  d_1_1_0|_{i} = (f_2_1_0  p_1_0_0| + y_ab * (d_2_0_0  p_1_0_0|
            BRA_d_d[1 * ncart_ket + iket] = B_2_1_0_1_0_0 + ( bAB_y * B_2_0_0_1_0_0 );

            // (d_2_0_0  d_1_0_1|_{i} = (f_2_0_1  p_1_0_0| + z_ab * (d_2_0_0  p_1_0_0|
            BRA_d_d[2 * ncart_ket + iket] = B_2_0_1_1_0_0 + ( bAB_z * B_2_0_0_1_0_0 );

            // (d_2_0_0  d_0_2_0|_{i} = (f_2_1_0  p_0_1_0| + y_ab * (d_2_0_0  p_0_1_0|
            BRA_d_d[3 * ncart_ket + iket] = B_2_1_0_0_1_0 + ( bAB_y * B_2_0_0_0_1_0 );

            // (d_2_0_0  d_0_1_1|_{i} = (f_2_0_1  p_0_1_0| + z_ab * (d_2_0_0  p_0_1_0|
            BRA_d_d[4 * ncart_ket + iket] = B_2_0_1_0_1_0 + ( bAB_z * B_2_0_0_0_1_0 );

            // (d_2_0_0  d_0_0_2|_{i} = (f_2_0_1  p_0_0_1| + z_ab * (d_2_0_0  p_0_0_1|
            BRA_d_d[5 * ncart_ket + iket] = B_2_0_1_0_0_1 + ( bAB_z * B_2_0_0_0_0_1 );

            // (d_1_1_0  d_2_0_0|_{i} = (f_2_1_0  p_1_0_0| + x_ab * (d_1_1_0  p_1_0_0|
            BRA_d_d[6 * ncart_ket + iket] = B_2_1_0_1_0_0 + ( bAB_x * B_1_1_0_1_0_0 );

            // (d_1_1_0  d_1_1_0|_{i} = (f_1_2_0  p_1_0_0| + y_ab * (d_1_1_0  p_1_0_0|
            BRA_d_d[7 * ncart_ket + iket] = B_1_2_0_1_0_0 + ( bAB_y * B_1_1_0_1_0_0 );

            // (d_1_1_0  d_1_0_1|_{i} = (f_1_1_1  p_1_0_0| + z_ab * (d_1_1_0  p_1_0_0|
            BRA_d_d[8 * ncart_ket + iket] = B_1_1_1_1_0_0 + ( bAB_z * B_1_1_0_1_0_0 );

            // (d_1_1_0  d_0_2_0|_{i} = (f_1_2_0  p_0_1_0| + y_ab * (d_1_1_0  p_0_1_0|
            BRA_d_d[9 * ncart_ket + iket] = B_1_2_0_0_1_0 + ( bAB_y * B_1_1_0_0_1_0 );

            // (d_1_1_0  d_0_1_1|_{i} = (f_1_1_1  p_0_1_0| + z_ab * (d_1_1_0  p_0_1_0|
            BRA_d_d[10 * ncart_ket + iket] = B_1_1_1_0_1_0 + ( bAB_z * B_1_1_0_0_1_0 );

            // (d_1_1_0  d_0_0_2|_{i} = (f_1_1_1  p_0_0_1| + z_ab * (d_1_1_0  p_0_0_1|
            BRA_d_d[11 * ncart_ket + iket] = B_1_1_1_0_0_1 + ( bAB_z * B_1_1_0_0_0_1 );

            // (d_1_0_1  d_2_0_0|_{i} = (f_2_0_1  p_1_0_0| + x_ab * (d_1_0_1  p_1_0_0|
            BRA_d_d[12 * ncart_ket + iket] = B_2_0_1_1_0_0 + ( bAB_x * B_1_0_1_1_0_0 );

            // (d_1_0_1  d_1_1_0|_{i} = (f_1_1_1  p_1_0_0| + y_ab * (d_1_0_1  p_1_0_0|
            BRA_d_d[13 * ncart_ket + iket] = B_1_1_1_1_0_0 + ( bAB_y * B_1_0_1_1_0_0 );

            // (d_1_0_1  d_1_0_1|_{i} = (f_1_0_2  p_1_0_0| + z_ab * (d_1_0_1  p_1_0_0|
            BRA_d_d[14 * ncart_ket + iket] = B_1_0_2_1_0_0 + ( bAB_z * B_1_0_1_1_0_0 );

            // (d_1_0_1  d_0_2_0|_{i} = (f_1_1_1  p_0_1_0| + y_ab * (d_1_0_1  p_0_1_0|
            BRA_d_d[15 * ncart_ket + iket] = B_1_1_1_0_1_0 + ( bAB_y * B_1_0_1_0_1_0 );

            // (d_1_0_1  d_0_1_1|_{i} = (f_1_0_2  p_0_1_0| + z_ab * (d_1_0_1  p_0_1_0|
            BRA_d_d[16 * ncart_ket + iket] = B_1_0_2_0_1_0 + ( bAB_z * B_1_0_1_0_1_0 );

            // (d_1_0_1  d_0_0_2|_{i} = (f_1_0_2  p_0_0_1| + z_ab * (d_1_0_1  p_0_0_1|
            BRA_d_d[17 * ncart_ket + iket] = B_1_0_2_0_0_1 + ( bAB_z * B_1_0_1_0_0_1 );

            // (d_0_2_0  d_2_0_0|_{i} = (f_1_2_0  p_1_0_0| + x_ab * (d_0_2_0  p_1_0_0|
            BRA_d_d[18 * ncart_ket + iket] = B_1_2_0_1_0_0 + ( bAB_x * B_0_2_0_1_0_0 );

            // (d_0_2_0  d_1_1_0|_{i} = (f_0_3_0  p_1_0_0| + y_ab * (d_0_2_0  p_1_0_0|
            BRA_d_d[19 * ncart_ket + iket] = B_0_3_0_1_0_0 + ( bAB_y * B_0_2_0_1_0_0 );

            // (d_0_2_0  d_1_0_1|_{i} = (f_0_2_1  p_1_0_0| + z_ab * (d_0_2_0  p_1_0_0|
            BRA_d_d[20 * ncart_ket + iket] = B_0_2_1_1_0_0 + ( bAB_z * B_0_2_0_1_0_0 );

            // (d_0_2_0  d_0_2_0|_{i} = (f_0_3_0  p_0_1_0| + y_ab * (d_0_2_0  p_0_1_0|
            BRA_d_d[21 * ncart_ket + iket] = B_0_3_0_0_1_0 + ( bAB_y * B_0_2_0_0_1_0 );

            // (d_0_2_0  d_0_1_1|_{i} = (f_0_2_1  p_0_1_0| + z_ab * (d_0_2_0  p_0_1_0|
            BRA_d_d[22 * ncart_ket + iket] = B_0_2_1_0_1_0 + ( bAB_z * B_0_2_0_0_1_0 );

            // (d_0_2_0  d_0_0_2|_{i} = (f_0_2_1  p_0_0_1| + z_ab * (d_0_2_0  p_0_0_1|
            BRA_d_d[23 * ncart_ket + iket] = B_0_2_1_0_0_1 + ( bAB_z * B_0_2_0_0_0_1 );

            // (d_0_1_1  d_2_0_0|_{i} = (f_1_1_1  p_1_0_0| + x_ab * (d_0_1_1  p_1_0_0|
            BRA_d_d[24 * ncart_ket + iket] = B_1_1_1_1_0_0 + ( bAB_x * B_0_1_1_1_0_0 );

            // (d_0_1_1  d_1_1_0|_{i} = (f_0_2_1  p_1_0_0| + y_ab * (d_0_1_1  p_1_0_0|
            BRA_d_d[25 * ncart_ket + iket] = B_0_2_1_1_0_0 + ( bAB_y * B_0_1_1_1_0_0 );

            // (d_0_1_1  d_1_0_1|_{i} = (f_0_1_2  p_1_0_0| + z_ab * (d_0_1_1  p_1_0_0|
            BRA_d_d[26 * ncart_ket + iket] = B_0_1_2_1_0_0 + ( bAB_z * B_0_1_1_1_0_0 );

            // (d_0_1_1  d_0_2_0|_{i} = (f_0_2_1  p_0_1_0| + y_ab * (d_0_1_1  p_0_1_0|
            BRA_d_d[27 * ncart_ket + iket] = B_0_2_1_0_1_0 + ( bAB_y * B_0_1_1_0_1_0 );

            // (d_0_1_1  d_0_1_1|_{i} = (f_0_1_2  p_0_1_0| + z_ab * (d_0_1_1  p_0_1_0|
            BRA_d_d[28 * ncart_ket + iket] = B_0_1_2_0_1_0 + ( bAB_z * B_0_1_1_0_1_0 );

            // (d_0_1_1  d_0_0_2|_{i} = (f_0_1_2  p_0_0_1| + z_ab * (d_0_1_1  p_0_0_1|
            BRA_d_d[29 * ncart_ket + iket] = B_0_1_2_0_0_1 + ( bAB_z * B_0_1_1_0_0_1 );

            // (d_0_0_2  d_2_0_0|_{i} = (f_1_0_2  p_1_0_0| + x_ab * (d_0_0_2  p_1_0_0|
            BRA_d_d[30 * ncart_ket + iket] = B_1_0_2_1_0_0 + ( bAB_x * B_0_0_2_1_0_0 );

            // (d_0_0_2  d_1_1_0|_{i} = (f_0_1_2  p_1_0_0| + y_ab * (d_0_0_2  p_1_0_0|
            BRA_d_d[31 * ncart_ket + iket] = B_0_1_2_1_0_0 + ( bAB_y * B_0_0_2_1_0_0 );

            // (d_0_0_2  d_1_0_1|_{i} = (f_0_0_3  p_1_0_0| + z_ab * (d_0_0_2  p_1_0_0|
            BRA_d_d[32 * ncart_ket + iket] = B_0_0_3_1_0_0 + ( bAB_z * B_0_0_2_1_0_0 );

            // (d_0_0_2  d_0_2_0|_{i} = (f_0_1_2  p_0_1_0| + y_ab * (d_0_0_2  p_0_1_0|
            BRA_d_d[33 * ncart_ket + iket] = B_0_1_2_0_1_0 + ( bAB_y * B_0_0_2_0_1_0 );

            // (d_0_0_2  d_0_1_1|_{i} = (f_0_0_3  p_0_1_0| + z_ab * (d_0_0_2  p_0_1_0|
            BRA_d_d[34 * ncart_ket + iket] = B_0_0_3_0_1_0 + ( bAB_z * B_0_0_2_0_1_0 );

            // (d_0_0_2  d_0_0_2|_{i} = (f_0_0_3  p_0_0_1| + z_ab * (d_0_0_2  p_0_0_1|
            BRA_d_d[35 * ncart_ket + iket] = B_0_0_3_0_0_1 + ( bAB_z * B_0_0_2_0_0_1 );

        }


}


    //////////////////////////////////////////////
    // BRA: ( f p |
    // Steps: 30
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_BRA_f_p(
                   double * const restrict BRA_f_s,
                   double * const restrict BRA_f_p,
                   double * const restrict BRA_g_s,
                   const double bAB_x, const double bAB_y, const double bAB_z, const int ncart_ket
                 )
{
    int iket;


        for(iket = 0; iket < ncart_ket; ++iket)
        {
            // (f_3_0_0  p_1_0_0|_{i} = (g_4_0_0  s_0_0_0|_{t} + x_ab * (f_3_0_0  s_0_0_0|_{t}
            BRA_f_p[0 * ncart_ket + iket] = BRA_g_s[0 * ncart_ket + iket] + ( bAB_x * BRA_f_s[0 * ncart_ket + iket] );

            // (f_3_0_0  p_0_1_0|_{i} = (g_3_1_0  s_0_0_0|_{t} + y_ab * (f_3_0_0  s_0_0_0|_{t}
            BRA_f_p[1 * ncart_ket + iket] = BRA_g_s[1 * ncart_ket + iket] + ( bAB_y * BRA_f_s[0 * ncart_ket + iket] );

            // (f_3_0_0  p_0_0_1|_{i} = (g_3_0_1  s_0_0_0|_{t} + z_ab * (f_3_0_0  s_0_0_0|_{t}
            BRA_f_p[2 * ncart_ket + iket] = BRA_g_s[2 * ncart_ket + iket] + ( bAB_z * BRA_f_s[0 * ncart_ket + iket] );

            // (f_2_1_0  p_1_0_0|_{i} = (g_3_1_0  s_0_0_0|_{t} + x_ab * (f_2_1_0  s_0_0_0|_{t}
            BRA_f_p[3 * ncart_ket + iket] = BRA_g_s[1 * ncart_ket + iket] + ( bAB_x * BRA_f_s[1 * ncart_ket + iket] );

            // (f_2_1_0  p_0_1_0|_{i} = (g_2_2_0  s_0_0_0|_{t} + y_ab * (f_2_1_0  s_0_0_0|_{t}
            BRA_f_p[4 * ncart_ket + iket] = BRA_g_s[3 * ncart_ket + iket] + ( bAB_y * BRA_f_s[1 * ncart_ket + iket] );

            // (f_2_1_0  p_0_0_1|_{i} = (g_2_1_1  s_0_0_0|_{t} + z_ab * (f_2_1_0  s_0_0_0|_{t}
            BRA_f_p[5 * ncart_ket + iket] = BRA_g_s[4 * ncart_ket + iket] + ( bAB_z * BRA_f_s[1 * ncart_ket + iket] );

            // (f_2_0_1  p_1_0_0|_{i} = (g_3_0_1  s_0_0_0|_{t} + x_ab * (f_2_0_1  s_0_0_0|_{t}
            BRA_f_p[6 * ncart_ket + iket] = BRA_g_s[2 * ncart_ket + iket] + ( bAB_x * BRA_f_s[2 * ncart_ket + iket] );

            // (f_2_0_1  p_0_1_0|_{i} = (g_2_1_1  s_0_0_0|_{t} + y_ab * (f_2_0_1  s_0_0_0|_{t}
            BRA_f_p[7 * ncart_ket + iket] = BRA_g_s[4 * ncart_ket + iket] + ( bAB_y * BRA_f_s[2 * ncart_ket + iket] );

            // (f_2_0_1  p_0_0_1|_{i} = (g_2_0_2  s_0_0_0|_{t} + z_ab * (f_2_0_1  s_0_0_0|_{t}
            BRA_f_p[8 * ncart_ket + iket] = BRA_g_s[5 * ncart_ket + iket] + ( bAB_z * BRA_f_s[2 * ncart_ket + iket] );

            // (f_1_2_0  p_1_0_0|_{i} = (g_2_2_0  s_0_0_0|_{t} + x_ab * (f_1_2_0  s_0_0_0|_{t}
            BRA_f_p[9 * ncart_ket + iket] = BRA_g_s[3 * ncart_ket + iket] + ( bAB_x * BRA_f_s[3 * ncart_ket + iket] );

            // (f_1_2_0  p_0_1_0|_{i} = (g_1_3_0  s_0_0_0|_{t} + y_ab * (f_1_2_0  s_0_0_0|_{t}
            BRA_f_p[10 * ncart_ket + iket] = BRA_g_s[6 * ncart_ket + iket] + ( bAB_y * BRA_f_s[3 * ncart_ket + iket] );

            // (f_1_2_0  p_0_0_1|_{i} = (g_1_2_1  s_0_0_0|_{t} + z_ab * (f_1_2_0  s_0_0_0|_{t}
            BRA_f_p[11 * ncart_ket + iket] = BRA_g_s[7 * ncart_ket + iket] + ( bAB_z * BRA_f_s[3 * ncart_ket + iket] );

            // (f_1_1_1  p_1_0_0|_{i} = (g_2_1_1  s_0_0_0|_{t} + x_ab * (f_1_1_1  s_0_0_0|_{t}
            BRA_f_p[12 * ncart_ket + iket] = BRA_g_s[4 * ncart_ket + iket] + ( bAB_x * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_1_1  p_0_1_0|_{i} = (g_1_2_1  s_0_0_0|_{t} + y_ab * (f_1_1_1  s_0_0_0|_{t}
            BRA_f_p[13 * ncart_ket + iket] = BRA_g_s[7 * ncart_ket + iket] + ( bAB_y * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_1_1  p_0_0_1|_{i} = (g_1_1_2  s_0_0_0|_{t} + z_ab * (f_1_1_1  s_0_0_0|_{t}
            BRA_f_p[14 * ncart_ket + iket] = BRA_g_s[8 * ncart_ket + iket] + ( bAB_z * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_0_2  p_1_0_0|_{i} = (g_2_0_2  s_0_0_0|_{t} + x_ab * (f_1_0_2  s_0_0_0|_{t}
            BRA_f_p[15 * ncart_ket + iket] = BRA_g_s[5 * ncart_ket + iket] + ( bAB_x * BRA_f_s[5 * ncart_ket + iket] );

            // (f_1_0_2  p_0_1_0|_{i} = (g_1_1_2  s_0_0_0|_{t} + y_ab * (f_1_0_2  s_0_0_0|_{t}
            BRA_f_p[16 * ncart_ket + iket] = BRA_g_s[8 * ncart_ket + iket] + ( bAB_y * BRA_f_s[5 * ncart_ket + iket] );

            // (f_1_0_2  p_0_0_1|_{i} = (g_1_0_3  s_0_0_0|_{t} + z_ab * (f_1_0_2  s_0_0_0|_{t}
            BRA_f_p[17 * ncart_ket + iket] = BRA_g_s[9 * ncart_ket + iket] + ( bAB_z * BRA_f_s[5 * ncart_ket + iket] );

            // (f_0_3_0  p_1_0_0|_{i} = (g_1_3_0  s_0_0_0|_{t} + x_ab * (f_0_3_0  s_0_0_0|_{t}
            BRA_f_p[18 * ncart_ket + iket] = BRA_g_s[6 * ncart_ket + iket] + ( bAB_x * BRA_f_s[6 * ncart_ket + iket] );

            // (f_0_3_0  p_0_1_0|_{i} = (g_0_4_0  s_0_0_0|_{t} + y_ab * (f_0_3_0  s_0_0_0|_{t}
            BRA_f_p[19 * ncart_ket + iket] = BRA_g_s[10 * ncart_ket + iket] + ( bAB_y * BRA_f_s[6 * ncart_ket + iket] );

            // (f_0_3_0  p_0_0_1|_{i} = (g_0_3_1  s_0_0_0|_{t} + z_ab * (f_0_3_0  s_0_0_0|_{t}
            BRA_f_p[20 * ncart_ket + iket] = BRA_g_s[11 * ncart_ket + iket] + ( bAB_z * BRA_f_s[6 * ncart_ket + iket] );

            // (f_0_2_1  p_1_0_0|_{i} = (g_1_2_1  s_0_0_0|_{t} + x_ab * (f_0_2_1  s_0_0_0|_{t}
            BRA_f_p[21 * ncart_ket + iket] = BRA_g_s[7 * ncart_ket + iket] + ( bAB_x * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_2_1  p_0_1_0|_{i} = (g_0_3_1  s_0_0_0|_{t} + y_ab * (f_0_2_1  s_0_0_0|_{t}
            BRA_f_p[22 * ncart_ket + iket] = BRA_g_s[11 * ncart_ket + iket] + ( bAB_y * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_2_1  p_0_0_1|_{i} = (g_0_2_2  s_0_0_0|_{t} + z_ab * (f_0_2_1  s_0_0_0|_{t}
            BRA_f_p[23 * ncart_ket + iket] = BRA_g_s[12 * ncart_ket + iket] + ( bAB_z * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_1_2  p_1_0_0|_{i} = (g_1_1_2  s_0_0_0|_{t} + x_ab * (f_0_1_2  s_0_0_0|_{t}
            BRA_f_p[24 * ncart_ket + iket] = BRA_g_s[8 * ncart_ket + iket] + ( bAB_x * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_1_2  p_0_1_0|_{i} = (g_0_2_2  s_0_0_0|_{t} + y_ab * (f_0_1_2  s_0_0_0|_{t}
            BRA_f_p[25 * ncart_ket + iket] = BRA_g_s[12 * ncart_ket + iket] + ( bAB_y * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_1_2  p_0_0_1|_{i} = (g_0_1_3  s_0_0_0|_{t} + z_ab * (f_0_1_2  s_0_0_0|_{t}
            BRA_f_p[26 * ncart_ket + iket] = BRA_g_s[13 * ncart_ket + iket] + ( bAB_z * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_0_3  p_1_0_0|_{i} = (g_1_0_3  s_0_0_0|_{t} + x_ab * (f_0_0_3  s_0_0_0|_{t}
            BRA_f_p[27 * ncart_ket + iket] = BRA_g_s[9 * ncart_ket + iket] + ( bAB_x * BRA_f_s[9 * ncart_ket + iket] );

            // (f_0_0_3  p_0_1_0|_{i} = (g_0_1_3  s_0_0_0|_{t} + y_ab * (f_0_0_3  s_0_0_0|_{t}
            BRA_f_p[28 * ncart_ket + iket] = BRA_g_s[13 * ncart_ket + iket] + ( bAB_y * BRA_f_s[9 * ncart_ket + iket] );

            // (f_0_0_3  p_0_0_1|_{i} = (g_0_0_4  s_0_0_0|_{t} + z_ab * (f_0_0_3  s_0_0_0|_{t}
            BRA_f_p[29 * ncart_ket + iket] = BRA_g_s[14 * ncart_ket + iket] + ( bAB_z * BRA_f_s[9 * ncart_ket + iket] );

        }


}


    //////////////////////////////////////////////
    // BRA: ( f d |
    // Steps: 129
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_BRA_f_d(
                   double * const restrict BRA_f_s,
                   double * const restrict BRA_f_d,
                   double * const restrict BRA_g_s,
                   double * const restrict BRA_h_s,
                   const double bAB_x, const double bAB_y, const double bAB_z, const int ncart_ket
                 )
{
    int iket;


        for(iket = 0; iket < ncart_ket; ++iket)
        {
            // (f_3_0_0  p_1_0_0| = (g_4_0_0  s_0_0_0|_{t} + x_ab * (f_3_0_0  s_0_0_0|_{t}
            const double B_3_0_0_1_0_0 = BRA_g_s[0 * ncart_ket + iket] + ( bAB_x * BRA_f_s[0 * ncart_ket + iket] );

            // (f_3_0_0  p_0_1_0| = (g_3_1_0  s_0_0_0|_{t} + y_ab * (f_3_0_0  s_0_0_0|_{t}
            const double B_3_0_0_0_1_0 = BRA_g_s[1 * ncart_ket + iket] + ( bAB_y * BRA_f_s[0 * ncart_ket + iket] );

            // (f_3_0_0  p_0_0_1| = (g_3_0_1  s_0_0_0|_{t} + z_ab * (f_3_0_0  s_0_0_0|_{t}
            const double B_3_0_0_0_0_1 = BRA_g_s[2 * ncart_ket + iket] + ( bAB_z * BRA_f_s[0 * ncart_ket + iket] );

            // (f_2_1_0  p_1_0_0| = (g_3_1_0  s_0_0_0|_{t} + x_ab * (f_2_1_0  s_0_0_0|_{t}
            const double B_2_1_0_1_0_0 = BRA_g_s[1 * ncart_ket + iket] + ( bAB_x * BRA_f_s[1 * ncart_ket + iket] );

            // (f_2_1_0  p_0_1_0| = (g_2_2_0  s_0_0_0|_{t} + y_ab * (f_2_1_0  s_0_0_0|_{t}
            const double B_2_1_0_0_1_0 = BRA_g_s[3 * ncart_ket + iket] + ( bAB_y * BRA_f_s[1 * ncart_ket + iket] );

            // (f_2_1_0  p_0_0_1| = (g_2_1_1  s_0_0_0|_{t} + z_ab * (f_2_1_0  s_0_0_0|_{t}
            const double B_2_1_0_0_0_1 = BRA_g_s[4 * ncart_ket + iket] + ( bAB_z * BRA_f_s[1 * ncart_ket + iket] );

            // (f_2_0_1  p_1_0_0| = (g_3_0_1  s_0_0_0|_{t} + x_ab * (f_2_0_1  s_0_0_0|_{t}
            const double B_2_0_1_1_0_0 = BRA_g_s[2 * ncart_ket + iket] + ( bAB_x * BRA_f_s[2 * ncart_ket + iket] );

            // (f_2_0_1  p_0_1_0| = (g_2_1_1  s_0_0_0|_{t} + y_ab * (f_2_0_1  s_0_0_0|_{t}
            const double B_2_0_1_0_1_0 = BRA_g_s[4 * ncart_ket + iket] + ( bAB_y * BRA_f_s[2 * ncart_ket + iket] );

            // (f_2_0_1  p_0_0_1| = (g_2_0_2  s_0_0_0|_{t} + z_ab * (f_2_0_1  s_0_0_0|_{t}
            const double B_2_0_1_0_0_1 = BRA_g_s[5 * ncart_ket + iket] + ( bAB_z * BRA_f_s[2 * ncart_ket + iket] );

            // (f_1_2_0  p_1_0_0| = (g_2_2_0  s_0_0_0|_{t} + x_ab * (f_1_2_0  s_0_0_0|_{t}
            const double B_1_2_0_1_0_0 = BRA_g_s[3 * ncart_ket + iket] + ( bAB_x * BRA_f_s[3 * ncart_ket + iket] );

            // (f_1_2_0  p_0_1_0| = (g_1_3_0  s_0_0_0|_{t} + y_ab * (f_1_2_0  s_0_0_0|_{t}
            const double B_1_2_0_0_1_0 = BRA_g_s[6 * ncart_ket + iket] + ( bAB_y * BRA_f_s[3 * ncart_ket + iket] );

            // (f_1_2_0  p_0_0_1| = (g_1_2_1  s_0_0_0|_{t} + z_ab * (f_1_2_0  s_0_0_0|_{t}
            const double B_1_2_0_0_0_1 = BRA_g_s[7 * ncart_ket + iket] + ( bAB_z * BRA_f_s[3 * ncart_ket + iket] );

            // (f_1_1_1  p_1_0_0| = (g_2_1_1  s_0_0_0|_{t} + x_ab * (f_1_1_1  s_0_0_0|_{t}
            const double B_1_1_1_1_0_0 = BRA_g_s[4 * ncart_ket + iket] + ( bAB_x * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_1_1  p_0_1_0| = (g_1_2_1  s_0_0_0|_{t} + y_ab * (f_1_1_1  s_0_0_0|_{t}
            const double B_1_1_1_0_1_0 = BRA_g_s[7 * ncart_ket + iket] + ( bAB_y * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_1_1  p_0_0_1| = (g_1_1_2  s_0_0_0|_{t} + z_ab * (f_1_1_1  s_0_0_0|_{t}
            const double B_1_1_1_0_0_1 = BRA_g_s[8 * ncart_ket + iket] + ( bAB_z * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_0_2  p_1_0_0| = (g_2_0_2  s_0_0_0|_{t} + x_ab * (f_1_0_2  s_0_0_0|_{t}
            const double B_1_0_2_1_0_0 = BRA_g_s[5 * ncart_ket + iket] + ( bAB_x * BRA_f_s[5 * ncart_ket + iket] );

            // (f_1_0_2  p_0_1_0| = (g_1_1_2  s_0_0_0|_{t} + y_ab * (f_1_0_2  s_0_0_0|_{t}
            const double B_1_0_2_0_1_0 = BRA_g_s[8 * ncart_ket + iket] + ( bAB_y * BRA_f_s[5 * ncart_ket + iket] );

            // (f_1_0_2  p_0_0_1| = (g_1_0_3  s_0_0_0|_{t} + z_ab * (f_1_0_2  s_0_0_0|_{t}
            const double B_1_0_2_0_0_1 = BRA_g_s[9 * ncart_ket + iket] + ( bAB_z * BRA_f_s[5 * ncart_ket + iket] );

            // (f_0_3_0  p_1_0_0| = (g_1_3_0  s_0_0_0|_{t} + x_ab * (f_0_3_0  s_0_0_0|_{t}
            const double B_0_3_0_1_0_0 = BRA_g_s[6 * ncart_ket + iket] + ( bAB_x * BRA_f_s[6 * ncart_ket + iket] );

            // (f_0_3_0  p_0_1_0| = (g_0_4_0  s_0_0_0|_{t} + y_ab * (f_0_3_0  s_0_0_0|_{t}
            const double B_0_3_0_0_1_0 = BRA_g_s[10 * ncart_ket + iket] + ( bAB_y * BRA_f_s[6 * ncart_ket + iket] );

            // (f_0_3_0  p_0_0_1| = (g_0_3_1  s_0_0_0|_{t} + z_ab * (f_0_3_0  s_0_0_0|_{t}
            const double B_0_3_0_0_0_1 = BRA_g_s[11 * ncart_ket + iket] + ( bAB_z * BRA_f_s[6 * ncart_ket + iket] );

            // (f_0_2_1  p_1_0_0| = (g_1_2_1  s_0_0_0|_{t} + x_ab * (f_0_2_1  s_0_0_0|_{t}
            const double B_0_2_1_1_0_0 = BRA_g_s[7 * ncart_ket + iket] + ( bAB_x * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_2_1  p_0_1_0| = (g_0_3_1  s_0_0_0|_{t} + y_ab * (f_0_2_1  s_0_0_0|_{t}
            const double B_0_2_1_0_1_0 = BRA_g_s[11 * ncart_ket + iket] + ( bAB_y * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_2_1  p_0_0_1| = (g_0_2_2  s_0_0_0|_{t} + z_ab * (f_0_2_1  s_0_0_0|_{t}
            const double B_0_2_1_0_0_1 = BRA_g_s[12 * ncart_ket + iket] + ( bAB_z * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_1_2  p_1_0_0| = (g_1_1_2  s_0_0_0|_{t} + x_ab * (f_0_1_2  s_0_0_0|_{t}
            const double B_0_1_2_1_0_0 = BRA_g_s[8 * ncart_ket + iket] + ( bAB_x * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_1_2  p_0_1_0| = (g_0_2_2  s_0_0_0|_{t} + y_ab * (f_0_1_2  s_0_0_0|_{t}
            const double B_0_1_2_0_1_0 = BRA_g_s[12 * ncart_ket + iket] + ( bAB_y * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_1_2  p_0_0_1| = (g_0_1_3  s_0_0_0|_{t} + z_ab * (f_0_1_2  s_0_0_0|_{t}
            const double B_0_1_2_0_0_1 = BRA_g_s[13 * ncart_ket + iket] + ( bAB_z * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_0_3  p_1_0_0| = (g_1_0_3  s_0_0_0|_{t} + x_ab * (f_0_0_3  s_0_0_0|_{t}
            const double B_0_0_3_1_0_0 = BRA_g_s[9 * ncart_ket + iket] + ( bAB_x * BRA_f_s[9 * ncart_ket + iket] );

            // (f_0_0_3  p_0_1_0| = (g_0_1_3  s_0_0_0|_{t} + y_ab * (f_0_0_3  s_0_0_0|_{t}
            const double B_0_0_3_0_1_0 = BRA_g_s[13 * ncart_ket + iket] + ( bAB_y * BRA_f_s[9 * ncart_ket + iket] );

            // (f_0_0_3  p_0_0_1| = (g_0_0_4  s_0_0_0|_{t} + z_ab * (f_0_0_3  s_0_0_0|_{t}
            const double B_0_0_3_0_0_1 = BRA_g_s[14 * ncart_ket + iket] + ( bAB_z * BRA_f_s[9 * ncart_ket + iket] );

            // (g_4_0_0  p_1_0_0| = (h_5_0_0  s_0_0_0|_{t} + x_ab * (g_4_0_0  s_0_0_0|_{t}
            const double B_4_0_0_1_0_0 = BRA_h_s[0 * ncart_ket + iket] + ( bAB_x * BRA_g_s[0 * ncart_ket + iket] );

            // (g_3_1_0  p_1_0_0| = (h_4_1_0  s_0_0_0|_{t} + x_ab * (g_3_1_0  s_0_0_0|_{t}
            const double B_3_1_0_1_0_0 = BRA_h_s[1 * ncart_ket + iket] + ( bAB_x * BRA_g_s[1 * ncart_ket + iket] );

            // (g_3_1_0  p_0_1_0| = (h_3_2_0  s_0_0_0|_{t} + y_ab * (g_3_1_0  s_0_0_0|_{t}
            const double B_3_1_0_0_1_0 = BRA_h_s[3 * ncart_ket + iket] + ( bAB_y * BRA_g_s[1 * ncart_ket + iket] );

            // (g_3_0_1  p_1_0_0| = (h_4_0_1  s_0_0_0|_{t} + x_ab * (g_3_0_1  s_0_0_0|_{t}
            const double B_3_0_1_1_0_0 = BRA_h_s[2 * ncart_ket + iket] + ( bAB_x * BRA_g_s[2 * ncart_ket + iket] );

            // (g_3_0_1  p_0_1_0| = (h_3_1_1  s_0_0_0|_{t} + y_ab * (g_3_0_1  s_0_0_0|_{t}
            const double B_3_0_1_0_1_0 = BRA_h_s[4 * ncart_ket + iket] + ( bAB_y * BRA_g_s[2 * ncart_ket + iket] );

            // (g_3_0_1  p_0_0_1| = (h_3_0_2  s_0_0_0|_{t} + z_ab * (g_3_0_1  s_0_0_0|_{t}
            const double B_3_0_1_0_0_1 = BRA_h_s[5 * ncart_ket + iket] + ( bAB_z * BRA_g_s[2 * ncart_ket + iket] );

            // (g_2_2_0  p_1_0_0| = (h_3_2_0  s_0_0_0|_{t} + x_ab * (g_2_2_0  s_0_0_0|_{t}
            const double B_2_2_0_1_0_0 = BRA_h_s[3 * ncart_ket + iket] + ( bAB_x * BRA_g_s[3 * ncart_ket + iket] );

            // (g_2_2_0  p_0_1_0| = (h_2_3_0  s_0_0_0|_{t} + y_ab * (g_2_2_0  s_0_0_0|_{t}
            const double B_2_2_0_0_1_0 = BRA_h_s[6 * ncart_ket + iket] + ( bAB_y * BRA_g_s[3 * ncart_ket + iket] );

            // (g_2_1_1  p_1_0_0| = (h_3_1_1  s_0_0_0|_{t} + x_ab * (g_2_1_1  s_0_0_0|_{t}
            const double B_2_1_1_1_0_0 = BRA_h_s[4 * ncart_ket + iket] + ( bAB_x * BRA_g_s[4 * ncart_ket + iket] );

            // (g_2_1_1  p_0_1_0| = (h_2_2_1  s_0_0_0|_{t} + y_ab * (g_2_1_1  s_0_0_0|_{t}
            const double B_2_1_1_0_1_0 = BRA_h_s[7 * ncart_ket + iket] + ( bAB_y * BRA_g_s[4 * ncart_ket + iket] );

            // (g_2_1_1  p_0_0_1| = (h_2_1_2  s_0_0_0|_{t} + z_ab * (g_2_1_1  s_0_0_0|_{t}
            const double B_2_1_1_0_0_1 = BRA_h_s[8 * ncart_ket + iket] + ( bAB_z * BRA_g_s[4 * ncart_ket + iket] );

            // (g_2_0_2  p_1_0_0| = (h_3_0_2  s_0_0_0|_{t} + x_ab * (g_2_0_2  s_0_0_0|_{t}
            const double B_2_0_2_1_0_0 = BRA_h_s[5 * ncart_ket + iket] + ( bAB_x * BRA_g_s[5 * ncart_ket + iket] );

            // (g_2_0_2  p_0_1_0| = (h_2_1_2  s_0_0_0|_{t} + y_ab * (g_2_0_2  s_0_0_0|_{t}
            const double B_2_0_2_0_1_0 = BRA_h_s[8 * ncart_ket + iket] + ( bAB_y * BRA_g_s[5 * ncart_ket + iket] );

            // (g_2_0_2  p_0_0_1| = (h_2_0_3  s_0_0_0|_{t} + z_ab * (g_2_0_2  s_0_0_0|_{t}
            const double B_2_0_2_0_0_1 = BRA_h_s[9 * ncart_ket + iket] + ( bAB_z * BRA_g_s[5 * ncart_ket + iket] );

            // (g_1_3_0  p_1_0_0| = (h_2_3_0  s_0_0_0|_{t} + x_ab * (g_1_3_0  s_0_0_0|_{t}
            const double B_1_3_0_1_0_0 = BRA_h_s[6 * ncart_ket + iket] + ( bAB_x * BRA_g_s[6 * ncart_ket + iket] );

            // (g_1_3_0  p_0_1_0| = (h_1_4_0  s_0_0_0|_{t} + y_ab * (g_1_3_0  s_0_0_0|_{t}
            const double B_1_3_0_0_1_0 = BRA_h_s[10 * ncart_ket + iket] + ( bAB_y * BRA_g_s[6 * ncart_ket + iket] );

            // (g_1_2_1  p_1_0_0| = (h_2_2_1  s_0_0_0|_{t} + x_ab * (g_1_2_1  s_0_0_0|_{t}
            const double B_1_2_1_1_0_0 = BRA_h_s[7 * ncart_ket + iket] + ( bAB_x * BRA_g_s[7 * ncart_ket + iket] );

            // (g_1_2_1  p_0_1_0| = (h_1_3_1  s_0_0_0|_{t} + y_ab * (g_1_2_1  s_0_0_0|_{t}
            const double B_1_2_1_0_1_0 = BRA_h_s[11 * ncart_ket + iket] + ( bAB_y * BRA_g_s[7 * ncart_ket + iket] );

            // (g_1_2_1  p_0_0_1| = (h_1_2_2  s_0_0_0|_{t} + z_ab * (g_1_2_1  s_0_0_0|_{t}
            const double B_1_2_1_0_0_1 = BRA_h_s[12 * ncart_ket + iket] + ( bAB_z * BRA_g_s[7 * ncart_ket + iket] );

            // (g_1_1_2  p_1_0_0| = (h_2_1_2  s_0_0_0|_{t} + x_ab * (g_1_1_2  s_0_0_0|_{t}
            const double B_1_1_2_1_0_0 = BRA_h_s[8 * ncart_ket + iket] + ( bAB_x * BRA_g_s[8 * ncart_ket + iket] );

            // (g_1_1_2  p_0_1_0| = (h_1_2_2  s_0_0_0|_{t} + y_ab * (g_1_1_2  s_0_0_0|_{t}
            const double B_1_1_2_0_1_0 = BRA_h_s[12 * ncart_ket + iket] + ( bAB_y * BRA_g_s[8 * ncart_ket + iket] );

            // (g_1_1_2  p_0_0_1| = (h_1_1_3  s_0_0_0|_{t} + z_ab * (g_1_1_2  s_0_0_0|_{t}
            const double B_1_1_2_0_0_1 = BRA_h_s[13 * ncart_ket + iket] + ( bAB_z * BRA_g_s[8 * ncart_ket + iket] );

            // (g_1_0_3  p_1_0_0| = (h_2_0_3  s_0_0_0|_{t} + x_ab * (g_1_0_3  s_0_0_0|_{t}
            const double B_1_0_3_1_0_0 = BRA_h_s[9 * ncart_ket + iket] + ( bAB_x * BRA_g_s[9 * ncart_ket + iket] );

            // (g_1_0_3  p_0_1_0| = (h_1_1_3  s_0_0_0|_{t} + y_ab * (g_1_0_3  s_0_0_0|_{t}
            const double B_1_0_3_0_1_0 = BRA_h_s[13 * ncart_ket + iket] + ( bAB_y * BRA_g_s[9 * ncart_ket + iket] );

            // (g_1_0_3  p_0_0_1| = (h_1_0_4  s_0_0_0|_{t} + z_ab * (g_1_0_3  s_0_0_0|_{t}
            const double B_1_0_3_0_0_1 = BRA_h_s[14 * ncart_ket + iket] + ( bAB_z * BRA_g_s[9 * ncart_ket + iket] );

            // (g_0_4_0  p_1_0_0| = (h_1_4_0  s_0_0_0|_{t} + x_ab * (g_0_4_0  s_0_0_0|_{t}
            const double B_0_4_0_1_0_0 = BRA_h_s[10 * ncart_ket + iket] + ( bAB_x * BRA_g_s[10 * ncart_ket + iket] );

            // (g_0_4_0  p_0_1_0| = (h_0_5_0  s_0_0_0|_{t} + y_ab * (g_0_4_0  s_0_0_0|_{t}
            const double B_0_4_0_0_1_0 = BRA_h_s[15 * ncart_ket + iket] + ( bAB_y * BRA_g_s[10 * ncart_ket + iket] );

            // (g_0_3_1  p_1_0_0| = (h_1_3_1  s_0_0_0|_{t} + x_ab * (g_0_3_1  s_0_0_0|_{t}
            const double B_0_3_1_1_0_0 = BRA_h_s[11 * ncart_ket + iket] + ( bAB_x * BRA_g_s[11 * ncart_ket + iket] );

            // (g_0_3_1  p_0_1_0| = (h_0_4_1  s_0_0_0|_{t} + y_ab * (g_0_3_1  s_0_0_0|_{t}
            const double B_0_3_1_0_1_0 = BRA_h_s[16 * ncart_ket + iket] + ( bAB_y * BRA_g_s[11 * ncart_ket + iket] );

            // (g_0_3_1  p_0_0_1| = (h_0_3_2  s_0_0_0|_{t} + z_ab * (g_0_3_1  s_0_0_0|_{t}
            const double B_0_3_1_0_0_1 = BRA_h_s[17 * ncart_ket + iket] + ( bAB_z * BRA_g_s[11 * ncart_ket + iket] );

            // (g_0_2_2  p_1_0_0| = (h_1_2_2  s_0_0_0|_{t} + x_ab * (g_0_2_2  s_0_0_0|_{t}
            const double B_0_2_2_1_0_0 = BRA_h_s[12 * ncart_ket + iket] + ( bAB_x * BRA_g_s[12 * ncart_ket + iket] );

            // (g_0_2_2  p_0_1_0| = (h_0_3_2  s_0_0_0|_{t} + y_ab * (g_0_2_2  s_0_0_0|_{t}
            const double B_0_2_2_0_1_0 = BRA_h_s[17 * ncart_ket + iket] + ( bAB_y * BRA_g_s[12 * ncart_ket + iket] );

            // (g_0_2_2  p_0_0_1| = (h_0_2_3  s_0_0_0|_{t} + z_ab * (g_0_2_2  s_0_0_0|_{t}
            const double B_0_2_2_0_0_1 = BRA_h_s[18 * ncart_ket + iket] + ( bAB_z * BRA_g_s[12 * ncart_ket + iket] );

            // (g_0_1_3  p_1_0_0| = (h_1_1_3  s_0_0_0|_{t} + x_ab * (g_0_1_3  s_0_0_0|_{t}
            const double B_0_1_3_1_0_0 = BRA_h_s[13 * ncart_ket + iket] + ( bAB_x * BRA_g_s[13 * ncart_ket + iket] );

            // (g_0_1_3  p_0_1_0| = (h_0_2_3  s_0_0_0|_{t} + y_ab * (g_0_1_3  s_0_0_0|_{t}
            const double B_0_1_3_0_1_0 = BRA_h_s[18 * ncart_ket + iket] + ( bAB_y * BRA_g_s[13 * ncart_ket + iket] );

            // (g_0_1_3  p_0_0_1| = (h_0_1_4  s_0_0_0|_{t} + z_ab * (g_0_1_3  s_0_0_0|_{t}
            const double B_0_1_3_0_0_1 = BRA_h_s[19 * ncart_ket + iket] + ( bAB_z * BRA_g_s[13 * ncart_ket + iket] );

            // (g_0_0_4  p_1_0_0| = (h_1_0_4  s_0_0_0|_{t} + x_ab * (g_0_0_4  s_0_0_0|_{t}
            const double B_0_0_4_1_0_0 = BRA_h_s[14 * ncart_ket + iket] + ( bAB_x * BRA_g_s[14 * ncart_ket + iket] );

            // (g_0_0_4  p_0_1_0| = (h_0_1_4  s_0_0_0|_{t} + y_ab * (g_0_0_4  s_0_0_0|_{t}
            const double B_0_0_4_0_1_0 = BRA_h_s[19 * ncart_ket + iket] + ( bAB_y * BRA_g_s[14 * ncart_ket + iket] );

            // (g_0_0_4  p_0_0_1| = (h_0_0_5  s_0_0_0|_{t} + z_ab * (g_0_0_4  s_0_0_0|_{t}
            const double B_0_0_4_0_0_1 = BRA_h_s[20 * ncart_ket + iket] + ( bAB_z * BRA_g_s[14 * ncart_ket + iket] );

            // (f_3_0_0  d_2_0_0|_{i} = (g_4_0_0  p_1_0_0| + x_ab * (f_3_0_0  p_1_0_0|
            BRA_f_d[0 * ncart_ket + iket] = B_4_0_0_1_0_0 + ( bAB_x * B_3_0_0_1_0_0 );

            // (f_3_0_0  d_1_1_0|_{i} = (g_3_1_0  p_1_0_0| + y_ab * (f_3_0_0  p_1_0_0|
            BRA_f_d[1 * ncart_ket + iket] = B_3_1_0_1_0_0 + ( bAB_y * B_3_0_0_1_0_0 );

            // (f_3_0_0  d_1_0_1|_{i} = (g_3_0_1  p_1_0_0| + z_ab * (f_3_0_0  p_1_0_0|
            BRA_f_d[2 * ncart_ket + iket] = B_3_0_1_1_0_0 + ( bAB_z * B_3_0_0_1_0_0 );

            // (f_3_0_0  d_0_2_0|_{i} = (g_3_1_0  p_0_1_0| + y_ab * (f_3_0_0  p_0_1_0|
            BRA_f_d[3 * ncart_ket + iket] = B_3_1_0_0_1_0 + ( bAB_y * B_3_0_0_0_1_0 );

            // (f_3_0_0  d_0_1_1|_{i} = (g_3_0_1  p_0_1_0| + z_ab * (f_3_0_0  p_0_1_0|
            BRA_f_d[4 * ncart_ket + iket] = B_3_0_1_0_1_0 + ( bAB_z * B_3_0_0_0_1_0 );

            // (f_3_0_0  d_0_0_2|_{i} = (g_3_0_1  p_0_0_1| + z_ab * (f_3_0_0  p_0_0_1|
            BRA_f_d[5 * ncart_ket + iket] = B_3_0_1_0_0_1 + ( bAB_z * B_3_0_0_0_0_1 );

            // (f_2_1_0  d_2_0_0|_{i} = (g_3_1_0  p_1_0_0| + x_ab * (f_2_1_0  p_1_0_0|
            BRA_f_d[6 * ncart_ket + iket] = B_3_1_0_1_0_0 + ( bAB_x * B_2_1_0_1_0_0 );

            // (f_2_1_0  d_1_1_0|_{i} = (g_2_2_0  p_1_0_0| + y_ab * (f_2_1_0  p_1_0_0|
            BRA_f_d[7 * ncart_ket + iket] = B_2_2_0_1_0_0 + ( bAB_y * B_2_1_0_1_0_0 );

            // (f_2_1_0  d_1_0_1|_{i} = (g_2_1_1  p_1_0_0| + z_ab * (f_2_1_0  p_1_0_0|
            BRA_f_d[8 * ncart_ket + iket] = B_2_1_1_1_0_0 + ( bAB_z * B_2_1_0_1_0_0 );

            // (f_2_1_0  d_0_2_0|_{i} = (g_2_2_0  p_0_1_0| + y_ab * (f_2_1_0  p_0_1_0|
            BRA_f_d[9 * ncart_ket + iket] = B_2_2_0_0_1_0 + ( bAB_y * B_2_1_0_0_1_0 );

            // (f_2_1_0  d_0_1_1|_{i} = (g_2_1_1  p_0_1_0| + z_ab * (f_2_1_0  p_0_1_0|
            BRA_f_d[10 * ncart_ket + iket] = B_2_1_1_0_1_0 + ( bAB_z * B_2_1_0_0_1_0 );

            // (f_2_1_0  d_0_0_2|_{i} = (g_2_1_1  p_0_0_1| + z_ab * (f_2_1_0  p_0_0_1|
            BRA_f_d[11 * ncart_ket + iket] = B_2_1_1_0_0_1 + ( bAB_z * B_2_1_0_0_0_1 );

            // (f_2_0_1  d_2_0_0|_{i} = (g_3_0_1  p_1_0_0| + x_ab * (f_2_0_1  p_1_0_0|
            BRA_f_d[12 * ncart_ket + iket] = B_3_0_1_1_0_0 + ( bAB_x * B_2_0_1_1_0_0 );

            // (f_2_0_1  d_1_1_0|_{i} = (g_2_1_1  p_1_0_0| + y_ab * (f_2_0_1  p_1_0_0|
            BRA_f_d[13 * ncart_ket + iket] = B_2_1_1_1_0_0 + ( bAB_y * B_2_0_1_1_0_0 );

            // (f_2_0_1  d_1_0_1|_{i} = (g_2_0_2  p_1_0_0| + z_ab * (f_2_0_1  p_1_0_0|
            BRA_f_d[14 * ncart_ket + iket] = B_2_0_2_1_0_0 + ( bAB_z * B_2_0_1_1_0_0 );

            // (f_2_0_1  d_0_2_0|_{i} = (g_2_1_1  p_0_1_0| + y_ab * (f_2_0_1  p_0_1_0|
            BRA_f_d[15 * ncart_ket + iket] = B_2_1_1_0_1_0 + ( bAB_y * B_2_0_1_0_1_0 );

            // (f_2_0_1  d_0_1_1|_{i} = (g_2_0_2  p_0_1_0| + z_ab * (f_2_0_1  p_0_1_0|
            BRA_f_d[16 * ncart_ket + iket] = B_2_0_2_0_1_0 + ( bAB_z * B_2_0_1_0_1_0 );

            // (f_2_0_1  d_0_0_2|_{i} = (g_2_0_2  p_0_0_1| + z_ab * (f_2_0_1  p_0_0_1|
            BRA_f_d[17 * ncart_ket + iket] = B_2_0_2_0_0_1 + ( bAB_z * B_2_0_1_0_0_1 );

            // (f_1_2_0  d_2_0_0|_{i} = (g_2_2_0  p_1_0_0| + x_ab * (f_1_2_0  p_1_0_0|
            BRA_f_d[18 * ncart_ket + iket] = B_2_2_0_1_0_0 + ( bAB_x * B_1_2_0_1_0_0 );

            // (f_1_2_0  d_1_1_0|_{i} = (g_1_3_0  p_1_0_0| + y_ab * (f_1_2_0  p_1_0_0|
            BRA_f_d[19 * ncart_ket + iket] = B_1_3_0_1_0_0 + ( bAB_y * B_1_2_0_1_0_0 );

            // (f_1_2_0  d_1_0_1|_{i} = (g_1_2_1  p_1_0_0| + z_ab * (f_1_2_0  p_1_0_0|
            BRA_f_d[20 * ncart_ket + iket] = B_1_2_1_1_0_0 + ( bAB_z * B_1_2_0_1_0_0 );

            // (f_1_2_0  d_0_2_0|_{i} = (g_1_3_0  p_0_1_0| + y_ab * (f_1_2_0  p_0_1_0|
            BRA_f_d[21 * ncart_ket + iket] = B_1_3_0_0_1_0 + ( bAB_y * B_1_2_0_0_1_0 );

            // (f_1_2_0  d_0_1_1|_{i} = (g_1_2_1  p_0_1_0| + z_ab * (f_1_2_0  p_0_1_0|
            BRA_f_d[22 * ncart_ket + iket] = B_1_2_1_0_1_0 + ( bAB_z * B_1_2_0_0_1_0 );

            // (f_1_2_0  d_0_0_2|_{i} = (g_1_2_1  p_0_0_1| + z_ab * (f_1_2_0  p_0_0_1|
            BRA_f_d[23 * ncart_ket + iket] = B_1_2_1_0_0_1 + ( bAB_z * B_1_2_0_0_0_1 );

            // (f_1_1_1  d_2_0_0|_{i} = (g_2_1_1  p_1_0_0| + x_ab * (f_1_1_1  p_1_0_0|
            BRA_f_d[24 * ncart_ket + iket] = B_2_1_1_1_0_0 + ( bAB_x * B_1_1_1_1_0_0 );

            // (f_1_1_1  d_1_1_0|_{i} = (g_1_2_1  p_1_0_0| + y_ab * (f_1_1_1  p_1_0_0|
            BRA_f_d[25 * ncart_ket + iket] = B_1_2_1_1_0_0 + ( bAB_y * B_1_1_1_1_0_0 );

            // (f_1_1_1  d_1_0_1|_{i} = (g_1_1_2  p_1_0_0| + z_ab * (f_1_1_1  p_1_0_0|
            BRA_f_d[26 * ncart_ket + iket] = B_1_1_2_1_0_0 + ( bAB_z * B_1_1_1_1_0_0 );

            // (f_1_1_1  d_0_2_0|_{i} = (g_1_2_1  p_0_1_0| + y_ab * (f_1_1_1  p_0_1_0|
            BRA_f_d[27 * ncart_ket + iket] = B_1_2_1_0_1_0 + ( bAB_y * B_1_1_1_0_1_0 );

            // (f_1_1_1  d_0_1_1|_{i} = (g_1_1_2  p_0_1_0| + z_ab * (f_1_1_1  p_0_1_0|
            BRA_f_d[28 * ncart_ket + iket] = B_1_1_2_0_1_0 + ( bAB_z * B_1_1_1_0_1_0 );

            // (f_1_1_1  d_0_0_2|_{i} = (g_1_1_2  p_0_0_1| + z_ab * (f_1_1_1  p_0_0_1|
            BRA_f_d[29 * ncart_ket + iket] = B_1_1_2_0_0_1 + ( bAB_z * B_1_1_1_0_0_1 );

            // (f_1_0_2  d_2_0_0|_{i} = (g_2_0_2  p_1_0_0| + x_ab * (f_1_0_2  p_1_0_0|
            BRA_f_d[30 * ncart_ket + iket] = B_2_0_2_1_0_0 + ( bAB_x * B_1_0_2_1_0_0 );

            // (f_1_0_2  d_1_1_0|_{i} = (g_1_1_2  p_1_0_0| + y_ab * (f_1_0_2  p_1_0_0|
            BRA_f_d[31 * ncart_ket + iket] = B_1_1_2_1_0_0 + ( bAB_y * B_1_0_2_1_0_0 );

            // (f_1_0_2  d_1_0_1|_{i} = (g_1_0_3  p_1_0_0| + z_ab * (f_1_0_2  p_1_0_0|
            BRA_f_d[32 * ncart_ket + iket] = B_1_0_3_1_0_0 + ( bAB_z * B_1_0_2_1_0_0 );

            // (f_1_0_2  d_0_2_0|_{i} = (g_1_1_2  p_0_1_0| + y_ab * (f_1_0_2  p_0_1_0|
            BRA_f_d[33 * ncart_ket + iket] = B_1_1_2_0_1_0 + ( bAB_y * B_1_0_2_0_1_0 );

            // (f_1_0_2  d_0_1_1|_{i} = (g_1_0_3  p_0_1_0| + z_ab * (f_1_0_2  p_0_1_0|
            BRA_f_d[34 * ncart_ket + iket] = B_1_0_3_0_1_0 + ( bAB_z * B_1_0_2_0_1_0 );

            // (f_1_0_2  d_0_0_2|_{i} = (g_1_0_3  p_0_0_1| + z_ab * (f_1_0_2  p_0_0_1|
            BRA_f_d[35 * ncart_ket + iket] = B_1_0_3_0_0_1 + ( bAB_z * B_1_0_2_0_0_1 );

            // (f_0_3_0  d_2_0_0|_{i} = (g_1_3_0  p_1_0_0| + x_ab * (f_0_3_0  p_1_0_0|
            BRA_f_d[36 * ncart_ket + iket] = B_1_3_0_1_0_0 + ( bAB_x * B_0_3_0_1_0_0 );

            // (f_0_3_0  d_1_1_0|_{i} = (g_0_4_0  p_1_0_0| + y_ab * (f_0_3_0  p_1_0_0|
            BRA_f_d[37 * ncart_ket + iket] = B_0_4_0_1_0_0 + ( bAB_y * B_0_3_0_1_0_0 );

            // (f_0_3_0  d_1_0_1|_{i} = (g_0_3_1  p_1_0_0| + z_ab * (f_0_3_0  p_1_0_0|
            BRA_f_d[38 * ncart_ket + iket] = B_0_3_1_1_0_0 + ( bAB_z * B_0_3_0_1_0_0 );

            // (f_0_3_0  d_0_2_0|_{i} = (g_0_4_0  p_0_1_0| + y_ab * (f_0_3_0  p_0_1_0|
            BRA_f_d[39 * ncart_ket + iket] = B_0_4_0_0_1_0 + ( bAB_y * B_0_3_0_0_1_0 );

            // (f_0_3_0  d_0_1_1|_{i} = (g_0_3_1  p_0_1_0| + z_ab * (f_0_3_0  p_0_1_0|
            BRA_f_d[40 * ncart_ket + iket] = B_0_3_1_0_1_0 + ( bAB_z * B_0_3_0_0_1_0 );

            // (f_0_3_0  d_0_0_2|_{i} = (g_0_3_1  p_0_0_1| + z_ab * (f_0_3_0  p_0_0_1|
            BRA_f_d[41 * ncart_ket + iket] = B_0_3_1_0_0_1 + ( bAB_z * B_0_3_0_0_0_1 );

            // (f_0_2_1  d_2_0_0|_{i} = (g_1_2_1  p_1_0_0| + x_ab * (f_0_2_1  p_1_0_0|
            BRA_f_d[42 * ncart_ket + iket] = B_1_2_1_1_0_0 + ( bAB_x * B_0_2_1_1_0_0 );

            // (f_0_2_1  d_1_1_0|_{i} = (g_0_3_1  p_1_0_0| + y_ab * (f_0_2_1  p_1_0_0|
            BRA_f_d[43 * ncart_ket + iket] = B_0_3_1_1_0_0 + ( bAB_y * B_0_2_1_1_0_0 );

            // (f_0_2_1  d_1_0_1|_{i} = (g_0_2_2  p_1_0_0| + z_ab * (f_0_2_1  p_1_0_0|
            BRA_f_d[44 * ncart_ket + iket] = B_0_2_2_1_0_0 + ( bAB_z * B_0_2_1_1_0_0 );

            // (f_0_2_1  d_0_2_0|_{i} = (g_0_3_1  p_0_1_0| + y_ab * (f_0_2_1  p_0_1_0|
            BRA_f_d[45 * ncart_ket + iket] = B_0_3_1_0_1_0 + ( bAB_y * B_0_2_1_0_1_0 );

            // (f_0_2_1  d_0_1_1|_{i} = (g_0_2_2  p_0_1_0| + z_ab * (f_0_2_1  p_0_1_0|
            BRA_f_d[46 * ncart_ket + iket] = B_0_2_2_0_1_0 + ( bAB_z * B_0_2_1_0_1_0 );

            // (f_0_2_1  d_0_0_2|_{i} = (g_0_2_2  p_0_0_1| + z_ab * (f_0_2_1  p_0_0_1|
            BRA_f_d[47 * ncart_ket + iket] = B_0_2_2_0_0_1 + ( bAB_z * B_0_2_1_0_0_1 );

            // (f_0_1_2  d_2_0_0|_{i} = (g_1_1_2  p_1_0_0| + x_ab * (f_0_1_2  p_1_0_0|
            BRA_f_d[48 * ncart_ket + iket] = B_1_1_2_1_0_0 + ( bAB_x * B_0_1_2_1_0_0 );

            // (f_0_1_2  d_1_1_0|_{i} = (g_0_2_2  p_1_0_0| + y_ab * (f_0_1_2  p_1_0_0|
            BRA_f_d[49 * ncart_ket + iket] = B_0_2_2_1_0_0 + ( bAB_y * B_0_1_2_1_0_0 );

            // (f_0_1_2  d_1_0_1|_{i} = (g_0_1_3  p_1_0_0| + z_ab * (f_0_1_2  p_1_0_0|
            BRA_f_d[50 * ncart_ket + iket] = B_0_1_3_1_0_0 + ( bAB_z * B_0_1_2_1_0_0 );

            // (f_0_1_2  d_0_2_0|_{i} = (g_0_2_2  p_0_1_0| + y_ab * (f_0_1_2  p_0_1_0|
            BRA_f_d[51 * ncart_ket + iket] = B_0_2_2_0_1_0 + ( bAB_y * B_0_1_2_0_1_0 );

            // (f_0_1_2  d_0_1_1|_{i} = (g_0_1_3  p_0_1_0| + z_ab * (f_0_1_2  p_0_1_0|
            BRA_f_d[52 * ncart_ket + iket] = B_0_1_3_0_1_0 + ( bAB_z * B_0_1_2_0_1_0 );

            // (f_0_1_2  d_0_0_2|_{i} = (g_0_1_3  p_0_0_1| + z_ab * (f_0_1_2  p_0_0_1|
            BRA_f_d[53 * ncart_ket + iket] = B_0_1_3_0_0_1 + ( bAB_z * B_0_1_2_0_0_1 );

            // (f_0_0_3  d_2_0_0|_{i} = (g_1_0_3  p_1_0_0| + x_ab * (f_0_0_3  p_1_0_0|
            BRA_f_d[54 * ncart_ket + iket] = B_1_0_3_1_0_0 + ( bAB_x * B_0_0_3_1_0_0 );

            // (f_0_0_3  d_1_1_0|_{i} = (g_0_1_3  p_1_0_0| + y_ab * (f_0_0_3  p_1_0_0|
            BRA_f_d[55 * ncart_ket + iket] = B_0_1_3_1_0_0 + ( bAB_y * B_0_0_3_1_0_0 );

            // (f_0_0_3  d_1_0_1|_{i} = (g_0_0_4  p_1_0_0| + z_ab * (f_0_0_3  p_1_0_0|
            BRA_f_d[56 * ncart_ket + iket] = B_0_0_4_1_0_0 + ( bAB_z * B_0_0_3_1_0_0 );

            // (f_0_0_3  d_0_2_0|_{i} = (g_0_1_3  p_0_1_0| + y_ab * (f_0_0_3  p_0_1_0|
            BRA_f_d[57 * ncart_ket + iket] = B_0_1_3_0_1_0 + ( bAB_y * B_0_0_3_0_1_0 );

            // (f_0_0_3  d_0_1_1|_{i} = (g_0_0_4  p_0_1_0| + z_ab * (f_0_0_3  p_0_1_0|
            BRA_f_d[58 * ncart_ket + iket] = B_0_0_4_0_1_0 + ( bAB_z * B_0_0_3_0_1_0 );

            // (f_0_0_3  d_0_0_2|_{i} = (g_0_0_4  p_0_0_1| + z_ab * (f_0_0_3  p_0_0_1|
            BRA_f_d[59 * ncart_ket + iket] = B_0_0_4_0_0_1 + ( bAB_z * B_0_0_3_0_0_1 );

        }


}


    //////////////////////////////////////////////
    // BRA: ( f f |
    // Steps: 319
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_BRA_f_f(
                   double * const restrict BRA_f_s,
                   double * const restrict BRA_f_f,
                   double * const restrict BRA_g_s,
                   double * const restrict BRA_h_s,
                   double * const restrict BRA_i_s,
                   const double bAB_x, const double bAB_y, const double bAB_z, const int ncart_ket
                 )
{
    int iket;


        for(iket = 0; iket < ncart_ket; ++iket)
        {
            // (f_3_0_0  p_1_0_0| = (g_4_0_0  s_0_0_0|_{t} + x_ab * (f_3_0_0  s_0_0_0|_{t}
            const double B_3_0_0_1_0_0 = BRA_g_s[0 * ncart_ket + iket] + ( bAB_x * BRA_f_s[0 * ncart_ket + iket] );

            // (f_3_0_0  p_0_1_0| = (g_3_1_0  s_0_0_0|_{t} + y_ab * (f_3_0_0  s_0_0_0|_{t}
            const double B_3_0_0_0_1_0 = BRA_g_s[1 * ncart_ket + iket] + ( bAB_y * BRA_f_s[0 * ncart_ket + iket] );

            // (f_3_0_0  p_0_0_1| = (g_3_0_1  s_0_0_0|_{t} + z_ab * (f_3_0_0  s_0_0_0|_{t}
            const double B_3_0_0_0_0_1 = BRA_g_s[2 * ncart_ket + iket] + ( bAB_z * BRA_f_s[0 * ncart_ket + iket] );

            // (f_2_1_0  p_1_0_0| = (g_3_1_0  s_0_0_0|_{t} + x_ab * (f_2_1_0  s_0_0_0|_{t}
            const double B_2_1_0_1_0_0 = BRA_g_s[1 * ncart_ket + iket] + ( bAB_x * BRA_f_s[1 * ncart_ket + iket] );

            // (f_2_1_0  p_0_1_0| = (g_2_2_0  s_0_0_0|_{t} + y_ab * (f_2_1_0  s_0_0_0|_{t}
            const double B_2_1_0_0_1_0 = BRA_g_s[3 * ncart_ket + iket] + ( bAB_y * BRA_f_s[1 * ncart_ket + iket] );

            // (f_2_1_0  p_0_0_1| = (g_2_1_1  s_0_0_0|_{t} + z_ab * (f_2_1_0  s_0_0_0|_{t}
            const double B_2_1_0_0_0_1 = BRA_g_s[4 * ncart_ket + iket] + ( bAB_z * BRA_f_s[1 * ncart_ket + iket] );

            // (f_2_0_1  p_1_0_0| = (g_3_0_1  s_0_0_0|_{t} + x_ab * (f_2_0_1  s_0_0_0|_{t}
            const double B_2_0_1_1_0_0 = BRA_g_s[2 * ncart_ket + iket] + ( bAB_x * BRA_f_s[2 * ncart_ket + iket] );

            // (f_2_0_1  p_0_1_0| = (g_2_1_1  s_0_0_0|_{t} + y_ab * (f_2_0_1  s_0_0_0|_{t}
            const double B_2_0_1_0_1_0 = BRA_g_s[4 * ncart_ket + iket] + ( bAB_y * BRA_f_s[2 * ncart_ket + iket] );

            // (f_2_0_1  p_0_0_1| = (g_2_0_2  s_0_0_0|_{t} + z_ab * (f_2_0_1  s_0_0_0|_{t}
            const double B_2_0_1_0_0_1 = BRA_g_s[5 * ncart_ket + iket] + ( bAB_z * BRA_f_s[2 * ncart_ket + iket] );

            // (f_1_2_0  p_1_0_0| = (g_2_2_0  s_0_0_0|_{t} + x_ab * (f_1_2_0  s_0_0_0|_{t}
            const double B_1_2_0_1_0_0 = BRA_g_s[3 * ncart_ket + iket] + ( bAB_x * BRA_f_s[3 * ncart_ket + iket] );

            // (f_1_2_0  p_0_1_0| = (g_1_3_0  s_0_0_0|_{t} + y_ab * (f_1_2_0  s_0_0_0|_{t}
            const double B_1_2_0_0_1_0 = BRA_g_s[6 * ncart_ket + iket] + ( bAB_y * BRA_f_s[3 * ncart_ket + iket] );

            // (f_1_2_0  p_0_0_1| = (g_1_2_1  s_0_0_0|_{t} + z_ab * (f_1_2_0  s_0_0_0|_{t}
            const double B_1_2_0_0_0_1 = BRA_g_s[7 * ncart_ket + iket] + ( bAB_z * BRA_f_s[3 * ncart_ket + iket] );

            // (f_1_1_1  p_1_0_0| = (g_2_1_1  s_0_0_0|_{t} + x_ab * (f_1_1_1  s_0_0_0|_{t}
            const double B_1_1_1_1_0_0 = BRA_g_s[4 * ncart_ket + iket] + ( bAB_x * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_1_1  p_0_1_0| = (g_1_2_1  s_0_0_0|_{t} + y_ab * (f_1_1_1  s_0_0_0|_{t}
            const double B_1_1_1_0_1_0 = BRA_g_s[7 * ncart_ket + iket] + ( bAB_y * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_1_1  p_0_0_1| = (g_1_1_2  s_0_0_0|_{t} + z_ab * (f_1_1_1  s_0_0_0|_{t}
            const double B_1_1_1_0_0_1 = BRA_g_s[8 * ncart_ket + iket] + ( bAB_z * BRA_f_s[4 * ncart_ket + iket] );

            // (f_1_0_2  p_1_0_0| = (g_2_0_2  s_0_0_0|_{t} + x_ab * (f_1_0_2  s_0_0_0|_{t}
            const double B_1_0_2_1_0_0 = BRA_g_s[5 * ncart_ket + iket] + ( bAB_x * BRA_f_s[5 * ncart_ket + iket] );

            // (f_1_0_2  p_0_1_0| = (g_1_1_2  s_0_0_0|_{t} + y_ab * (f_1_0_2  s_0_0_0|_{t}
            const double B_1_0_2_0_1_0 = BRA_g_s[8 * ncart_ket + iket] + ( bAB_y * BRA_f_s[5 * ncart_ket + iket] );

            // (f_1_0_2  p_0_0_1| = (g_1_0_3  s_0_0_0|_{t} + z_ab * (f_1_0_2  s_0_0_0|_{t}
            const double B_1_0_2_0_0_1 = BRA_g_s[9 * ncart_ket + iket] + ( bAB_z * BRA_f_s[5 * ncart_ket + iket] );

            // (f_0_3_0  p_1_0_0| = (g_1_3_0  s_0_0_0|_{t} + x_ab * (f_0_3_0  s_0_0_0|_{t}
            const double B_0_3_0_1_0_0 = BRA_g_s[6 * ncart_ket + iket] + ( bAB_x * BRA_f_s[6 * ncart_ket + iket] );

            // (f_0_3_0  p_0_1_0| = (g_0_4_0  s_0_0_0|_{t} + y_ab * (f_0_3_0  s_0_0_0|_{t}
            const double B_0_3_0_0_1_0 = BRA_g_s[10 * ncart_ket + iket] + ( bAB_y * BRA_f_s[6 * ncart_ket + iket] );

            // (f_0_3_0  p_0_0_1| = (g_0_3_1  s_0_0_0|_{t} + z_ab * (f_0_3_0  s_0_0_0|_{t}
            const double B_0_3_0_0_0_1 = BRA_g_s[11 * ncart_ket + iket] + ( bAB_z * BRA_f_s[6 * ncart_ket + iket] );

            // (f_0_2_1  p_1_0_0| = (g_1_2_1  s_0_0_0|_{t} + x_ab * (f_0_2_1  s_0_0_0|_{t}
            const double B_0_2_1_1_0_0 = BRA_g_s[7 * ncart_ket + iket] + ( bAB_x * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_2_1  p_0_1_0| = (g_0_3_1  s_0_0_0|_{t} + y_ab * (f_0_2_1  s_0_0_0|_{t}
            const double B_0_2_1_0_1_0 = BRA_g_s[11 * ncart_ket + iket] + ( bAB_y * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_2_1  p_0_0_1| = (g_0_2_2  s_0_0_0|_{t} + z_ab * (f_0_2_1  s_0_0_0|_{t}
            const double B_0_2_1_0_0_1 = BRA_g_s[12 * ncart_ket + iket] + ( bAB_z * BRA_f_s[7 * ncart_ket + iket] );

            // (f_0_1_2  p_1_0_0| = (g_1_1_2  s_0_0_0|_{t} + x_ab * (f_0_1_2  s_0_0_0|_{t}
            const double B_0_1_2_1_0_0 = BRA_g_s[8 * ncart_ket + iket] + ( bAB_x * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_1_2  p_0_1_0| = (g_0_2_2  s_0_0_0|_{t} + y_ab * (f_0_1_2  s_0_0_0|_{t}
            const double B_0_1_2_0_1_0 = BRA_g_s[12 * ncart_ket + iket] + ( bAB_y * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_1_2  p_0_0_1| = (g_0_1_3  s_0_0_0|_{t} + z_ab * (f_0_1_2  s_0_0_0|_{t}
            const double B_0_1_2_0_0_1 = BRA_g_s[13 * ncart_ket + iket] + ( bAB_z * BRA_f_s[8 * ncart_ket + iket] );

            // (f_0_0_3  p_1_0_0| = (g_1_0_3  s_0_0_0|_{t} + x_ab * (f_0_0_3  s_0_0_0|_{t}
            const double B_0_0_3_1_0_0 = BRA_g_s[9 * ncart_ket + iket] + ( bAB_x * BRA_f_s[9 * ncart_ket + iket] );

            // (f_0_0_3  p_0_1_0| = (g_0_1_3  s_0_0_0|_{t} + y_ab * (f_0_0_3  s_0_0_0|_{t}
            const double B_0_0_3_0_1_0 = BRA_g_s[13 * ncart_ket + iket] + ( bAB_y * BRA_f_s[9 * ncart_ket + iket] );

            // (f_0_0_3  p_0_0_1| = (g_0_0_4  s_0_0_0|_{t} + z_ab * (f_0_0_3  s_0_0_0|_{t}
            const double B_0_0_3_0_0_1 = BRA_g_s[14 * ncart_ket + iket] + ( bAB_z * BRA_f_s[9 * ncart_ket + iket] );

            // (g_4_0_0  p_1_0_0| = (h_5_0_0  s_0_0_0|_{t} + x_ab * (g_4_0_0  s_0_0_0|_{t}
            const double B_4_0_0_1_0_0 = BRA_h_s[0 * ncart_ket + iket] + ( bAB_x * BRA_g_s[0 * ncart_ket + iket] );

            // (g_4_0_0  p_0_1_0| = (h_4_1_0  s_0_0_0|_{t} + y_ab * (g_4_0_0  s_0_0_0|_{t}
            const double B_4_0_0_0_1_0 = BRA_h_s[1 * ncart_ket + iket] + ( bAB_y * BRA_g_s[0 * ncart_ket + iket] );

            // (g_4_0_0  p_0_0_1| = (h_4_0_1  s_0_0_0|_{t} + z_ab * (g_4_0_0  s_0_0_0|_{t}
            const double B_4_0_0_0_0_1 = BRA_h_s[2 * ncart_ket + iket] + ( bAB_z * BRA_g_s[0 * ncart_ket + iket] );

            // (g_3_1_0  p_1_0_0| = (h_4_1_0  s_0_0_0|_{t} + x_ab * (g_3_1_0  s_0_0_0|_{t}
            const double B_3_1_0_1_0_0 = BRA_h_s[1 * ncart_ket + iket] + ( bAB_x * BRA_g_s[1 * ncart_ket + iket] );

            // (g_3_1_0  p_0_1_0| = (h_3_2_0  s_0_0_0|_{t} + y_ab * (g_3_1_0  s_0_0_0|_{t}
            const double B_3_1_0_0_1_0 = BRA_h_s[3 * ncart_ket + iket] + ( bAB_y * BRA_g_s[1 * ncart_ket + iket] );

            // (g_3_1_0  p_0_0_1| = (h_3_1_1  s_0_0_0|_{t} + z_ab * (g_3_1_0  s_0_0_0|_{t}
            const double B_3_1_0_0_0_1 = BRA_h_s[4 * ncart_ket + iket] + ( bAB_z * BRA_g_s[1 * ncart_ket + iket] );

            // (g_3_0_1  p_1_0_0| = (h_4_0_1  s_0_0_0|_{t} + x_ab * (g_3_0_1  s_0_0_0|_{t}
            const double B_3_0_1_1_0_0 = BRA_h_s[2 * ncart_ket + iket] + ( bAB_x * BRA_g_s[2 * ncart_ket + iket] );

            // (g_3_0_1  p_0_1_0| = (h_3_1_1  s_0_0_0|_{t} + y_ab * (g_3_0_1  s_0_0_0|_{t}
            const double B_3_0_1_0_1_0 = BRA_h_s[4 * ncart_ket + iket] + ( bAB_y * BRA_g_s[2 * ncart_ket + iket] );

            // (g_3_0_1  p_0_0_1| = (h_3_0_2  s_0_0_0|_{t} + z_ab * (g_3_0_1  s_0_0_0|_{t}
            const double B_3_0_1_0_0_1 = BRA_h_s[5 * ncart_ket + iket] + ( bAB_z * BRA_g_s[2 * ncart_ket + iket] );

            // (g_2_2_0  p_1_0_0| = (h_3_2_0  s_0_0_0|_{t} + x_ab * (g_2_2_0  s_0_0_0|_{t}
            const double B_2_2_0_1_0_0 = BRA_h_s[3 * ncart_ket + iket] + ( bAB_x * BRA_g_s[3 * ncart_ket + iket] );

            // (g_2_2_0  p_0_1_0| = (h_2_3_0  s_0_0_0|_{t} + y_ab * (g_2_2_0  s_0_0_0|_{t}
            const double B_2_2_0_0_1_0 = BRA_h_s[6 * ncart_ket + iket] + ( bAB_y * BRA_g_s[3 * ncart_ket + iket] );

            // (g_2_2_0  p_0_0_1| = (h_2_2_1  s_0_0_0|_{t} + z_ab * (g_2_2_0  s_0_0_0|_{t}
            const double B_2_2_0_0_0_1 = BRA_h_s[7 * ncart_ket + iket] + ( bAB_z * BRA_g_s[3 * ncart_ket + iket] );

            // (g_2_1_1  p_1_0_0| = (h_3_1_1  s_0_0_0|_{t} + x_ab * (g_2_1_1  s_0_0_0|_{t}
            const double B_2_1_1_1_0_0 = BRA_h_s[4 * ncart_ket + iket] + ( bAB_x * BRA_g_s[4 * ncart_ket + iket] );

            // (g_2_1_1  p_0_1_0| = (h_2_2_1  s_0_0_0|_{t} + y_ab * (g_2_1_1  s_0_0_0|_{t}
            const double B_2_1_1_0_1_0 = BRA_h_s[7 * ncart_ket + iket] + ( bAB_y * BRA_g_s[4 * ncart_ket + iket] );

            // (g_2_1_1  p_0_0_1| = (h_2_1_2  s_0_0_0|_{t} + z_ab * (g_2_1_1  s_0_0_0|_{t}
            const double B_2_1_1_0_0_1 = BRA_h_s[8 * ncart_ket + iket] + ( bAB_z * BRA_g_s[4 * ncart_ket + iket] );

            // (g_2_0_2  p_1_0_0| = (h_3_0_2  s_0_0_0|_{t} + x_ab * (g_2_0_2  s_0_0_0|_{t}
            const double B_2_0_2_1_0_0 = BRA_h_s[5 * ncart_ket + iket] + ( bAB_x * BRA_g_s[5 * ncart_ket + iket] );

            // (g_2_0_2  p_0_1_0| = (h_2_1_2  s_0_0_0|_{t} + y_ab * (g_2_0_2  s_0_0_0|_{t}
            const double B_2_0_2_0_1_0 = BRA_h_s[8 * ncart_ket + iket] + ( bAB_y * BRA_g_s[5 * ncart_ket + iket] );

            // (g_2_0_2  p_0_0_1| = (h_2_0_3  s_0_0_0|_{t} + z_ab * (g_2_0_2  s_0_0_0|_{t}
            const double B_2_0_2_0_0_1 = BRA_h_s[9 * ncart_ket + iket] + ( bAB_z * BRA_g_s[5 * ncart_ket + iket] );

            // (g_1_3_0  p_1_0_0| = (h_2_3_0  s_0_0_0|_{t} + x_ab * (g_1_3_0  s_0_0_0|_{t}
            const double B_1_3_0_1_0_0 = BRA_h_s[6 * ncart_ket + iket] + ( bAB_x * BRA_g_s[6 * ncart_ket + iket] );

            // (g_1_3_0  p_0_1_0| = (h_1_4_0  s_0_0_0|_{t} + y_ab * (g_1_3_0  s_0_0_0|_{t}
            const double B_1_3_0_0_1_0 = BRA_h_s[10 * ncart_ket + iket] + ( bAB_y * BRA_g_s[6 * ncart_ket + iket] );

            // (g_1_3_0  p_0_0_1| = (h_1_3_1  s_0_0_0|_{t} + z_ab * (g_1_3_0  s_0_0_0|_{t}
            const double B_1_3_0_0_0_1 = BRA_h_s[11 * ncart_ket + iket] + ( bAB_z * BRA_g_s[6 * ncart_ket + iket] );

            // (g_1_2_1  p_1_0_0| = (h_2_2_1  s_0_0_0|_{t} + x_ab * (g_1_2_1  s_0_0_0|_{t}
            const double B_1_2_1_1_0_0 = BRA_h_s[7 * ncart_ket + iket] + ( bAB_x * BRA_g_s[7 * ncart_ket + iket] );

            // (g_1_2_1  p_0_1_0| = (h_1_3_1  s_0_0_0|_{t} + y_ab * (g_1_2_1  s_0_0_0|_{t}
            const double B_1_2_1_0_1_0 = BRA_h_s[11 * ncart_ket + iket] + ( bAB_y * BRA_g_s[7 * ncart_ket + iket] );

            // (g_1_2_1  p_0_0_1| = (h_1_2_2  s_0_0_0|_{t} + z_ab * (g_1_2_1  s_0_0_0|_{t}
            const double B_1_2_1_0_0_1 = BRA_h_s[12 * ncart_ket + iket] + ( bAB_z * BRA_g_s[7 * ncart_ket + iket] );

            // (g_1_1_2  p_1_0_0| = (h_2_1_2  s_0_0_0|_{t} + x_ab * (g_1_1_2  s_0_0_0|_{t}
            const double B_1_1_2_1_0_0 = BRA_h_s[8 * ncart_ket + iket] + ( bAB_x * BRA_g_s[8 * ncart_ket + iket] );

            // (g_1_1_2  p_0_1_0| = (h_1_2_2  s_0_0_0|_{t} + y_ab * (g_1_1_2  s_0_0_0|_{t}
            const double B_1_1_2_0_1_0 = BRA_h_s[12 * ncart_ket + iket] + ( bAB_y * BRA_g_s[8 * ncart_ket + iket] );

            // (g_1_1_2  p_0_0_1| = (h_1_1_3  s_0_0_0|_{t} + z_ab * (g_1_1_2  s_0_0_0|_{t}
            const double B_1_1_2_0_0_1 = BRA_h_s[13 * ncart_ket + iket] + ( bAB_z * BRA_g_s[8 * ncart_ket + iket] );

            // (g_1_0_3  p_1_0_0| = (h_2_0_3  s_0_0_0|_{t} + x_ab * (g_1_0_3  s_0_0_0|_{t}
            const double B_1_0_3_1_0_0 = BRA_h_s[9 * ncart_ket + iket] + ( bAB_x * BRA_g_s[9 * ncart_ket + iket] );

            // (g_1_0_3  p_0_1_0| = (h_1_1_3  s_0_0_0|_{t} + y_ab * (g_1_0_3  s_0_0_0|_{t}
            const double B_1_0_3_0_1_0 = BRA_h_s[13 * ncart_ket + iket] + ( bAB_y * BRA_g_s[9 * ncart_ket + iket] );

            // (g_1_0_3  p_0_0_1| = (h_1_0_4  s_0_0_0|_{t} + z_ab * (g_1_0_3  s_0_0_0|_{t}
            const double B_1_0_3_0_0_1 = BRA_h_s[14 * ncart_ket + iket] + ( bAB_z * BRA_g_s[9 * ncart_ket + iket] );

            // (g_0_4_0  p_1_0_0| = (h_1_4_0  s_0_0_0|_{t} + x_ab * (g_0_4_0  s_0_0_0|_{t}
            const double B_0_4_0_1_0_0 = BRA_h_s[10 * ncart_ket + iket] + ( bAB_x * BRA_g_s[10 * ncart_ket + iket] );

            // (g_0_4_0  p_0_1_0| = (h_0_5_0  s_0_0_0|_{t} + y_ab * (g_0_4_0  s_0_0_0|_{t}
            const double B_0_4_0_0_1_0 = BRA_h_s[15 * ncart_ket + iket] + ( bAB_y * BRA_g_s[10 * ncart_ket + iket] );

            // (g_0_4_0  p_0_0_1| = (h_0_4_1  s_0_0_0|_{t} + z_ab * (g_0_4_0  s_0_0_0|_{t}
            const double B_0_4_0_0_0_1 = BRA_h_s[16 * ncart_ket + iket] + ( bAB_z * BRA_g_s[10 * ncart_ket + iket] );

            // (g_0_3_1  p_1_0_0| = (h_1_3_1  s_0_0_0|_{t} + x_ab * (g_0_3_1  s_0_0_0|_{t}
            const double B_0_3_1_1_0_0 = BRA_h_s[11 * ncart_ket + iket] + ( bAB_x * BRA_g_s[11 * ncart_ket + iket] );

            // (g_0_3_1  p_0_1_0| = (h_0_4_1  s_0_0_0|_{t} + y_ab * (g_0_3_1  s_0_0_0|_{t}
            const double B_0_3_1_0_1_0 = BRA_h_s[16 * ncart_ket + iket] + ( bAB_y * BRA_g_s[11 * ncart_ket + iket] );

            // (g_0_3_1  p_0_0_1| = (h_0_3_2  s_0_0_0|_{t} + z_ab * (g_0_3_1  s_0_0_0|_{t}
            const double B_0_3_1_0_0_1 = BRA_h_s[17 * ncart_ket + iket] + ( bAB_z * BRA_g_s[11 * ncart_ket + iket] );

            // (g_0_2_2  p_1_0_0| = (h_1_2_2  s_0_0_0|_{t} + x_ab * (g_0_2_2  s_0_0_0|_{t}
            const double B_0_2_2_1_0_0 = BRA_h_s[12 * ncart_ket + iket] + ( bAB_x * BRA_g_s[12 * ncart_ket + iket] );

            // (g_0_2_2  p_0_1_0| = (h_0_3_2  s_0_0_0|_{t} + y_ab * (g_0_2_2  s_0_0_0|_{t}
            const double B_0_2_2_0_1_0 = BRA_h_s[17 * ncart_ket + iket] + ( bAB_y * BRA_g_s[12 * ncart_ket + iket] );

            // (g_0_2_2  p_0_0_1| = (h_0_2_3  s_0_0_0|_{t} + z_ab * (g_0_2_2  s_0_0_0|_{t}
            const double B_0_2_2_0_0_1 = BRA_h_s[18 * ncart_ket + iket] + ( bAB_z * BRA_g_s[12 * ncart_ket + iket] );

            // (g_0_1_3  p_1_0_0| = (h_1_1_3  s_0_0_0|_{t} + x_ab * (g_0_1_3  s_0_0_0|_{t}
            const double B_0_1_3_1_0_0 = BRA_h_s[13 * ncart_ket + iket] + ( bAB_x * BRA_g_s[13 * ncart_ket + iket] );

            // (g_0_1_3  p_0_1_0| = (h_0_2_3  s_0_0_0|_{t} + y_ab * (g_0_1_3  s_0_0_0|_{t}
            const double B_0_1_3_0_1_0 = BRA_h_s[18 * ncart_ket + iket] + ( bAB_y * BRA_g_s[13 * ncart_ket + iket] );

            // (g_0_1_3  p_0_0_1| = (h_0_1_4  s_0_0_0|_{t} + z_ab * (g_0_1_3  s_0_0_0|_{t}
            const double B_0_1_3_0_0_1 = BRA_h_s[19 * ncart_ket + iket] + ( bAB_z * BRA_g_s[13 * ncart_ket + iket] );

            // (g_0_0_4  p_1_0_0| = (h_1_0_4  s_0_0_0|_{t} + x_ab * (g_0_0_4  s_0_0_0|_{t}
            const double B_0_0_4_1_0_0 = BRA_h_s[14 * ncart_ket + iket] + ( bAB_x * BRA_g_s[14 * ncart_ket + iket] );

            // (g_0_0_4  p_0_1_0| = (h_0_1_4  s_0_0_0|_{t} + y_ab * (g_0_0_4  s_0_0_0|_{t}
            const double B_0_0_4_0_1_0 = BRA_h_s[19 * ncart_ket + iket] + ( bAB_y * BRA_g_s[14 * ncart_ket + iket] );

            // (g_0_0_4  p_0_0_1| = (h_0_0_5  s_0_0_0|_{t} + z_ab * (g_0_0_4  s_0_0_0|_{t}
            const double B_0_0_4_0_0_1 = BRA_h_s[20 * ncart_ket + iket] + ( bAB_z * BRA_g_s[14 * ncart_ket + iket] );

            // (h_5_0_0  p_1_0_0| = (i_6_0_0  s_0_0_0|_{t} + x_ab * (h_5_0_0  s_0_0_0|_{t}
            const double B_5_0_0_1_0_0 = BRA_i_s[0 * ncart_ket + iket] + ( bAB_x * BRA_h_s[0 * ncart_ket + iket] );

            // (h_4_1_0  p_1_0_0| = (i_5_1_0  s_0_0_0|_{t} + x_ab * (h_4_1_0  s_0_0_0|_{t}
            const double B_4_1_0_1_0_0 = BRA_i_s[1 * ncart_ket + iket] + ( bAB_x * BRA_h_s[1 * ncart_ket + iket] );

            // (h_4_1_0  p_0_1_0| = (i_4_2_0  s_0_0_0|_{t} + y_ab * (h_4_1_0  s_0_0_0|_{t}
            const double B_4_1_0_0_1_0 = BRA_i_s[3 * ncart_ket + iket] + ( bAB_y * BRA_h_s[1 * ncart_ket + iket] );

            // (h_4_0_1  p_1_0_0| = (i_5_0_1  s_0_0_0|_{t} + x_ab * (h_4_0_1  s_0_0_0|_{t}
            const double B_4_0_1_1_0_0 = BRA_i_s[2 * ncart_ket + iket] + ( bAB_x * BRA_h_s[2 * ncart_ket + iket] );

            // (h_4_0_1  p_0_0_1| = (i_4_0_2  s_0_0_0|_{t} + z_ab * (h_4_0_1  s_0_0_0|_{t}
            const double B_4_0_1_0_0_1 = BRA_i_s[5 * ncart_ket + iket] + ( bAB_z * BRA_h_s[2 * ncart_ket + iket] );

            // (h_3_2_0  p_1_0_0| = (i_4_2_0  s_0_0_0|_{t} + x_ab * (h_3_2_0  s_0_0_0|_{t}
            const double B_3_2_0_1_0_0 = BRA_i_s[3 * ncart_ket + iket] + ( bAB_x * BRA_h_s[3 * ncart_ket + iket] );

            // (h_3_2_0  p_0_1_0| = (i_3_3_0  s_0_0_0|_{t} + y_ab * (h_3_2_0  s_0_0_0|_{t}
            const double B_3_2_0_0_1_0 = BRA_i_s[6 * ncart_ket + iket] + ( bAB_y * BRA_h_s[3 * ncart_ket + iket] );

            // (h_3_1_1  p_1_0_0| = (i_4_1_1  s_0_0_0|_{t} + x_ab * (h_3_1_1  s_0_0_0|_{t}
            const double B_3_1_1_1_0_0 = BRA_i_s[4 * ncart_ket + iket] + ( bAB_x * BRA_h_s[4 * ncart_ket + iket] );

            // (h_3_1_1  p_0_1_0| = (i_3_2_1  s_0_0_0|_{t} + y_ab * (h_3_1_1  s_0_0_0|_{t}
            const double B_3_1_1_0_1_0 = BRA_i_s[7 * ncart_ket + iket] + ( bAB_y * BRA_h_s[4 * ncart_ket + iket] );

            // (h_3_1_1  p_0_0_1| = (i_3_1_2  s_0_0_0|_{t} + z_ab * (h_3_1_1  s_0_0_0|_{t}
            const double B_3_1_1_0_0_1 = BRA_i_s[8 * ncart_ket + iket] + ( bAB_z * BRA_h_s[4 * ncart_ket + iket] );

            // (h_3_0_2  p_1_0_0| = (i_4_0_2  s_0_0_0|_{t} + x_ab * (h_3_0_2  s_0_0_0|_{t}
            const double B_3_0_2_1_0_0 = BRA_i_s[5 * ncart_ket + iket] + ( bAB_x * BRA_h_s[5 * ncart_ket + iket] );

            // (h_3_0_2  p_0_0_1| = (i_3_0_3  s_0_0_0|_{t} + z_ab * (h_3_0_2  s_0_0_0|_{t}
            const double B_3_0_2_0_0_1 = BRA_i_s[9 * ncart_ket + iket] + ( bAB_z * BRA_h_s[5 * ncart_ket + iket] );

            // (h_2_3_0  p_1_0_0| = (i_3_3_0  s_0_0_0|_{t} + x_ab * (h_2_3_0  s_0_0_0|_{t}
            const double B_2_3_0_1_0_0 = BRA_i_s[6 * ncart_ket + iket] + ( bAB_x * BRA_h_s[6 * ncart_ket + iket] );

            // (h_2_3_0  p_0_1_0| = (i_2_4_0  s_0_0_0|_{t} + y_ab * (h_2_3_0  s_0_0_0|_{t}
            const double B_2_3_0_0_1_0 = BRA_i_s[10 * ncart_ket + iket] + ( bAB_y * BRA_h_s[6 * ncart_ket + iket] );

            // (h_2_2_1  p_1_0_0| = (i_3_2_1  s_0_0_0|_{t} + x_ab * (h_2_2_1  s_0_0_0|_{t}
            const double B_2_2_1_1_0_0 = BRA_i_s[7 * ncart_ket + iket] + ( bAB_x * BRA_h_s[7 * ncart_ket + iket] );

            // (h_2_2_1  p_0_1_0| = (i_2_3_1  s_0_0_0|_{t} + y_ab * (h_2_2_1  s_0_0_0|_{t}
            const double B_2_2_1_0_1_0 = BRA_i_s[11 * ncart_ket + iket] + ( bAB_y * BRA_h_s[7 * ncart_ket + iket] );

            // (h_2_2_1  p_0_0_1| = (i_2_2_2  s_0_0_0|_{t} + z_ab * (h_2_2_1  s_0_0_0|_{t}
            const double B_2_2_1_0_0_1 = BRA_i_s[12 * ncart_ket + iket] + ( bAB_z * BRA_h_s[7 * ncart_ket + iket] );

            // (h_2_1_2  p_1_0_0| = (i_3_1_2  s_0_0_0|_{t} + x_ab * (h_2_1_2  s_0_0_0|_{t}
            const double B_2_1_2_1_0_0 = BRA_i_s[8 * ncart_ket + iket] + ( bAB_x * BRA_h_s[8 * ncart_ket + iket] );

            // (h_2_1_2  p_0_1_0| = (i_2_2_2  s_0_0_0|_{t} + y_ab * (h_2_1_2  s_0_0_0|_{t}
            const double B_2_1_2_0_1_0 = BRA_i_s[12 * ncart_ket + iket] + ( bAB_y * BRA_h_s[8 * ncart_ket + iket] );

            // (h_2_1_2  p_0_0_1| = (i_2_1_3  s_0_0_0|_{t} + z_ab * (h_2_1_2  s_0_0_0|_{t}
            const double B_2_1_2_0_0_1 = BRA_i_s[13 * ncart_ket + iket] + ( bAB_z * BRA_h_s[8 * ncart_ket + iket] );

            // (h_2_0_3  p_1_0_0| = (i_3_0_3  s_0_0_0|_{t} + x_ab * (h_2_0_3  s_0_0_0|_{t}
            const double B_2_0_3_1_0_0 = BRA_i_s[9 * ncart_ket + iket] + ( bAB_x * BRA_h_s[9 * ncart_ket + iket] );

            // (h_2_0_3  p_0_0_1| = (i_2_0_4  s_0_0_0|_{t} + z_ab * (h_2_0_3  s_0_0_0|_{t}
            const double B_2_0_3_0_0_1 = BRA_i_s[14 * ncart_ket + iket] + ( bAB_z * BRA_h_s[9 * ncart_ket + iket] );

            // (h_1_4_0  p_1_0_0| = (i_2_4_0  s_0_0_0|_{t} + x_ab * (h_1_4_0  s_0_0_0|_{t}
            const double B_1_4_0_1_0_0 = BRA_i_s[10 * ncart_ket + iket] + ( bAB_x * BRA_h_s[10 * ncart_ket + iket] );

            // (h_1_4_0  p_0_1_0| = (i_1_5_0  s_0_0_0|_{t} + y_ab * (h_1_4_0  s_0_0_0|_{t}
            const double B_1_4_0_0_1_0 = BRA_i_s[15 * ncart_ket + iket] + ( bAB_y * BRA_h_s[10 * ncart_ket + iket] );

            // (h_1_3_1  p_1_0_0| = (i_2_3_1  s_0_0_0|_{t} + x_ab * (h_1_3_1  s_0_0_0|_{t}
            const double B_1_3_1_1_0_0 = BRA_i_s[11 * ncart_ket + iket] + ( bAB_x * BRA_h_s[11 * ncart_ket + iket] );

            // (h_1_3_1  p_0_1_0| = (i_1_4_1  s_0_0_0|_{t} + y_ab * (h_1_3_1  s_0_0_0|_{t}
            const double B_1_3_1_0_1_0 = BRA_i_s[16 * ncart_ket + iket] + ( bAB_y * BRA_h_s[11 * ncart_ket + iket] );

            // (h_1_3_1  p_0_0_1| = (i_1_3_2  s_0_0_0|_{t} + z_ab * (h_1_3_1  s_0_0_0|_{t}
            const double B_1_3_1_0_0_1 = BRA_i_s[17 * ncart_ket + iket] + ( bAB_z * BRA_h_s[11 * ncart_ket + iket] );

            // (h_1_2_2  p_1_0_0| = (i_2_2_2  s_0_0_0|_{t} + x_ab * (h_1_2_2  s_0_0_0|_{t}
            const double B_1_2_2_1_0_0 = BRA_i_s[12 * ncart_ket + iket] + ( bAB_x * BRA_h_s[12 * ncart_ket + iket] );

            // (h_1_2_2  p_0_1_0| = (i_1_3_2  s_0_0_0|_{t} + y_ab * (h_1_2_2  s_0_0_0|_{t}
            const double B_1_2_2_0_1_0 = BRA_i_s[17 * ncart_ket + iket] + ( bAB_y * BRA_h_s[12 * ncart_ket + iket] );

            // (h_1_2_2  p_0_0_1| = (i_1_2_3  s_0_0_0|_{t} + z_ab * (h_1_2_2  s_0_0_0|_{t}
            const double B_1_2_2_0_0_1 = BRA_i_s[18 * ncart_ket + iket] + ( bAB_z * BRA_h_s[12 * ncart_ket + iket] );

            // (h_1_1_3  p_1_0_0| = (i_2_1_3  s_0_0_0|_{t} + x_ab * (h_1_1_3  s_0_0_0|_{t}
            const double B_1_1_3_1_0_0 = BRA_i_s[13 * ncart_ket + iket] + ( bAB_x * BRA_h_s[13 * ncart_ket + iket] );

            // (h_1_1_3  p_0_1_0| = (i_1_2_3  s_0_0_0|_{t} + y_ab * (h_1_1_3  s_0_0_0|_{t}
            const double B_1_1_3_0_1_0 = BRA_i_s[18 * ncart_ket + iket] + ( bAB_y * BRA_h_s[13 * ncart_ket + iket] );

            // (h_1_1_3  p_0_0_1| = (i_1_1_4  s_0_0_0|_{t} + z_ab * (h_1_1_3  s_0_0_0|_{t}
            const double B_1_1_3_0_0_1 = BRA_i_s[19 * ncart_ket + iket] + ( bAB_z * BRA_h_s[13 * ncart_ket + iket] );

            // (h_1_0_4  p_1_0_0| = (i_2_0_4  s_0_0_0|_{t} + x_ab * (h_1_0_4  s_0_0_0|_{t}
            const double B_1_0_4_1_0_0 = BRA_i_s[14 * ncart_ket + iket] + ( bAB_x * BRA_h_s[14 * ncart_ket + iket] );

            // (h_1_0_4  p_0_0_1| = (i_1_0_5  s_0_0_0|_{t} + z_ab * (h_1_0_4  s_0_0_0|_{t}
            const double B_1_0_4_0_0_1 = BRA_i_s[20 * ncart_ket + iket] + ( bAB_z * BRA_h_s[14 * ncart_ket + iket] );

            // (h_0_5_0  p_0_1_0| = (i_0_6_0  s_0_0_0|_{t} + y_ab * (h_0_5_0  s_0_0_0|_{t}
            const double B_0_5_0_0_1_0 = BRA_i_s[21 * ncart_ket + iket] + ( bAB_y * BRA_h_s[15 * ncart_ket + iket] );

            // (h_0_4_1  p_1_0_0| = (i_1_4_1  s_0_0_0|_{t} + x_ab * (h_0_4_1  s_0_0_0|_{t}
            const double B_0_4_1_1_0_0 = BRA_i_s[16 * ncart_ket + iket] + ( bAB_x * BRA_h_s[16 * ncart_ket + iket] );

            // (h_0_4_1  p_0_1_0| = (i_0_5_1  s_0_0_0|_{t} + y_ab * (h_0_4_1  s_0_0_0|_{t}
            const double B_0_4_1_0_1_0 = BRA_i_s[22 * ncart_ket + iket] + ( bAB_y * BRA_h_s[16 * ncart_ket + iket] );

            // (h_0_4_1  p_0_0_1| = (i_0_4_2  s_0_0_0|_{t} + z_ab * (h_0_4_1  s_0_0_0|_{t}
            const double B_0_4_1_0_0_1 = BRA_i_s[23 * ncart_ket + iket] + ( bAB_z * BRA_h_s[16 * ncart_ket + iket] );

            // (h_0_3_2  p_1_0_0| = (i_1_3_2  s_0_0_0|_{t} + x_ab * (h_0_3_2  s_0_0_0|_{t}
            const double B_0_3_2_1_0_0 = BRA_i_s[17 * ncart_ket + iket] + ( bAB_x * BRA_h_s[17 * ncart_ket + iket] );

            // (h_0_3_2  p_0_1_0| = (i_0_4_2  s_0_0_0|_{t} + y_ab * (h_0_3_2  s_0_0_0|_{t}
            const double B_0_3_2_0_1_0 = BRA_i_s[23 * ncart_ket + iket] + ( bAB_y * BRA_h_s[17 * ncart_ket + iket] );

            // (h_0_3_2  p_0_0_1| = (i_0_3_3  s_0_0_0|_{t} + z_ab * (h_0_3_2  s_0_0_0|_{t}
            const double B_0_3_2_0_0_1 = BRA_i_s[24 * ncart_ket + iket] + ( bAB_z * BRA_h_s[17 * ncart_ket + iket] );

            // (h_0_2_3  p_1_0_0| = (i_1_2_3  s_0_0_0|_{t} + x_ab * (h_0_2_3  s_0_0_0|_{t}
            const double B_0_2_3_1_0_0 = BRA_i_s[18 * ncart_ket + iket] + ( bAB_x * BRA_h_s[18 * ncart_ket + iket] );

            // (h_0_2_3  p_0_1_0| = (i_0_3_3  s_0_0_0|_{t} + y_ab * (h_0_2_3  s_0_0_0|_{t}
            const double B_0_2_3_0_1_0 = BRA_i_s[24 * ncart_ket + iket] + ( bAB_y * BRA_h_s[18 * ncart_ket + iket] );

            // (h_0_2_3  p_0_0_1| = (i_0_2_4  s_0_0_0|_{t} + z_ab * (h_0_2_3  s_0_0_0|_{t}
            const double B_0_2_3_0_0_1 = BRA_i_s[25 * ncart_ket + iket] + ( bAB_z * BRA_h_s[18 * ncart_ket + iket] );

            // (h_0_1_4  p_1_0_0| = (i_1_1_4  s_0_0_0|_{t} + x_ab * (h_0_1_4  s_0_0_0|_{t}
            const double B_0_1_4_1_0_0 = BRA_i_s[19 * ncart_ket + iket] + ( bAB_x * BRA_h_s[19 * ncart_ket + iket] );

            // (h_0_1_4  p_0_1_0| = (i_0_2_4  s_0_0_0|_{t} + y_ab * (h_0_1_4  s_0_0_0|_{t}
            const double B_0_1_4_0_1_0 = BRA_i_s[25 * ncart_ket + iket] + ( bAB_y * BRA_h_s[19 * ncart_ket + iket] );

            // (h_0_1_4  p_0_0_1| = (i_0_1_5  s_0_0_0|_{t} + z_ab * (h_0_1_4  s_0_0_0|_{t}
            const double B_0_1_4_0_0_1 = BRA_i_s[26 * ncart_ket + iket] + ( bAB_z * BRA_h_s[19 * ncart_ket + iket] );

            // (h_0_0_5  p_0_0_1| = (i_0_0_6  s_0_0_0|_{t} + z_ab * (h_0_0_5  s_0_0_0|_{t}
            const double B_0_0_5_0_0_1 = BRA_i_s[27 * ncart_ket + iket] + ( bAB_z * BRA_h_s[20 * ncart_ket + iket] );

            // (f_3_0_0  d_2_0_0| = (g_4_0_0  p_1_0_0| + x_ab * (f_3_0_0  p_1_0_0|
            const double B_3_0_0_2_0_0 = B_4_0_0_1_0_0 + ( bAB_x * B_3_0_0_1_0_0 );

            // (f_3_0_0  d_1_1_0| = (g_3_1_0  p_1_0_0| + y_ab * (f_3_0_0  p_1_0_0|
            const double B_3_0_0_1_1_0 = B_3_1_0_1_0_0 + ( bAB_y * B_3_0_0_1_0_0 );

            // (f_3_0_0  d_0_2_0| = (g_3_1_0  p_0_1_0| + y_ab * (f_3_0_0  p_0_1_0|
            const double B_3_0_0_0_2_0 = B_3_1_0_0_1_0 + ( bAB_y * B_3_0_0_0_1_0 );

            // (f_3_0_0  d_0_0_2| = (g_3_0_1  p_0_0_1| + z_ab * (f_3_0_0  p_0_0_1|
            const double B_3_0_0_0_0_2 = B_3_0_1_0_0_1 + ( bAB_z * B_3_0_0_0_0_1 );

            // (f_2_1_0  d_2_0_0| = (g_3_1_0  p_1_0_0| + x_ab * (f_2_1_0  p_1_0_0|
            const double B_2_1_0_2_0_0 = B_3_1_0_1_0_0 + ( bAB_x * B_2_1_0_1_0_0 );

            // (f_2_1_0  d_1_1_0| = (g_2_2_0  p_1_0_0| + y_ab * (f_2_1_0  p_1_0_0|
            const double B_2_1_0_1_1_0 = B_2_2_0_1_0_0 + ( bAB_y * B_2_1_0_1_0_0 );

            // (f_2_1_0  d_0_2_0| = (g_2_2_0  p_0_1_0| + y_ab * (f_2_1_0  p_0_1_0|
            const double B_2_1_0_0_2_0 = B_2_2_0_0_1_0 + ( bAB_y * B_2_1_0_0_1_0 );

            // (f_2_1_0  d_0_0_2| = (g_2_1_1  p_0_0_1| + z_ab * (f_2_1_0  p_0_0_1|
            const double B_2_1_0_0_0_2 = B_2_1_1_0_0_1 + ( bAB_z * B_2_1_0_0_0_1 );

            // (f_2_0_1  d_2_0_0| = (g_3_0_1  p_1_0_0| + x_ab * (f_2_0_1  p_1_0_0|
            const double B_2_0_1_2_0_0 = B_3_0_1_1_0_0 + ( bAB_x * B_2_0_1_1_0_0 );

            // (f_2_0_1  d_1_1_0| = (g_2_1_1  p_1_0_0| + y_ab * (f_2_0_1  p_1_0_0|
            const double B_2_0_1_1_1_0 = B_2_1_1_1_0_0 + ( bAB_y * B_2_0_1_1_0_0 );

            // (f_2_0_1  d_0_2_0| = (g_2_1_1  p_0_1_0| + y_ab * (f_2_0_1  p_0_1_0|
            const double B_2_0_1_0_2_0 = B_2_1_1_0_1_0 + ( bAB_y * B_2_0_1_0_1_0 );

            // (f_2_0_1  d_0_0_2| = (g_2_0_2  p_0_0_1| + z_ab * (f_2_0_1  p_0_0_1|
            const double B_2_0_1_0_0_2 = B_2_0_2_0_0_1 + ( bAB_z * B_2_0_1_0_0_1 );

            // (f_1_2_0  d_2_0_0| = (g_2_2_0  p_1_0_0| + x_ab * (f_1_2_0  p_1_0_0|
            const double B_1_2_0_2_0_0 = B_2_2_0_1_0_0 + ( bAB_x * B_1_2_0_1_0_0 );

            // (f_1_2_0  d_1_1_0| = (g_1_3_0  p_1_0_0| + y_ab * (f_1_2_0  p_1_0_0|
            const double B_1_2_0_1_1_0 = B_1_3_0_1_0_0 + ( bAB_y * B_1_2_0_1_0_0 );

            // (f_1_2_0  d_0_2_0| = (g_1_3_0  p_0_1_0| + y_ab * (f_1_2_0  p_0_1_0|
            const double B_1_2_0_0_2_0 = B_1_3_0_0_1_0 + ( bAB_y * B_1_2_0_0_1_0 );

            // (f_1_2_0  d_0_0_2| = (g_1_2_1  p_0_0_1| + z_ab * (f_1_2_0  p_0_0_1|
            const double B_1_2_0_0_0_2 = B_1_2_1_0_0_1 + ( bAB_z * B_1_2_0_0_0_1 );

            // (f_1_1_1  d_2_0_0| = (g_2_1_1  p_1_0_0| + x_ab * (f_1_1_1  p_1_0_0|
            const double B_1_1_1_2_0_0 = B_2_1_1_1_0_0 + ( bAB_x * B_1_1_1_1_0_0 );

            // (f_1_1_1  d_1_1_0| = (g_1_2_1  p_1_0_0| + y_ab * (f_1_1_1  p_1_0_0|
            const double B_1_1_1_1_1_0 = B_1_2_1_1_0_0 + ( bAB_y * B_1_1_1_1_0_0 );

            // (f_1_1_1  d_0_2_0| = (g_1_2_1  p_0_1_0| + y_ab * (f_1_1_1  p_0_1_0|
            const double B_1_1_1_0_2_0 = B_1_2_1_0_1_0 + ( bAB_y * B_1_1_1_0_1_0 );

            // (f_1_1_1  d_0_0_2| = (g_1_1_2  p_0_0_1| + z_ab * (f_1_1_1  p_0_0_1|
            const double B_1_1_1_0_0_2 = B_1_1_2_0_0_1 + ( bAB_z * B_1_1_1_0_0_1 );

            // (f_1_0_2  d_2_0_0| = (g_2_0_2  p_1_0_0| + x_ab * (f_1_0_2  p_1_0_0|
            const double B_1_0_2_2_0_0 = B_2_0_2_1_0_0 + ( bAB_x * B_1_0_2_1_0_0 );

            // (f_1_0_2  d_1_1_0| = (g_1_1_2  p_1_0_0| + y_ab * (f_1_0_2  p_1_0_0|
            const double B_1_0_2_1_1_0 = B_1_1_2_1_0_0 + ( bAB_y * B_1_0_2_1_0_0 );

            // (f_1_0_2  d_0_2_0| = (g_1_1_2  p_0_1_0| + y_ab * (f_1_0_2  p_0_1_0|
            const double B_1_0_2_0_2_0 = B_1_1_2_0_1_0 + ( bAB_y * B_1_0_2_0_1_0 );

            // (f_1_0_2  d_0_0_2| = (g_1_0_3  p_0_0_1| + z_ab * (f_1_0_2  p_0_0_1|
            const double B_1_0_2_0_0_2 = B_1_0_3_0_0_1 + ( bAB_z * B_1_0_2_0_0_1 );

            // (f_0_3_0  d_2_0_0| = (g_1_3_0  p_1_0_0| + x_ab * (f_0_3_0  p_1_0_0|
            const double B_0_3_0_2_0_0 = B_1_3_0_1_0_0 + ( bAB_x * B_0_3_0_1_0_0 );

            // (f_0_3_0  d_1_1_0| = (g_0_4_0  p_1_0_0| + y_ab * (f_0_3_0  p_1_0_0|
            const double B_0_3_0_1_1_0 = B_0_4_0_1_0_0 + ( bAB_y * B_0_3_0_1_0_0 );

            // (f_0_3_0  d_0_2_0| = (g_0_4_0  p_0_1_0| + y_ab * (f_0_3_0  p_0_1_0|
            const double B_0_3_0_0_2_0 = B_0_4_0_0_1_0 + ( bAB_y * B_0_3_0_0_1_0 );

            // (f_0_3_0  d_0_0_2| = (g_0_3_1  p_0_0_1| + z_ab * (f_0_3_0  p_0_0_1|
            const double B_0_3_0_0_0_2 = B_0_3_1_0_0_1 + ( bAB_z * B_0_3_0_0_0_1 );

            // (f_0_2_1  d_2_0_0| = (g_1_2_1  p_1_0_0| + x_ab * (f_0_2_1  p_1_0_0|
            const double B_0_2_1_2_0_0 = B_1_2_1_1_0_0 + ( bAB_x * B_0_2_1_1_0_0 );

            // (f_0_2_1  d_1_1_0| = (g_0_3_1  p_1_0_0| + y_ab * (f_0_2_1  p_1_0_0|
            const double B_0_2_1_1_1_0 = B_0_3_1_1_0_0 + ( bAB_y * B_0_2_1_1_0_0 );

            // (f_0_2_1  d_0_2_0| = (g_0_3_1  p_0_1_0| + y_ab * (f_0_2_1  p_0_1_0|
            const double B_0_2_1_0_2_0 = B_0_3_1_0_1_0 + ( bAB_y * B_0_2_1_0_1_0 );

            // (f_0_2_1  d_0_0_2| = (g_0_2_2  p_0_0_1| + z_ab * (f_0_2_1  p_0_0_1|
            const double B_0_2_1_0_0_2 = B_0_2_2_0_0_1 + ( bAB_z * B_0_2_1_0_0_1 );

            // (f_0_1_2  d_2_0_0| = (g_1_1_2  p_1_0_0| + x_ab * (f_0_1_2  p_1_0_0|
            const double B_0_1_2_2_0_0 = B_1_1_2_1_0_0 + ( bAB_x * B_0_1_2_1_0_0 );

            // (f_0_1_2  d_1_1_0| = (g_0_2_2  p_1_0_0| + y_ab * (f_0_1_2  p_1_0_0|
            const double B_0_1_2_1_1_0 = B_0_2_2_1_0_0 + ( bAB_y * B_0_1_2_1_0_0 );

            // (f_0_1_2  d_0_2_0| = (g_0_2_2  p_0_1_0| + y_ab * (f_0_1_2  p_0_1_0|
            const double B_0_1_2_0_2_0 = B_0_2_2_0_1_0 + ( bAB_y * B_0_1_2_0_1_0 );

            // (f_0_1_2  d_0_0_2| = (g_0_1_3  p_0_0_1| + z_ab * (f_0_1_2  p_0_0_1|
            const double B_0_1_2_0_0_2 = B_0_1_3_0_0_1 + ( bAB_z * B_0_1_2_0_0_1 );

            // (f_0_0_3  d_2_0_0| = (g_1_0_3  p_1_0_0| + x_ab * (f_0_0_3  p_1_0_0|
            const double B_0_0_3_2_0_0 = B_1_0_3_1_0_0 + ( bAB_x * B_0_0_3_1_0_0 );

            // (f_0_0_3  d_1_1_0| = (g_0_1_3  p_1_0_0| + y_ab * (f_0_0_3  p_1_0_0|
            const double B_0_0_3_1_1_0 = B_0_1_3_1_0_0 + ( bAB_y * B_0_0_3_1_0_0 );

            // (f_0_0_3  d_0_2_0| = (g_0_1_3  p_0_1_0| + y_ab * (f_0_0_3  p_0_1_0|
            const double B_0_0_3_0_2_0 = B_0_1_3_0_1_0 + ( bAB_y * B_0_0_3_0_1_0 );

            // (f_0_0_3  d_0_0_2| = (g_0_0_4  p_0_0_1| + z_ab * (f_0_0_3  p_0_0_1|
            const double B_0_0_3_0_0_2 = B_0_0_4_0_0_1 + ( bAB_z * B_0_0_3_0_0_1 );

            // (g_4_0_0  d_2_0_0| = (h_5_0_0  p_1_0_0| + x_ab * (g_4_0_0  p_1_0_0|
            const double B_4_0_0_2_0_0 = B_5_0_0_1_0_0 + ( bAB_x * B_4_0_0_1_0_0 );

            // (g_4_0_0  d_0_2_0| = (h_4_1_0  p_0_1_0| + y_ab * (g_4_0_0  p_0_1_0|
            const double B_4_0_0_0_2_0 = B_4_1_0_0_1_0 + ( bAB_y * B_4_0_0_0_1_0 );

            // (g_4_0_0  d_0_0_2| = (h_4_0_1  p_0_0_1| + z_ab * (g_4_0_0  p_0_0_1|
            const double B_4_0_0_0_0_2 = B_4_0_1_0_0_1 + ( bAB_z * B_4_0_0_0_0_1 );

            // (g_3_1_0  d_2_0_0| = (h_4_1_0  p_1_0_0| + x_ab * (g_3_1_0  p_1_0_0|
            const double B_3_1_0_2_0_0 = B_4_1_0_1_0_0 + ( bAB_x * B_3_1_0_1_0_0 );

            // (g_3_1_0  d_0_2_0| = (h_3_2_0  p_0_1_0| + y_ab * (g_3_1_0  p_0_1_0|
            const double B_3_1_0_0_2_0 = B_3_2_0_0_1_0 + ( bAB_y * B_3_1_0_0_1_0 );

            // (g_3_1_0  d_0_0_2| = (h_3_1_1  p_0_0_1| + z_ab * (g_3_1_0  p_0_0_1|
            const double B_3_1_0_0_0_2 = B_3_1_1_0_0_1 + ( bAB_z * B_3_1_0_0_0_1 );

            // (g_3_0_1  d_2_0_0| = (h_4_0_1  p_1_0_0| + x_ab * (g_3_0_1  p_1_0_0|
            const double B_3_0_1_2_0_0 = B_4_0_1_1_0_0 + ( bAB_x * B_3_0_1_1_0_0 );

            // (g_3_0_1  d_1_1_0| = (h_3_1_1  p_1_0_0| + y_ab * (g_3_0_1  p_1_0_0|
            const double B_3_0_1_1_1_0 = B_3_1_1_1_0_0 + ( bAB_y * B_3_0_1_1_0_0 );

            // (g_3_0_1  d_0_2_0| = (h_3_1_1  p_0_1_0| + y_ab * (g_3_0_1  p_0_1_0|
            const double B_3_0_1_0_2_0 = B_3_1_1_0_1_0 + ( bAB_y * B_3_0_1_0_1_0 );

            // (g_3_0_1  d_0_0_2| = (h_3_0_2  p_0_0_1| + z_ab * (g_3_0_1  p_0_0_1|
            const double B_3_0_1_0_0_2 = B_3_0_2_0_0_1 + ( bAB_z * B_3_0_1_0_0_1 );

            // (g_2_2_0  d_2_0_0| = (h_3_2_0  p_1_0_0| + x_ab * (g_2_2_0  p_1_0_0|
            const double B_2_2_0_2_0_0 = B_3_2_0_1_0_0 + ( bAB_x * B_2_2_0_1_0_0 );

            // (g_2_2_0  d_0_2_0| = (h_2_3_0  p_0_1_0| + y_ab * (g_2_2_0  p_0_1_0|
            const double B_2_2_0_0_2_0 = B_2_3_0_0_1_0 + ( bAB_y * B_2_2_0_0_1_0 );

            // (g_2_2_0  d_0_0_2| = (h_2_2_1  p_0_0_1| + z_ab * (g_2_2_0  p_0_0_1|
            const double B_2_2_0_0_0_2 = B_2_2_1_0_0_1 + ( bAB_z * B_2_2_0_0_0_1 );

            // (g_2_1_1  d_2_0_0| = (h_3_1_1  p_1_0_0| + x_ab * (g_2_1_1  p_1_0_0|
            const double B_2_1_1_2_0_0 = B_3_1_1_1_0_0 + ( bAB_x * B_2_1_1_1_0_0 );

            // (g_2_1_1  d_1_1_0| = (h_2_2_1  p_1_0_0| + y_ab * (g_2_1_1  p_1_0_0|
            const double B_2_1_1_1_1_0 = B_2_2_1_1_0_0 + ( bAB_y * B_2_1_1_1_0_0 );

            // (g_2_1_1  d_0_2_0| = (h_2_2_1  p_0_1_0| + y_ab * (g_2_1_1  p_0_1_0|
            const double B_2_1_1_0_2_0 = B_2_2_1_0_1_0 + ( bAB_y * B_2_1_1_0_1_0 );

            // (g_2_1_1  d_0_0_2| = (h_2_1_2  p_0_0_1| + z_ab * (g_2_1_1  p_0_0_1|
            const double B_2_1_1_0_0_2 = B_2_1_2_0_0_1 + ( bAB_z * B_2_1_1_0_0_1 );

            // (g_2_0_2  d_2_0_0| = (h_3_0_2  p_1_0_0| + x_ab * (g_2_0_2  p_1_0_0|
            const double B_2_0_2_2_0_0 = B_3_0_2_1_0_0 + ( bAB_x * B_2_0_2_1_0_0 );

            // (g_2_0_2  d_1_1_0| = (h_2_1_2  p_1_0_0| + y_ab * (g_2_0_2  p_1_0_0|
            const double B_2_0_2_1_1_0 = B_2_1_2_1_0_0 + ( bAB_y * B_2_0_2_1_0_0 );

            // (g_2_0_2  d_0_2_0| = (h_2_1_2  p_0_1_0| + y_ab * (g_2_0_2  p_0_1_0|
            const double B_2_0_2_0_2_0 = B_2_1_2_0_1_0 + ( bAB_y * B_2_0_2_0_1_0 );

            // (g_2_0_2  d_0_0_2| = (h_2_0_3  p_0_0_1| + z_ab * (g_2_0_2  p_0_0_1|
            const double B_2_0_2_0_0_2 = B_2_0_3_0_0_1 + ( bAB_z * B_2_0_2_0_0_1 );

            // (g_1_3_0  d_2_0_0| = (h_2_3_0  p_1_0_0| + x_ab * (g_1_3_0  p_1_0_0|
            const double B_1_3_0_2_0_0 = B_2_3_0_1_0_0 + ( bAB_x * B_1_3_0_1_0_0 );

            // (g_1_3_0  d_0_2_0| = (h_1_4_0  p_0_1_0| + y_ab * (g_1_3_0  p_0_1_0|
            const double B_1_3_0_0_2_0 = B_1_4_0_0_1_0 + ( bAB_y * B_1_3_0_0_1_0 );

            // (g_1_3_0  d_0_0_2| = (h_1_3_1  p_0_0_1| + z_ab * (g_1_3_0  p_0_0_1|
            const double B_1_3_0_0_0_2 = B_1_3_1_0_0_1 + ( bAB_z * B_1_3_0_0_0_1 );

            // (g_1_2_1  d_2_0_0| = (h_2_2_1  p_1_0_0| + x_ab * (g_1_2_1  p_1_0_0|
            const double B_1_2_1_2_0_0 = B_2_2_1_1_0_0 + ( bAB_x * B_1_2_1_1_0_0 );

            // (g_1_2_1  d_1_1_0| = (h_1_3_1  p_1_0_0| + y_ab * (g_1_2_1  p_1_0_0|
            const double B_1_2_1_1_1_0 = B_1_3_1_1_0_0 + ( bAB_y * B_1_2_1_1_0_0 );

            // (g_1_2_1  d_0_2_0| = (h_1_3_1  p_0_1_0| + y_ab * (g_1_2_1  p_0_1_0|
            const double B_1_2_1_0_2_0 = B_1_3_1_0_1_0 + ( bAB_y * B_1_2_1_0_1_0 );

            // (g_1_2_1  d_0_0_2| = (h_1_2_2  p_0_0_1| + z_ab * (g_1_2_1  p_0_0_1|
            const double B_1_2_1_0_0_2 = B_1_2_2_0_0_1 + ( bAB_z * B_1_2_1_0_0_1 );

            // (g_1_1_2  d_2_0_0| = (h_2_1_2  p_1_0_0| + x_ab * (g_1_1_2  p_1_0_0|
            const double B_1_1_2_2_0_0 = B_2_1_2_1_0_0 + ( bAB_x * B_1_1_2_1_0_0 );

            // (g_1_1_2  d_1_1_0| = (h_1_2_2  p_1_0_0| + y_ab * (g_1_1_2  p_1_0_0|
            const double B_1_1_2_1_1_0 = B_1_2_2_1_0_0 + ( bAB_y * B_1_1_2_1_0_0 );

            // (g_1_1_2  d_0_2_0| = (h_1_2_2  p_0_1_0| + y_ab * (g_1_1_2  p_0_1_0|
            const double B_1_1_2_0_2_0 = B_1_2_2_0_1_0 + ( bAB_y * B_1_1_2_0_1_0 );

            // (g_1_1_2  d_0_0_2| = (h_1_1_3  p_0_0_1| + z_ab * (g_1_1_2  p_0_0_1|
            const double B_1_1_2_0_0_2 = B_1_1_3_0_0_1 + ( bAB_z * B_1_1_2_0_0_1 );

            // (g_1_0_3  d_2_0_0| = (h_2_0_3  p_1_0_0| + x_ab * (g_1_0_3  p_1_0_0|
            const double B_1_0_3_2_0_0 = B_2_0_3_1_0_0 + ( bAB_x * B_1_0_3_1_0_0 );

            // (g_1_0_3  d_1_1_0| = (h_1_1_3  p_1_0_0| + y_ab * (g_1_0_3  p_1_0_0|
            const double B_1_0_3_1_1_0 = B_1_1_3_1_0_0 + ( bAB_y * B_1_0_3_1_0_0 );

            // (g_1_0_3  d_0_2_0| = (h_1_1_3  p_0_1_0| + y_ab * (g_1_0_3  p_0_1_0|
            const double B_1_0_3_0_2_0 = B_1_1_3_0_1_0 + ( bAB_y * B_1_0_3_0_1_0 );

            // (g_1_0_3  d_0_0_2| = (h_1_0_4  p_0_0_1| + z_ab * (g_1_0_3  p_0_0_1|
            const double B_1_0_3_0_0_2 = B_1_0_4_0_0_1 + ( bAB_z * B_1_0_3_0_0_1 );

            // (g_0_4_0  d_2_0_0| = (h_1_4_0  p_1_0_0| + x_ab * (g_0_4_0  p_1_0_0|
            const double B_0_4_0_2_0_0 = B_1_4_0_1_0_0 + ( bAB_x * B_0_4_0_1_0_0 );

            // (g_0_4_0  d_0_2_0| = (h_0_5_0  p_0_1_0| + y_ab * (g_0_4_0  p_0_1_0|
            const double B_0_4_0_0_2_0 = B_0_5_0_0_1_0 + ( bAB_y * B_0_4_0_0_1_0 );

            // (g_0_4_0  d_0_0_2| = (h_0_4_1  p_0_0_1| + z_ab * (g_0_4_0  p_0_0_1|
            const double B_0_4_0_0_0_2 = B_0_4_1_0_0_1 + ( bAB_z * B_0_4_0_0_0_1 );

            // (g_0_3_1  d_2_0_0| = (h_1_3_1  p_1_0_0| + x_ab * (g_0_3_1  p_1_0_0|
            const double B_0_3_1_2_0_0 = B_1_3_1_1_0_0 + ( bAB_x * B_0_3_1_1_0_0 );

            // (g_0_3_1  d_1_1_0| = (h_0_4_1  p_1_0_0| + y_ab * (g_0_3_1  p_1_0_0|
            const double B_0_3_1_1_1_0 = B_0_4_1_1_0_0 + ( bAB_y * B_0_3_1_1_0_0 );

            // (g_0_3_1  d_0_2_0| = (h_0_4_1  p_0_1_0| + y_ab * (g_0_3_1  p_0_1_0|
            const double B_0_3_1_0_2_0 = B_0_4_1_0_1_0 + ( bAB_y * B_0_3_1_0_1_0 );

            // (g_0_3_1  d_0_0_2| = (h_0_3_2  p_0_0_1| + z_ab * (g_0_3_1  p_0_0_1|
            const double B_0_3_1_0_0_2 = B_0_3_2_0_0_1 + ( bAB_z * B_0_3_1_0_0_1 );

            // (g_0_2_2  d_2_0_0| = (h_1_2_2  p_1_0_0| + x_ab * (g_0_2_2  p_1_0_0|
            const double B_0_2_2_2_0_0 = B_1_2_2_1_0_0 + ( bAB_x * B_0_2_2_1_0_0 );

            // (g_0_2_2  d_1_1_0| = (h_0_3_2  p_1_0_0| + y_ab * (g_0_2_2  p_1_0_0|
            const double B_0_2_2_1_1_0 = B_0_3_2_1_0_0 + ( bAB_y * B_0_2_2_1_0_0 );

            // (g_0_2_2  d_0_2_0| = (h_0_3_2  p_0_1_0| + y_ab * (g_0_2_2  p_0_1_0|
            const double B_0_2_2_0_2_0 = B_0_3_2_0_1_0 + ( bAB_y * B_0_2_2_0_1_0 );

            // (g_0_2_2  d_0_0_2| = (h_0_2_3  p_0_0_1| + z_ab * (g_0_2_2  p_0_0_1|
            const double B_0_2_2_0_0_2 = B_0_2_3_0_0_1 + ( bAB_z * B_0_2_2_0_0_1 );

            // (g_0_1_3  d_2_0_0| = (h_1_1_3  p_1_0_0| + x_ab * (g_0_1_3  p_1_0_0|
            const double B_0_1_3_2_0_0 = B_1_1_3_1_0_0 + ( bAB_x * B_0_1_3_1_0_0 );

            // (g_0_1_3  d_1_1_0| = (h_0_2_3  p_1_0_0| + y_ab * (g_0_1_3  p_1_0_0|
            const double B_0_1_3_1_1_0 = B_0_2_3_1_0_0 + ( bAB_y * B_0_1_3_1_0_0 );

            // (g_0_1_3  d_0_2_0| = (h_0_2_3  p_0_1_0| + y_ab * (g_0_1_3  p_0_1_0|
            const double B_0_1_3_0_2_0 = B_0_2_3_0_1_0 + ( bAB_y * B_0_1_3_0_1_0 );

            // (g_0_1_3  d_0_0_2| = (h_0_1_4  p_0_0_1| + z_ab * (g_0_1_3  p_0_0_1|
            const double B_0_1_3_0_0_2 = B_0_1_4_0_0_1 + ( bAB_z * B_0_1_3_0_0_1 );

            // (g_0_0_4  d_2_0_0| = (h_1_0_4  p_1_0_0| + x_ab * (g_0_0_4  p_1_0_0|
            const double B_0_0_4_2_0_0 = B_1_0_4_1_0_0 + ( bAB_x * B_0_0_4_1_0_0 );

            // (g_0_0_4  d_1_1_0| = (h_0_1_4  p_1_0_0| + y_ab * (g_0_0_4  p_1_0_0|
            const double B_0_0_4_1_1_0 = B_0_1_4_1_0_0 + ( bAB_y * B_0_0_4_1_0_0 );

            // (g_0_0_4  d_0_2_0| = (h_0_1_4  p_0_1_0| + y_ab * (g_0_0_4  p_0_1_0|
            const double B_0_0_4_0_2_0 = B_0_1_4_0_1_0 + ( bAB_y * B_0_0_4_0_1_0 );

            // (g_0_0_4  d_0_0_2| = (h_0_0_5  p_0_0_1| + z_ab * (g_0_0_4  p_0_0_1|
            const double B_0_0_4_0_0_2 = B_0_0_5_0_0_1 + ( bAB_z * B_0_0_4_0_0_1 );

            // (f_3_0_0  f_3_0_0|_{i} = (g_4_0_0  d_2_0_0| + x_ab * (f_3_0_0  d_2_0_0|
            BRA_f_f[0 * ncart_ket + iket] = B_4_0_0_2_0_0 + ( bAB_x * B_3_0_0_2_0_0 );

            // (f_3_0_0  f_2_1_0|_{i} = (g_3_1_0  d_2_0_0| + y_ab * (f_3_0_0  d_2_0_0|
            BRA_f_f[1 * ncart_ket + iket] = B_3_1_0_2_0_0 + ( bAB_y * B_3_0_0_2_0_0 );

            // (f_3_0_0  f_2_0_1|_{i} = (g_3_0_1  d_2_0_0| + z_ab * (f_3_0_0  d_2_0_0|
            BRA_f_f[2 * ncart_ket + iket] = B_3_0_1_2_0_0 + ( bAB_z * B_3_0_0_2_0_0 );

            // (f_3_0_0  f_1_2_0|_{i} = (g_4_0_0  d_0_2_0| + x_ab * (f_3_0_0  d_0_2_0|
            BRA_f_f[3 * ncart_ket + iket] = B_4_0_0_0_2_0 + ( bAB_x * B_3_0_0_0_2_0 );

            // (f_3_0_0  f_1_1_1|_{i} = (g_3_0_1  d_1_1_0| + z_ab * (f_3_0_0  d_1_1_0|
            BRA_f_f[4 * ncart_ket + iket] = B_3_0_1_1_1_0 + ( bAB_z * B_3_0_0_1_1_0 );

            // (f_3_0_0  f_1_0_2|_{i} = (g_4_0_0  d_0_0_2| + x_ab * (f_3_0_0  d_0_0_2|
            BRA_f_f[5 * ncart_ket + iket] = B_4_0_0_0_0_2 + ( bAB_x * B_3_0_0_0_0_2 );

            // (f_3_0_0  f_0_3_0|_{i} = (g_3_1_0  d_0_2_0| + y_ab * (f_3_0_0  d_0_2_0|
            BRA_f_f[6 * ncart_ket + iket] = B_3_1_0_0_2_0 + ( bAB_y * B_3_0_0_0_2_0 );

            // (f_3_0_0  f_0_2_1|_{i} = (g_3_0_1  d_0_2_0| + z_ab * (f_3_0_0  d_0_2_0|
            BRA_f_f[7 * ncart_ket + iket] = B_3_0_1_0_2_0 + ( bAB_z * B_3_0_0_0_2_0 );

            // (f_3_0_0  f_0_1_2|_{i} = (g_3_1_0  d_0_0_2| + y_ab * (f_3_0_0  d_0_0_2|
            BRA_f_f[8 * ncart_ket + iket] = B_3_1_0_0_0_2 + ( bAB_y * B_3_0_0_0_0_2 );

            // (f_3_0_0  f_0_0_3|_{i} = (g_3_0_1  d_0_0_2| + z_ab * (f_3_0_0  d_0_0_2|
            BRA_f_f[9 * ncart_ket + iket] = B_3_0_1_0_0_2 + ( bAB_z * B_3_0_0_0_0_2 );

            // (f_2_1_0  f_3_0_0|_{i} = (g_3_1_0  d_2_0_0| + x_ab * (f_2_1_0  d_2_0_0|
            BRA_f_f[10 * ncart_ket + iket] = B_3_1_0_2_0_0 + ( bAB_x * B_2_1_0_2_0_0 );

            // (f_2_1_0  f_2_1_0|_{i} = (g_2_2_0  d_2_0_0| + y_ab * (f_2_1_0  d_2_0_0|
            BRA_f_f[11 * ncart_ket + iket] = B_2_2_0_2_0_0 + ( bAB_y * B_2_1_0_2_0_0 );

            // (f_2_1_0  f_2_0_1|_{i} = (g_2_1_1  d_2_0_0| + z_ab * (f_2_1_0  d_2_0_0|
            BRA_f_f[12 * ncart_ket + iket] = B_2_1_1_2_0_0 + ( bAB_z * B_2_1_0_2_0_0 );

            // (f_2_1_0  f_1_2_0|_{i} = (g_3_1_0  d_0_2_0| + x_ab * (f_2_1_0  d_0_2_0|
            BRA_f_f[13 * ncart_ket + iket] = B_3_1_0_0_2_0 + ( bAB_x * B_2_1_0_0_2_0 );

            // (f_2_1_0  f_1_1_1|_{i} = (g_2_1_1  d_1_1_0| + z_ab * (f_2_1_0  d_1_1_0|
            BRA_f_f[14 * ncart_ket + iket] = B_2_1_1_1_1_0 + ( bAB_z * B_2_1_0_1_1_0 );

            // (f_2_1_0  f_1_0_2|_{i} = (g_3_1_0  d_0_0_2| + x_ab * (f_2_1_0  d_0_0_2|
            BRA_f_f[15 * ncart_ket + iket] = B_3_1_0_0_0_2 + ( bAB_x * B_2_1_0_0_0_2 );

            // (f_2_1_0  f_0_3_0|_{i} = (g_2_2_0  d_0_2_0| + y_ab * (f_2_1_0  d_0_2_0|
            BRA_f_f[16 * ncart_ket + iket] = B_2_2_0_0_2_0 + ( bAB_y * B_2_1_0_0_2_0 );

            // (f_2_1_0  f_0_2_1|_{i} = (g_2_1_1  d_0_2_0| + z_ab * (f_2_1_0  d_0_2_0|
            BRA_f_f[17 * ncart_ket + iket] = B_2_1_1_0_2_0 + ( bAB_z * B_2_1_0_0_2_0 );

            // (f_2_1_0  f_0_1_2|_{i} = (g_2_2_0  d_0_0_2| + y_ab * (f_2_1_0  d_0_0_2|
            BRA_f_f[18 * ncart_ket + iket] = B_2_2_0_0_0_2 + ( bAB_y * B_2_1_0_0_0_2 );

            // (f_2_1_0  f_0_0_3|_{i} = (g_2_1_1  d_0_0_2| + z_ab * (f_2_1_0  d_0_0_2|
            BRA_f_f[19 * ncart_ket + iket] = B_2_1_1_0_0_2 + ( bAB_z * B_2_1_0_0_0_2 );

            // (f_2_0_1  f_3_0_0|_{i} = (g_3_0_1  d_2_0_0| + x_ab * (f_2_0_1  d_2_0_0|
            BRA_f_f[20 * ncart_ket + iket] = B_3_0_1_2_0_0 + ( bAB_x * B_2_0_1_2_0_0 );

            // (f_2_0_1  f_2_1_0|_{i} = (g_2_1_1  d_2_0_0| + y_ab * (f_2_0_1  d_2_0_0|
            BRA_f_f[21 * ncart_ket + iket] = B_2_1_1_2_0_0 + ( bAB_y * B_2_0_1_2_0_0 );

            // (f_2_0_1  f_2_0_1|_{i} = (g_2_0_2  d_2_0_0| + z_ab * (f_2_0_1  d_2_0_0|
            BRA_f_f[22 * ncart_ket + iket] = B_2_0_2_2_0_0 + ( bAB_z * B_2_0_1_2_0_0 );

            // (f_2_0_1  f_1_2_0|_{i} = (g_3_0_1  d_0_2_0| + x_ab * (f_2_0_1  d_0_2_0|
            BRA_f_f[23 * ncart_ket + iket] = B_3_0_1_0_2_0 + ( bAB_x * B_2_0_1_0_2_0 );

            // (f_2_0_1  f_1_1_1|_{i} = (g_2_0_2  d_1_1_0| + z_ab * (f_2_0_1  d_1_1_0|
            BRA_f_f[24 * ncart_ket + iket] = B_2_0_2_1_1_0 + ( bAB_z * B_2_0_1_1_1_0 );

            // (f_2_0_1  f_1_0_2|_{i} = (g_3_0_1  d_0_0_2| + x_ab * (f_2_0_1  d_0_0_2|
            BRA_f_f[25 * ncart_ket + iket] = B_3_0_1_0_0_2 + ( bAB_x * B_2_0_1_0_0_2 );

            // (f_2_0_1  f_0_3_0|_{i} = (g_2_1_1  d_0_2_0| + y_ab * (f_2_0_1  d_0_2_0|
            BRA_f_f[26 * ncart_ket + iket] = B_2_1_1_0_2_0 + ( bAB_y * B_2_0_1_0_2_0 );

            // (f_2_0_1  f_0_2_1|_{i} = (g_2_0_2  d_0_2_0| + z_ab * (f_2_0_1  d_0_2_0|
            BRA_f_f[27 * ncart_ket + iket] = B_2_0_2_0_2_0 + ( bAB_z * B_2_0_1_0_2_0 );

            // (f_2_0_1  f_0_1_2|_{i} = (g_2_1_1  d_0_0_2| + y_ab * (f_2_0_1  d_0_0_2|
            BRA_f_f[28 * ncart_ket + iket] = B_2_1_1_0_0_2 + ( bAB_y * B_2_0_1_0_0_2 );

            // (f_2_0_1  f_0_0_3|_{i} = (g_2_0_2  d_0_0_2| + z_ab * (f_2_0_1  d_0_0_2|
            BRA_f_f[29 * ncart_ket + iket] = B_2_0_2_0_0_2 + ( bAB_z * B_2_0_1_0_0_2 );

            // (f_1_2_0  f_3_0_0|_{i} = (g_2_2_0  d_2_0_0| + x_ab * (f_1_2_0  d_2_0_0|
            BRA_f_f[30 * ncart_ket + iket] = B_2_2_0_2_0_0 + ( bAB_x * B_1_2_0_2_0_0 );

            // (f_1_2_0  f_2_1_0|_{i} = (g_1_3_0  d_2_0_0| + y_ab * (f_1_2_0  d_2_0_0|
            BRA_f_f[31 * ncart_ket + iket] = B_1_3_0_2_0_0 + ( bAB_y * B_1_2_0_2_0_0 );

            // (f_1_2_0  f_2_0_1|_{i} = (g_1_2_1  d_2_0_0| + z_ab * (f_1_2_0  d_2_0_0|
            BRA_f_f[32 * ncart_ket + iket] = B_1_2_1_2_0_0 + ( bAB_z * B_1_2_0_2_0_0 );

            // (f_1_2_0  f_1_2_0|_{i} = (g_2_2_0  d_0_2_0| + x_ab * (f_1_2_0  d_0_2_0|
            BRA_f_f[33 * ncart_ket + iket] = B_2_2_0_0_2_0 + ( bAB_x * B_1_2_0_0_2_0 );

            // (f_1_2_0  f_1_1_1|_{i} = (g_1_2_1  d_1_1_0| + z_ab * (f_1_2_0  d_1_1_0|
            BRA_f_f[34 * ncart_ket + iket] = B_1_2_1_1_1_0 + ( bAB_z * B_1_2_0_1_1_0 );

            // (f_1_2_0  f_1_0_2|_{i} = (g_2_2_0  d_0_0_2| + x_ab * (f_1_2_0  d_0_0_2|
            BRA_f_f[35 * ncart_ket + iket] = B_2_2_0_0_0_2 + ( bAB_x * B_1_2_0_0_0_2 );

            // (f_1_2_0  f_0_3_0|_{i} = (g_1_3_0  d_0_2_0| + y_ab * (f_1_2_0  d_0_2_0|
            BRA_f_f[36 * ncart_ket + iket] = B_1_3_0_0_2_0 + ( bAB_y * B_1_2_0_0_2_0 );

            // (f_1_2_0  f_0_2_1|_{i} = (g_1_2_1  d_0_2_0| + z_ab * (f_1_2_0  d_0_2_0|
            BRA_f_f[37 * ncart_ket + iket] = B_1_2_1_0_2_0 + ( bAB_z * B_1_2_0_0_2_0 );

            // (f_1_2_0  f_0_1_2|_{i} = (g_1_3_0  d_0_0_2| + y_ab * (f_1_2_0  d_0_0_2|
            BRA_f_f[38 * ncart_ket + iket] = B_1_3_0_0_0_2 + ( bAB_y * B_1_2_0_0_0_2 );

            // (f_1_2_0  f_0_0_3|_{i} = (g_1_2_1  d_0_0_2| + z_ab * (f_1_2_0  d_0_0_2|
            BRA_f_f[39 * ncart_ket + iket] = B_1_2_1_0_0_2 + ( bAB_z * B_1_2_0_0_0_2 );

            // (f_1_1_1  f_3_0_0|_{i} = (g_2_1_1  d_2_0_0| + x_ab * (f_1_1_1  d_2_0_0|
            BRA_f_f[40 * ncart_ket + iket] = B_2_1_1_2_0_0 + ( bAB_x * B_1_1_1_2_0_0 );

            // (f_1_1_1  f_2_1_0|_{i} = (g_1_2_1  d_2_0_0| + y_ab * (f_1_1_1  d_2_0_0|
            BRA_f_f[41 * ncart_ket + iket] = B_1_2_1_2_0_0 + ( bAB_y * B_1_1_1_2_0_0 );

            // (f_1_1_1  f_2_0_1|_{i} = (g_1_1_2  d_2_0_0| + z_ab * (f_1_1_1  d_2_0_0|
            BRA_f_f[42 * ncart_ket + iket] = B_1_1_2_2_0_0 + ( bAB_z * B_1_1_1_2_0_0 );

            // (f_1_1_1  f_1_2_0|_{i} = (g_2_1_1  d_0_2_0| + x_ab * (f_1_1_1  d_0_2_0|
            BRA_f_f[43 * ncart_ket + iket] = B_2_1_1_0_2_0 + ( bAB_x * B_1_1_1_0_2_0 );

            // (f_1_1_1  f_1_1_1|_{i} = (g_1_1_2  d_1_1_0| + z_ab * (f_1_1_1  d_1_1_0|
            BRA_f_f[44 * ncart_ket + iket] = B_1_1_2_1_1_0 + ( bAB_z * B_1_1_1_1_1_0 );

            // (f_1_1_1  f_1_0_2|_{i} = (g_2_1_1  d_0_0_2| + x_ab * (f_1_1_1  d_0_0_2|
            BRA_f_f[45 * ncart_ket + iket] = B_2_1_1_0_0_2 + ( bAB_x * B_1_1_1_0_0_2 );

            // (f_1_1_1  f_0_3_0|_{i} = (g_1_2_1  d_0_2_0| + y_ab * (f_1_1_1  d_0_2_0|
            BRA_f_f[46 * ncart_ket + iket] = B_1_2_1_0_2_0 + ( bAB_y * B_1_1_1_0_2_0 );

            // (f_1_1_1  f_0_2_1|_{i} = (g_1_1_2  d_0_2_0| + z_ab * (f_1_1_1  d_0_2_0|
            BRA_f_f[47 * ncart_ket + iket] = B_1_1_2_0_2_0 + ( bAB_z * B_1_1_1_0_2_0 );

            // (f_1_1_1  f_0_1_2|_{i} = (g_1_2_1  d_0_0_2| + y_ab * (f_1_1_1  d_0_0_2|
            BRA_f_f[48 * ncart_ket + iket] = B_1_2_1_0_0_2 + ( bAB_y * B_1_1_1_0_0_2 );

            // (f_1_1_1  f_0_0_3|_{i} = (g_1_1_2  d_0_0_2| + z_ab * (f_1_1_1  d_0_0_2|
            BRA_f_f[49 * ncart_ket + iket] = B_1_1_2_0_0_2 + ( bAB_z * B_1_1_1_0_0_2 );

            // (f_1_0_2  f_3_0_0|_{i} = (g_2_0_2  d_2_0_0| + x_ab * (f_1_0_2  d_2_0_0|
            BRA_f_f[50 * ncart_ket + iket] = B_2_0_2_2_0_0 + ( bAB_x * B_1_0_2_2_0_0 );

            // (f_1_0_2  f_2_1_0|_{i} = (g_1_1_2  d_2_0_0| + y_ab * (f_1_0_2  d_2_0_0|
            BRA_f_f[51 * ncart_ket + iket] = B_1_1_2_2_0_0 + ( bAB_y * B_1_0_2_2_0_0 );

            // (f_1_0_2  f_2_0_1|_{i} = (g_1_0_3  d_2_0_0| + z_ab * (f_1_0_2  d_2_0_0|
            BRA_f_f[52 * ncart_ket + iket] = B_1_0_3_2_0_0 + ( bAB_z * B_1_0_2_2_0_0 );

            // (f_1_0_2  f_1_2_0|_{i} = (g_2_0_2  d_0_2_0| + x_ab * (f_1_0_2  d_0_2_0|
            BRA_f_f[53 * ncart_ket + iket] = B_2_0_2_0_2_0 + ( bAB_x * B_1_0_2_0_2_0 );

            // (f_1_0_2  f_1_1_1|_{i} = (g_1_0_3  d_1_1_0| + z_ab * (f_1_0_2  d_1_1_0|
            BRA_f_f[54 * ncart_ket + iket] = B_1_0_3_1_1_0 + ( bAB_z * B_1_0_2_1_1_0 );

            // (f_1_0_2  f_1_0_2|_{i} = (g_2_0_2  d_0_0_2| + x_ab * (f_1_0_2  d_0_0_2|
            BRA_f_f[55 * ncart_ket + iket] = B_2_0_2_0_0_2 + ( bAB_x * B_1_0_2_0_0_2 );

            // (f_1_0_2  f_0_3_0|_{i} = (g_1_1_2  d_0_2_0| + y_ab * (f_1_0_2  d_0_2_0|
            BRA_f_f[56 * ncart_ket + iket] = B_1_1_2_0_2_0 + ( bAB_y * B_1_0_2_0_2_0 );

            // (f_1_0_2  f_0_2_1|_{i} = (g_1_0_3  d_0_2_0| + z_ab * (f_1_0_2  d_0_2_0|
            BRA_f_f[57 * ncart_ket + iket] = B_1_0_3_0_2_0 + ( bAB_z * B_1_0_2_0_2_0 );

            // (f_1_0_2  f_0_1_2|_{i} = (g_1_1_2  d_0_0_2| + y_ab * (f_1_0_2  d_0_0_2|
            BRA_f_f[58 * ncart_ket + iket] = B_1_1_2_0_0_2 + ( bAB_y * B_1_0_2_0_0_2 );

            // (f_1_0_2  f_0_0_3|_{i} = (g_1_0_3  d_0_0_2| + z_ab * (f_1_0_2  d_0_0_2|
            BRA_f_f[59 * ncart_ket + iket] = B_1_0_3_0_0_2 + ( bAB_z * B_1_0_2_0_0_2 );

            // (f_0_3_0  f_3_0_0|_{i} = (g_1_3_0  d_2_0_0| + x_ab * (f_0_3_0  d_2_0_0|
            BRA_f_f[60 * ncart_ket + iket] = B_1_3_0_2_0_0 + ( bAB_x * B_0_3_0_2_0_0 );

            // (f_0_3_0  f_2_1_0|_{i} = (g_0_4_0  d_2_0_0| + y_ab * (f_0_3_0  d_2_0_0|
            BRA_f_f[61 * ncart_ket + iket] = B_0_4_0_2_0_0 + ( bAB_y * B_0_3_0_2_0_0 );

            // (f_0_3_0  f_2_0_1|_{i} = (g_0_3_1  d_2_0_0| + z_ab * (f_0_3_0  d_2_0_0|
            BRA_f_f[62 * ncart_ket + iket] = B_0_3_1_2_0_0 + ( bAB_z * B_0_3_0_2_0_0 );

            // (f_0_3_0  f_1_2_0|_{i} = (g_1_3_0  d_0_2_0| + x_ab * (f_0_3_0  d_0_2_0|
            BRA_f_f[63 * ncart_ket + iket] = B_1_3_0_0_2_0 + ( bAB_x * B_0_3_0_0_2_0 );

            // (f_0_3_0  f_1_1_1|_{i} = (g_0_3_1  d_1_1_0| + z_ab * (f_0_3_0  d_1_1_0|
            BRA_f_f[64 * ncart_ket + iket] = B_0_3_1_1_1_0 + ( bAB_z * B_0_3_0_1_1_0 );

            // (f_0_3_0  f_1_0_2|_{i} = (g_1_3_0  d_0_0_2| + x_ab * (f_0_3_0  d_0_0_2|
            BRA_f_f[65 * ncart_ket + iket] = B_1_3_0_0_0_2 + ( bAB_x * B_0_3_0_0_0_2 );

            // (f_0_3_0  f_0_3_0|_{i} = (g_0_4_0  d_0_2_0| + y_ab * (f_0_3_0  d_0_2_0|
            BRA_f_f[66 * ncart_ket + iket] = B_0_4_0_0_2_0 + ( bAB_y * B_0_3_0_0_2_0 );

            // (f_0_3_0  f_0_2_1|_{i} = (g_0_3_1  d_0_2_0| + z_ab * (f_0_3_0  d_0_2_0|
            BRA_f_f[67 * ncart_ket + iket] = B_0_3_1_0_2_0 + ( bAB_z * B_0_3_0_0_2_0 );

            // (f_0_3_0  f_0_1_2|_{i} = (g_0_4_0  d_0_0_2| + y_ab * (f_0_3_0  d_0_0_2|
            BRA_f_f[68 * ncart_ket + iket] = B_0_4_0_0_0_2 + ( bAB_y * B_0_3_0_0_0_2 );

            // (f_0_3_0  f_0_0_3|_{i} = (g_0_3_1  d_0_0_2| + z_ab * (f_0_3_0  d_0_0_2|
            BRA_f_f[69 * ncart_ket + iket] = B_0_3_1_0_0_2 + ( bAB_z * B_0_3_0_0_0_2 );

            // (f_0_2_1  f_3_0_0|_{i} = (g_1_2_1  d_2_0_0| + x_ab * (f_0_2_1  d_2_0_0|
            BRA_f_f[70 * ncart_ket + iket] = B_1_2_1_2_0_0 + ( bAB_x * B_0_2_1_2_0_0 );

            // (f_0_2_1  f_2_1_0|_{i} = (g_0_3_1  d_2_0_0| + y_ab * (f_0_2_1  d_2_0_0|
            BRA_f_f[71 * ncart_ket + iket] = B_0_3_1_2_0_0 + ( bAB_y * B_0_2_1_2_0_0 );

            // (f_0_2_1  f_2_0_1|_{i} = (g_0_2_2  d_2_0_0| + z_ab * (f_0_2_1  d_2_0_0|
            BRA_f_f[72 * ncart_ket + iket] = B_0_2_2_2_0_0 + ( bAB_z * B_0_2_1_2_0_0 );

            // (f_0_2_1  f_1_2_0|_{i} = (g_1_2_1  d_0_2_0| + x_ab * (f_0_2_1  d_0_2_0|
            BRA_f_f[73 * ncart_ket + iket] = B_1_2_1_0_2_0 + ( bAB_x * B_0_2_1_0_2_0 );

            // (f_0_2_1  f_1_1_1|_{i} = (g_0_2_2  d_1_1_0| + z_ab * (f_0_2_1  d_1_1_0|
            BRA_f_f[74 * ncart_ket + iket] = B_0_2_2_1_1_0 + ( bAB_z * B_0_2_1_1_1_0 );

            // (f_0_2_1  f_1_0_2|_{i} = (g_1_2_1  d_0_0_2| + x_ab * (f_0_2_1  d_0_0_2|
            BRA_f_f[75 * ncart_ket + iket] = B_1_2_1_0_0_2 + ( bAB_x * B_0_2_1_0_0_2 );

            // (f_0_2_1  f_0_3_0|_{i} = (g_0_3_1  d_0_2_0| + y_ab * (f_0_2_1  d_0_2_0|
            BRA_f_f[76 * ncart_ket + iket] = B_0_3_1_0_2_0 + ( bAB_y * B_0_2_1_0_2_0 );

            // (f_0_2_1  f_0_2_1|_{i} = (g_0_2_2  d_0_2_0| + z_ab * (f_0_2_1  d_0_2_0|
            BRA_f_f[77 * ncart_ket + iket] = B_0_2_2_0_2_0 + ( bAB_z * B_0_2_1_0_2_0 );

            // (f_0_2_1  f_0_1_2|_{i} = (g_0_3_1  d_0_0_2| + y_ab * (f_0_2_1  d_0_0_2|
            BRA_f_f[78 * ncart_ket + iket] = B_0_3_1_0_0_2 + ( bAB_y * B_0_2_1_0_0_2 );

            // (f_0_2_1  f_0_0_3|_{i} = (g_0_2_2  d_0_0_2| + z_ab * (f_0_2_1  d_0_0_2|
            BRA_f_f[79 * ncart_ket + iket] = B_0_2_2_0_0_2 + ( bAB_z * B_0_2_1_0_0_2 );

            // (f_0_1_2  f_3_0_0|_{i} = (g_1_1_2  d_2_0_0| + x_ab * (f_0_1_2  d_2_0_0|
            BRA_f_f[80 * ncart_ket + iket] = B_1_1_2_2_0_0 + ( bAB_x * B_0_1_2_2_0_0 );

            // (f_0_1_2  f_2_1_0|_{i} = (g_0_2_2  d_2_0_0| + y_ab * (f_0_1_2  d_2_0_0|
            BRA_f_f[81 * ncart_ket + iket] = B_0_2_2_2_0_0 + ( bAB_y * B_0_1_2_2_0_0 );

            // (f_0_1_2  f_2_0_1|_{i} = (g_0_1_3  d_2_0_0| + z_ab * (f_0_1_2  d_2_0_0|
            BRA_f_f[82 * ncart_ket + iket] = B_0_1_3_2_0_0 + ( bAB_z * B_0_1_2_2_0_0 );

            // (f_0_1_2  f_1_2_0|_{i} = (g_1_1_2  d_0_2_0| + x_ab * (f_0_1_2  d_0_2_0|
            BRA_f_f[83 * ncart_ket + iket] = B_1_1_2_0_2_0 + ( bAB_x * B_0_1_2_0_2_0 );

            // (f_0_1_2  f_1_1_1|_{i} = (g_0_1_3  d_1_1_0| + z_ab * (f_0_1_2  d_1_1_0|
            BRA_f_f[84 * ncart_ket + iket] = B_0_1_3_1_1_0 + ( bAB_z * B_0_1_2_1_1_0 );

            // (f_0_1_2  f_1_0_2|_{i} = (g_1_1_2  d_0_0_2| + x_ab * (f_0_1_2  d_0_0_2|
            BRA_f_f[85 * ncart_ket + iket] = B_1_1_2_0_0_2 + ( bAB_x * B_0_1_2_0_0_2 );

            // (f_0_1_2  f_0_3_0|_{i} = (g_0_2_2  d_0_2_0| + y_ab * (f_0_1_2  d_0_2_0|
            BRA_f_f[86 * ncart_ket + iket] = B_0_2_2_0_2_0 + ( bAB_y * B_0_1_2_0_2_0 );

            // (f_0_1_2  f_0_2_1|_{i} = (g_0_1_3  d_0_2_0| + z_ab * (f_0_1_2  d_0_2_0|
            BRA_f_f[87 * ncart_ket + iket] = B_0_1_3_0_2_0 + ( bAB_z * B_0_1_2_0_2_0 );

            // (f_0_1_2  f_0_1_2|_{i} = (g_0_2_2  d_0_0_2| + y_ab * (f_0_1_2  d_0_0_2|
            BRA_f_f[88 * ncart_ket + iket] = B_0_2_2_0_0_2 + ( bAB_y * B_0_1_2_0_0_2 );

            // (f_0_1_2  f_0_0_3|_{i} = (g_0_1_3  d_0_0_2| + z_ab * (f_0_1_2  d_0_0_2|
            BRA_f_f[89 * ncart_ket + iket] = B_0_1_3_0_0_2 + ( bAB_z * B_0_1_2_0_0_2 );

            // (f_0_0_3  f_3_0_0|_{i} = (g_1_0_3  d_2_0_0| + x_ab * (f_0_0_3  d_2_0_0|
            BRA_f_f[90 * ncart_ket + iket] = B_1_0_3_2_0_0 + ( bAB_x * B_0_0_3_2_0_0 );

            // (f_0_0_3  f_2_1_0|_{i} = (g_0_1_3  d_2_0_0| + y_ab * (f_0_0_3  d_2_0_0|
            BRA_f_f[91 * ncart_ket + iket] = B_0_1_3_2_0_0 + ( bAB_y * B_0_0_3_2_0_0 );

            // (f_0_0_3  f_2_0_1|_{i} = (g_0_0_4  d_2_0_0| + z_ab * (f_0_0_3  d_2_0_0|
            BRA_f_f[92 * ncart_ket + iket] = B_0_0_4_2_0_0 + ( bAB_z * B_0_0_3_2_0_0 );

            // (f_0_0_3  f_1_2_0|_{i} = (g_1_0_3  d_0_2_0| + x_ab * (f_0_0_3  d_0_2_0|
            BRA_f_f[93 * ncart_ket + iket] = B_1_0_3_0_2_0 + ( bAB_x * B_0_0_3_0_2_0 );

            // (f_0_0_3  f_1_1_1|_{i} = (g_0_0_4  d_1_1_0| + z_ab * (f_0_0_3  d_1_1_0|
            BRA_f_f[94 * ncart_ket + iket] = B_0_0_4_1_1_0 + ( bAB_z * B_0_0_3_1_1_0 );

            // (f_0_0_3  f_1_0_2|_{i} = (g_1_0_3  d_0_0_2| + x_ab * (f_0_0_3  d_0_0_2|
            BRA_f_f[95 * ncart_ket + iket] = B_1_0_3_0_0_2 + ( bAB_x * B_0_0_3_0_0_2 );

            // (f_0_0_3  f_0_3_0|_{i} = (g_0_1_3  d_0_2_0| + y_ab * (f_0_0_3  d_0_2_0|
            BRA_f_f[96 * ncart_ket + iket] = B_0_1_3_0_2_0 + ( bAB_y * B_0_0_3_0_2_0 );

            // (f_0_0_3  f_0_2_1|_{i} = (g_0_0_4  d_0_2_0| + z_ab * (f_0_0_3  d_0_2_0|
            BRA_f_f[97 * ncart_ket + iket] = B_0_0_4_0_2_0 + ( bAB_z * B_0_0_3_0_2_0 );

            // (f_0_0_3  f_0_1_2|_{i} = (g_0_1_3  d_0_0_2| + y_ab * (f_0_0_3  d_0_0_2|
            BRA_f_f[98 * ncart_ket + iket] = B_0_1_3_0_0_2 + ( bAB_y * B_0_0_3_0_0_2 );

            // (f_0_0_3  f_0_0_3|_{i} = (g_0_0_4  d_0_0_2| + z_ab * (f_0_0_3  d_0_0_2|
            BRA_f_f[99 * ncart_ket + iket] = B_0_0_4_0_0_2 + ( bAB_z * B_0_0_3_0_0_2 );

        }


}


