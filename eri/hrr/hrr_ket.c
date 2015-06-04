
#include "vectorization.h"


    //////////////////////////////////////////////
    // KET: ( p p |
    // Steps: 9
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET__p_p(
                  double * const restrict KET_p_s,
                  double * const restrict KET_p_p,
                  double * const restrict KET_d_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 )
{
    int ibra;

        for(ibra = 0; ibra < ncart_bra; ++ibra)
        {
            // |p_1_0_0  p_1_0_0)_{i} = |d_2_0_0  s_0_0_0)_{t} + x_cd * |p_1_0_0  s_0_0_0)_{t}
            KET_p_p[ibra * 9 + 0] = KET_d_s[ibra * 6 + 0] + ( kCD_x * KET_p_s[ibra * 3 + 0] );

            // |p_1_0_0  p_0_1_0)_{i} = |d_1_1_0  s_0_0_0)_{t} + y_cd * |p_1_0_0  s_0_0_0)_{t}
            KET_p_p[ibra * 9 + 1] = KET_d_s[ibra * 6 + 1] + ( kCD_y * KET_p_s[ibra * 3 + 0] );

            // |p_1_0_0  p_0_0_1)_{i} = |d_1_0_1  s_0_0_0)_{t} + z_cd * |p_1_0_0  s_0_0_0)_{t}
            KET_p_p[ibra * 9 + 2] = KET_d_s[ibra * 6 + 2] + ( kCD_z * KET_p_s[ibra * 3 + 0] );

            // |p_0_1_0  p_1_0_0)_{i} = |d_1_1_0  s_0_0_0)_{t} + x_cd * |p_0_1_0  s_0_0_0)_{t}
            KET_p_p[ibra * 9 + 3] = KET_d_s[ibra * 6 + 1] + ( kCD_x * KET_p_s[ibra * 3 + 1] );

            // |p_0_1_0  p_0_1_0)_{i} = |d_0_2_0  s_0_0_0)_{t} + y_cd * |p_0_1_0  s_0_0_0)_{t}
            KET_p_p[ibra * 9 + 4] = KET_d_s[ibra * 6 + 3] + ( kCD_y * KET_p_s[ibra * 3 + 1] );

            // |p_0_1_0  p_0_0_1)_{i} = |d_0_1_1  s_0_0_0)_{t} + z_cd * |p_0_1_0  s_0_0_0)_{t}
            KET_p_p[ibra * 9 + 5] = KET_d_s[ibra * 6 + 4] + ( kCD_z * KET_p_s[ibra * 3 + 1] );

            // |p_0_0_1  p_1_0_0)_{i} = |d_1_0_1  s_0_0_0)_{t} + x_cd * |p_0_0_1  s_0_0_0)_{t}
            KET_p_p[ibra * 9 + 6] = KET_d_s[ibra * 6 + 2] + ( kCD_x * KET_p_s[ibra * 3 + 2] );

            // |p_0_0_1  p_0_1_0)_{i} = |d_0_1_1  s_0_0_0)_{t} + y_cd * |p_0_0_1  s_0_0_0)_{t}
            KET_p_p[ibra * 9 + 7] = KET_d_s[ibra * 6 + 4] + ( kCD_y * KET_p_s[ibra * 3 + 2] );

            // |p_0_0_1  p_0_0_1)_{i} = |d_0_0_2  s_0_0_0)_{t} + z_cd * |p_0_0_1  s_0_0_0)_{t}
            KET_p_p[ibra * 9 + 8] = KET_d_s[ibra * 6 + 5] + ( kCD_z * KET_p_s[ibra * 3 + 2] );

        }

}


    //////////////////////////////////////////////
    // KET: ( d p |
    // Steps: 18
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET__d_p(
                  double * const restrict KET_d_s,
                  double * const restrict KET_d_p,
                  double * const restrict KET_f_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 )
{
    int ibra;

        for(ibra = 0; ibra < ncart_bra; ++ibra)
        {
            // |d_2_0_0  p_1_0_0)_{i} = |f_3_0_0  s_0_0_0)_{t} + x_cd * |d_2_0_0  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 0] = KET_f_s[ibra * 10 + 0] + ( kCD_x * KET_d_s[ibra * 6 + 0] );

            // |d_2_0_0  p_0_1_0)_{i} = |f_2_1_0  s_0_0_0)_{t} + y_cd * |d_2_0_0  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 1] = KET_f_s[ibra * 10 + 1] + ( kCD_y * KET_d_s[ibra * 6 + 0] );

            // |d_2_0_0  p_0_0_1)_{i} = |f_2_0_1  s_0_0_0)_{t} + z_cd * |d_2_0_0  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 2] = KET_f_s[ibra * 10 + 2] + ( kCD_z * KET_d_s[ibra * 6 + 0] );

            // |d_1_1_0  p_1_0_0)_{i} = |f_2_1_0  s_0_0_0)_{t} + x_cd * |d_1_1_0  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 3] = KET_f_s[ibra * 10 + 1] + ( kCD_x * KET_d_s[ibra * 6 + 1] );

            // |d_1_1_0  p_0_1_0)_{i} = |f_1_2_0  s_0_0_0)_{t} + y_cd * |d_1_1_0  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 4] = KET_f_s[ibra * 10 + 3] + ( kCD_y * KET_d_s[ibra * 6 + 1] );

            // |d_1_1_0  p_0_0_1)_{i} = |f_1_1_1  s_0_0_0)_{t} + z_cd * |d_1_1_0  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 5] = KET_f_s[ibra * 10 + 4] + ( kCD_z * KET_d_s[ibra * 6 + 1] );

            // |d_1_0_1  p_1_0_0)_{i} = |f_2_0_1  s_0_0_0)_{t} + x_cd * |d_1_0_1  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 6] = KET_f_s[ibra * 10 + 2] + ( kCD_x * KET_d_s[ibra * 6 + 2] );

            // |d_1_0_1  p_0_1_0)_{i} = |f_1_1_1  s_0_0_0)_{t} + y_cd * |d_1_0_1  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 7] = KET_f_s[ibra * 10 + 4] + ( kCD_y * KET_d_s[ibra * 6 + 2] );

            // |d_1_0_1  p_0_0_1)_{i} = |f_1_0_2  s_0_0_0)_{t} + z_cd * |d_1_0_1  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 8] = KET_f_s[ibra * 10 + 5] + ( kCD_z * KET_d_s[ibra * 6 + 2] );

            // |d_0_2_0  p_1_0_0)_{i} = |f_1_2_0  s_0_0_0)_{t} + x_cd * |d_0_2_0  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 9] = KET_f_s[ibra * 10 + 3] + ( kCD_x * KET_d_s[ibra * 6 + 3] );

            // |d_0_2_0  p_0_1_0)_{i} = |f_0_3_0  s_0_0_0)_{t} + y_cd * |d_0_2_0  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 10] = KET_f_s[ibra * 10 + 6] + ( kCD_y * KET_d_s[ibra * 6 + 3] );

            // |d_0_2_0  p_0_0_1)_{i} = |f_0_2_1  s_0_0_0)_{t} + z_cd * |d_0_2_0  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 11] = KET_f_s[ibra * 10 + 7] + ( kCD_z * KET_d_s[ibra * 6 + 3] );

            // |d_0_1_1  p_1_0_0)_{i} = |f_1_1_1  s_0_0_0)_{t} + x_cd * |d_0_1_1  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 12] = KET_f_s[ibra * 10 + 4] + ( kCD_x * KET_d_s[ibra * 6 + 4] );

            // |d_0_1_1  p_0_1_0)_{i} = |f_0_2_1  s_0_0_0)_{t} + y_cd * |d_0_1_1  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 13] = KET_f_s[ibra * 10 + 7] + ( kCD_y * KET_d_s[ibra * 6 + 4] );

            // |d_0_1_1  p_0_0_1)_{i} = |f_0_1_2  s_0_0_0)_{t} + z_cd * |d_0_1_1  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 14] = KET_f_s[ibra * 10 + 8] + ( kCD_z * KET_d_s[ibra * 6 + 4] );

            // |d_0_0_2  p_1_0_0)_{i} = |f_1_0_2  s_0_0_0)_{t} + x_cd * |d_0_0_2  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 15] = KET_f_s[ibra * 10 + 5] + ( kCD_x * KET_d_s[ibra * 6 + 5] );

            // |d_0_0_2  p_0_1_0)_{i} = |f_0_1_2  s_0_0_0)_{t} + y_cd * |d_0_0_2  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 16] = KET_f_s[ibra * 10 + 8] + ( kCD_y * KET_d_s[ibra * 6 + 5] );

            // |d_0_0_2  p_0_0_1)_{i} = |f_0_0_3  s_0_0_0)_{t} + z_cd * |d_0_0_2  s_0_0_0)_{t}
            KET_d_p[ibra * 18 + 17] = KET_f_s[ibra * 10 + 9] + ( kCD_z * KET_d_s[ibra * 6 + 5] );

        }

}


    //////////////////////////////////////////////
    // KET: ( d d |
    // Steps: 79
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET__d_d(
                  double * const restrict KET_d_s,
                  double * const restrict KET_d_d,
                  double * const restrict KET_f_s,
                  double * const restrict KET_g_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 )
{
    int ibra;

        for(ibra = 0; ibra < ncart_bra; ++ibra)
        {
            // |d_2_0_0  p_1_0_0) = |f_3_0_0  s_0_0_0)_{t} + x_cd * |d_2_0_0  s_0_0_0)_{t}
            const double K_2_0_0_1_0_0 = KET_f_s[ibra * 10 + 0] + ( kCD_x * KET_d_s[ibra * 6 + 0] );

            // |d_2_0_0  p_0_1_0) = |f_2_1_0  s_0_0_0)_{t} + y_cd * |d_2_0_0  s_0_0_0)_{t}
            const double K_2_0_0_0_1_0 = KET_f_s[ibra * 10 + 1] + ( kCD_y * KET_d_s[ibra * 6 + 0] );

            // |d_2_0_0  p_0_0_1) = |f_2_0_1  s_0_0_0)_{t} + z_cd * |d_2_0_0  s_0_0_0)_{t}
            const double K_2_0_0_0_0_1 = KET_f_s[ibra * 10 + 2] + ( kCD_z * KET_d_s[ibra * 6 + 0] );

            // |d_1_1_0  p_1_0_0) = |f_2_1_0  s_0_0_0)_{t} + x_cd * |d_1_1_0  s_0_0_0)_{t}
            const double K_1_1_0_1_0_0 = KET_f_s[ibra * 10 + 1] + ( kCD_x * KET_d_s[ibra * 6 + 1] );

            // |d_1_1_0  p_0_1_0) = |f_1_2_0  s_0_0_0)_{t} + y_cd * |d_1_1_0  s_0_0_0)_{t}
            const double K_1_1_0_0_1_0 = KET_f_s[ibra * 10 + 3] + ( kCD_y * KET_d_s[ibra * 6 + 1] );

            // |d_1_1_0  p_0_0_1) = |f_1_1_1  s_0_0_0)_{t} + z_cd * |d_1_1_0  s_0_0_0)_{t}
            const double K_1_1_0_0_0_1 = KET_f_s[ibra * 10 + 4] + ( kCD_z * KET_d_s[ibra * 6 + 1] );

            // |d_1_0_1  p_1_0_0) = |f_2_0_1  s_0_0_0)_{t} + x_cd * |d_1_0_1  s_0_0_0)_{t}
            const double K_1_0_1_1_0_0 = KET_f_s[ibra * 10 + 2] + ( kCD_x * KET_d_s[ibra * 6 + 2] );

            // |d_1_0_1  p_0_1_0) = |f_1_1_1  s_0_0_0)_{t} + y_cd * |d_1_0_1  s_0_0_0)_{t}
            const double K_1_0_1_0_1_0 = KET_f_s[ibra * 10 + 4] + ( kCD_y * KET_d_s[ibra * 6 + 2] );

            // |d_1_0_1  p_0_0_1) = |f_1_0_2  s_0_0_0)_{t} + z_cd * |d_1_0_1  s_0_0_0)_{t}
            const double K_1_0_1_0_0_1 = KET_f_s[ibra * 10 + 5] + ( kCD_z * KET_d_s[ibra * 6 + 2] );

            // |d_0_2_0  p_1_0_0) = |f_1_2_0  s_0_0_0)_{t} + x_cd * |d_0_2_0  s_0_0_0)_{t}
            const double K_0_2_0_1_0_0 = KET_f_s[ibra * 10 + 3] + ( kCD_x * KET_d_s[ibra * 6 + 3] );

            // |d_0_2_0  p_0_1_0) = |f_0_3_0  s_0_0_0)_{t} + y_cd * |d_0_2_0  s_0_0_0)_{t}
            const double K_0_2_0_0_1_0 = KET_f_s[ibra * 10 + 6] + ( kCD_y * KET_d_s[ibra * 6 + 3] );

            // |d_0_2_0  p_0_0_1) = |f_0_2_1  s_0_0_0)_{t} + z_cd * |d_0_2_0  s_0_0_0)_{t}
            const double K_0_2_0_0_0_1 = KET_f_s[ibra * 10 + 7] + ( kCD_z * KET_d_s[ibra * 6 + 3] );

            // |d_0_1_1  p_1_0_0) = |f_1_1_1  s_0_0_0)_{t} + x_cd * |d_0_1_1  s_0_0_0)_{t}
            const double K_0_1_1_1_0_0 = KET_f_s[ibra * 10 + 4] + ( kCD_x * KET_d_s[ibra * 6 + 4] );

            // |d_0_1_1  p_0_1_0) = |f_0_2_1  s_0_0_0)_{t} + y_cd * |d_0_1_1  s_0_0_0)_{t}
            const double K_0_1_1_0_1_0 = KET_f_s[ibra * 10 + 7] + ( kCD_y * KET_d_s[ibra * 6 + 4] );

            // |d_0_1_1  p_0_0_1) = |f_0_1_2  s_0_0_0)_{t} + z_cd * |d_0_1_1  s_0_0_0)_{t}
            const double K_0_1_1_0_0_1 = KET_f_s[ibra * 10 + 8] + ( kCD_z * KET_d_s[ibra * 6 + 4] );

            // |d_0_0_2  p_1_0_0) = |f_1_0_2  s_0_0_0)_{t} + x_cd * |d_0_0_2  s_0_0_0)_{t}
            const double K_0_0_2_1_0_0 = KET_f_s[ibra * 10 + 5] + ( kCD_x * KET_d_s[ibra * 6 + 5] );

            // |d_0_0_2  p_0_1_0) = |f_0_1_2  s_0_0_0)_{t} + y_cd * |d_0_0_2  s_0_0_0)_{t}
            const double K_0_0_2_0_1_0 = KET_f_s[ibra * 10 + 8] + ( kCD_y * KET_d_s[ibra * 6 + 5] );

            // |d_0_0_2  p_0_0_1) = |f_0_0_3  s_0_0_0)_{t} + z_cd * |d_0_0_2  s_0_0_0)_{t}
            const double K_0_0_2_0_0_1 = KET_f_s[ibra * 10 + 9] + ( kCD_z * KET_d_s[ibra * 6 + 5] );

            // |f_3_0_0  p_1_0_0) = |g_4_0_0  s_0_0_0)_{t} + x_cd * |f_3_0_0  s_0_0_0)_{t}
            const double K_3_0_0_1_0_0 = KET_g_s[ibra * 15 + 0] + ( kCD_x * KET_f_s[ibra * 10 + 0] );

            // |f_2_1_0  p_1_0_0) = |g_3_1_0  s_0_0_0)_{t} + x_cd * |f_2_1_0  s_0_0_0)_{t}
            const double K_2_1_0_1_0_0 = KET_g_s[ibra * 15 + 1] + ( kCD_x * KET_f_s[ibra * 10 + 1] );

            // |f_2_1_0  p_0_1_0) = |g_2_2_0  s_0_0_0)_{t} + y_cd * |f_2_1_0  s_0_0_0)_{t}
            const double K_2_1_0_0_1_0 = KET_g_s[ibra * 15 + 3] + ( kCD_y * KET_f_s[ibra * 10 + 1] );

            // |f_2_0_1  p_1_0_0) = |g_3_0_1  s_0_0_0)_{t} + x_cd * |f_2_0_1  s_0_0_0)_{t}
            const double K_2_0_1_1_0_0 = KET_g_s[ibra * 15 + 2] + ( kCD_x * KET_f_s[ibra * 10 + 2] );

            // |f_2_0_1  p_0_1_0) = |g_2_1_1  s_0_0_0)_{t} + y_cd * |f_2_0_1  s_0_0_0)_{t}
            const double K_2_0_1_0_1_0 = KET_g_s[ibra * 15 + 4] + ( kCD_y * KET_f_s[ibra * 10 + 2] );

            // |f_2_0_1  p_0_0_1) = |g_2_0_2  s_0_0_0)_{t} + z_cd * |f_2_0_1  s_0_0_0)_{t}
            const double K_2_0_1_0_0_1 = KET_g_s[ibra * 15 + 5] + ( kCD_z * KET_f_s[ibra * 10 + 2] );

            // |f_1_2_0  p_1_0_0) = |g_2_2_0  s_0_0_0)_{t} + x_cd * |f_1_2_0  s_0_0_0)_{t}
            const double K_1_2_0_1_0_0 = KET_g_s[ibra * 15 + 3] + ( kCD_x * KET_f_s[ibra * 10 + 3] );

            // |f_1_2_0  p_0_1_0) = |g_1_3_0  s_0_0_0)_{t} + y_cd * |f_1_2_0  s_0_0_0)_{t}
            const double K_1_2_0_0_1_0 = KET_g_s[ibra * 15 + 6] + ( kCD_y * KET_f_s[ibra * 10 + 3] );

            // |f_1_1_1  p_1_0_0) = |g_2_1_1  s_0_0_0)_{t} + x_cd * |f_1_1_1  s_0_0_0)_{t}
            const double K_1_1_1_1_0_0 = KET_g_s[ibra * 15 + 4] + ( kCD_x * KET_f_s[ibra * 10 + 4] );

            // |f_1_1_1  p_0_1_0) = |g_1_2_1  s_0_0_0)_{t} + y_cd * |f_1_1_1  s_0_0_0)_{t}
            const double K_1_1_1_0_1_0 = KET_g_s[ibra * 15 + 7] + ( kCD_y * KET_f_s[ibra * 10 + 4] );

            // |f_1_1_1  p_0_0_1) = |g_1_1_2  s_0_0_0)_{t} + z_cd * |f_1_1_1  s_0_0_0)_{t}
            const double K_1_1_1_0_0_1 = KET_g_s[ibra * 15 + 8] + ( kCD_z * KET_f_s[ibra * 10 + 4] );

            // |f_1_0_2  p_1_0_0) = |g_2_0_2  s_0_0_0)_{t} + x_cd * |f_1_0_2  s_0_0_0)_{t}
            const double K_1_0_2_1_0_0 = KET_g_s[ibra * 15 + 5] + ( kCD_x * KET_f_s[ibra * 10 + 5] );

            // |f_1_0_2  p_0_1_0) = |g_1_1_2  s_0_0_0)_{t} + y_cd * |f_1_0_2  s_0_0_0)_{t}
            const double K_1_0_2_0_1_0 = KET_g_s[ibra * 15 + 8] + ( kCD_y * KET_f_s[ibra * 10 + 5] );

            // |f_1_0_2  p_0_0_1) = |g_1_0_3  s_0_0_0)_{t} + z_cd * |f_1_0_2  s_0_0_0)_{t}
            const double K_1_0_2_0_0_1 = KET_g_s[ibra * 15 + 9] + ( kCD_z * KET_f_s[ibra * 10 + 5] );

            // |f_0_3_0  p_1_0_0) = |g_1_3_0  s_0_0_0)_{t} + x_cd * |f_0_3_0  s_0_0_0)_{t}
            const double K_0_3_0_1_0_0 = KET_g_s[ibra * 15 + 6] + ( kCD_x * KET_f_s[ibra * 10 + 6] );

            // |f_0_3_0  p_0_1_0) = |g_0_4_0  s_0_0_0)_{t} + y_cd * |f_0_3_0  s_0_0_0)_{t}
            const double K_0_3_0_0_1_0 = KET_g_s[ibra * 15 + 10] + ( kCD_y * KET_f_s[ibra * 10 + 6] );

            // |f_0_2_1  p_1_0_0) = |g_1_2_1  s_0_0_0)_{t} + x_cd * |f_0_2_1  s_0_0_0)_{t}
            const double K_0_2_1_1_0_0 = KET_g_s[ibra * 15 + 7] + ( kCD_x * KET_f_s[ibra * 10 + 7] );

            // |f_0_2_1  p_0_1_0) = |g_0_3_1  s_0_0_0)_{t} + y_cd * |f_0_2_1  s_0_0_0)_{t}
            const double K_0_2_1_0_1_0 = KET_g_s[ibra * 15 + 11] + ( kCD_y * KET_f_s[ibra * 10 + 7] );

            // |f_0_2_1  p_0_0_1) = |g_0_2_2  s_0_0_0)_{t} + z_cd * |f_0_2_1  s_0_0_0)_{t}
            const double K_0_2_1_0_0_1 = KET_g_s[ibra * 15 + 12] + ( kCD_z * KET_f_s[ibra * 10 + 7] );

            // |f_0_1_2  p_1_0_0) = |g_1_1_2  s_0_0_0)_{t} + x_cd * |f_0_1_2  s_0_0_0)_{t}
            const double K_0_1_2_1_0_0 = KET_g_s[ibra * 15 + 8] + ( kCD_x * KET_f_s[ibra * 10 + 8] );

            // |f_0_1_2  p_0_1_0) = |g_0_2_2  s_0_0_0)_{t} + y_cd * |f_0_1_2  s_0_0_0)_{t}
            const double K_0_1_2_0_1_0 = KET_g_s[ibra * 15 + 12] + ( kCD_y * KET_f_s[ibra * 10 + 8] );

            // |f_0_1_2  p_0_0_1) = |g_0_1_3  s_0_0_0)_{t} + z_cd * |f_0_1_2  s_0_0_0)_{t}
            const double K_0_1_2_0_0_1 = KET_g_s[ibra * 15 + 13] + ( kCD_z * KET_f_s[ibra * 10 + 8] );

            // |f_0_0_3  p_1_0_0) = |g_1_0_3  s_0_0_0)_{t} + x_cd * |f_0_0_3  s_0_0_0)_{t}
            const double K_0_0_3_1_0_0 = KET_g_s[ibra * 15 + 9] + ( kCD_x * KET_f_s[ibra * 10 + 9] );

            // |f_0_0_3  p_0_1_0) = |g_0_1_3  s_0_0_0)_{t} + y_cd * |f_0_0_3  s_0_0_0)_{t}
            const double K_0_0_3_0_1_0 = KET_g_s[ibra * 15 + 13] + ( kCD_y * KET_f_s[ibra * 10 + 9] );

            // |f_0_0_3  p_0_0_1) = |g_0_0_4  s_0_0_0)_{t} + z_cd * |f_0_0_3  s_0_0_0)_{t}
            const double K_0_0_3_0_0_1 = KET_g_s[ibra * 15 + 14] + ( kCD_z * KET_f_s[ibra * 10 + 9] );

            // |d_2_0_0  d_2_0_0)_{i} = |f_3_0_0  p_1_0_0) + x_cd * |d_2_0_0  p_1_0_0)
            KET_d_d[ibra * 36 + 0] = K_3_0_0_1_0_0 + ( kCD_x * K_2_0_0_1_0_0 );

            // |d_2_0_0  d_1_1_0)_{i} = |f_2_1_0  p_1_0_0) + y_cd * |d_2_0_0  p_1_0_0)
            KET_d_d[ibra * 36 + 1] = K_2_1_0_1_0_0 + ( kCD_y * K_2_0_0_1_0_0 );

            // |d_2_0_0  d_1_0_1)_{i} = |f_2_0_1  p_1_0_0) + z_cd * |d_2_0_0  p_1_0_0)
            KET_d_d[ibra * 36 + 2] = K_2_0_1_1_0_0 + ( kCD_z * K_2_0_0_1_0_0 );

            // |d_2_0_0  d_0_2_0)_{i} = |f_2_1_0  p_0_1_0) + y_cd * |d_2_0_0  p_0_1_0)
            KET_d_d[ibra * 36 + 3] = K_2_1_0_0_1_0 + ( kCD_y * K_2_0_0_0_1_0 );

            // |d_2_0_0  d_0_1_1)_{i} = |f_2_0_1  p_0_1_0) + z_cd * |d_2_0_0  p_0_1_0)
            KET_d_d[ibra * 36 + 4] = K_2_0_1_0_1_0 + ( kCD_z * K_2_0_0_0_1_0 );

            // |d_2_0_0  d_0_0_2)_{i} = |f_2_0_1  p_0_0_1) + z_cd * |d_2_0_0  p_0_0_1)
            KET_d_d[ibra * 36 + 5] = K_2_0_1_0_0_1 + ( kCD_z * K_2_0_0_0_0_1 );

            // |d_1_1_0  d_2_0_0)_{i} = |f_2_1_0  p_1_0_0) + x_cd * |d_1_1_0  p_1_0_0)
            KET_d_d[ibra * 36 + 6] = K_2_1_0_1_0_0 + ( kCD_x * K_1_1_0_1_0_0 );

            // |d_1_1_0  d_1_1_0)_{i} = |f_1_2_0  p_1_0_0) + y_cd * |d_1_1_0  p_1_0_0)
            KET_d_d[ibra * 36 + 7] = K_1_2_0_1_0_0 + ( kCD_y * K_1_1_0_1_0_0 );

            // |d_1_1_0  d_1_0_1)_{i} = |f_1_1_1  p_1_0_0) + z_cd * |d_1_1_0  p_1_0_0)
            KET_d_d[ibra * 36 + 8] = K_1_1_1_1_0_0 + ( kCD_z * K_1_1_0_1_0_0 );

            // |d_1_1_0  d_0_2_0)_{i} = |f_1_2_0  p_0_1_0) + y_cd * |d_1_1_0  p_0_1_0)
            KET_d_d[ibra * 36 + 9] = K_1_2_0_0_1_0 + ( kCD_y * K_1_1_0_0_1_0 );

            // |d_1_1_0  d_0_1_1)_{i} = |f_1_1_1  p_0_1_0) + z_cd * |d_1_1_0  p_0_1_0)
            KET_d_d[ibra * 36 + 10] = K_1_1_1_0_1_0 + ( kCD_z * K_1_1_0_0_1_0 );

            // |d_1_1_0  d_0_0_2)_{i} = |f_1_1_1  p_0_0_1) + z_cd * |d_1_1_0  p_0_0_1)
            KET_d_d[ibra * 36 + 11] = K_1_1_1_0_0_1 + ( kCD_z * K_1_1_0_0_0_1 );

            // |d_1_0_1  d_2_0_0)_{i} = |f_2_0_1  p_1_0_0) + x_cd * |d_1_0_1  p_1_0_0)
            KET_d_d[ibra * 36 + 12] = K_2_0_1_1_0_0 + ( kCD_x * K_1_0_1_1_0_0 );

            // |d_1_0_1  d_1_1_0)_{i} = |f_1_1_1  p_1_0_0) + y_cd * |d_1_0_1  p_1_0_0)
            KET_d_d[ibra * 36 + 13] = K_1_1_1_1_0_0 + ( kCD_y * K_1_0_1_1_0_0 );

            // |d_1_0_1  d_1_0_1)_{i} = |f_1_0_2  p_1_0_0) + z_cd * |d_1_0_1  p_1_0_0)
            KET_d_d[ibra * 36 + 14] = K_1_0_2_1_0_0 + ( kCD_z * K_1_0_1_1_0_0 );

            // |d_1_0_1  d_0_2_0)_{i} = |f_1_1_1  p_0_1_0) + y_cd * |d_1_0_1  p_0_1_0)
            KET_d_d[ibra * 36 + 15] = K_1_1_1_0_1_0 + ( kCD_y * K_1_0_1_0_1_0 );

            // |d_1_0_1  d_0_1_1)_{i} = |f_1_0_2  p_0_1_0) + z_cd * |d_1_0_1  p_0_1_0)
            KET_d_d[ibra * 36 + 16] = K_1_0_2_0_1_0 + ( kCD_z * K_1_0_1_0_1_0 );

            // |d_1_0_1  d_0_0_2)_{i} = |f_1_0_2  p_0_0_1) + z_cd * |d_1_0_1  p_0_0_1)
            KET_d_d[ibra * 36 + 17] = K_1_0_2_0_0_1 + ( kCD_z * K_1_0_1_0_0_1 );

            // |d_0_2_0  d_2_0_0)_{i} = |f_1_2_0  p_1_0_0) + x_cd * |d_0_2_0  p_1_0_0)
            KET_d_d[ibra * 36 + 18] = K_1_2_0_1_0_0 + ( kCD_x * K_0_2_0_1_0_0 );

            // |d_0_2_0  d_1_1_0)_{i} = |f_0_3_0  p_1_0_0) + y_cd * |d_0_2_0  p_1_0_0)
            KET_d_d[ibra * 36 + 19] = K_0_3_0_1_0_0 + ( kCD_y * K_0_2_0_1_0_0 );

            // |d_0_2_0  d_1_0_1)_{i} = |f_0_2_1  p_1_0_0) + z_cd * |d_0_2_0  p_1_0_0)
            KET_d_d[ibra * 36 + 20] = K_0_2_1_1_0_0 + ( kCD_z * K_0_2_0_1_0_0 );

            // |d_0_2_0  d_0_2_0)_{i} = |f_0_3_0  p_0_1_0) + y_cd * |d_0_2_0  p_0_1_0)
            KET_d_d[ibra * 36 + 21] = K_0_3_0_0_1_0 + ( kCD_y * K_0_2_0_0_1_0 );

            // |d_0_2_0  d_0_1_1)_{i} = |f_0_2_1  p_0_1_0) + z_cd * |d_0_2_0  p_0_1_0)
            KET_d_d[ibra * 36 + 22] = K_0_2_1_0_1_0 + ( kCD_z * K_0_2_0_0_1_0 );

            // |d_0_2_0  d_0_0_2)_{i} = |f_0_2_1  p_0_0_1) + z_cd * |d_0_2_0  p_0_0_1)
            KET_d_d[ibra * 36 + 23] = K_0_2_1_0_0_1 + ( kCD_z * K_0_2_0_0_0_1 );

            // |d_0_1_1  d_2_0_0)_{i} = |f_1_1_1  p_1_0_0) + x_cd * |d_0_1_1  p_1_0_0)
            KET_d_d[ibra * 36 + 24] = K_1_1_1_1_0_0 + ( kCD_x * K_0_1_1_1_0_0 );

            // |d_0_1_1  d_1_1_0)_{i} = |f_0_2_1  p_1_0_0) + y_cd * |d_0_1_1  p_1_0_0)
            KET_d_d[ibra * 36 + 25] = K_0_2_1_1_0_0 + ( kCD_y * K_0_1_1_1_0_0 );

            // |d_0_1_1  d_1_0_1)_{i} = |f_0_1_2  p_1_0_0) + z_cd * |d_0_1_1  p_1_0_0)
            KET_d_d[ibra * 36 + 26] = K_0_1_2_1_0_0 + ( kCD_z * K_0_1_1_1_0_0 );

            // |d_0_1_1  d_0_2_0)_{i} = |f_0_2_1  p_0_1_0) + y_cd * |d_0_1_1  p_0_1_0)
            KET_d_d[ibra * 36 + 27] = K_0_2_1_0_1_0 + ( kCD_y * K_0_1_1_0_1_0 );

            // |d_0_1_1  d_0_1_1)_{i} = |f_0_1_2  p_0_1_0) + z_cd * |d_0_1_1  p_0_1_0)
            KET_d_d[ibra * 36 + 28] = K_0_1_2_0_1_0 + ( kCD_z * K_0_1_1_0_1_0 );

            // |d_0_1_1  d_0_0_2)_{i} = |f_0_1_2  p_0_0_1) + z_cd * |d_0_1_1  p_0_0_1)
            KET_d_d[ibra * 36 + 29] = K_0_1_2_0_0_1 + ( kCD_z * K_0_1_1_0_0_1 );

            // |d_0_0_2  d_2_0_0)_{i} = |f_1_0_2  p_1_0_0) + x_cd * |d_0_0_2  p_1_0_0)
            KET_d_d[ibra * 36 + 30] = K_1_0_2_1_0_0 + ( kCD_x * K_0_0_2_1_0_0 );

            // |d_0_0_2  d_1_1_0)_{i} = |f_0_1_2  p_1_0_0) + y_cd * |d_0_0_2  p_1_0_0)
            KET_d_d[ibra * 36 + 31] = K_0_1_2_1_0_0 + ( kCD_y * K_0_0_2_1_0_0 );

            // |d_0_0_2  d_1_0_1)_{i} = |f_0_0_3  p_1_0_0) + z_cd * |d_0_0_2  p_1_0_0)
            KET_d_d[ibra * 36 + 32] = K_0_0_3_1_0_0 + ( kCD_z * K_0_0_2_1_0_0 );

            // |d_0_0_2  d_0_2_0)_{i} = |f_0_1_2  p_0_1_0) + y_cd * |d_0_0_2  p_0_1_0)
            KET_d_d[ibra * 36 + 33] = K_0_1_2_0_1_0 + ( kCD_y * K_0_0_2_0_1_0 );

            // |d_0_0_2  d_0_1_1)_{i} = |f_0_0_3  p_0_1_0) + z_cd * |d_0_0_2  p_0_1_0)
            KET_d_d[ibra * 36 + 34] = K_0_0_3_0_1_0 + ( kCD_z * K_0_0_2_0_1_0 );

            // |d_0_0_2  d_0_0_2)_{i} = |f_0_0_3  p_0_0_1) + z_cd * |d_0_0_2  p_0_0_1)
            KET_d_d[ibra * 36 + 35] = K_0_0_3_0_0_1 + ( kCD_z * K_0_0_2_0_0_1 );

        }

}


    //////////////////////////////////////////////
    // KET: ( f p |
    // Steps: 30
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET__f_p(
                  double * const restrict KET_f_s,
                  double * const restrict KET_f_p,
                  double * const restrict KET_g_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 )
{
    int ibra;

        for(ibra = 0; ibra < ncart_bra; ++ibra)
        {
            // |f_3_0_0  p_1_0_0)_{i} = |g_4_0_0  s_0_0_0)_{t} + x_cd * |f_3_0_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 0] = KET_g_s[ibra * 15 + 0] + ( kCD_x * KET_f_s[ibra * 10 + 0] );

            // |f_3_0_0  p_0_1_0)_{i} = |g_3_1_0  s_0_0_0)_{t} + y_cd * |f_3_0_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 1] = KET_g_s[ibra * 15 + 1] + ( kCD_y * KET_f_s[ibra * 10 + 0] );

            // |f_3_0_0  p_0_0_1)_{i} = |g_3_0_1  s_0_0_0)_{t} + z_cd * |f_3_0_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 2] = KET_g_s[ibra * 15 + 2] + ( kCD_z * KET_f_s[ibra * 10 + 0] );

            // |f_2_1_0  p_1_0_0)_{i} = |g_3_1_0  s_0_0_0)_{t} + x_cd * |f_2_1_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 3] = KET_g_s[ibra * 15 + 1] + ( kCD_x * KET_f_s[ibra * 10 + 1] );

            // |f_2_1_0  p_0_1_0)_{i} = |g_2_2_0  s_0_0_0)_{t} + y_cd * |f_2_1_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 4] = KET_g_s[ibra * 15 + 3] + ( kCD_y * KET_f_s[ibra * 10 + 1] );

            // |f_2_1_0  p_0_0_1)_{i} = |g_2_1_1  s_0_0_0)_{t} + z_cd * |f_2_1_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 5] = KET_g_s[ibra * 15 + 4] + ( kCD_z * KET_f_s[ibra * 10 + 1] );

            // |f_2_0_1  p_1_0_0)_{i} = |g_3_0_1  s_0_0_0)_{t} + x_cd * |f_2_0_1  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 6] = KET_g_s[ibra * 15 + 2] + ( kCD_x * KET_f_s[ibra * 10 + 2] );

            // |f_2_0_1  p_0_1_0)_{i} = |g_2_1_1  s_0_0_0)_{t} + y_cd * |f_2_0_1  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 7] = KET_g_s[ibra * 15 + 4] + ( kCD_y * KET_f_s[ibra * 10 + 2] );

            // |f_2_0_1  p_0_0_1)_{i} = |g_2_0_2  s_0_0_0)_{t} + z_cd * |f_2_0_1  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 8] = KET_g_s[ibra * 15 + 5] + ( kCD_z * KET_f_s[ibra * 10 + 2] );

            // |f_1_2_0  p_1_0_0)_{i} = |g_2_2_0  s_0_0_0)_{t} + x_cd * |f_1_2_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 9] = KET_g_s[ibra * 15 + 3] + ( kCD_x * KET_f_s[ibra * 10 + 3] );

            // |f_1_2_0  p_0_1_0)_{i} = |g_1_3_0  s_0_0_0)_{t} + y_cd * |f_1_2_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 10] = KET_g_s[ibra * 15 + 6] + ( kCD_y * KET_f_s[ibra * 10 + 3] );

            // |f_1_2_0  p_0_0_1)_{i} = |g_1_2_1  s_0_0_0)_{t} + z_cd * |f_1_2_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 11] = KET_g_s[ibra * 15 + 7] + ( kCD_z * KET_f_s[ibra * 10 + 3] );

            // |f_1_1_1  p_1_0_0)_{i} = |g_2_1_1  s_0_0_0)_{t} + x_cd * |f_1_1_1  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 12] = KET_g_s[ibra * 15 + 4] + ( kCD_x * KET_f_s[ibra * 10 + 4] );

            // |f_1_1_1  p_0_1_0)_{i} = |g_1_2_1  s_0_0_0)_{t} + y_cd * |f_1_1_1  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 13] = KET_g_s[ibra * 15 + 7] + ( kCD_y * KET_f_s[ibra * 10 + 4] );

            // |f_1_1_1  p_0_0_1)_{i} = |g_1_1_2  s_0_0_0)_{t} + z_cd * |f_1_1_1  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 14] = KET_g_s[ibra * 15 + 8] + ( kCD_z * KET_f_s[ibra * 10 + 4] );

            // |f_1_0_2  p_1_0_0)_{i} = |g_2_0_2  s_0_0_0)_{t} + x_cd * |f_1_0_2  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 15] = KET_g_s[ibra * 15 + 5] + ( kCD_x * KET_f_s[ibra * 10 + 5] );

            // |f_1_0_2  p_0_1_0)_{i} = |g_1_1_2  s_0_0_0)_{t} + y_cd * |f_1_0_2  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 16] = KET_g_s[ibra * 15 + 8] + ( kCD_y * KET_f_s[ibra * 10 + 5] );

            // |f_1_0_2  p_0_0_1)_{i} = |g_1_0_3  s_0_0_0)_{t} + z_cd * |f_1_0_2  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 17] = KET_g_s[ibra * 15 + 9] + ( kCD_z * KET_f_s[ibra * 10 + 5] );

            // |f_0_3_0  p_1_0_0)_{i} = |g_1_3_0  s_0_0_0)_{t} + x_cd * |f_0_3_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 18] = KET_g_s[ibra * 15 + 6] + ( kCD_x * KET_f_s[ibra * 10 + 6] );

            // |f_0_3_0  p_0_1_0)_{i} = |g_0_4_0  s_0_0_0)_{t} + y_cd * |f_0_3_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 19] = KET_g_s[ibra * 15 + 10] + ( kCD_y * KET_f_s[ibra * 10 + 6] );

            // |f_0_3_0  p_0_0_1)_{i} = |g_0_3_1  s_0_0_0)_{t} + z_cd * |f_0_3_0  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 20] = KET_g_s[ibra * 15 + 11] + ( kCD_z * KET_f_s[ibra * 10 + 6] );

            // |f_0_2_1  p_1_0_0)_{i} = |g_1_2_1  s_0_0_0)_{t} + x_cd * |f_0_2_1  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 21] = KET_g_s[ibra * 15 + 7] + ( kCD_x * KET_f_s[ibra * 10 + 7] );

            // |f_0_2_1  p_0_1_0)_{i} = |g_0_3_1  s_0_0_0)_{t} + y_cd * |f_0_2_1  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 22] = KET_g_s[ibra * 15 + 11] + ( kCD_y * KET_f_s[ibra * 10 + 7] );

            // |f_0_2_1  p_0_0_1)_{i} = |g_0_2_2  s_0_0_0)_{t} + z_cd * |f_0_2_1  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 23] = KET_g_s[ibra * 15 + 12] + ( kCD_z * KET_f_s[ibra * 10 + 7] );

            // |f_0_1_2  p_1_0_0)_{i} = |g_1_1_2  s_0_0_0)_{t} + x_cd * |f_0_1_2  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 24] = KET_g_s[ibra * 15 + 8] + ( kCD_x * KET_f_s[ibra * 10 + 8] );

            // |f_0_1_2  p_0_1_0)_{i} = |g_0_2_2  s_0_0_0)_{t} + y_cd * |f_0_1_2  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 25] = KET_g_s[ibra * 15 + 12] + ( kCD_y * KET_f_s[ibra * 10 + 8] );

            // |f_0_1_2  p_0_0_1)_{i} = |g_0_1_3  s_0_0_0)_{t} + z_cd * |f_0_1_2  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 26] = KET_g_s[ibra * 15 + 13] + ( kCD_z * KET_f_s[ibra * 10 + 8] );

            // |f_0_0_3  p_1_0_0)_{i} = |g_1_0_3  s_0_0_0)_{t} + x_cd * |f_0_0_3  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 27] = KET_g_s[ibra * 15 + 9] + ( kCD_x * KET_f_s[ibra * 10 + 9] );

            // |f_0_0_3  p_0_1_0)_{i} = |g_0_1_3  s_0_0_0)_{t} + y_cd * |f_0_0_3  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 28] = KET_g_s[ibra * 15 + 13] + ( kCD_y * KET_f_s[ibra * 10 + 9] );

            // |f_0_0_3  p_0_0_1)_{i} = |g_0_0_4  s_0_0_0)_{t} + z_cd * |f_0_0_3  s_0_0_0)_{t}
            KET_f_p[ibra * 30 + 29] = KET_g_s[ibra * 15 + 14] + ( kCD_z * KET_f_s[ibra * 10 + 9] );

        }

}


    //////////////////////////////////////////////
    // KET: ( f d |
    // Steps: 129
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET__f_d(
                  double * const restrict KET_f_s,
                  double * const restrict KET_f_d,
                  double * const restrict KET_g_s,
                  double * const restrict KET_h_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 )
{
    int ibra;

        for(ibra = 0; ibra < ncart_bra; ++ibra)
        {
            // |f_3_0_0  p_1_0_0) = |g_4_0_0  s_0_0_0)_{t} + x_cd * |f_3_0_0  s_0_0_0)_{t}
            const double K_3_0_0_1_0_0 = KET_g_s[ibra * 15 + 0] + ( kCD_x * KET_f_s[ibra * 10 + 0] );

            // |f_3_0_0  p_0_1_0) = |g_3_1_0  s_0_0_0)_{t} + y_cd * |f_3_0_0  s_0_0_0)_{t}
            const double K_3_0_0_0_1_0 = KET_g_s[ibra * 15 + 1] + ( kCD_y * KET_f_s[ibra * 10 + 0] );

            // |f_3_0_0  p_0_0_1) = |g_3_0_1  s_0_0_0)_{t} + z_cd * |f_3_0_0  s_0_0_0)_{t}
            const double K_3_0_0_0_0_1 = KET_g_s[ibra * 15 + 2] + ( kCD_z * KET_f_s[ibra * 10 + 0] );

            // |f_2_1_0  p_1_0_0) = |g_3_1_0  s_0_0_0)_{t} + x_cd * |f_2_1_0  s_0_0_0)_{t}
            const double K_2_1_0_1_0_0 = KET_g_s[ibra * 15 + 1] + ( kCD_x * KET_f_s[ibra * 10 + 1] );

            // |f_2_1_0  p_0_1_0) = |g_2_2_0  s_0_0_0)_{t} + y_cd * |f_2_1_0  s_0_0_0)_{t}
            const double K_2_1_0_0_1_0 = KET_g_s[ibra * 15 + 3] + ( kCD_y * KET_f_s[ibra * 10 + 1] );

            // |f_2_1_0  p_0_0_1) = |g_2_1_1  s_0_0_0)_{t} + z_cd * |f_2_1_0  s_0_0_0)_{t}
            const double K_2_1_0_0_0_1 = KET_g_s[ibra * 15 + 4] + ( kCD_z * KET_f_s[ibra * 10 + 1] );

            // |f_2_0_1  p_1_0_0) = |g_3_0_1  s_0_0_0)_{t} + x_cd * |f_2_0_1  s_0_0_0)_{t}
            const double K_2_0_1_1_0_0 = KET_g_s[ibra * 15 + 2] + ( kCD_x * KET_f_s[ibra * 10 + 2] );

            // |f_2_0_1  p_0_1_0) = |g_2_1_1  s_0_0_0)_{t} + y_cd * |f_2_0_1  s_0_0_0)_{t}
            const double K_2_0_1_0_1_0 = KET_g_s[ibra * 15 + 4] + ( kCD_y * KET_f_s[ibra * 10 + 2] );

            // |f_2_0_1  p_0_0_1) = |g_2_0_2  s_0_0_0)_{t} + z_cd * |f_2_0_1  s_0_0_0)_{t}
            const double K_2_0_1_0_0_1 = KET_g_s[ibra * 15 + 5] + ( kCD_z * KET_f_s[ibra * 10 + 2] );

            // |f_1_2_0  p_1_0_0) = |g_2_2_0  s_0_0_0)_{t} + x_cd * |f_1_2_0  s_0_0_0)_{t}
            const double K_1_2_0_1_0_0 = KET_g_s[ibra * 15 + 3] + ( kCD_x * KET_f_s[ibra * 10 + 3] );

            // |f_1_2_0  p_0_1_0) = |g_1_3_0  s_0_0_0)_{t} + y_cd * |f_1_2_0  s_0_0_0)_{t}
            const double K_1_2_0_0_1_0 = KET_g_s[ibra * 15 + 6] + ( kCD_y * KET_f_s[ibra * 10 + 3] );

            // |f_1_2_0  p_0_0_1) = |g_1_2_1  s_0_0_0)_{t} + z_cd * |f_1_2_0  s_0_0_0)_{t}
            const double K_1_2_0_0_0_1 = KET_g_s[ibra * 15 + 7] + ( kCD_z * KET_f_s[ibra * 10 + 3] );

            // |f_1_1_1  p_1_0_0) = |g_2_1_1  s_0_0_0)_{t} + x_cd * |f_1_1_1  s_0_0_0)_{t}
            const double K_1_1_1_1_0_0 = KET_g_s[ibra * 15 + 4] + ( kCD_x * KET_f_s[ibra * 10 + 4] );

            // |f_1_1_1  p_0_1_0) = |g_1_2_1  s_0_0_0)_{t} + y_cd * |f_1_1_1  s_0_0_0)_{t}
            const double K_1_1_1_0_1_0 = KET_g_s[ibra * 15 + 7] + ( kCD_y * KET_f_s[ibra * 10 + 4] );

            // |f_1_1_1  p_0_0_1) = |g_1_1_2  s_0_0_0)_{t} + z_cd * |f_1_1_1  s_0_0_0)_{t}
            const double K_1_1_1_0_0_1 = KET_g_s[ibra * 15 + 8] + ( kCD_z * KET_f_s[ibra * 10 + 4] );

            // |f_1_0_2  p_1_0_0) = |g_2_0_2  s_0_0_0)_{t} + x_cd * |f_1_0_2  s_0_0_0)_{t}
            const double K_1_0_2_1_0_0 = KET_g_s[ibra * 15 + 5] + ( kCD_x * KET_f_s[ibra * 10 + 5] );

            // |f_1_0_2  p_0_1_0) = |g_1_1_2  s_0_0_0)_{t} + y_cd * |f_1_0_2  s_0_0_0)_{t}
            const double K_1_0_2_0_1_0 = KET_g_s[ibra * 15 + 8] + ( kCD_y * KET_f_s[ibra * 10 + 5] );

            // |f_1_0_2  p_0_0_1) = |g_1_0_3  s_0_0_0)_{t} + z_cd * |f_1_0_2  s_0_0_0)_{t}
            const double K_1_0_2_0_0_1 = KET_g_s[ibra * 15 + 9] + ( kCD_z * KET_f_s[ibra * 10 + 5] );

            // |f_0_3_0  p_1_0_0) = |g_1_3_0  s_0_0_0)_{t} + x_cd * |f_0_3_0  s_0_0_0)_{t}
            const double K_0_3_0_1_0_0 = KET_g_s[ibra * 15 + 6] + ( kCD_x * KET_f_s[ibra * 10 + 6] );

            // |f_0_3_0  p_0_1_0) = |g_0_4_0  s_0_0_0)_{t} + y_cd * |f_0_3_0  s_0_0_0)_{t}
            const double K_0_3_0_0_1_0 = KET_g_s[ibra * 15 + 10] + ( kCD_y * KET_f_s[ibra * 10 + 6] );

            // |f_0_3_0  p_0_0_1) = |g_0_3_1  s_0_0_0)_{t} + z_cd * |f_0_3_0  s_0_0_0)_{t}
            const double K_0_3_0_0_0_1 = KET_g_s[ibra * 15 + 11] + ( kCD_z * KET_f_s[ibra * 10 + 6] );

            // |f_0_2_1  p_1_0_0) = |g_1_2_1  s_0_0_0)_{t} + x_cd * |f_0_2_1  s_0_0_0)_{t}
            const double K_0_2_1_1_0_0 = KET_g_s[ibra * 15 + 7] + ( kCD_x * KET_f_s[ibra * 10 + 7] );

            // |f_0_2_1  p_0_1_0) = |g_0_3_1  s_0_0_0)_{t} + y_cd * |f_0_2_1  s_0_0_0)_{t}
            const double K_0_2_1_0_1_0 = KET_g_s[ibra * 15 + 11] + ( kCD_y * KET_f_s[ibra * 10 + 7] );

            // |f_0_2_1  p_0_0_1) = |g_0_2_2  s_0_0_0)_{t} + z_cd * |f_0_2_1  s_0_0_0)_{t}
            const double K_0_2_1_0_0_1 = KET_g_s[ibra * 15 + 12] + ( kCD_z * KET_f_s[ibra * 10 + 7] );

            // |f_0_1_2  p_1_0_0) = |g_1_1_2  s_0_0_0)_{t} + x_cd * |f_0_1_2  s_0_0_0)_{t}
            const double K_0_1_2_1_0_0 = KET_g_s[ibra * 15 + 8] + ( kCD_x * KET_f_s[ibra * 10 + 8] );

            // |f_0_1_2  p_0_1_0) = |g_0_2_2  s_0_0_0)_{t} + y_cd * |f_0_1_2  s_0_0_0)_{t}
            const double K_0_1_2_0_1_0 = KET_g_s[ibra * 15 + 12] + ( kCD_y * KET_f_s[ibra * 10 + 8] );

            // |f_0_1_2  p_0_0_1) = |g_0_1_3  s_0_0_0)_{t} + z_cd * |f_0_1_2  s_0_0_0)_{t}
            const double K_0_1_2_0_0_1 = KET_g_s[ibra * 15 + 13] + ( kCD_z * KET_f_s[ibra * 10 + 8] );

            // |f_0_0_3  p_1_0_0) = |g_1_0_3  s_0_0_0)_{t} + x_cd * |f_0_0_3  s_0_0_0)_{t}
            const double K_0_0_3_1_0_0 = KET_g_s[ibra * 15 + 9] + ( kCD_x * KET_f_s[ibra * 10 + 9] );

            // |f_0_0_3  p_0_1_0) = |g_0_1_3  s_0_0_0)_{t} + y_cd * |f_0_0_3  s_0_0_0)_{t}
            const double K_0_0_3_0_1_0 = KET_g_s[ibra * 15 + 13] + ( kCD_y * KET_f_s[ibra * 10 + 9] );

            // |f_0_0_3  p_0_0_1) = |g_0_0_4  s_0_0_0)_{t} + z_cd * |f_0_0_3  s_0_0_0)_{t}
            const double K_0_0_3_0_0_1 = KET_g_s[ibra * 15 + 14] + ( kCD_z * KET_f_s[ibra * 10 + 9] );

            // |g_4_0_0  p_1_0_0) = |h_5_0_0  s_0_0_0)_{t} + x_cd * |g_4_0_0  s_0_0_0)_{t}
            const double K_4_0_0_1_0_0 = KET_h_s[ibra * 21 + 0] + ( kCD_x * KET_g_s[ibra * 15 + 0] );

            // |g_3_1_0  p_1_0_0) = |h_4_1_0  s_0_0_0)_{t} + x_cd * |g_3_1_0  s_0_0_0)_{t}
            const double K_3_1_0_1_0_0 = KET_h_s[ibra * 21 + 1] + ( kCD_x * KET_g_s[ibra * 15 + 1] );

            // |g_3_1_0  p_0_1_0) = |h_3_2_0  s_0_0_0)_{t} + y_cd * |g_3_1_0  s_0_0_0)_{t}
            const double K_3_1_0_0_1_0 = KET_h_s[ibra * 21 + 3] + ( kCD_y * KET_g_s[ibra * 15 + 1] );

            // |g_3_0_1  p_1_0_0) = |h_4_0_1  s_0_0_0)_{t} + x_cd * |g_3_0_1  s_0_0_0)_{t}
            const double K_3_0_1_1_0_0 = KET_h_s[ibra * 21 + 2] + ( kCD_x * KET_g_s[ibra * 15 + 2] );

            // |g_3_0_1  p_0_1_0) = |h_3_1_1  s_0_0_0)_{t} + y_cd * |g_3_0_1  s_0_0_0)_{t}
            const double K_3_0_1_0_1_0 = KET_h_s[ibra * 21 + 4] + ( kCD_y * KET_g_s[ibra * 15 + 2] );

            // |g_3_0_1  p_0_0_1) = |h_3_0_2  s_0_0_0)_{t} + z_cd * |g_3_0_1  s_0_0_0)_{t}
            const double K_3_0_1_0_0_1 = KET_h_s[ibra * 21 + 5] + ( kCD_z * KET_g_s[ibra * 15 + 2] );

            // |g_2_2_0  p_1_0_0) = |h_3_2_0  s_0_0_0)_{t} + x_cd * |g_2_2_0  s_0_0_0)_{t}
            const double K_2_2_0_1_0_0 = KET_h_s[ibra * 21 + 3] + ( kCD_x * KET_g_s[ibra * 15 + 3] );

            // |g_2_2_0  p_0_1_0) = |h_2_3_0  s_0_0_0)_{t} + y_cd * |g_2_2_0  s_0_0_0)_{t}
            const double K_2_2_0_0_1_0 = KET_h_s[ibra * 21 + 6] + ( kCD_y * KET_g_s[ibra * 15 + 3] );

            // |g_2_1_1  p_1_0_0) = |h_3_1_1  s_0_0_0)_{t} + x_cd * |g_2_1_1  s_0_0_0)_{t}
            const double K_2_1_1_1_0_0 = KET_h_s[ibra * 21 + 4] + ( kCD_x * KET_g_s[ibra * 15 + 4] );

            // |g_2_1_1  p_0_1_0) = |h_2_2_1  s_0_0_0)_{t} + y_cd * |g_2_1_1  s_0_0_0)_{t}
            const double K_2_1_1_0_1_0 = KET_h_s[ibra * 21 + 7] + ( kCD_y * KET_g_s[ibra * 15 + 4] );

            // |g_2_1_1  p_0_0_1) = |h_2_1_2  s_0_0_0)_{t} + z_cd * |g_2_1_1  s_0_0_0)_{t}
            const double K_2_1_1_0_0_1 = KET_h_s[ibra * 21 + 8] + ( kCD_z * KET_g_s[ibra * 15 + 4] );

            // |g_2_0_2  p_1_0_0) = |h_3_0_2  s_0_0_0)_{t} + x_cd * |g_2_0_2  s_0_0_0)_{t}
            const double K_2_0_2_1_0_0 = KET_h_s[ibra * 21 + 5] + ( kCD_x * KET_g_s[ibra * 15 + 5] );

            // |g_2_0_2  p_0_1_0) = |h_2_1_2  s_0_0_0)_{t} + y_cd * |g_2_0_2  s_0_0_0)_{t}
            const double K_2_0_2_0_1_0 = KET_h_s[ibra * 21 + 8] + ( kCD_y * KET_g_s[ibra * 15 + 5] );

            // |g_2_0_2  p_0_0_1) = |h_2_0_3  s_0_0_0)_{t} + z_cd * |g_2_0_2  s_0_0_0)_{t}
            const double K_2_0_2_0_0_1 = KET_h_s[ibra * 21 + 9] + ( kCD_z * KET_g_s[ibra * 15 + 5] );

            // |g_1_3_0  p_1_0_0) = |h_2_3_0  s_0_0_0)_{t} + x_cd * |g_1_3_0  s_0_0_0)_{t}
            const double K_1_3_0_1_0_0 = KET_h_s[ibra * 21 + 6] + ( kCD_x * KET_g_s[ibra * 15 + 6] );

            // |g_1_3_0  p_0_1_0) = |h_1_4_0  s_0_0_0)_{t} + y_cd * |g_1_3_0  s_0_0_0)_{t}
            const double K_1_3_0_0_1_0 = KET_h_s[ibra * 21 + 10] + ( kCD_y * KET_g_s[ibra * 15 + 6] );

            // |g_1_2_1  p_1_0_0) = |h_2_2_1  s_0_0_0)_{t} + x_cd * |g_1_2_1  s_0_0_0)_{t}
            const double K_1_2_1_1_0_0 = KET_h_s[ibra * 21 + 7] + ( kCD_x * KET_g_s[ibra * 15 + 7] );

            // |g_1_2_1  p_0_1_0) = |h_1_3_1  s_0_0_0)_{t} + y_cd * |g_1_2_1  s_0_0_0)_{t}
            const double K_1_2_1_0_1_0 = KET_h_s[ibra * 21 + 11] + ( kCD_y * KET_g_s[ibra * 15 + 7] );

            // |g_1_2_1  p_0_0_1) = |h_1_2_2  s_0_0_0)_{t} + z_cd * |g_1_2_1  s_0_0_0)_{t}
            const double K_1_2_1_0_0_1 = KET_h_s[ibra * 21 + 12] + ( kCD_z * KET_g_s[ibra * 15 + 7] );

            // |g_1_1_2  p_1_0_0) = |h_2_1_2  s_0_0_0)_{t} + x_cd * |g_1_1_2  s_0_0_0)_{t}
            const double K_1_1_2_1_0_0 = KET_h_s[ibra * 21 + 8] + ( kCD_x * KET_g_s[ibra * 15 + 8] );

            // |g_1_1_2  p_0_1_0) = |h_1_2_2  s_0_0_0)_{t} + y_cd * |g_1_1_2  s_0_0_0)_{t}
            const double K_1_1_2_0_1_0 = KET_h_s[ibra * 21 + 12] + ( kCD_y * KET_g_s[ibra * 15 + 8] );

            // |g_1_1_2  p_0_0_1) = |h_1_1_3  s_0_0_0)_{t} + z_cd * |g_1_1_2  s_0_0_0)_{t}
            const double K_1_1_2_0_0_1 = KET_h_s[ibra * 21 + 13] + ( kCD_z * KET_g_s[ibra * 15 + 8] );

            // |g_1_0_3  p_1_0_0) = |h_2_0_3  s_0_0_0)_{t} + x_cd * |g_1_0_3  s_0_0_0)_{t}
            const double K_1_0_3_1_0_0 = KET_h_s[ibra * 21 + 9] + ( kCD_x * KET_g_s[ibra * 15 + 9] );

            // |g_1_0_3  p_0_1_0) = |h_1_1_3  s_0_0_0)_{t} + y_cd * |g_1_0_3  s_0_0_0)_{t}
            const double K_1_0_3_0_1_0 = KET_h_s[ibra * 21 + 13] + ( kCD_y * KET_g_s[ibra * 15 + 9] );

            // |g_1_0_3  p_0_0_1) = |h_1_0_4  s_0_0_0)_{t} + z_cd * |g_1_0_3  s_0_0_0)_{t}
            const double K_1_0_3_0_0_1 = KET_h_s[ibra * 21 + 14] + ( kCD_z * KET_g_s[ibra * 15 + 9] );

            // |g_0_4_0  p_1_0_0) = |h_1_4_0  s_0_0_0)_{t} + x_cd * |g_0_4_0  s_0_0_0)_{t}
            const double K_0_4_0_1_0_0 = KET_h_s[ibra * 21 + 10] + ( kCD_x * KET_g_s[ibra * 15 + 10] );

            // |g_0_4_0  p_0_1_0) = |h_0_5_0  s_0_0_0)_{t} + y_cd * |g_0_4_0  s_0_0_0)_{t}
            const double K_0_4_0_0_1_0 = KET_h_s[ibra * 21 + 15] + ( kCD_y * KET_g_s[ibra * 15 + 10] );

            // |g_0_3_1  p_1_0_0) = |h_1_3_1  s_0_0_0)_{t} + x_cd * |g_0_3_1  s_0_0_0)_{t}
            const double K_0_3_1_1_0_0 = KET_h_s[ibra * 21 + 11] + ( kCD_x * KET_g_s[ibra * 15 + 11] );

            // |g_0_3_1  p_0_1_0) = |h_0_4_1  s_0_0_0)_{t} + y_cd * |g_0_3_1  s_0_0_0)_{t}
            const double K_0_3_1_0_1_0 = KET_h_s[ibra * 21 + 16] + ( kCD_y * KET_g_s[ibra * 15 + 11] );

            // |g_0_3_1  p_0_0_1) = |h_0_3_2  s_0_0_0)_{t} + z_cd * |g_0_3_1  s_0_0_0)_{t}
            const double K_0_3_1_0_0_1 = KET_h_s[ibra * 21 + 17] + ( kCD_z * KET_g_s[ibra * 15 + 11] );

            // |g_0_2_2  p_1_0_0) = |h_1_2_2  s_0_0_0)_{t} + x_cd * |g_0_2_2  s_0_0_0)_{t}
            const double K_0_2_2_1_0_0 = KET_h_s[ibra * 21 + 12] + ( kCD_x * KET_g_s[ibra * 15 + 12] );

            // |g_0_2_2  p_0_1_0) = |h_0_3_2  s_0_0_0)_{t} + y_cd * |g_0_2_2  s_0_0_0)_{t}
            const double K_0_2_2_0_1_0 = KET_h_s[ibra * 21 + 17] + ( kCD_y * KET_g_s[ibra * 15 + 12] );

            // |g_0_2_2  p_0_0_1) = |h_0_2_3  s_0_0_0)_{t} + z_cd * |g_0_2_2  s_0_0_0)_{t}
            const double K_0_2_2_0_0_1 = KET_h_s[ibra * 21 + 18] + ( kCD_z * KET_g_s[ibra * 15 + 12] );

            // |g_0_1_3  p_1_0_0) = |h_1_1_3  s_0_0_0)_{t} + x_cd * |g_0_1_3  s_0_0_0)_{t}
            const double K_0_1_3_1_0_0 = KET_h_s[ibra * 21 + 13] + ( kCD_x * KET_g_s[ibra * 15 + 13] );

            // |g_0_1_3  p_0_1_0) = |h_0_2_3  s_0_0_0)_{t} + y_cd * |g_0_1_3  s_0_0_0)_{t}
            const double K_0_1_3_0_1_0 = KET_h_s[ibra * 21 + 18] + ( kCD_y * KET_g_s[ibra * 15 + 13] );

            // |g_0_1_3  p_0_0_1) = |h_0_1_4  s_0_0_0)_{t} + z_cd * |g_0_1_3  s_0_0_0)_{t}
            const double K_0_1_3_0_0_1 = KET_h_s[ibra * 21 + 19] + ( kCD_z * KET_g_s[ibra * 15 + 13] );

            // |g_0_0_4  p_1_0_0) = |h_1_0_4  s_0_0_0)_{t} + x_cd * |g_0_0_4  s_0_0_0)_{t}
            const double K_0_0_4_1_0_0 = KET_h_s[ibra * 21 + 14] + ( kCD_x * KET_g_s[ibra * 15 + 14] );

            // |g_0_0_4  p_0_1_0) = |h_0_1_4  s_0_0_0)_{t} + y_cd * |g_0_0_4  s_0_0_0)_{t}
            const double K_0_0_4_0_1_0 = KET_h_s[ibra * 21 + 19] + ( kCD_y * KET_g_s[ibra * 15 + 14] );

            // |g_0_0_4  p_0_0_1) = |h_0_0_5  s_0_0_0)_{t} + z_cd * |g_0_0_4  s_0_0_0)_{t}
            const double K_0_0_4_0_0_1 = KET_h_s[ibra * 21 + 20] + ( kCD_z * KET_g_s[ibra * 15 + 14] );

            // |f_3_0_0  d_2_0_0)_{i} = |g_4_0_0  p_1_0_0) + x_cd * |f_3_0_0  p_1_0_0)
            KET_f_d[ibra * 60 + 0] = K_4_0_0_1_0_0 + ( kCD_x * K_3_0_0_1_0_0 );

            // |f_3_0_0  d_1_1_0)_{i} = |g_3_1_0  p_1_0_0) + y_cd * |f_3_0_0  p_1_0_0)
            KET_f_d[ibra * 60 + 1] = K_3_1_0_1_0_0 + ( kCD_y * K_3_0_0_1_0_0 );

            // |f_3_0_0  d_1_0_1)_{i} = |g_3_0_1  p_1_0_0) + z_cd * |f_3_0_0  p_1_0_0)
            KET_f_d[ibra * 60 + 2] = K_3_0_1_1_0_0 + ( kCD_z * K_3_0_0_1_0_0 );

            // |f_3_0_0  d_0_2_0)_{i} = |g_3_1_0  p_0_1_0) + y_cd * |f_3_0_0  p_0_1_0)
            KET_f_d[ibra * 60 + 3] = K_3_1_0_0_1_0 + ( kCD_y * K_3_0_0_0_1_0 );

            // |f_3_0_0  d_0_1_1)_{i} = |g_3_0_1  p_0_1_0) + z_cd * |f_3_0_0  p_0_1_0)
            KET_f_d[ibra * 60 + 4] = K_3_0_1_0_1_0 + ( kCD_z * K_3_0_0_0_1_0 );

            // |f_3_0_0  d_0_0_2)_{i} = |g_3_0_1  p_0_0_1) + z_cd * |f_3_0_0  p_0_0_1)
            KET_f_d[ibra * 60 + 5] = K_3_0_1_0_0_1 + ( kCD_z * K_3_0_0_0_0_1 );

            // |f_2_1_0  d_2_0_0)_{i} = |g_3_1_0  p_1_0_0) + x_cd * |f_2_1_0  p_1_0_0)
            KET_f_d[ibra * 60 + 6] = K_3_1_0_1_0_0 + ( kCD_x * K_2_1_0_1_0_0 );

            // |f_2_1_0  d_1_1_0)_{i} = |g_2_2_0  p_1_0_0) + y_cd * |f_2_1_0  p_1_0_0)
            KET_f_d[ibra * 60 + 7] = K_2_2_0_1_0_0 + ( kCD_y * K_2_1_0_1_0_0 );

            // |f_2_1_0  d_1_0_1)_{i} = |g_2_1_1  p_1_0_0) + z_cd * |f_2_1_0  p_1_0_0)
            KET_f_d[ibra * 60 + 8] = K_2_1_1_1_0_0 + ( kCD_z * K_2_1_0_1_0_0 );

            // |f_2_1_0  d_0_2_0)_{i} = |g_2_2_0  p_0_1_0) + y_cd * |f_2_1_0  p_0_1_0)
            KET_f_d[ibra * 60 + 9] = K_2_2_0_0_1_0 + ( kCD_y * K_2_1_0_0_1_0 );

            // |f_2_1_0  d_0_1_1)_{i} = |g_2_1_1  p_0_1_0) + z_cd * |f_2_1_0  p_0_1_0)
            KET_f_d[ibra * 60 + 10] = K_2_1_1_0_1_0 + ( kCD_z * K_2_1_0_0_1_0 );

            // |f_2_1_0  d_0_0_2)_{i} = |g_2_1_1  p_0_0_1) + z_cd * |f_2_1_0  p_0_0_1)
            KET_f_d[ibra * 60 + 11] = K_2_1_1_0_0_1 + ( kCD_z * K_2_1_0_0_0_1 );

            // |f_2_0_1  d_2_0_0)_{i} = |g_3_0_1  p_1_0_0) + x_cd * |f_2_0_1  p_1_0_0)
            KET_f_d[ibra * 60 + 12] = K_3_0_1_1_0_0 + ( kCD_x * K_2_0_1_1_0_0 );

            // |f_2_0_1  d_1_1_0)_{i} = |g_2_1_1  p_1_0_0) + y_cd * |f_2_0_1  p_1_0_0)
            KET_f_d[ibra * 60 + 13] = K_2_1_1_1_0_0 + ( kCD_y * K_2_0_1_1_0_0 );

            // |f_2_0_1  d_1_0_1)_{i} = |g_2_0_2  p_1_0_0) + z_cd * |f_2_0_1  p_1_0_0)
            KET_f_d[ibra * 60 + 14] = K_2_0_2_1_0_0 + ( kCD_z * K_2_0_1_1_0_0 );

            // |f_2_0_1  d_0_2_0)_{i} = |g_2_1_1  p_0_1_0) + y_cd * |f_2_0_1  p_0_1_0)
            KET_f_d[ibra * 60 + 15] = K_2_1_1_0_1_0 + ( kCD_y * K_2_0_1_0_1_0 );

            // |f_2_0_1  d_0_1_1)_{i} = |g_2_0_2  p_0_1_0) + z_cd * |f_2_0_1  p_0_1_0)
            KET_f_d[ibra * 60 + 16] = K_2_0_2_0_1_0 + ( kCD_z * K_2_0_1_0_1_0 );

            // |f_2_0_1  d_0_0_2)_{i} = |g_2_0_2  p_0_0_1) + z_cd * |f_2_0_1  p_0_0_1)
            KET_f_d[ibra * 60 + 17] = K_2_0_2_0_0_1 + ( kCD_z * K_2_0_1_0_0_1 );

            // |f_1_2_0  d_2_0_0)_{i} = |g_2_2_0  p_1_0_0) + x_cd * |f_1_2_0  p_1_0_0)
            KET_f_d[ibra * 60 + 18] = K_2_2_0_1_0_0 + ( kCD_x * K_1_2_0_1_0_0 );

            // |f_1_2_0  d_1_1_0)_{i} = |g_1_3_0  p_1_0_0) + y_cd * |f_1_2_0  p_1_0_0)
            KET_f_d[ibra * 60 + 19] = K_1_3_0_1_0_0 + ( kCD_y * K_1_2_0_1_0_0 );

            // |f_1_2_0  d_1_0_1)_{i} = |g_1_2_1  p_1_0_0) + z_cd * |f_1_2_0  p_1_0_0)
            KET_f_d[ibra * 60 + 20] = K_1_2_1_1_0_0 + ( kCD_z * K_1_2_0_1_0_0 );

            // |f_1_2_0  d_0_2_0)_{i} = |g_1_3_0  p_0_1_0) + y_cd * |f_1_2_0  p_0_1_0)
            KET_f_d[ibra * 60 + 21] = K_1_3_0_0_1_0 + ( kCD_y * K_1_2_0_0_1_0 );

            // |f_1_2_0  d_0_1_1)_{i} = |g_1_2_1  p_0_1_0) + z_cd * |f_1_2_0  p_0_1_0)
            KET_f_d[ibra * 60 + 22] = K_1_2_1_0_1_0 + ( kCD_z * K_1_2_0_0_1_0 );

            // |f_1_2_0  d_0_0_2)_{i} = |g_1_2_1  p_0_0_1) + z_cd * |f_1_2_0  p_0_0_1)
            KET_f_d[ibra * 60 + 23] = K_1_2_1_0_0_1 + ( kCD_z * K_1_2_0_0_0_1 );

            // |f_1_1_1  d_2_0_0)_{i} = |g_2_1_1  p_1_0_0) + x_cd * |f_1_1_1  p_1_0_0)
            KET_f_d[ibra * 60 + 24] = K_2_1_1_1_0_0 + ( kCD_x * K_1_1_1_1_0_0 );

            // |f_1_1_1  d_1_1_0)_{i} = |g_1_2_1  p_1_0_0) + y_cd * |f_1_1_1  p_1_0_0)
            KET_f_d[ibra * 60 + 25] = K_1_2_1_1_0_0 + ( kCD_y * K_1_1_1_1_0_0 );

            // |f_1_1_1  d_1_0_1)_{i} = |g_1_1_2  p_1_0_0) + z_cd * |f_1_1_1  p_1_0_0)
            KET_f_d[ibra * 60 + 26] = K_1_1_2_1_0_0 + ( kCD_z * K_1_1_1_1_0_0 );

            // |f_1_1_1  d_0_2_0)_{i} = |g_1_2_1  p_0_1_0) + y_cd * |f_1_1_1  p_0_1_0)
            KET_f_d[ibra * 60 + 27] = K_1_2_1_0_1_0 + ( kCD_y * K_1_1_1_0_1_0 );

            // |f_1_1_1  d_0_1_1)_{i} = |g_1_1_2  p_0_1_0) + z_cd * |f_1_1_1  p_0_1_0)
            KET_f_d[ibra * 60 + 28] = K_1_1_2_0_1_0 + ( kCD_z * K_1_1_1_0_1_0 );

            // |f_1_1_1  d_0_0_2)_{i} = |g_1_1_2  p_0_0_1) + z_cd * |f_1_1_1  p_0_0_1)
            KET_f_d[ibra * 60 + 29] = K_1_1_2_0_0_1 + ( kCD_z * K_1_1_1_0_0_1 );

            // |f_1_0_2  d_2_0_0)_{i} = |g_2_0_2  p_1_0_0) + x_cd * |f_1_0_2  p_1_0_0)
            KET_f_d[ibra * 60 + 30] = K_2_0_2_1_0_0 + ( kCD_x * K_1_0_2_1_0_0 );

            // |f_1_0_2  d_1_1_0)_{i} = |g_1_1_2  p_1_0_0) + y_cd * |f_1_0_2  p_1_0_0)
            KET_f_d[ibra * 60 + 31] = K_1_1_2_1_0_0 + ( kCD_y * K_1_0_2_1_0_0 );

            // |f_1_0_2  d_1_0_1)_{i} = |g_1_0_3  p_1_0_0) + z_cd * |f_1_0_2  p_1_0_0)
            KET_f_d[ibra * 60 + 32] = K_1_0_3_1_0_0 + ( kCD_z * K_1_0_2_1_0_0 );

            // |f_1_0_2  d_0_2_0)_{i} = |g_1_1_2  p_0_1_0) + y_cd * |f_1_0_2  p_0_1_0)
            KET_f_d[ibra * 60 + 33] = K_1_1_2_0_1_0 + ( kCD_y * K_1_0_2_0_1_0 );

            // |f_1_0_2  d_0_1_1)_{i} = |g_1_0_3  p_0_1_0) + z_cd * |f_1_0_2  p_0_1_0)
            KET_f_d[ibra * 60 + 34] = K_1_0_3_0_1_0 + ( kCD_z * K_1_0_2_0_1_0 );

            // |f_1_0_2  d_0_0_2)_{i} = |g_1_0_3  p_0_0_1) + z_cd * |f_1_0_2  p_0_0_1)
            KET_f_d[ibra * 60 + 35] = K_1_0_3_0_0_1 + ( kCD_z * K_1_0_2_0_0_1 );

            // |f_0_3_0  d_2_0_0)_{i} = |g_1_3_0  p_1_0_0) + x_cd * |f_0_3_0  p_1_0_0)
            KET_f_d[ibra * 60 + 36] = K_1_3_0_1_0_0 + ( kCD_x * K_0_3_0_1_0_0 );

            // |f_0_3_0  d_1_1_0)_{i} = |g_0_4_0  p_1_0_0) + y_cd * |f_0_3_0  p_1_0_0)
            KET_f_d[ibra * 60 + 37] = K_0_4_0_1_0_0 + ( kCD_y * K_0_3_0_1_0_0 );

            // |f_0_3_0  d_1_0_1)_{i} = |g_0_3_1  p_1_0_0) + z_cd * |f_0_3_0  p_1_0_0)
            KET_f_d[ibra * 60 + 38] = K_0_3_1_1_0_0 + ( kCD_z * K_0_3_0_1_0_0 );

            // |f_0_3_0  d_0_2_0)_{i} = |g_0_4_0  p_0_1_0) + y_cd * |f_0_3_0  p_0_1_0)
            KET_f_d[ibra * 60 + 39] = K_0_4_0_0_1_0 + ( kCD_y * K_0_3_0_0_1_0 );

            // |f_0_3_0  d_0_1_1)_{i} = |g_0_3_1  p_0_1_0) + z_cd * |f_0_3_0  p_0_1_0)
            KET_f_d[ibra * 60 + 40] = K_0_3_1_0_1_0 + ( kCD_z * K_0_3_0_0_1_0 );

            // |f_0_3_0  d_0_0_2)_{i} = |g_0_3_1  p_0_0_1) + z_cd * |f_0_3_0  p_0_0_1)
            KET_f_d[ibra * 60 + 41] = K_0_3_1_0_0_1 + ( kCD_z * K_0_3_0_0_0_1 );

            // |f_0_2_1  d_2_0_0)_{i} = |g_1_2_1  p_1_0_0) + x_cd * |f_0_2_1  p_1_0_0)
            KET_f_d[ibra * 60 + 42] = K_1_2_1_1_0_0 + ( kCD_x * K_0_2_1_1_0_0 );

            // |f_0_2_1  d_1_1_0)_{i} = |g_0_3_1  p_1_0_0) + y_cd * |f_0_2_1  p_1_0_0)
            KET_f_d[ibra * 60 + 43] = K_0_3_1_1_0_0 + ( kCD_y * K_0_2_1_1_0_0 );

            // |f_0_2_1  d_1_0_1)_{i} = |g_0_2_2  p_1_0_0) + z_cd * |f_0_2_1  p_1_0_0)
            KET_f_d[ibra * 60 + 44] = K_0_2_2_1_0_0 + ( kCD_z * K_0_2_1_1_0_0 );

            // |f_0_2_1  d_0_2_0)_{i} = |g_0_3_1  p_0_1_0) + y_cd * |f_0_2_1  p_0_1_0)
            KET_f_d[ibra * 60 + 45] = K_0_3_1_0_1_0 + ( kCD_y * K_0_2_1_0_1_0 );

            // |f_0_2_1  d_0_1_1)_{i} = |g_0_2_2  p_0_1_0) + z_cd * |f_0_2_1  p_0_1_0)
            KET_f_d[ibra * 60 + 46] = K_0_2_2_0_1_0 + ( kCD_z * K_0_2_1_0_1_0 );

            // |f_0_2_1  d_0_0_2)_{i} = |g_0_2_2  p_0_0_1) + z_cd * |f_0_2_1  p_0_0_1)
            KET_f_d[ibra * 60 + 47] = K_0_2_2_0_0_1 + ( kCD_z * K_0_2_1_0_0_1 );

            // |f_0_1_2  d_2_0_0)_{i} = |g_1_1_2  p_1_0_0) + x_cd * |f_0_1_2  p_1_0_0)
            KET_f_d[ibra * 60 + 48] = K_1_1_2_1_0_0 + ( kCD_x * K_0_1_2_1_0_0 );

            // |f_0_1_2  d_1_1_0)_{i} = |g_0_2_2  p_1_0_0) + y_cd * |f_0_1_2  p_1_0_0)
            KET_f_d[ibra * 60 + 49] = K_0_2_2_1_0_0 + ( kCD_y * K_0_1_2_1_0_0 );

            // |f_0_1_2  d_1_0_1)_{i} = |g_0_1_3  p_1_0_0) + z_cd * |f_0_1_2  p_1_0_0)
            KET_f_d[ibra * 60 + 50] = K_0_1_3_1_0_0 + ( kCD_z * K_0_1_2_1_0_0 );

            // |f_0_1_2  d_0_2_0)_{i} = |g_0_2_2  p_0_1_0) + y_cd * |f_0_1_2  p_0_1_0)
            KET_f_d[ibra * 60 + 51] = K_0_2_2_0_1_0 + ( kCD_y * K_0_1_2_0_1_0 );

            // |f_0_1_2  d_0_1_1)_{i} = |g_0_1_3  p_0_1_0) + z_cd * |f_0_1_2  p_0_1_0)
            KET_f_d[ibra * 60 + 52] = K_0_1_3_0_1_0 + ( kCD_z * K_0_1_2_0_1_0 );

            // |f_0_1_2  d_0_0_2)_{i} = |g_0_1_3  p_0_0_1) + z_cd * |f_0_1_2  p_0_0_1)
            KET_f_d[ibra * 60 + 53] = K_0_1_3_0_0_1 + ( kCD_z * K_0_1_2_0_0_1 );

            // |f_0_0_3  d_2_0_0)_{i} = |g_1_0_3  p_1_0_0) + x_cd * |f_0_0_3  p_1_0_0)
            KET_f_d[ibra * 60 + 54] = K_1_0_3_1_0_0 + ( kCD_x * K_0_0_3_1_0_0 );

            // |f_0_0_3  d_1_1_0)_{i} = |g_0_1_3  p_1_0_0) + y_cd * |f_0_0_3  p_1_0_0)
            KET_f_d[ibra * 60 + 55] = K_0_1_3_1_0_0 + ( kCD_y * K_0_0_3_1_0_0 );

            // |f_0_0_3  d_1_0_1)_{i} = |g_0_0_4  p_1_0_0) + z_cd * |f_0_0_3  p_1_0_0)
            KET_f_d[ibra * 60 + 56] = K_0_0_4_1_0_0 + ( kCD_z * K_0_0_3_1_0_0 );

            // |f_0_0_3  d_0_2_0)_{i} = |g_0_1_3  p_0_1_0) + y_cd * |f_0_0_3  p_0_1_0)
            KET_f_d[ibra * 60 + 57] = K_0_1_3_0_1_0 + ( kCD_y * K_0_0_3_0_1_0 );

            // |f_0_0_3  d_0_1_1)_{i} = |g_0_0_4  p_0_1_0) + z_cd * |f_0_0_3  p_0_1_0)
            KET_f_d[ibra * 60 + 58] = K_0_0_4_0_1_0 + ( kCD_z * K_0_0_3_0_1_0 );

            // |f_0_0_3  d_0_0_2)_{i} = |g_0_0_4  p_0_0_1) + z_cd * |f_0_0_3  p_0_0_1)
            KET_f_d[ibra * 60 + 59] = K_0_0_4_0_0_1 + ( kCD_z * K_0_0_3_0_0_1 );

        }

}


    //////////////////////////////////////////////
    // KET: ( f f |
    // Steps: 319
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET__f_f(
                  double * const restrict KET_f_s,
                  double * const restrict KET_f_f,
                  double * const restrict KET_g_s,
                  double * const restrict KET_h_s,
                  double * const restrict KET_i_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 )
{
    int ibra;

        for(ibra = 0; ibra < ncart_bra; ++ibra)
        {
            // |f_3_0_0  p_1_0_0) = |g_4_0_0  s_0_0_0)_{t} + x_cd * |f_3_0_0  s_0_0_0)_{t}
            const double K_3_0_0_1_0_0 = KET_g_s[ibra * 15 + 0] + ( kCD_x * KET_f_s[ibra * 10 + 0] );

            // |f_3_0_0  p_0_1_0) = |g_3_1_0  s_0_0_0)_{t} + y_cd * |f_3_0_0  s_0_0_0)_{t}
            const double K_3_0_0_0_1_0 = KET_g_s[ibra * 15 + 1] + ( kCD_y * KET_f_s[ibra * 10 + 0] );

            // |f_3_0_0  p_0_0_1) = |g_3_0_1  s_0_0_0)_{t} + z_cd * |f_3_0_0  s_0_0_0)_{t}
            const double K_3_0_0_0_0_1 = KET_g_s[ibra * 15 + 2] + ( kCD_z * KET_f_s[ibra * 10 + 0] );

            // |f_2_1_0  p_1_0_0) = |g_3_1_0  s_0_0_0)_{t} + x_cd * |f_2_1_0  s_0_0_0)_{t}
            const double K_2_1_0_1_0_0 = KET_g_s[ibra * 15 + 1] + ( kCD_x * KET_f_s[ibra * 10 + 1] );

            // |f_2_1_0  p_0_1_0) = |g_2_2_0  s_0_0_0)_{t} + y_cd * |f_2_1_0  s_0_0_0)_{t}
            const double K_2_1_0_0_1_0 = KET_g_s[ibra * 15 + 3] + ( kCD_y * KET_f_s[ibra * 10 + 1] );

            // |f_2_1_0  p_0_0_1) = |g_2_1_1  s_0_0_0)_{t} + z_cd * |f_2_1_0  s_0_0_0)_{t}
            const double K_2_1_0_0_0_1 = KET_g_s[ibra * 15 + 4] + ( kCD_z * KET_f_s[ibra * 10 + 1] );

            // |f_2_0_1  p_1_0_0) = |g_3_0_1  s_0_0_0)_{t} + x_cd * |f_2_0_1  s_0_0_0)_{t}
            const double K_2_0_1_1_0_0 = KET_g_s[ibra * 15 + 2] + ( kCD_x * KET_f_s[ibra * 10 + 2] );

            // |f_2_0_1  p_0_1_0) = |g_2_1_1  s_0_0_0)_{t} + y_cd * |f_2_0_1  s_0_0_0)_{t}
            const double K_2_0_1_0_1_0 = KET_g_s[ibra * 15 + 4] + ( kCD_y * KET_f_s[ibra * 10 + 2] );

            // |f_2_0_1  p_0_0_1) = |g_2_0_2  s_0_0_0)_{t} + z_cd * |f_2_0_1  s_0_0_0)_{t}
            const double K_2_0_1_0_0_1 = KET_g_s[ibra * 15 + 5] + ( kCD_z * KET_f_s[ibra * 10 + 2] );

            // |f_1_2_0  p_1_0_0) = |g_2_2_0  s_0_0_0)_{t} + x_cd * |f_1_2_0  s_0_0_0)_{t}
            const double K_1_2_0_1_0_0 = KET_g_s[ibra * 15 + 3] + ( kCD_x * KET_f_s[ibra * 10 + 3] );

            // |f_1_2_0  p_0_1_0) = |g_1_3_0  s_0_0_0)_{t} + y_cd * |f_1_2_0  s_0_0_0)_{t}
            const double K_1_2_0_0_1_0 = KET_g_s[ibra * 15 + 6] + ( kCD_y * KET_f_s[ibra * 10 + 3] );

            // |f_1_2_0  p_0_0_1) = |g_1_2_1  s_0_0_0)_{t} + z_cd * |f_1_2_0  s_0_0_0)_{t}
            const double K_1_2_0_0_0_1 = KET_g_s[ibra * 15 + 7] + ( kCD_z * KET_f_s[ibra * 10 + 3] );

            // |f_1_1_1  p_1_0_0) = |g_2_1_1  s_0_0_0)_{t} + x_cd * |f_1_1_1  s_0_0_0)_{t}
            const double K_1_1_1_1_0_0 = KET_g_s[ibra * 15 + 4] + ( kCD_x * KET_f_s[ibra * 10 + 4] );

            // |f_1_1_1  p_0_1_0) = |g_1_2_1  s_0_0_0)_{t} + y_cd * |f_1_1_1  s_0_0_0)_{t}
            const double K_1_1_1_0_1_0 = KET_g_s[ibra * 15 + 7] + ( kCD_y * KET_f_s[ibra * 10 + 4] );

            // |f_1_1_1  p_0_0_1) = |g_1_1_2  s_0_0_0)_{t} + z_cd * |f_1_1_1  s_0_0_0)_{t}
            const double K_1_1_1_0_0_1 = KET_g_s[ibra * 15 + 8] + ( kCD_z * KET_f_s[ibra * 10 + 4] );

            // |f_1_0_2  p_1_0_0) = |g_2_0_2  s_0_0_0)_{t} + x_cd * |f_1_0_2  s_0_0_0)_{t}
            const double K_1_0_2_1_0_0 = KET_g_s[ibra * 15 + 5] + ( kCD_x * KET_f_s[ibra * 10 + 5] );

            // |f_1_0_2  p_0_1_0) = |g_1_1_2  s_0_0_0)_{t} + y_cd * |f_1_0_2  s_0_0_0)_{t}
            const double K_1_0_2_0_1_0 = KET_g_s[ibra * 15 + 8] + ( kCD_y * KET_f_s[ibra * 10 + 5] );

            // |f_1_0_2  p_0_0_1) = |g_1_0_3  s_0_0_0)_{t} + z_cd * |f_1_0_2  s_0_0_0)_{t}
            const double K_1_0_2_0_0_1 = KET_g_s[ibra * 15 + 9] + ( kCD_z * KET_f_s[ibra * 10 + 5] );

            // |f_0_3_0  p_1_0_0) = |g_1_3_0  s_0_0_0)_{t} + x_cd * |f_0_3_0  s_0_0_0)_{t}
            const double K_0_3_0_1_0_0 = KET_g_s[ibra * 15 + 6] + ( kCD_x * KET_f_s[ibra * 10 + 6] );

            // |f_0_3_0  p_0_1_0) = |g_0_4_0  s_0_0_0)_{t} + y_cd * |f_0_3_0  s_0_0_0)_{t}
            const double K_0_3_0_0_1_0 = KET_g_s[ibra * 15 + 10] + ( kCD_y * KET_f_s[ibra * 10 + 6] );

            // |f_0_3_0  p_0_0_1) = |g_0_3_1  s_0_0_0)_{t} + z_cd * |f_0_3_0  s_0_0_0)_{t}
            const double K_0_3_0_0_0_1 = KET_g_s[ibra * 15 + 11] + ( kCD_z * KET_f_s[ibra * 10 + 6] );

            // |f_0_2_1  p_1_0_0) = |g_1_2_1  s_0_0_0)_{t} + x_cd * |f_0_2_1  s_0_0_0)_{t}
            const double K_0_2_1_1_0_0 = KET_g_s[ibra * 15 + 7] + ( kCD_x * KET_f_s[ibra * 10 + 7] );

            // |f_0_2_1  p_0_1_0) = |g_0_3_1  s_0_0_0)_{t} + y_cd * |f_0_2_1  s_0_0_0)_{t}
            const double K_0_2_1_0_1_0 = KET_g_s[ibra * 15 + 11] + ( kCD_y * KET_f_s[ibra * 10 + 7] );

            // |f_0_2_1  p_0_0_1) = |g_0_2_2  s_0_0_0)_{t} + z_cd * |f_0_2_1  s_0_0_0)_{t}
            const double K_0_2_1_0_0_1 = KET_g_s[ibra * 15 + 12] + ( kCD_z * KET_f_s[ibra * 10 + 7] );

            // |f_0_1_2  p_1_0_0) = |g_1_1_2  s_0_0_0)_{t} + x_cd * |f_0_1_2  s_0_0_0)_{t}
            const double K_0_1_2_1_0_0 = KET_g_s[ibra * 15 + 8] + ( kCD_x * KET_f_s[ibra * 10 + 8] );

            // |f_0_1_2  p_0_1_0) = |g_0_2_2  s_0_0_0)_{t} + y_cd * |f_0_1_2  s_0_0_0)_{t}
            const double K_0_1_2_0_1_0 = KET_g_s[ibra * 15 + 12] + ( kCD_y * KET_f_s[ibra * 10 + 8] );

            // |f_0_1_2  p_0_0_1) = |g_0_1_3  s_0_0_0)_{t} + z_cd * |f_0_1_2  s_0_0_0)_{t}
            const double K_0_1_2_0_0_1 = KET_g_s[ibra * 15 + 13] + ( kCD_z * KET_f_s[ibra * 10 + 8] );

            // |f_0_0_3  p_1_0_0) = |g_1_0_3  s_0_0_0)_{t} + x_cd * |f_0_0_3  s_0_0_0)_{t}
            const double K_0_0_3_1_0_0 = KET_g_s[ibra * 15 + 9] + ( kCD_x * KET_f_s[ibra * 10 + 9] );

            // |f_0_0_3  p_0_1_0) = |g_0_1_3  s_0_0_0)_{t} + y_cd * |f_0_0_3  s_0_0_0)_{t}
            const double K_0_0_3_0_1_0 = KET_g_s[ibra * 15 + 13] + ( kCD_y * KET_f_s[ibra * 10 + 9] );

            // |f_0_0_3  p_0_0_1) = |g_0_0_4  s_0_0_0)_{t} + z_cd * |f_0_0_3  s_0_0_0)_{t}
            const double K_0_0_3_0_0_1 = KET_g_s[ibra * 15 + 14] + ( kCD_z * KET_f_s[ibra * 10 + 9] );

            // |g_4_0_0  p_1_0_0) = |h_5_0_0  s_0_0_0)_{t} + x_cd * |g_4_0_0  s_0_0_0)_{t}
            const double K_4_0_0_1_0_0 = KET_h_s[ibra * 21 + 0] + ( kCD_x * KET_g_s[ibra * 15 + 0] );

            // |g_4_0_0  p_0_1_0) = |h_4_1_0  s_0_0_0)_{t} + y_cd * |g_4_0_0  s_0_0_0)_{t}
            const double K_4_0_0_0_1_0 = KET_h_s[ibra * 21 + 1] + ( kCD_y * KET_g_s[ibra * 15 + 0] );

            // |g_4_0_0  p_0_0_1) = |h_4_0_1  s_0_0_0)_{t} + z_cd * |g_4_0_0  s_0_0_0)_{t}
            const double K_4_0_0_0_0_1 = KET_h_s[ibra * 21 + 2] + ( kCD_z * KET_g_s[ibra * 15 + 0] );

            // |g_3_1_0  p_1_0_0) = |h_4_1_0  s_0_0_0)_{t} + x_cd * |g_3_1_0  s_0_0_0)_{t}
            const double K_3_1_0_1_0_0 = KET_h_s[ibra * 21 + 1] + ( kCD_x * KET_g_s[ibra * 15 + 1] );

            // |g_3_1_0  p_0_1_0) = |h_3_2_0  s_0_0_0)_{t} + y_cd * |g_3_1_0  s_0_0_0)_{t}
            const double K_3_1_0_0_1_0 = KET_h_s[ibra * 21 + 3] + ( kCD_y * KET_g_s[ibra * 15 + 1] );

            // |g_3_1_0  p_0_0_1) = |h_3_1_1  s_0_0_0)_{t} + z_cd * |g_3_1_0  s_0_0_0)_{t}
            const double K_3_1_0_0_0_1 = KET_h_s[ibra * 21 + 4] + ( kCD_z * KET_g_s[ibra * 15 + 1] );

            // |g_3_0_1  p_1_0_0) = |h_4_0_1  s_0_0_0)_{t} + x_cd * |g_3_0_1  s_0_0_0)_{t}
            const double K_3_0_1_1_0_0 = KET_h_s[ibra * 21 + 2] + ( kCD_x * KET_g_s[ibra * 15 + 2] );

            // |g_3_0_1  p_0_1_0) = |h_3_1_1  s_0_0_0)_{t} + y_cd * |g_3_0_1  s_0_0_0)_{t}
            const double K_3_0_1_0_1_0 = KET_h_s[ibra * 21 + 4] + ( kCD_y * KET_g_s[ibra * 15 + 2] );

            // |g_3_0_1  p_0_0_1) = |h_3_0_2  s_0_0_0)_{t} + z_cd * |g_3_0_1  s_0_0_0)_{t}
            const double K_3_0_1_0_0_1 = KET_h_s[ibra * 21 + 5] + ( kCD_z * KET_g_s[ibra * 15 + 2] );

            // |g_2_2_0  p_1_0_0) = |h_3_2_0  s_0_0_0)_{t} + x_cd * |g_2_2_0  s_0_0_0)_{t}
            const double K_2_2_0_1_0_0 = KET_h_s[ibra * 21 + 3] + ( kCD_x * KET_g_s[ibra * 15 + 3] );

            // |g_2_2_0  p_0_1_0) = |h_2_3_0  s_0_0_0)_{t} + y_cd * |g_2_2_0  s_0_0_0)_{t}
            const double K_2_2_0_0_1_0 = KET_h_s[ibra * 21 + 6] + ( kCD_y * KET_g_s[ibra * 15 + 3] );

            // |g_2_2_0  p_0_0_1) = |h_2_2_1  s_0_0_0)_{t} + z_cd * |g_2_2_0  s_0_0_0)_{t}
            const double K_2_2_0_0_0_1 = KET_h_s[ibra * 21 + 7] + ( kCD_z * KET_g_s[ibra * 15 + 3] );

            // |g_2_1_1  p_1_0_0) = |h_3_1_1  s_0_0_0)_{t} + x_cd * |g_2_1_1  s_0_0_0)_{t}
            const double K_2_1_1_1_0_0 = KET_h_s[ibra * 21 + 4] + ( kCD_x * KET_g_s[ibra * 15 + 4] );

            // |g_2_1_1  p_0_1_0) = |h_2_2_1  s_0_0_0)_{t} + y_cd * |g_2_1_1  s_0_0_0)_{t}
            const double K_2_1_1_0_1_0 = KET_h_s[ibra * 21 + 7] + ( kCD_y * KET_g_s[ibra * 15 + 4] );

            // |g_2_1_1  p_0_0_1) = |h_2_1_2  s_0_0_0)_{t} + z_cd * |g_2_1_1  s_0_0_0)_{t}
            const double K_2_1_1_0_0_1 = KET_h_s[ibra * 21 + 8] + ( kCD_z * KET_g_s[ibra * 15 + 4] );

            // |g_2_0_2  p_1_0_0) = |h_3_0_2  s_0_0_0)_{t} + x_cd * |g_2_0_2  s_0_0_0)_{t}
            const double K_2_0_2_1_0_0 = KET_h_s[ibra * 21 + 5] + ( kCD_x * KET_g_s[ibra * 15 + 5] );

            // |g_2_0_2  p_0_1_0) = |h_2_1_2  s_0_0_0)_{t} + y_cd * |g_2_0_2  s_0_0_0)_{t}
            const double K_2_0_2_0_1_0 = KET_h_s[ibra * 21 + 8] + ( kCD_y * KET_g_s[ibra * 15 + 5] );

            // |g_2_0_2  p_0_0_1) = |h_2_0_3  s_0_0_0)_{t} + z_cd * |g_2_0_2  s_0_0_0)_{t}
            const double K_2_0_2_0_0_1 = KET_h_s[ibra * 21 + 9] + ( kCD_z * KET_g_s[ibra * 15 + 5] );

            // |g_1_3_0  p_1_0_0) = |h_2_3_0  s_0_0_0)_{t} + x_cd * |g_1_3_0  s_0_0_0)_{t}
            const double K_1_3_0_1_0_0 = KET_h_s[ibra * 21 + 6] + ( kCD_x * KET_g_s[ibra * 15 + 6] );

            // |g_1_3_0  p_0_1_0) = |h_1_4_0  s_0_0_0)_{t} + y_cd * |g_1_3_0  s_0_0_0)_{t}
            const double K_1_3_0_0_1_0 = KET_h_s[ibra * 21 + 10] + ( kCD_y * KET_g_s[ibra * 15 + 6] );

            // |g_1_3_0  p_0_0_1) = |h_1_3_1  s_0_0_0)_{t} + z_cd * |g_1_3_0  s_0_0_0)_{t}
            const double K_1_3_0_0_0_1 = KET_h_s[ibra * 21 + 11] + ( kCD_z * KET_g_s[ibra * 15 + 6] );

            // |g_1_2_1  p_1_0_0) = |h_2_2_1  s_0_0_0)_{t} + x_cd * |g_1_2_1  s_0_0_0)_{t}
            const double K_1_2_1_1_0_0 = KET_h_s[ibra * 21 + 7] + ( kCD_x * KET_g_s[ibra * 15 + 7] );

            // |g_1_2_1  p_0_1_0) = |h_1_3_1  s_0_0_0)_{t} + y_cd * |g_1_2_1  s_0_0_0)_{t}
            const double K_1_2_1_0_1_0 = KET_h_s[ibra * 21 + 11] + ( kCD_y * KET_g_s[ibra * 15 + 7] );

            // |g_1_2_1  p_0_0_1) = |h_1_2_2  s_0_0_0)_{t} + z_cd * |g_1_2_1  s_0_0_0)_{t}
            const double K_1_2_1_0_0_1 = KET_h_s[ibra * 21 + 12] + ( kCD_z * KET_g_s[ibra * 15 + 7] );

            // |g_1_1_2  p_1_0_0) = |h_2_1_2  s_0_0_0)_{t} + x_cd * |g_1_1_2  s_0_0_0)_{t}
            const double K_1_1_2_1_0_0 = KET_h_s[ibra * 21 + 8] + ( kCD_x * KET_g_s[ibra * 15 + 8] );

            // |g_1_1_2  p_0_1_0) = |h_1_2_2  s_0_0_0)_{t} + y_cd * |g_1_1_2  s_0_0_0)_{t}
            const double K_1_1_2_0_1_0 = KET_h_s[ibra * 21 + 12] + ( kCD_y * KET_g_s[ibra * 15 + 8] );

            // |g_1_1_2  p_0_0_1) = |h_1_1_3  s_0_0_0)_{t} + z_cd * |g_1_1_2  s_0_0_0)_{t}
            const double K_1_1_2_0_0_1 = KET_h_s[ibra * 21 + 13] + ( kCD_z * KET_g_s[ibra * 15 + 8] );

            // |g_1_0_3  p_1_0_0) = |h_2_0_3  s_0_0_0)_{t} + x_cd * |g_1_0_3  s_0_0_0)_{t}
            const double K_1_0_3_1_0_0 = KET_h_s[ibra * 21 + 9] + ( kCD_x * KET_g_s[ibra * 15 + 9] );

            // |g_1_0_3  p_0_1_0) = |h_1_1_3  s_0_0_0)_{t} + y_cd * |g_1_0_3  s_0_0_0)_{t}
            const double K_1_0_3_0_1_0 = KET_h_s[ibra * 21 + 13] + ( kCD_y * KET_g_s[ibra * 15 + 9] );

            // |g_1_0_3  p_0_0_1) = |h_1_0_4  s_0_0_0)_{t} + z_cd * |g_1_0_3  s_0_0_0)_{t}
            const double K_1_0_3_0_0_1 = KET_h_s[ibra * 21 + 14] + ( kCD_z * KET_g_s[ibra * 15 + 9] );

            // |g_0_4_0  p_1_0_0) = |h_1_4_0  s_0_0_0)_{t} + x_cd * |g_0_4_0  s_0_0_0)_{t}
            const double K_0_4_0_1_0_0 = KET_h_s[ibra * 21 + 10] + ( kCD_x * KET_g_s[ibra * 15 + 10] );

            // |g_0_4_0  p_0_1_0) = |h_0_5_0  s_0_0_0)_{t} + y_cd * |g_0_4_0  s_0_0_0)_{t}
            const double K_0_4_0_0_1_0 = KET_h_s[ibra * 21 + 15] + ( kCD_y * KET_g_s[ibra * 15 + 10] );

            // |g_0_4_0  p_0_0_1) = |h_0_4_1  s_0_0_0)_{t} + z_cd * |g_0_4_0  s_0_0_0)_{t}
            const double K_0_4_0_0_0_1 = KET_h_s[ibra * 21 + 16] + ( kCD_z * KET_g_s[ibra * 15 + 10] );

            // |g_0_3_1  p_1_0_0) = |h_1_3_1  s_0_0_0)_{t} + x_cd * |g_0_3_1  s_0_0_0)_{t}
            const double K_0_3_1_1_0_0 = KET_h_s[ibra * 21 + 11] + ( kCD_x * KET_g_s[ibra * 15 + 11] );

            // |g_0_3_1  p_0_1_0) = |h_0_4_1  s_0_0_0)_{t} + y_cd * |g_0_3_1  s_0_0_0)_{t}
            const double K_0_3_1_0_1_0 = KET_h_s[ibra * 21 + 16] + ( kCD_y * KET_g_s[ibra * 15 + 11] );

            // |g_0_3_1  p_0_0_1) = |h_0_3_2  s_0_0_0)_{t} + z_cd * |g_0_3_1  s_0_0_0)_{t}
            const double K_0_3_1_0_0_1 = KET_h_s[ibra * 21 + 17] + ( kCD_z * KET_g_s[ibra * 15 + 11] );

            // |g_0_2_2  p_1_0_0) = |h_1_2_2  s_0_0_0)_{t} + x_cd * |g_0_2_2  s_0_0_0)_{t}
            const double K_0_2_2_1_0_0 = KET_h_s[ibra * 21 + 12] + ( kCD_x * KET_g_s[ibra * 15 + 12] );

            // |g_0_2_2  p_0_1_0) = |h_0_3_2  s_0_0_0)_{t} + y_cd * |g_0_2_2  s_0_0_0)_{t}
            const double K_0_2_2_0_1_0 = KET_h_s[ibra * 21 + 17] + ( kCD_y * KET_g_s[ibra * 15 + 12] );

            // |g_0_2_2  p_0_0_1) = |h_0_2_3  s_0_0_0)_{t} + z_cd * |g_0_2_2  s_0_0_0)_{t}
            const double K_0_2_2_0_0_1 = KET_h_s[ibra * 21 + 18] + ( kCD_z * KET_g_s[ibra * 15 + 12] );

            // |g_0_1_3  p_1_0_0) = |h_1_1_3  s_0_0_0)_{t} + x_cd * |g_0_1_3  s_0_0_0)_{t}
            const double K_0_1_3_1_0_0 = KET_h_s[ibra * 21 + 13] + ( kCD_x * KET_g_s[ibra * 15 + 13] );

            // |g_0_1_3  p_0_1_0) = |h_0_2_3  s_0_0_0)_{t} + y_cd * |g_0_1_3  s_0_0_0)_{t}
            const double K_0_1_3_0_1_0 = KET_h_s[ibra * 21 + 18] + ( kCD_y * KET_g_s[ibra * 15 + 13] );

            // |g_0_1_3  p_0_0_1) = |h_0_1_4  s_0_0_0)_{t} + z_cd * |g_0_1_3  s_0_0_0)_{t}
            const double K_0_1_3_0_0_1 = KET_h_s[ibra * 21 + 19] + ( kCD_z * KET_g_s[ibra * 15 + 13] );

            // |g_0_0_4  p_1_0_0) = |h_1_0_4  s_0_0_0)_{t} + x_cd * |g_0_0_4  s_0_0_0)_{t}
            const double K_0_0_4_1_0_0 = KET_h_s[ibra * 21 + 14] + ( kCD_x * KET_g_s[ibra * 15 + 14] );

            // |g_0_0_4  p_0_1_0) = |h_0_1_4  s_0_0_0)_{t} + y_cd * |g_0_0_4  s_0_0_0)_{t}
            const double K_0_0_4_0_1_0 = KET_h_s[ibra * 21 + 19] + ( kCD_y * KET_g_s[ibra * 15 + 14] );

            // |g_0_0_4  p_0_0_1) = |h_0_0_5  s_0_0_0)_{t} + z_cd * |g_0_0_4  s_0_0_0)_{t}
            const double K_0_0_4_0_0_1 = KET_h_s[ibra * 21 + 20] + ( kCD_z * KET_g_s[ibra * 15 + 14] );

            // |h_5_0_0  p_1_0_0) = |i_6_0_0  s_0_0_0)_{t} + x_cd * |h_5_0_0  s_0_0_0)_{t}
            const double K_5_0_0_1_0_0 = KET_i_s[ibra * 28 + 0] + ( kCD_x * KET_h_s[ibra * 21 + 0] );

            // |h_4_1_0  p_1_0_0) = |i_5_1_0  s_0_0_0)_{t} + x_cd * |h_4_1_0  s_0_0_0)_{t}
            const double K_4_1_0_1_0_0 = KET_i_s[ibra * 28 + 1] + ( kCD_x * KET_h_s[ibra * 21 + 1] );

            // |h_4_1_0  p_0_1_0) = |i_4_2_0  s_0_0_0)_{t} + y_cd * |h_4_1_0  s_0_0_0)_{t}
            const double K_4_1_0_0_1_0 = KET_i_s[ibra * 28 + 3] + ( kCD_y * KET_h_s[ibra * 21 + 1] );

            // |h_4_0_1  p_1_0_0) = |i_5_0_1  s_0_0_0)_{t} + x_cd * |h_4_0_1  s_0_0_0)_{t}
            const double K_4_0_1_1_0_0 = KET_i_s[ibra * 28 + 2] + ( kCD_x * KET_h_s[ibra * 21 + 2] );

            // |h_4_0_1  p_0_0_1) = |i_4_0_2  s_0_0_0)_{t} + z_cd * |h_4_0_1  s_0_0_0)_{t}
            const double K_4_0_1_0_0_1 = KET_i_s[ibra * 28 + 5] + ( kCD_z * KET_h_s[ibra * 21 + 2] );

            // |h_3_2_0  p_1_0_0) = |i_4_2_0  s_0_0_0)_{t} + x_cd * |h_3_2_0  s_0_0_0)_{t}
            const double K_3_2_0_1_0_0 = KET_i_s[ibra * 28 + 3] + ( kCD_x * KET_h_s[ibra * 21 + 3] );

            // |h_3_2_0  p_0_1_0) = |i_3_3_0  s_0_0_0)_{t} + y_cd * |h_3_2_0  s_0_0_0)_{t}
            const double K_3_2_0_0_1_0 = KET_i_s[ibra * 28 + 6] + ( kCD_y * KET_h_s[ibra * 21 + 3] );

            // |h_3_1_1  p_1_0_0) = |i_4_1_1  s_0_0_0)_{t} + x_cd * |h_3_1_1  s_0_0_0)_{t}
            const double K_3_1_1_1_0_0 = KET_i_s[ibra * 28 + 4] + ( kCD_x * KET_h_s[ibra * 21 + 4] );

            // |h_3_1_1  p_0_1_0) = |i_3_2_1  s_0_0_0)_{t} + y_cd * |h_3_1_1  s_0_0_0)_{t}
            const double K_3_1_1_0_1_0 = KET_i_s[ibra * 28 + 7] + ( kCD_y * KET_h_s[ibra * 21 + 4] );

            // |h_3_1_1  p_0_0_1) = |i_3_1_2  s_0_0_0)_{t} + z_cd * |h_3_1_1  s_0_0_0)_{t}
            const double K_3_1_1_0_0_1 = KET_i_s[ibra * 28 + 8] + ( kCD_z * KET_h_s[ibra * 21 + 4] );

            // |h_3_0_2  p_1_0_0) = |i_4_0_2  s_0_0_0)_{t} + x_cd * |h_3_0_2  s_0_0_0)_{t}
            const double K_3_0_2_1_0_0 = KET_i_s[ibra * 28 + 5] + ( kCD_x * KET_h_s[ibra * 21 + 5] );

            // |h_3_0_2  p_0_0_1) = |i_3_0_3  s_0_0_0)_{t} + z_cd * |h_3_0_2  s_0_0_0)_{t}
            const double K_3_0_2_0_0_1 = KET_i_s[ibra * 28 + 9] + ( kCD_z * KET_h_s[ibra * 21 + 5] );

            // |h_2_3_0  p_1_0_0) = |i_3_3_0  s_0_0_0)_{t} + x_cd * |h_2_3_0  s_0_0_0)_{t}
            const double K_2_3_0_1_0_0 = KET_i_s[ibra * 28 + 6] + ( kCD_x * KET_h_s[ibra * 21 + 6] );

            // |h_2_3_0  p_0_1_0) = |i_2_4_0  s_0_0_0)_{t} + y_cd * |h_2_3_0  s_0_0_0)_{t}
            const double K_2_3_0_0_1_0 = KET_i_s[ibra * 28 + 10] + ( kCD_y * KET_h_s[ibra * 21 + 6] );

            // |h_2_2_1  p_1_0_0) = |i_3_2_1  s_0_0_0)_{t} + x_cd * |h_2_2_1  s_0_0_0)_{t}
            const double K_2_2_1_1_0_0 = KET_i_s[ibra * 28 + 7] + ( kCD_x * KET_h_s[ibra * 21 + 7] );

            // |h_2_2_1  p_0_1_0) = |i_2_3_1  s_0_0_0)_{t} + y_cd * |h_2_2_1  s_0_0_0)_{t}
            const double K_2_2_1_0_1_0 = KET_i_s[ibra * 28 + 11] + ( kCD_y * KET_h_s[ibra * 21 + 7] );

            // |h_2_2_1  p_0_0_1) = |i_2_2_2  s_0_0_0)_{t} + z_cd * |h_2_2_1  s_0_0_0)_{t}
            const double K_2_2_1_0_0_1 = KET_i_s[ibra * 28 + 12] + ( kCD_z * KET_h_s[ibra * 21 + 7] );

            // |h_2_1_2  p_1_0_0) = |i_3_1_2  s_0_0_0)_{t} + x_cd * |h_2_1_2  s_0_0_0)_{t}
            const double K_2_1_2_1_0_0 = KET_i_s[ibra * 28 + 8] + ( kCD_x * KET_h_s[ibra * 21 + 8] );

            // |h_2_1_2  p_0_1_0) = |i_2_2_2  s_0_0_0)_{t} + y_cd * |h_2_1_2  s_0_0_0)_{t}
            const double K_2_1_2_0_1_0 = KET_i_s[ibra * 28 + 12] + ( kCD_y * KET_h_s[ibra * 21 + 8] );

            // |h_2_1_2  p_0_0_1) = |i_2_1_3  s_0_0_0)_{t} + z_cd * |h_2_1_2  s_0_0_0)_{t}
            const double K_2_1_2_0_0_1 = KET_i_s[ibra * 28 + 13] + ( kCD_z * KET_h_s[ibra * 21 + 8] );

            // |h_2_0_3  p_1_0_0) = |i_3_0_3  s_0_0_0)_{t} + x_cd * |h_2_0_3  s_0_0_0)_{t}
            const double K_2_0_3_1_0_0 = KET_i_s[ibra * 28 + 9] + ( kCD_x * KET_h_s[ibra * 21 + 9] );

            // |h_2_0_3  p_0_0_1) = |i_2_0_4  s_0_0_0)_{t} + z_cd * |h_2_0_3  s_0_0_0)_{t}
            const double K_2_0_3_0_0_1 = KET_i_s[ibra * 28 + 14] + ( kCD_z * KET_h_s[ibra * 21 + 9] );

            // |h_1_4_0  p_1_0_0) = |i_2_4_0  s_0_0_0)_{t} + x_cd * |h_1_4_0  s_0_0_0)_{t}
            const double K_1_4_0_1_0_0 = KET_i_s[ibra * 28 + 10] + ( kCD_x * KET_h_s[ibra * 21 + 10] );

            // |h_1_4_0  p_0_1_0) = |i_1_5_0  s_0_0_0)_{t} + y_cd * |h_1_4_0  s_0_0_0)_{t}
            const double K_1_4_0_0_1_0 = KET_i_s[ibra * 28 + 15] + ( kCD_y * KET_h_s[ibra * 21 + 10] );

            // |h_1_3_1  p_1_0_0) = |i_2_3_1  s_0_0_0)_{t} + x_cd * |h_1_3_1  s_0_0_0)_{t}
            const double K_1_3_1_1_0_0 = KET_i_s[ibra * 28 + 11] + ( kCD_x * KET_h_s[ibra * 21 + 11] );

            // |h_1_3_1  p_0_1_0) = |i_1_4_1  s_0_0_0)_{t} + y_cd * |h_1_3_1  s_0_0_0)_{t}
            const double K_1_3_1_0_1_0 = KET_i_s[ibra * 28 + 16] + ( kCD_y * KET_h_s[ibra * 21 + 11] );

            // |h_1_3_1  p_0_0_1) = |i_1_3_2  s_0_0_0)_{t} + z_cd * |h_1_3_1  s_0_0_0)_{t}
            const double K_1_3_1_0_0_1 = KET_i_s[ibra * 28 + 17] + ( kCD_z * KET_h_s[ibra * 21 + 11] );

            // |h_1_2_2  p_1_0_0) = |i_2_2_2  s_0_0_0)_{t} + x_cd * |h_1_2_2  s_0_0_0)_{t}
            const double K_1_2_2_1_0_0 = KET_i_s[ibra * 28 + 12] + ( kCD_x * KET_h_s[ibra * 21 + 12] );

            // |h_1_2_2  p_0_1_0) = |i_1_3_2  s_0_0_0)_{t} + y_cd * |h_1_2_2  s_0_0_0)_{t}
            const double K_1_2_2_0_1_0 = KET_i_s[ibra * 28 + 17] + ( kCD_y * KET_h_s[ibra * 21 + 12] );

            // |h_1_2_2  p_0_0_1) = |i_1_2_3  s_0_0_0)_{t} + z_cd * |h_1_2_2  s_0_0_0)_{t}
            const double K_1_2_2_0_0_1 = KET_i_s[ibra * 28 + 18] + ( kCD_z * KET_h_s[ibra * 21 + 12] );

            // |h_1_1_3  p_1_0_0) = |i_2_1_3  s_0_0_0)_{t} + x_cd * |h_1_1_3  s_0_0_0)_{t}
            const double K_1_1_3_1_0_0 = KET_i_s[ibra * 28 + 13] + ( kCD_x * KET_h_s[ibra * 21 + 13] );

            // |h_1_1_3  p_0_1_0) = |i_1_2_3  s_0_0_0)_{t} + y_cd * |h_1_1_3  s_0_0_0)_{t}
            const double K_1_1_3_0_1_0 = KET_i_s[ibra * 28 + 18] + ( kCD_y * KET_h_s[ibra * 21 + 13] );

            // |h_1_1_3  p_0_0_1) = |i_1_1_4  s_0_0_0)_{t} + z_cd * |h_1_1_3  s_0_0_0)_{t}
            const double K_1_1_3_0_0_1 = KET_i_s[ibra * 28 + 19] + ( kCD_z * KET_h_s[ibra * 21 + 13] );

            // |h_1_0_4  p_1_0_0) = |i_2_0_4  s_0_0_0)_{t} + x_cd * |h_1_0_4  s_0_0_0)_{t}
            const double K_1_0_4_1_0_0 = KET_i_s[ibra * 28 + 14] + ( kCD_x * KET_h_s[ibra * 21 + 14] );

            // |h_1_0_4  p_0_0_1) = |i_1_0_5  s_0_0_0)_{t} + z_cd * |h_1_0_4  s_0_0_0)_{t}
            const double K_1_0_4_0_0_1 = KET_i_s[ibra * 28 + 20] + ( kCD_z * KET_h_s[ibra * 21 + 14] );

            // |h_0_5_0  p_0_1_0) = |i_0_6_0  s_0_0_0)_{t} + y_cd * |h_0_5_0  s_0_0_0)_{t}
            const double K_0_5_0_0_1_0 = KET_i_s[ibra * 28 + 21] + ( kCD_y * KET_h_s[ibra * 21 + 15] );

            // |h_0_4_1  p_1_0_0) = |i_1_4_1  s_0_0_0)_{t} + x_cd * |h_0_4_1  s_0_0_0)_{t}
            const double K_0_4_1_1_0_0 = KET_i_s[ibra * 28 + 16] + ( kCD_x * KET_h_s[ibra * 21 + 16] );

            // |h_0_4_1  p_0_1_0) = |i_0_5_1  s_0_0_0)_{t} + y_cd * |h_0_4_1  s_0_0_0)_{t}
            const double K_0_4_1_0_1_0 = KET_i_s[ibra * 28 + 22] + ( kCD_y * KET_h_s[ibra * 21 + 16] );

            // |h_0_4_1  p_0_0_1) = |i_0_4_2  s_0_0_0)_{t} + z_cd * |h_0_4_1  s_0_0_0)_{t}
            const double K_0_4_1_0_0_1 = KET_i_s[ibra * 28 + 23] + ( kCD_z * KET_h_s[ibra * 21 + 16] );

            // |h_0_3_2  p_1_0_0) = |i_1_3_2  s_0_0_0)_{t} + x_cd * |h_0_3_2  s_0_0_0)_{t}
            const double K_0_3_2_1_0_0 = KET_i_s[ibra * 28 + 17] + ( kCD_x * KET_h_s[ibra * 21 + 17] );

            // |h_0_3_2  p_0_1_0) = |i_0_4_2  s_0_0_0)_{t} + y_cd * |h_0_3_2  s_0_0_0)_{t}
            const double K_0_3_2_0_1_0 = KET_i_s[ibra * 28 + 23] + ( kCD_y * KET_h_s[ibra * 21 + 17] );

            // |h_0_3_2  p_0_0_1) = |i_0_3_3  s_0_0_0)_{t} + z_cd * |h_0_3_2  s_0_0_0)_{t}
            const double K_0_3_2_0_0_1 = KET_i_s[ibra * 28 + 24] + ( kCD_z * KET_h_s[ibra * 21 + 17] );

            // |h_0_2_3  p_1_0_0) = |i_1_2_3  s_0_0_0)_{t} + x_cd * |h_0_2_3  s_0_0_0)_{t}
            const double K_0_2_3_1_0_0 = KET_i_s[ibra * 28 + 18] + ( kCD_x * KET_h_s[ibra * 21 + 18] );

            // |h_0_2_3  p_0_1_0) = |i_0_3_3  s_0_0_0)_{t} + y_cd * |h_0_2_3  s_0_0_0)_{t}
            const double K_0_2_3_0_1_0 = KET_i_s[ibra * 28 + 24] + ( kCD_y * KET_h_s[ibra * 21 + 18] );

            // |h_0_2_3  p_0_0_1) = |i_0_2_4  s_0_0_0)_{t} + z_cd * |h_0_2_3  s_0_0_0)_{t}
            const double K_0_2_3_0_0_1 = KET_i_s[ibra * 28 + 25] + ( kCD_z * KET_h_s[ibra * 21 + 18] );

            // |h_0_1_4  p_1_0_0) = |i_1_1_4  s_0_0_0)_{t} + x_cd * |h_0_1_4  s_0_0_0)_{t}
            const double K_0_1_4_1_0_0 = KET_i_s[ibra * 28 + 19] + ( kCD_x * KET_h_s[ibra * 21 + 19] );

            // |h_0_1_4  p_0_1_0) = |i_0_2_4  s_0_0_0)_{t} + y_cd * |h_0_1_4  s_0_0_0)_{t}
            const double K_0_1_4_0_1_0 = KET_i_s[ibra * 28 + 25] + ( kCD_y * KET_h_s[ibra * 21 + 19] );

            // |h_0_1_4  p_0_0_1) = |i_0_1_5  s_0_0_0)_{t} + z_cd * |h_0_1_4  s_0_0_0)_{t}
            const double K_0_1_4_0_0_1 = KET_i_s[ibra * 28 + 26] + ( kCD_z * KET_h_s[ibra * 21 + 19] );

            // |h_0_0_5  p_0_0_1) = |i_0_0_6  s_0_0_0)_{t} + z_cd * |h_0_0_5  s_0_0_0)_{t}
            const double K_0_0_5_0_0_1 = KET_i_s[ibra * 28 + 27] + ( kCD_z * KET_h_s[ibra * 21 + 20] );

            // |f_3_0_0  d_2_0_0) = |g_4_0_0  p_1_0_0) + x_cd * |f_3_0_0  p_1_0_0)
            const double K_3_0_0_2_0_0 = K_4_0_0_1_0_0 + ( kCD_x * K_3_0_0_1_0_0 );

            // |f_3_0_0  d_1_1_0) = |g_3_1_0  p_1_0_0) + y_cd * |f_3_0_0  p_1_0_0)
            const double K_3_0_0_1_1_0 = K_3_1_0_1_0_0 + ( kCD_y * K_3_0_0_1_0_0 );

            // |f_3_0_0  d_0_2_0) = |g_3_1_0  p_0_1_0) + y_cd * |f_3_0_0  p_0_1_0)
            const double K_3_0_0_0_2_0 = K_3_1_0_0_1_0 + ( kCD_y * K_3_0_0_0_1_0 );

            // |f_3_0_0  d_0_0_2) = |g_3_0_1  p_0_0_1) + z_cd * |f_3_0_0  p_0_0_1)
            const double K_3_0_0_0_0_2 = K_3_0_1_0_0_1 + ( kCD_z * K_3_0_0_0_0_1 );

            // |f_2_1_0  d_2_0_0) = |g_3_1_0  p_1_0_0) + x_cd * |f_2_1_0  p_1_0_0)
            const double K_2_1_0_2_0_0 = K_3_1_0_1_0_0 + ( kCD_x * K_2_1_0_1_0_0 );

            // |f_2_1_0  d_1_1_0) = |g_2_2_0  p_1_0_0) + y_cd * |f_2_1_0  p_1_0_0)
            const double K_2_1_0_1_1_0 = K_2_2_0_1_0_0 + ( kCD_y * K_2_1_0_1_0_0 );

            // |f_2_1_0  d_0_2_0) = |g_2_2_0  p_0_1_0) + y_cd * |f_2_1_0  p_0_1_0)
            const double K_2_1_0_0_2_0 = K_2_2_0_0_1_0 + ( kCD_y * K_2_1_0_0_1_0 );

            // |f_2_1_0  d_0_0_2) = |g_2_1_1  p_0_0_1) + z_cd * |f_2_1_0  p_0_0_1)
            const double K_2_1_0_0_0_2 = K_2_1_1_0_0_1 + ( kCD_z * K_2_1_0_0_0_1 );

            // |f_2_0_1  d_2_0_0) = |g_3_0_1  p_1_0_0) + x_cd * |f_2_0_1  p_1_0_0)
            const double K_2_0_1_2_0_0 = K_3_0_1_1_0_0 + ( kCD_x * K_2_0_1_1_0_0 );

            // |f_2_0_1  d_1_1_0) = |g_2_1_1  p_1_0_0) + y_cd * |f_2_0_1  p_1_0_0)
            const double K_2_0_1_1_1_0 = K_2_1_1_1_0_0 + ( kCD_y * K_2_0_1_1_0_0 );

            // |f_2_0_1  d_0_2_0) = |g_2_1_1  p_0_1_0) + y_cd * |f_2_0_1  p_0_1_0)
            const double K_2_0_1_0_2_0 = K_2_1_1_0_1_0 + ( kCD_y * K_2_0_1_0_1_0 );

            // |f_2_0_1  d_0_0_2) = |g_2_0_2  p_0_0_1) + z_cd * |f_2_0_1  p_0_0_1)
            const double K_2_0_1_0_0_2 = K_2_0_2_0_0_1 + ( kCD_z * K_2_0_1_0_0_1 );

            // |f_1_2_0  d_2_0_0) = |g_2_2_0  p_1_0_0) + x_cd * |f_1_2_0  p_1_0_0)
            const double K_1_2_0_2_0_0 = K_2_2_0_1_0_0 + ( kCD_x * K_1_2_0_1_0_0 );

            // |f_1_2_0  d_1_1_0) = |g_1_3_0  p_1_0_0) + y_cd * |f_1_2_0  p_1_0_0)
            const double K_1_2_0_1_1_0 = K_1_3_0_1_0_0 + ( kCD_y * K_1_2_0_1_0_0 );

            // |f_1_2_0  d_0_2_0) = |g_1_3_0  p_0_1_0) + y_cd * |f_1_2_0  p_0_1_0)
            const double K_1_2_0_0_2_0 = K_1_3_0_0_1_0 + ( kCD_y * K_1_2_0_0_1_0 );

            // |f_1_2_0  d_0_0_2) = |g_1_2_1  p_0_0_1) + z_cd * |f_1_2_0  p_0_0_1)
            const double K_1_2_0_0_0_2 = K_1_2_1_0_0_1 + ( kCD_z * K_1_2_0_0_0_1 );

            // |f_1_1_1  d_2_0_0) = |g_2_1_1  p_1_0_0) + x_cd * |f_1_1_1  p_1_0_0)
            const double K_1_1_1_2_0_0 = K_2_1_1_1_0_0 + ( kCD_x * K_1_1_1_1_0_0 );

            // |f_1_1_1  d_1_1_0) = |g_1_2_1  p_1_0_0) + y_cd * |f_1_1_1  p_1_0_0)
            const double K_1_1_1_1_1_0 = K_1_2_1_1_0_0 + ( kCD_y * K_1_1_1_1_0_0 );

            // |f_1_1_1  d_0_2_0) = |g_1_2_1  p_0_1_0) + y_cd * |f_1_1_1  p_0_1_0)
            const double K_1_1_1_0_2_0 = K_1_2_1_0_1_0 + ( kCD_y * K_1_1_1_0_1_0 );

            // |f_1_1_1  d_0_0_2) = |g_1_1_2  p_0_0_1) + z_cd * |f_1_1_1  p_0_0_1)
            const double K_1_1_1_0_0_2 = K_1_1_2_0_0_1 + ( kCD_z * K_1_1_1_0_0_1 );

            // |f_1_0_2  d_2_0_0) = |g_2_0_2  p_1_0_0) + x_cd * |f_1_0_2  p_1_0_0)
            const double K_1_0_2_2_0_0 = K_2_0_2_1_0_0 + ( kCD_x * K_1_0_2_1_0_0 );

            // |f_1_0_2  d_1_1_0) = |g_1_1_2  p_1_0_0) + y_cd * |f_1_0_2  p_1_0_0)
            const double K_1_0_2_1_1_0 = K_1_1_2_1_0_0 + ( kCD_y * K_1_0_2_1_0_0 );

            // |f_1_0_2  d_0_2_0) = |g_1_1_2  p_0_1_0) + y_cd * |f_1_0_2  p_0_1_0)
            const double K_1_0_2_0_2_0 = K_1_1_2_0_1_0 + ( kCD_y * K_1_0_2_0_1_0 );

            // |f_1_0_2  d_0_0_2) = |g_1_0_3  p_0_0_1) + z_cd * |f_1_0_2  p_0_0_1)
            const double K_1_0_2_0_0_2 = K_1_0_3_0_0_1 + ( kCD_z * K_1_0_2_0_0_1 );

            // |f_0_3_0  d_2_0_0) = |g_1_3_0  p_1_0_0) + x_cd * |f_0_3_0  p_1_0_0)
            const double K_0_3_0_2_0_0 = K_1_3_0_1_0_0 + ( kCD_x * K_0_3_0_1_0_0 );

            // |f_0_3_0  d_1_1_0) = |g_0_4_0  p_1_0_0) + y_cd * |f_0_3_0  p_1_0_0)
            const double K_0_3_0_1_1_0 = K_0_4_0_1_0_0 + ( kCD_y * K_0_3_0_1_0_0 );

            // |f_0_3_0  d_0_2_0) = |g_0_4_0  p_0_1_0) + y_cd * |f_0_3_0  p_0_1_0)
            const double K_0_3_0_0_2_0 = K_0_4_0_0_1_0 + ( kCD_y * K_0_3_0_0_1_0 );

            // |f_0_3_0  d_0_0_2) = |g_0_3_1  p_0_0_1) + z_cd * |f_0_3_0  p_0_0_1)
            const double K_0_3_0_0_0_2 = K_0_3_1_0_0_1 + ( kCD_z * K_0_3_0_0_0_1 );

            // |f_0_2_1  d_2_0_0) = |g_1_2_1  p_1_0_0) + x_cd * |f_0_2_1  p_1_0_0)
            const double K_0_2_1_2_0_0 = K_1_2_1_1_0_0 + ( kCD_x * K_0_2_1_1_0_0 );

            // |f_0_2_1  d_1_1_0) = |g_0_3_1  p_1_0_0) + y_cd * |f_0_2_1  p_1_0_0)
            const double K_0_2_1_1_1_0 = K_0_3_1_1_0_0 + ( kCD_y * K_0_2_1_1_0_0 );

            // |f_0_2_1  d_0_2_0) = |g_0_3_1  p_0_1_0) + y_cd * |f_0_2_1  p_0_1_0)
            const double K_0_2_1_0_2_0 = K_0_3_1_0_1_0 + ( kCD_y * K_0_2_1_0_1_0 );

            // |f_0_2_1  d_0_0_2) = |g_0_2_2  p_0_0_1) + z_cd * |f_0_2_1  p_0_0_1)
            const double K_0_2_1_0_0_2 = K_0_2_2_0_0_1 + ( kCD_z * K_0_2_1_0_0_1 );

            // |f_0_1_2  d_2_0_0) = |g_1_1_2  p_1_0_0) + x_cd * |f_0_1_2  p_1_0_0)
            const double K_0_1_2_2_0_0 = K_1_1_2_1_0_0 + ( kCD_x * K_0_1_2_1_0_0 );

            // |f_0_1_2  d_1_1_0) = |g_0_2_2  p_1_0_0) + y_cd * |f_0_1_2  p_1_0_0)
            const double K_0_1_2_1_1_0 = K_0_2_2_1_0_0 + ( kCD_y * K_0_1_2_1_0_0 );

            // |f_0_1_2  d_0_2_0) = |g_0_2_2  p_0_1_0) + y_cd * |f_0_1_2  p_0_1_0)
            const double K_0_1_2_0_2_0 = K_0_2_2_0_1_0 + ( kCD_y * K_0_1_2_0_1_0 );

            // |f_0_1_2  d_0_0_2) = |g_0_1_3  p_0_0_1) + z_cd * |f_0_1_2  p_0_0_1)
            const double K_0_1_2_0_0_2 = K_0_1_3_0_0_1 + ( kCD_z * K_0_1_2_0_0_1 );

            // |f_0_0_3  d_2_0_0) = |g_1_0_3  p_1_0_0) + x_cd * |f_0_0_3  p_1_0_0)
            const double K_0_0_3_2_0_0 = K_1_0_3_1_0_0 + ( kCD_x * K_0_0_3_1_0_0 );

            // |f_0_0_3  d_1_1_0) = |g_0_1_3  p_1_0_0) + y_cd * |f_0_0_3  p_1_0_0)
            const double K_0_0_3_1_1_0 = K_0_1_3_1_0_0 + ( kCD_y * K_0_0_3_1_0_0 );

            // |f_0_0_3  d_0_2_0) = |g_0_1_3  p_0_1_0) + y_cd * |f_0_0_3  p_0_1_0)
            const double K_0_0_3_0_2_0 = K_0_1_3_0_1_0 + ( kCD_y * K_0_0_3_0_1_0 );

            // |f_0_0_3  d_0_0_2) = |g_0_0_4  p_0_0_1) + z_cd * |f_0_0_3  p_0_0_1)
            const double K_0_0_3_0_0_2 = K_0_0_4_0_0_1 + ( kCD_z * K_0_0_3_0_0_1 );

            // |g_4_0_0  d_2_0_0) = |h_5_0_0  p_1_0_0) + x_cd * |g_4_0_0  p_1_0_0)
            const double K_4_0_0_2_0_0 = K_5_0_0_1_0_0 + ( kCD_x * K_4_0_0_1_0_0 );

            // |g_4_0_0  d_0_2_0) = |h_4_1_0  p_0_1_0) + y_cd * |g_4_0_0  p_0_1_0)
            const double K_4_0_0_0_2_0 = K_4_1_0_0_1_0 + ( kCD_y * K_4_0_0_0_1_0 );

            // |g_4_0_0  d_0_0_2) = |h_4_0_1  p_0_0_1) + z_cd * |g_4_0_0  p_0_0_1)
            const double K_4_0_0_0_0_2 = K_4_0_1_0_0_1 + ( kCD_z * K_4_0_0_0_0_1 );

            // |g_3_1_0  d_2_0_0) = |h_4_1_0  p_1_0_0) + x_cd * |g_3_1_0  p_1_0_0)
            const double K_3_1_0_2_0_0 = K_4_1_0_1_0_0 + ( kCD_x * K_3_1_0_1_0_0 );

            // |g_3_1_0  d_0_2_0) = |h_3_2_0  p_0_1_0) + y_cd * |g_3_1_0  p_0_1_0)
            const double K_3_1_0_0_2_0 = K_3_2_0_0_1_0 + ( kCD_y * K_3_1_0_0_1_0 );

            // |g_3_1_0  d_0_0_2) = |h_3_1_1  p_0_0_1) + z_cd * |g_3_1_0  p_0_0_1)
            const double K_3_1_0_0_0_2 = K_3_1_1_0_0_1 + ( kCD_z * K_3_1_0_0_0_1 );

            // |g_3_0_1  d_2_0_0) = |h_4_0_1  p_1_0_0) + x_cd * |g_3_0_1  p_1_0_0)
            const double K_3_0_1_2_0_0 = K_4_0_1_1_0_0 + ( kCD_x * K_3_0_1_1_0_0 );

            // |g_3_0_1  d_1_1_0) = |h_3_1_1  p_1_0_0) + y_cd * |g_3_0_1  p_1_0_0)
            const double K_3_0_1_1_1_0 = K_3_1_1_1_0_0 + ( kCD_y * K_3_0_1_1_0_0 );

            // |g_3_0_1  d_0_2_0) = |h_3_1_1  p_0_1_0) + y_cd * |g_3_0_1  p_0_1_0)
            const double K_3_0_1_0_2_0 = K_3_1_1_0_1_0 + ( kCD_y * K_3_0_1_0_1_0 );

            // |g_3_0_1  d_0_0_2) = |h_3_0_2  p_0_0_1) + z_cd * |g_3_0_1  p_0_0_1)
            const double K_3_0_1_0_0_2 = K_3_0_2_0_0_1 + ( kCD_z * K_3_0_1_0_0_1 );

            // |g_2_2_0  d_2_0_0) = |h_3_2_0  p_1_0_0) + x_cd * |g_2_2_0  p_1_0_0)
            const double K_2_2_0_2_0_0 = K_3_2_0_1_0_0 + ( kCD_x * K_2_2_0_1_0_0 );

            // |g_2_2_0  d_0_2_0) = |h_2_3_0  p_0_1_0) + y_cd * |g_2_2_0  p_0_1_0)
            const double K_2_2_0_0_2_0 = K_2_3_0_0_1_0 + ( kCD_y * K_2_2_0_0_1_0 );

            // |g_2_2_0  d_0_0_2) = |h_2_2_1  p_0_0_1) + z_cd * |g_2_2_0  p_0_0_1)
            const double K_2_2_0_0_0_2 = K_2_2_1_0_0_1 + ( kCD_z * K_2_2_0_0_0_1 );

            // |g_2_1_1  d_2_0_0) = |h_3_1_1  p_1_0_0) + x_cd * |g_2_1_1  p_1_0_0)
            const double K_2_1_1_2_0_0 = K_3_1_1_1_0_0 + ( kCD_x * K_2_1_1_1_0_0 );

            // |g_2_1_1  d_1_1_0) = |h_2_2_1  p_1_0_0) + y_cd * |g_2_1_1  p_1_0_0)
            const double K_2_1_1_1_1_0 = K_2_2_1_1_0_0 + ( kCD_y * K_2_1_1_1_0_0 );

            // |g_2_1_1  d_0_2_0) = |h_2_2_1  p_0_1_0) + y_cd * |g_2_1_1  p_0_1_0)
            const double K_2_1_1_0_2_0 = K_2_2_1_0_1_0 + ( kCD_y * K_2_1_1_0_1_0 );

            // |g_2_1_1  d_0_0_2) = |h_2_1_2  p_0_0_1) + z_cd * |g_2_1_1  p_0_0_1)
            const double K_2_1_1_0_0_2 = K_2_1_2_0_0_1 + ( kCD_z * K_2_1_1_0_0_1 );

            // |g_2_0_2  d_2_0_0) = |h_3_0_2  p_1_0_0) + x_cd * |g_2_0_2  p_1_0_0)
            const double K_2_0_2_2_0_0 = K_3_0_2_1_0_0 + ( kCD_x * K_2_0_2_1_0_0 );

            // |g_2_0_2  d_1_1_0) = |h_2_1_2  p_1_0_0) + y_cd * |g_2_0_2  p_1_0_0)
            const double K_2_0_2_1_1_0 = K_2_1_2_1_0_0 + ( kCD_y * K_2_0_2_1_0_0 );

            // |g_2_0_2  d_0_2_0) = |h_2_1_2  p_0_1_0) + y_cd * |g_2_0_2  p_0_1_0)
            const double K_2_0_2_0_2_0 = K_2_1_2_0_1_0 + ( kCD_y * K_2_0_2_0_1_0 );

            // |g_2_0_2  d_0_0_2) = |h_2_0_3  p_0_0_1) + z_cd * |g_2_0_2  p_0_0_1)
            const double K_2_0_2_0_0_2 = K_2_0_3_0_0_1 + ( kCD_z * K_2_0_2_0_0_1 );

            // |g_1_3_0  d_2_0_0) = |h_2_3_0  p_1_0_0) + x_cd * |g_1_3_0  p_1_0_0)
            const double K_1_3_0_2_0_0 = K_2_3_0_1_0_0 + ( kCD_x * K_1_3_0_1_0_0 );

            // |g_1_3_0  d_0_2_0) = |h_1_4_0  p_0_1_0) + y_cd * |g_1_3_0  p_0_1_0)
            const double K_1_3_0_0_2_0 = K_1_4_0_0_1_0 + ( kCD_y * K_1_3_0_0_1_0 );

            // |g_1_3_0  d_0_0_2) = |h_1_3_1  p_0_0_1) + z_cd * |g_1_3_0  p_0_0_1)
            const double K_1_3_0_0_0_2 = K_1_3_1_0_0_1 + ( kCD_z * K_1_3_0_0_0_1 );

            // |g_1_2_1  d_2_0_0) = |h_2_2_1  p_1_0_0) + x_cd * |g_1_2_1  p_1_0_0)
            const double K_1_2_1_2_0_0 = K_2_2_1_1_0_0 + ( kCD_x * K_1_2_1_1_0_0 );

            // |g_1_2_1  d_1_1_0) = |h_1_3_1  p_1_0_0) + y_cd * |g_1_2_1  p_1_0_0)
            const double K_1_2_1_1_1_0 = K_1_3_1_1_0_0 + ( kCD_y * K_1_2_1_1_0_0 );

            // |g_1_2_1  d_0_2_0) = |h_1_3_1  p_0_1_0) + y_cd * |g_1_2_1  p_0_1_0)
            const double K_1_2_1_0_2_0 = K_1_3_1_0_1_0 + ( kCD_y * K_1_2_1_0_1_0 );

            // |g_1_2_1  d_0_0_2) = |h_1_2_2  p_0_0_1) + z_cd * |g_1_2_1  p_0_0_1)
            const double K_1_2_1_0_0_2 = K_1_2_2_0_0_1 + ( kCD_z * K_1_2_1_0_0_1 );

            // |g_1_1_2  d_2_0_0) = |h_2_1_2  p_1_0_0) + x_cd * |g_1_1_2  p_1_0_0)
            const double K_1_1_2_2_0_0 = K_2_1_2_1_0_0 + ( kCD_x * K_1_1_2_1_0_0 );

            // |g_1_1_2  d_1_1_0) = |h_1_2_2  p_1_0_0) + y_cd * |g_1_1_2  p_1_0_0)
            const double K_1_1_2_1_1_0 = K_1_2_2_1_0_0 + ( kCD_y * K_1_1_2_1_0_0 );

            // |g_1_1_2  d_0_2_0) = |h_1_2_2  p_0_1_0) + y_cd * |g_1_1_2  p_0_1_0)
            const double K_1_1_2_0_2_0 = K_1_2_2_0_1_0 + ( kCD_y * K_1_1_2_0_1_0 );

            // |g_1_1_2  d_0_0_2) = |h_1_1_3  p_0_0_1) + z_cd * |g_1_1_2  p_0_0_1)
            const double K_1_1_2_0_0_2 = K_1_1_3_0_0_1 + ( kCD_z * K_1_1_2_0_0_1 );

            // |g_1_0_3  d_2_0_0) = |h_2_0_3  p_1_0_0) + x_cd * |g_1_0_3  p_1_0_0)
            const double K_1_0_3_2_0_0 = K_2_0_3_1_0_0 + ( kCD_x * K_1_0_3_1_0_0 );

            // |g_1_0_3  d_1_1_0) = |h_1_1_3  p_1_0_0) + y_cd * |g_1_0_3  p_1_0_0)
            const double K_1_0_3_1_1_0 = K_1_1_3_1_0_0 + ( kCD_y * K_1_0_3_1_0_0 );

            // |g_1_0_3  d_0_2_0) = |h_1_1_3  p_0_1_0) + y_cd * |g_1_0_3  p_0_1_0)
            const double K_1_0_3_0_2_0 = K_1_1_3_0_1_0 + ( kCD_y * K_1_0_3_0_1_0 );

            // |g_1_0_3  d_0_0_2) = |h_1_0_4  p_0_0_1) + z_cd * |g_1_0_3  p_0_0_1)
            const double K_1_0_3_0_0_2 = K_1_0_4_0_0_1 + ( kCD_z * K_1_0_3_0_0_1 );

            // |g_0_4_0  d_2_0_0) = |h_1_4_0  p_1_0_0) + x_cd * |g_0_4_0  p_1_0_0)
            const double K_0_4_0_2_0_0 = K_1_4_0_1_0_0 + ( kCD_x * K_0_4_0_1_0_0 );

            // |g_0_4_0  d_0_2_0) = |h_0_5_0  p_0_1_0) + y_cd * |g_0_4_0  p_0_1_0)
            const double K_0_4_0_0_2_0 = K_0_5_0_0_1_0 + ( kCD_y * K_0_4_0_0_1_0 );

            // |g_0_4_0  d_0_0_2) = |h_0_4_1  p_0_0_1) + z_cd * |g_0_4_0  p_0_0_1)
            const double K_0_4_0_0_0_2 = K_0_4_1_0_0_1 + ( kCD_z * K_0_4_0_0_0_1 );

            // |g_0_3_1  d_2_0_0) = |h_1_3_1  p_1_0_0) + x_cd * |g_0_3_1  p_1_0_0)
            const double K_0_3_1_2_0_0 = K_1_3_1_1_0_0 + ( kCD_x * K_0_3_1_1_0_0 );

            // |g_0_3_1  d_1_1_0) = |h_0_4_1  p_1_0_0) + y_cd * |g_0_3_1  p_1_0_0)
            const double K_0_3_1_1_1_0 = K_0_4_1_1_0_0 + ( kCD_y * K_0_3_1_1_0_0 );

            // |g_0_3_1  d_0_2_0) = |h_0_4_1  p_0_1_0) + y_cd * |g_0_3_1  p_0_1_0)
            const double K_0_3_1_0_2_0 = K_0_4_1_0_1_0 + ( kCD_y * K_0_3_1_0_1_0 );

            // |g_0_3_1  d_0_0_2) = |h_0_3_2  p_0_0_1) + z_cd * |g_0_3_1  p_0_0_1)
            const double K_0_3_1_0_0_2 = K_0_3_2_0_0_1 + ( kCD_z * K_0_3_1_0_0_1 );

            // |g_0_2_2  d_2_0_0) = |h_1_2_2  p_1_0_0) + x_cd * |g_0_2_2  p_1_0_0)
            const double K_0_2_2_2_0_0 = K_1_2_2_1_0_0 + ( kCD_x * K_0_2_2_1_0_0 );

            // |g_0_2_2  d_1_1_0) = |h_0_3_2  p_1_0_0) + y_cd * |g_0_2_2  p_1_0_0)
            const double K_0_2_2_1_1_0 = K_0_3_2_1_0_0 + ( kCD_y * K_0_2_2_1_0_0 );

            // |g_0_2_2  d_0_2_0) = |h_0_3_2  p_0_1_0) + y_cd * |g_0_2_2  p_0_1_0)
            const double K_0_2_2_0_2_0 = K_0_3_2_0_1_0 + ( kCD_y * K_0_2_2_0_1_0 );

            // |g_0_2_2  d_0_0_2) = |h_0_2_3  p_0_0_1) + z_cd * |g_0_2_2  p_0_0_1)
            const double K_0_2_2_0_0_2 = K_0_2_3_0_0_1 + ( kCD_z * K_0_2_2_0_0_1 );

            // |g_0_1_3  d_2_0_0) = |h_1_1_3  p_1_0_0) + x_cd * |g_0_1_3  p_1_0_0)
            const double K_0_1_3_2_0_0 = K_1_1_3_1_0_0 + ( kCD_x * K_0_1_3_1_0_0 );

            // |g_0_1_3  d_1_1_0) = |h_0_2_3  p_1_0_0) + y_cd * |g_0_1_3  p_1_0_0)
            const double K_0_1_3_1_1_0 = K_0_2_3_1_0_0 + ( kCD_y * K_0_1_3_1_0_0 );

            // |g_0_1_3  d_0_2_0) = |h_0_2_3  p_0_1_0) + y_cd * |g_0_1_3  p_0_1_0)
            const double K_0_1_3_0_2_0 = K_0_2_3_0_1_0 + ( kCD_y * K_0_1_3_0_1_0 );

            // |g_0_1_3  d_0_0_2) = |h_0_1_4  p_0_0_1) + z_cd * |g_0_1_3  p_0_0_1)
            const double K_0_1_3_0_0_2 = K_0_1_4_0_0_1 + ( kCD_z * K_0_1_3_0_0_1 );

            // |g_0_0_4  d_2_0_0) = |h_1_0_4  p_1_0_0) + x_cd * |g_0_0_4  p_1_0_0)
            const double K_0_0_4_2_0_0 = K_1_0_4_1_0_0 + ( kCD_x * K_0_0_4_1_0_0 );

            // |g_0_0_4  d_1_1_0) = |h_0_1_4  p_1_0_0) + y_cd * |g_0_0_4  p_1_0_0)
            const double K_0_0_4_1_1_0 = K_0_1_4_1_0_0 + ( kCD_y * K_0_0_4_1_0_0 );

            // |g_0_0_4  d_0_2_0) = |h_0_1_4  p_0_1_0) + y_cd * |g_0_0_4  p_0_1_0)
            const double K_0_0_4_0_2_0 = K_0_1_4_0_1_0 + ( kCD_y * K_0_0_4_0_1_0 );

            // |g_0_0_4  d_0_0_2) = |h_0_0_5  p_0_0_1) + z_cd * |g_0_0_4  p_0_0_1)
            const double K_0_0_4_0_0_2 = K_0_0_5_0_0_1 + ( kCD_z * K_0_0_4_0_0_1 );

            // |f_3_0_0  f_3_0_0)_{i} = |g_4_0_0  d_2_0_0) + x_cd * |f_3_0_0  d_2_0_0)
            KET_f_f[ibra * 100 + 0] = K_4_0_0_2_0_0 + ( kCD_x * K_3_0_0_2_0_0 );

            // |f_3_0_0  f_2_1_0)_{i} = |g_3_1_0  d_2_0_0) + y_cd * |f_3_0_0  d_2_0_0)
            KET_f_f[ibra * 100 + 1] = K_3_1_0_2_0_0 + ( kCD_y * K_3_0_0_2_0_0 );

            // |f_3_0_0  f_2_0_1)_{i} = |g_3_0_1  d_2_0_0) + z_cd * |f_3_0_0  d_2_0_0)
            KET_f_f[ibra * 100 + 2] = K_3_0_1_2_0_0 + ( kCD_z * K_3_0_0_2_0_0 );

            // |f_3_0_0  f_1_2_0)_{i} = |g_4_0_0  d_0_2_0) + x_cd * |f_3_0_0  d_0_2_0)
            KET_f_f[ibra * 100 + 3] = K_4_0_0_0_2_0 + ( kCD_x * K_3_0_0_0_2_0 );

            // |f_3_0_0  f_1_1_1)_{i} = |g_3_0_1  d_1_1_0) + z_cd * |f_3_0_0  d_1_1_0)
            KET_f_f[ibra * 100 + 4] = K_3_0_1_1_1_0 + ( kCD_z * K_3_0_0_1_1_0 );

            // |f_3_0_0  f_1_0_2)_{i} = |g_4_0_0  d_0_0_2) + x_cd * |f_3_0_0  d_0_0_2)
            KET_f_f[ibra * 100 + 5] = K_4_0_0_0_0_2 + ( kCD_x * K_3_0_0_0_0_2 );

            // |f_3_0_0  f_0_3_0)_{i} = |g_3_1_0  d_0_2_0) + y_cd * |f_3_0_0  d_0_2_0)
            KET_f_f[ibra * 100 + 6] = K_3_1_0_0_2_0 + ( kCD_y * K_3_0_0_0_2_0 );

            // |f_3_0_0  f_0_2_1)_{i} = |g_3_0_1  d_0_2_0) + z_cd * |f_3_0_0  d_0_2_0)
            KET_f_f[ibra * 100 + 7] = K_3_0_1_0_2_0 + ( kCD_z * K_3_0_0_0_2_0 );

            // |f_3_0_0  f_0_1_2)_{i} = |g_3_1_0  d_0_0_2) + y_cd * |f_3_0_0  d_0_0_2)
            KET_f_f[ibra * 100 + 8] = K_3_1_0_0_0_2 + ( kCD_y * K_3_0_0_0_0_2 );

            // |f_3_0_0  f_0_0_3)_{i} = |g_3_0_1  d_0_0_2) + z_cd * |f_3_0_0  d_0_0_2)
            KET_f_f[ibra * 100 + 9] = K_3_0_1_0_0_2 + ( kCD_z * K_3_0_0_0_0_2 );

            // |f_2_1_0  f_3_0_0)_{i} = |g_3_1_0  d_2_0_0) + x_cd * |f_2_1_0  d_2_0_0)
            KET_f_f[ibra * 100 + 10] = K_3_1_0_2_0_0 + ( kCD_x * K_2_1_0_2_0_0 );

            // |f_2_1_0  f_2_1_0)_{i} = |g_2_2_0  d_2_0_0) + y_cd * |f_2_1_0  d_2_0_0)
            KET_f_f[ibra * 100 + 11] = K_2_2_0_2_0_0 + ( kCD_y * K_2_1_0_2_0_0 );

            // |f_2_1_0  f_2_0_1)_{i} = |g_2_1_1  d_2_0_0) + z_cd * |f_2_1_0  d_2_0_0)
            KET_f_f[ibra * 100 + 12] = K_2_1_1_2_0_0 + ( kCD_z * K_2_1_0_2_0_0 );

            // |f_2_1_0  f_1_2_0)_{i} = |g_3_1_0  d_0_2_0) + x_cd * |f_2_1_0  d_0_2_0)
            KET_f_f[ibra * 100 + 13] = K_3_1_0_0_2_0 + ( kCD_x * K_2_1_0_0_2_0 );

            // |f_2_1_0  f_1_1_1)_{i} = |g_2_1_1  d_1_1_0) + z_cd * |f_2_1_0  d_1_1_0)
            KET_f_f[ibra * 100 + 14] = K_2_1_1_1_1_0 + ( kCD_z * K_2_1_0_1_1_0 );

            // |f_2_1_0  f_1_0_2)_{i} = |g_3_1_0  d_0_0_2) + x_cd * |f_2_1_0  d_0_0_2)
            KET_f_f[ibra * 100 + 15] = K_3_1_0_0_0_2 + ( kCD_x * K_2_1_0_0_0_2 );

            // |f_2_1_0  f_0_3_0)_{i} = |g_2_2_0  d_0_2_0) + y_cd * |f_2_1_0  d_0_2_0)
            KET_f_f[ibra * 100 + 16] = K_2_2_0_0_2_0 + ( kCD_y * K_2_1_0_0_2_0 );

            // |f_2_1_0  f_0_2_1)_{i} = |g_2_1_1  d_0_2_0) + z_cd * |f_2_1_0  d_0_2_0)
            KET_f_f[ibra * 100 + 17] = K_2_1_1_0_2_0 + ( kCD_z * K_2_1_0_0_2_0 );

            // |f_2_1_0  f_0_1_2)_{i} = |g_2_2_0  d_0_0_2) + y_cd * |f_2_1_0  d_0_0_2)
            KET_f_f[ibra * 100 + 18] = K_2_2_0_0_0_2 + ( kCD_y * K_2_1_0_0_0_2 );

            // |f_2_1_0  f_0_0_3)_{i} = |g_2_1_1  d_0_0_2) + z_cd * |f_2_1_0  d_0_0_2)
            KET_f_f[ibra * 100 + 19] = K_2_1_1_0_0_2 + ( kCD_z * K_2_1_0_0_0_2 );

            // |f_2_0_1  f_3_0_0)_{i} = |g_3_0_1  d_2_0_0) + x_cd * |f_2_0_1  d_2_0_0)
            KET_f_f[ibra * 100 + 20] = K_3_0_1_2_0_0 + ( kCD_x * K_2_0_1_2_0_0 );

            // |f_2_0_1  f_2_1_0)_{i} = |g_2_1_1  d_2_0_0) + y_cd * |f_2_0_1  d_2_0_0)
            KET_f_f[ibra * 100 + 21] = K_2_1_1_2_0_0 + ( kCD_y * K_2_0_1_2_0_0 );

            // |f_2_0_1  f_2_0_1)_{i} = |g_2_0_2  d_2_0_0) + z_cd * |f_2_0_1  d_2_0_0)
            KET_f_f[ibra * 100 + 22] = K_2_0_2_2_0_0 + ( kCD_z * K_2_0_1_2_0_0 );

            // |f_2_0_1  f_1_2_0)_{i} = |g_3_0_1  d_0_2_0) + x_cd * |f_2_0_1  d_0_2_0)
            KET_f_f[ibra * 100 + 23] = K_3_0_1_0_2_0 + ( kCD_x * K_2_0_1_0_2_0 );

            // |f_2_0_1  f_1_1_1)_{i} = |g_2_0_2  d_1_1_0) + z_cd * |f_2_0_1  d_1_1_0)
            KET_f_f[ibra * 100 + 24] = K_2_0_2_1_1_0 + ( kCD_z * K_2_0_1_1_1_0 );

            // |f_2_0_1  f_1_0_2)_{i} = |g_3_0_1  d_0_0_2) + x_cd * |f_2_0_1  d_0_0_2)
            KET_f_f[ibra * 100 + 25] = K_3_0_1_0_0_2 + ( kCD_x * K_2_0_1_0_0_2 );

            // |f_2_0_1  f_0_3_0)_{i} = |g_2_1_1  d_0_2_0) + y_cd * |f_2_0_1  d_0_2_0)
            KET_f_f[ibra * 100 + 26] = K_2_1_1_0_2_0 + ( kCD_y * K_2_0_1_0_2_0 );

            // |f_2_0_1  f_0_2_1)_{i} = |g_2_0_2  d_0_2_0) + z_cd * |f_2_0_1  d_0_2_0)
            KET_f_f[ibra * 100 + 27] = K_2_0_2_0_2_0 + ( kCD_z * K_2_0_1_0_2_0 );

            // |f_2_0_1  f_0_1_2)_{i} = |g_2_1_1  d_0_0_2) + y_cd * |f_2_0_1  d_0_0_2)
            KET_f_f[ibra * 100 + 28] = K_2_1_1_0_0_2 + ( kCD_y * K_2_0_1_0_0_2 );

            // |f_2_0_1  f_0_0_3)_{i} = |g_2_0_2  d_0_0_2) + z_cd * |f_2_0_1  d_0_0_2)
            KET_f_f[ibra * 100 + 29] = K_2_0_2_0_0_2 + ( kCD_z * K_2_0_1_0_0_2 );

            // |f_1_2_0  f_3_0_0)_{i} = |g_2_2_0  d_2_0_0) + x_cd * |f_1_2_0  d_2_0_0)
            KET_f_f[ibra * 100 + 30] = K_2_2_0_2_0_0 + ( kCD_x * K_1_2_0_2_0_0 );

            // |f_1_2_0  f_2_1_0)_{i} = |g_1_3_0  d_2_0_0) + y_cd * |f_1_2_0  d_2_0_0)
            KET_f_f[ibra * 100 + 31] = K_1_3_0_2_0_0 + ( kCD_y * K_1_2_0_2_0_0 );

            // |f_1_2_0  f_2_0_1)_{i} = |g_1_2_1  d_2_0_0) + z_cd * |f_1_2_0  d_2_0_0)
            KET_f_f[ibra * 100 + 32] = K_1_2_1_2_0_0 + ( kCD_z * K_1_2_0_2_0_0 );

            // |f_1_2_0  f_1_2_0)_{i} = |g_2_2_0  d_0_2_0) + x_cd * |f_1_2_0  d_0_2_0)
            KET_f_f[ibra * 100 + 33] = K_2_2_0_0_2_0 + ( kCD_x * K_1_2_0_0_2_0 );

            // |f_1_2_0  f_1_1_1)_{i} = |g_1_2_1  d_1_1_0) + z_cd * |f_1_2_0  d_1_1_0)
            KET_f_f[ibra * 100 + 34] = K_1_2_1_1_1_0 + ( kCD_z * K_1_2_0_1_1_0 );

            // |f_1_2_0  f_1_0_2)_{i} = |g_2_2_0  d_0_0_2) + x_cd * |f_1_2_0  d_0_0_2)
            KET_f_f[ibra * 100 + 35] = K_2_2_0_0_0_2 + ( kCD_x * K_1_2_0_0_0_2 );

            // |f_1_2_0  f_0_3_0)_{i} = |g_1_3_0  d_0_2_0) + y_cd * |f_1_2_0  d_0_2_0)
            KET_f_f[ibra * 100 + 36] = K_1_3_0_0_2_0 + ( kCD_y * K_1_2_0_0_2_0 );

            // |f_1_2_0  f_0_2_1)_{i} = |g_1_2_1  d_0_2_0) + z_cd * |f_1_2_0  d_0_2_0)
            KET_f_f[ibra * 100 + 37] = K_1_2_1_0_2_0 + ( kCD_z * K_1_2_0_0_2_0 );

            // |f_1_2_0  f_0_1_2)_{i} = |g_1_3_0  d_0_0_2) + y_cd * |f_1_2_0  d_0_0_2)
            KET_f_f[ibra * 100 + 38] = K_1_3_0_0_0_2 + ( kCD_y * K_1_2_0_0_0_2 );

            // |f_1_2_0  f_0_0_3)_{i} = |g_1_2_1  d_0_0_2) + z_cd * |f_1_2_0  d_0_0_2)
            KET_f_f[ibra * 100 + 39] = K_1_2_1_0_0_2 + ( kCD_z * K_1_2_0_0_0_2 );

            // |f_1_1_1  f_3_0_0)_{i} = |g_2_1_1  d_2_0_0) + x_cd * |f_1_1_1  d_2_0_0)
            KET_f_f[ibra * 100 + 40] = K_2_1_1_2_0_0 + ( kCD_x * K_1_1_1_2_0_0 );

            // |f_1_1_1  f_2_1_0)_{i} = |g_1_2_1  d_2_0_0) + y_cd * |f_1_1_1  d_2_0_0)
            KET_f_f[ibra * 100 + 41] = K_1_2_1_2_0_0 + ( kCD_y * K_1_1_1_2_0_0 );

            // |f_1_1_1  f_2_0_1)_{i} = |g_1_1_2  d_2_0_0) + z_cd * |f_1_1_1  d_2_0_0)
            KET_f_f[ibra * 100 + 42] = K_1_1_2_2_0_0 + ( kCD_z * K_1_1_1_2_0_0 );

            // |f_1_1_1  f_1_2_0)_{i} = |g_2_1_1  d_0_2_0) + x_cd * |f_1_1_1  d_0_2_0)
            KET_f_f[ibra * 100 + 43] = K_2_1_1_0_2_0 + ( kCD_x * K_1_1_1_0_2_0 );

            // |f_1_1_1  f_1_1_1)_{i} = |g_1_1_2  d_1_1_0) + z_cd * |f_1_1_1  d_1_1_0)
            KET_f_f[ibra * 100 + 44] = K_1_1_2_1_1_0 + ( kCD_z * K_1_1_1_1_1_0 );

            // |f_1_1_1  f_1_0_2)_{i} = |g_2_1_1  d_0_0_2) + x_cd * |f_1_1_1  d_0_0_2)
            KET_f_f[ibra * 100 + 45] = K_2_1_1_0_0_2 + ( kCD_x * K_1_1_1_0_0_2 );

            // |f_1_1_1  f_0_3_0)_{i} = |g_1_2_1  d_0_2_0) + y_cd * |f_1_1_1  d_0_2_0)
            KET_f_f[ibra * 100 + 46] = K_1_2_1_0_2_0 + ( kCD_y * K_1_1_1_0_2_0 );

            // |f_1_1_1  f_0_2_1)_{i} = |g_1_1_2  d_0_2_0) + z_cd * |f_1_1_1  d_0_2_0)
            KET_f_f[ibra * 100 + 47] = K_1_1_2_0_2_0 + ( kCD_z * K_1_1_1_0_2_0 );

            // |f_1_1_1  f_0_1_2)_{i} = |g_1_2_1  d_0_0_2) + y_cd * |f_1_1_1  d_0_0_2)
            KET_f_f[ibra * 100 + 48] = K_1_2_1_0_0_2 + ( kCD_y * K_1_1_1_0_0_2 );

            // |f_1_1_1  f_0_0_3)_{i} = |g_1_1_2  d_0_0_2) + z_cd * |f_1_1_1  d_0_0_2)
            KET_f_f[ibra * 100 + 49] = K_1_1_2_0_0_2 + ( kCD_z * K_1_1_1_0_0_2 );

            // |f_1_0_2  f_3_0_0)_{i} = |g_2_0_2  d_2_0_0) + x_cd * |f_1_0_2  d_2_0_0)
            KET_f_f[ibra * 100 + 50] = K_2_0_2_2_0_0 + ( kCD_x * K_1_0_2_2_0_0 );

            // |f_1_0_2  f_2_1_0)_{i} = |g_1_1_2  d_2_0_0) + y_cd * |f_1_0_2  d_2_0_0)
            KET_f_f[ibra * 100 + 51] = K_1_1_2_2_0_0 + ( kCD_y * K_1_0_2_2_0_0 );

            // |f_1_0_2  f_2_0_1)_{i} = |g_1_0_3  d_2_0_0) + z_cd * |f_1_0_2  d_2_0_0)
            KET_f_f[ibra * 100 + 52] = K_1_0_3_2_0_0 + ( kCD_z * K_1_0_2_2_0_0 );

            // |f_1_0_2  f_1_2_0)_{i} = |g_2_0_2  d_0_2_0) + x_cd * |f_1_0_2  d_0_2_0)
            KET_f_f[ibra * 100 + 53] = K_2_0_2_0_2_0 + ( kCD_x * K_1_0_2_0_2_0 );

            // |f_1_0_2  f_1_1_1)_{i} = |g_1_0_3  d_1_1_0) + z_cd * |f_1_0_2  d_1_1_0)
            KET_f_f[ibra * 100 + 54] = K_1_0_3_1_1_0 + ( kCD_z * K_1_0_2_1_1_0 );

            // |f_1_0_2  f_1_0_2)_{i} = |g_2_0_2  d_0_0_2) + x_cd * |f_1_0_2  d_0_0_2)
            KET_f_f[ibra * 100 + 55] = K_2_0_2_0_0_2 + ( kCD_x * K_1_0_2_0_0_2 );

            // |f_1_0_2  f_0_3_0)_{i} = |g_1_1_2  d_0_2_0) + y_cd * |f_1_0_2  d_0_2_0)
            KET_f_f[ibra * 100 + 56] = K_1_1_2_0_2_0 + ( kCD_y * K_1_0_2_0_2_0 );

            // |f_1_0_2  f_0_2_1)_{i} = |g_1_0_3  d_0_2_0) + z_cd * |f_1_0_2  d_0_2_0)
            KET_f_f[ibra * 100 + 57] = K_1_0_3_0_2_0 + ( kCD_z * K_1_0_2_0_2_0 );

            // |f_1_0_2  f_0_1_2)_{i} = |g_1_1_2  d_0_0_2) + y_cd * |f_1_0_2  d_0_0_2)
            KET_f_f[ibra * 100 + 58] = K_1_1_2_0_0_2 + ( kCD_y * K_1_0_2_0_0_2 );

            // |f_1_0_2  f_0_0_3)_{i} = |g_1_0_3  d_0_0_2) + z_cd * |f_1_0_2  d_0_0_2)
            KET_f_f[ibra * 100 + 59] = K_1_0_3_0_0_2 + ( kCD_z * K_1_0_2_0_0_2 );

            // |f_0_3_0  f_3_0_0)_{i} = |g_1_3_0  d_2_0_0) + x_cd * |f_0_3_0  d_2_0_0)
            KET_f_f[ibra * 100 + 60] = K_1_3_0_2_0_0 + ( kCD_x * K_0_3_0_2_0_0 );

            // |f_0_3_0  f_2_1_0)_{i} = |g_0_4_0  d_2_0_0) + y_cd * |f_0_3_0  d_2_0_0)
            KET_f_f[ibra * 100 + 61] = K_0_4_0_2_0_0 + ( kCD_y * K_0_3_0_2_0_0 );

            // |f_0_3_0  f_2_0_1)_{i} = |g_0_3_1  d_2_0_0) + z_cd * |f_0_3_0  d_2_0_0)
            KET_f_f[ibra * 100 + 62] = K_0_3_1_2_0_0 + ( kCD_z * K_0_3_0_2_0_0 );

            // |f_0_3_0  f_1_2_0)_{i} = |g_1_3_0  d_0_2_0) + x_cd * |f_0_3_0  d_0_2_0)
            KET_f_f[ibra * 100 + 63] = K_1_3_0_0_2_0 + ( kCD_x * K_0_3_0_0_2_0 );

            // |f_0_3_0  f_1_1_1)_{i} = |g_0_3_1  d_1_1_0) + z_cd * |f_0_3_0  d_1_1_0)
            KET_f_f[ibra * 100 + 64] = K_0_3_1_1_1_0 + ( kCD_z * K_0_3_0_1_1_0 );

            // |f_0_3_0  f_1_0_2)_{i} = |g_1_3_0  d_0_0_2) + x_cd * |f_0_3_0  d_0_0_2)
            KET_f_f[ibra * 100 + 65] = K_1_3_0_0_0_2 + ( kCD_x * K_0_3_0_0_0_2 );

            // |f_0_3_0  f_0_3_0)_{i} = |g_0_4_0  d_0_2_0) + y_cd * |f_0_3_0  d_0_2_0)
            KET_f_f[ibra * 100 + 66] = K_0_4_0_0_2_0 + ( kCD_y * K_0_3_0_0_2_0 );

            // |f_0_3_0  f_0_2_1)_{i} = |g_0_3_1  d_0_2_0) + z_cd * |f_0_3_0  d_0_2_0)
            KET_f_f[ibra * 100 + 67] = K_0_3_1_0_2_0 + ( kCD_z * K_0_3_0_0_2_0 );

            // |f_0_3_0  f_0_1_2)_{i} = |g_0_4_0  d_0_0_2) + y_cd * |f_0_3_0  d_0_0_2)
            KET_f_f[ibra * 100 + 68] = K_0_4_0_0_0_2 + ( kCD_y * K_0_3_0_0_0_2 );

            // |f_0_3_0  f_0_0_3)_{i} = |g_0_3_1  d_0_0_2) + z_cd * |f_0_3_0  d_0_0_2)
            KET_f_f[ibra * 100 + 69] = K_0_3_1_0_0_2 + ( kCD_z * K_0_3_0_0_0_2 );

            // |f_0_2_1  f_3_0_0)_{i} = |g_1_2_1  d_2_0_0) + x_cd * |f_0_2_1  d_2_0_0)
            KET_f_f[ibra * 100 + 70] = K_1_2_1_2_0_0 + ( kCD_x * K_0_2_1_2_0_0 );

            // |f_0_2_1  f_2_1_0)_{i} = |g_0_3_1  d_2_0_0) + y_cd * |f_0_2_1  d_2_0_0)
            KET_f_f[ibra * 100 + 71] = K_0_3_1_2_0_0 + ( kCD_y * K_0_2_1_2_0_0 );

            // |f_0_2_1  f_2_0_1)_{i} = |g_0_2_2  d_2_0_0) + z_cd * |f_0_2_1  d_2_0_0)
            KET_f_f[ibra * 100 + 72] = K_0_2_2_2_0_0 + ( kCD_z * K_0_2_1_2_0_0 );

            // |f_0_2_1  f_1_2_0)_{i} = |g_1_2_1  d_0_2_0) + x_cd * |f_0_2_1  d_0_2_0)
            KET_f_f[ibra * 100 + 73] = K_1_2_1_0_2_0 + ( kCD_x * K_0_2_1_0_2_0 );

            // |f_0_2_1  f_1_1_1)_{i} = |g_0_2_2  d_1_1_0) + z_cd * |f_0_2_1  d_1_1_0)
            KET_f_f[ibra * 100 + 74] = K_0_2_2_1_1_0 + ( kCD_z * K_0_2_1_1_1_0 );

            // |f_0_2_1  f_1_0_2)_{i} = |g_1_2_1  d_0_0_2) + x_cd * |f_0_2_1  d_0_0_2)
            KET_f_f[ibra * 100 + 75] = K_1_2_1_0_0_2 + ( kCD_x * K_0_2_1_0_0_2 );

            // |f_0_2_1  f_0_3_0)_{i} = |g_0_3_1  d_0_2_0) + y_cd * |f_0_2_1  d_0_2_0)
            KET_f_f[ibra * 100 + 76] = K_0_3_1_0_2_0 + ( kCD_y * K_0_2_1_0_2_0 );

            // |f_0_2_1  f_0_2_1)_{i} = |g_0_2_2  d_0_2_0) + z_cd * |f_0_2_1  d_0_2_0)
            KET_f_f[ibra * 100 + 77] = K_0_2_2_0_2_0 + ( kCD_z * K_0_2_1_0_2_0 );

            // |f_0_2_1  f_0_1_2)_{i} = |g_0_3_1  d_0_0_2) + y_cd * |f_0_2_1  d_0_0_2)
            KET_f_f[ibra * 100 + 78] = K_0_3_1_0_0_2 + ( kCD_y * K_0_2_1_0_0_2 );

            // |f_0_2_1  f_0_0_3)_{i} = |g_0_2_2  d_0_0_2) + z_cd * |f_0_2_1  d_0_0_2)
            KET_f_f[ibra * 100 + 79] = K_0_2_2_0_0_2 + ( kCD_z * K_0_2_1_0_0_2 );

            // |f_0_1_2  f_3_0_0)_{i} = |g_1_1_2  d_2_0_0) + x_cd * |f_0_1_2  d_2_0_0)
            KET_f_f[ibra * 100 + 80] = K_1_1_2_2_0_0 + ( kCD_x * K_0_1_2_2_0_0 );

            // |f_0_1_2  f_2_1_0)_{i} = |g_0_2_2  d_2_0_0) + y_cd * |f_0_1_2  d_2_0_0)
            KET_f_f[ibra * 100 + 81] = K_0_2_2_2_0_0 + ( kCD_y * K_0_1_2_2_0_0 );

            // |f_0_1_2  f_2_0_1)_{i} = |g_0_1_3  d_2_0_0) + z_cd * |f_0_1_2  d_2_0_0)
            KET_f_f[ibra * 100 + 82] = K_0_1_3_2_0_0 + ( kCD_z * K_0_1_2_2_0_0 );

            // |f_0_1_2  f_1_2_0)_{i} = |g_1_1_2  d_0_2_0) + x_cd * |f_0_1_2  d_0_2_0)
            KET_f_f[ibra * 100 + 83] = K_1_1_2_0_2_0 + ( kCD_x * K_0_1_2_0_2_0 );

            // |f_0_1_2  f_1_1_1)_{i} = |g_0_1_3  d_1_1_0) + z_cd * |f_0_1_2  d_1_1_0)
            KET_f_f[ibra * 100 + 84] = K_0_1_3_1_1_0 + ( kCD_z * K_0_1_2_1_1_0 );

            // |f_0_1_2  f_1_0_2)_{i} = |g_1_1_2  d_0_0_2) + x_cd * |f_0_1_2  d_0_0_2)
            KET_f_f[ibra * 100 + 85] = K_1_1_2_0_0_2 + ( kCD_x * K_0_1_2_0_0_2 );

            // |f_0_1_2  f_0_3_0)_{i} = |g_0_2_2  d_0_2_0) + y_cd * |f_0_1_2  d_0_2_0)
            KET_f_f[ibra * 100 + 86] = K_0_2_2_0_2_0 + ( kCD_y * K_0_1_2_0_2_0 );

            // |f_0_1_2  f_0_2_1)_{i} = |g_0_1_3  d_0_2_0) + z_cd * |f_0_1_2  d_0_2_0)
            KET_f_f[ibra * 100 + 87] = K_0_1_3_0_2_0 + ( kCD_z * K_0_1_2_0_2_0 );

            // |f_0_1_2  f_0_1_2)_{i} = |g_0_2_2  d_0_0_2) + y_cd * |f_0_1_2  d_0_0_2)
            KET_f_f[ibra * 100 + 88] = K_0_2_2_0_0_2 + ( kCD_y * K_0_1_2_0_0_2 );

            // |f_0_1_2  f_0_0_3)_{i} = |g_0_1_3  d_0_0_2) + z_cd * |f_0_1_2  d_0_0_2)
            KET_f_f[ibra * 100 + 89] = K_0_1_3_0_0_2 + ( kCD_z * K_0_1_2_0_0_2 );

            // |f_0_0_3  f_3_0_0)_{i} = |g_1_0_3  d_2_0_0) + x_cd * |f_0_0_3  d_2_0_0)
            KET_f_f[ibra * 100 + 90] = K_1_0_3_2_0_0 + ( kCD_x * K_0_0_3_2_0_0 );

            // |f_0_0_3  f_2_1_0)_{i} = |g_0_1_3  d_2_0_0) + y_cd * |f_0_0_3  d_2_0_0)
            KET_f_f[ibra * 100 + 91] = K_0_1_3_2_0_0 + ( kCD_y * K_0_0_3_2_0_0 );

            // |f_0_0_3  f_2_0_1)_{i} = |g_0_0_4  d_2_0_0) + z_cd * |f_0_0_3  d_2_0_0)
            KET_f_f[ibra * 100 + 92] = K_0_0_4_2_0_0 + ( kCD_z * K_0_0_3_2_0_0 );

            // |f_0_0_3  f_1_2_0)_{i} = |g_1_0_3  d_0_2_0) + x_cd * |f_0_0_3  d_0_2_0)
            KET_f_f[ibra * 100 + 93] = K_1_0_3_0_2_0 + ( kCD_x * K_0_0_3_0_2_0 );

            // |f_0_0_3  f_1_1_1)_{i} = |g_0_0_4  d_1_1_0) + z_cd * |f_0_0_3  d_1_1_0)
            KET_f_f[ibra * 100 + 94] = K_0_0_4_1_1_0 + ( kCD_z * K_0_0_3_1_1_0 );

            // |f_0_0_3  f_1_0_2)_{i} = |g_1_0_3  d_0_0_2) + x_cd * |f_0_0_3  d_0_0_2)
            KET_f_f[ibra * 100 + 95] = K_1_0_3_0_0_2 + ( kCD_x * K_0_0_3_0_0_2 );

            // |f_0_0_3  f_0_3_0)_{i} = |g_0_1_3  d_0_2_0) + y_cd * |f_0_0_3  d_0_2_0)
            KET_f_f[ibra * 100 + 96] = K_0_1_3_0_2_0 + ( kCD_y * K_0_0_3_0_2_0 );

            // |f_0_0_3  f_0_2_1)_{i} = |g_0_0_4  d_0_2_0) + z_cd * |f_0_0_3  d_0_2_0)
            KET_f_f[ibra * 100 + 97] = K_0_0_4_0_2_0 + ( kCD_z * K_0_0_3_0_2_0 );

            // |f_0_0_3  f_0_1_2)_{i} = |g_0_1_3  d_0_0_2) + y_cd * |f_0_0_3  d_0_0_2)
            KET_f_f[ibra * 100 + 98] = K_0_1_3_0_0_2 + ( kCD_y * K_0_0_3_0_0_2 );

            // |f_0_0_3  f_0_0_3)_{i} = |g_0_0_4  d_0_0_2) + z_cd * |f_0_0_3  d_0_0_2)
            KET_f_f[ibra * 100 + 99] = K_0_0_4_0_0_2 + ( kCD_z * K_0_0_3_0_0_2 );

        }

}


