#ifndef HRR_H
#define HRR_H

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
                 );


    //////////////////////////////////////////////
    // KET: ( p p |
    // Steps: 9
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET_p_p(
                  double * const restrict KET_p_s,
                  double * const restrict KET_p_p,
                  double * const restrict KET_d_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 );


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
                 );


    //////////////////////////////////////////////
    // KET: ( d p |
    // Steps: 18
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET_d_p(
                  double * const restrict KET_d_s,
                  double * const restrict KET_d_p,
                  double * const restrict KET_f_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 );


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
                 );


    //////////////////////////////////////////////
    // KET: ( d d |
    // Steps: 79
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET_d_d(
                  double * const restrict KET_d_s,
                  double * const restrict KET_d_d,
                  double * const restrict KET_f_s,
                  double * const restrict KET_g_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 );


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
                 );


    //////////////////////////////////////////////
    // KET: ( f p |
    // Steps: 30
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET_f_p(
                  double * const restrict KET_f_s,
                  double * const restrict KET_f_p,
                  double * const restrict KET_g_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 );


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
                 );


    //////////////////////////////////////////////
    // KET: ( f d |
    // Steps: 129
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET_f_d(
                  double * const restrict KET_f_s,
                  double * const restrict KET_f_d,
                  double * const restrict KET_g_s,
                  double * const restrict KET_h_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 );


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
                 );


    //////////////////////////////////////////////
    // KET: ( f f |
    // Steps: 319
    //////////////////////////////////////////////

#pragma omp declare simd simdlen(SIMD_LEN)
void HRR_KET_f_f(
                  double * const restrict KET_f_s,
                  double * const restrict KET_f_f,
                  double * const restrict KET_g_s,
                  double * const restrict KET_h_s,
                  double * const restrict KET_i_s,
                  const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra
                 );


#endif
