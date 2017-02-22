/*
 Generated with:
   ../../../python/gen_fill.py -l 2 -d 1 ostei_deriv1_init.c
*/



#include "simint/ostei/ostei.h"
#include "simint/ostei/ostei_init.h"
#include "simint/ostei/ostei_config.h"
#include "simint/ostei/gen/ostei_deriv1_generated.h"


// Stores pointers to the ostei functions
#define AMSIZE   SIMINT_OSTEI_MAXAM+1
#define DERSIZE  SIMINT_OSTEI_MAXDER+1
extern simint_osteifunc simint_osteifunc_array[DERSIZE][AMSIZE][AMSIZE][AMSIZE][AMSIZE];



void simint_ostei_deriv1_finalize(void)
{
    // nothing to do
}


void simint_ostei_deriv1_init(void)
{
    #if SIMINT_OSTEI_DERIV1_MAXAM >= 0
    simint_osteifunc_array[1][0][0][0][0] = ostei_deriv1_s_s_s_s;
    #endif

    #if SIMINT_OSTEI_DERIV1_MAXAM >= 1
    simint_osteifunc_array[1][0][0][0][1] = ostei_deriv1_s_s_s_p;
    simint_osteifunc_array[1][0][0][1][0] = ostei_deriv1_s_s_p_s;
    simint_osteifunc_array[1][0][0][1][1] = ostei_deriv1_s_s_p_p;
    simint_osteifunc_array[1][0][1][0][0] = ostei_deriv1_s_p_s_s;
    simint_osteifunc_array[1][0][1][0][1] = ostei_deriv1_s_p_s_p;
    simint_osteifunc_array[1][0][1][1][0] = ostei_deriv1_s_p_p_s;
    simint_osteifunc_array[1][0][1][1][1] = ostei_deriv1_s_p_p_p;
    simint_osteifunc_array[1][1][0][0][0] = ostei_deriv1_p_s_s_s;
    simint_osteifunc_array[1][1][0][0][1] = ostei_deriv1_p_s_s_p;
    simint_osteifunc_array[1][1][0][1][0] = ostei_deriv1_p_s_p_s;
    simint_osteifunc_array[1][1][0][1][1] = ostei_deriv1_p_s_p_p;
    simint_osteifunc_array[1][1][1][0][0] = ostei_deriv1_p_p_s_s;
    simint_osteifunc_array[1][1][1][0][1] = ostei_deriv1_p_p_s_p;
    simint_osteifunc_array[1][1][1][1][0] = ostei_deriv1_p_p_p_s;
    simint_osteifunc_array[1][1][1][1][1] = ostei_deriv1_p_p_p_p;
    #endif

    #if SIMINT_OSTEI_DERIV1_MAXAM >= 2
    simint_osteifunc_array[1][0][0][0][2] = ostei_deriv1_s_s_s_d;
    simint_osteifunc_array[1][0][0][1][2] = ostei_deriv1_s_s_p_d;
    simint_osteifunc_array[1][0][0][2][0] = ostei_deriv1_s_s_d_s;
    simint_osteifunc_array[1][0][0][2][1] = ostei_deriv1_s_s_d_p;
    simint_osteifunc_array[1][0][0][2][2] = ostei_deriv1_s_s_d_d;
    simint_osteifunc_array[1][0][1][0][2] = ostei_deriv1_s_p_s_d;
    simint_osteifunc_array[1][0][1][1][2] = ostei_deriv1_s_p_p_d;
    simint_osteifunc_array[1][0][1][2][0] = ostei_deriv1_s_p_d_s;
    simint_osteifunc_array[1][0][1][2][1] = ostei_deriv1_s_p_d_p;
    simint_osteifunc_array[1][0][1][2][2] = ostei_deriv1_s_p_d_d;
    simint_osteifunc_array[1][0][2][0][0] = ostei_deriv1_s_d_s_s;
    simint_osteifunc_array[1][0][2][0][1] = ostei_deriv1_s_d_s_p;
    simint_osteifunc_array[1][0][2][0][2] = ostei_deriv1_s_d_s_d;
    simint_osteifunc_array[1][0][2][1][0] = ostei_deriv1_s_d_p_s;
    simint_osteifunc_array[1][0][2][1][1] = ostei_deriv1_s_d_p_p;
    simint_osteifunc_array[1][0][2][1][2] = ostei_deriv1_s_d_p_d;
    simint_osteifunc_array[1][0][2][2][0] = ostei_deriv1_s_d_d_s;
    simint_osteifunc_array[1][0][2][2][1] = ostei_deriv1_s_d_d_p;
    simint_osteifunc_array[1][0][2][2][2] = ostei_deriv1_s_d_d_d;
    simint_osteifunc_array[1][1][0][0][2] = ostei_deriv1_p_s_s_d;
    simint_osteifunc_array[1][1][0][1][2] = ostei_deriv1_p_s_p_d;
    simint_osteifunc_array[1][1][0][2][0] = ostei_deriv1_p_s_d_s;
    simint_osteifunc_array[1][1][0][2][1] = ostei_deriv1_p_s_d_p;
    simint_osteifunc_array[1][1][0][2][2] = ostei_deriv1_p_s_d_d;
    simint_osteifunc_array[1][1][1][0][2] = ostei_deriv1_p_p_s_d;
    simint_osteifunc_array[1][1][1][1][2] = ostei_deriv1_p_p_p_d;
    simint_osteifunc_array[1][1][1][2][0] = ostei_deriv1_p_p_d_s;
    simint_osteifunc_array[1][1][1][2][1] = ostei_deriv1_p_p_d_p;
    simint_osteifunc_array[1][1][1][2][2] = ostei_deriv1_p_p_d_d;
    simint_osteifunc_array[1][1][2][0][0] = ostei_deriv1_p_d_s_s;
    simint_osteifunc_array[1][1][2][0][1] = ostei_deriv1_p_d_s_p;
    simint_osteifunc_array[1][1][2][0][2] = ostei_deriv1_p_d_s_d;
    simint_osteifunc_array[1][1][2][1][0] = ostei_deriv1_p_d_p_s;
    simint_osteifunc_array[1][1][2][1][1] = ostei_deriv1_p_d_p_p;
    simint_osteifunc_array[1][1][2][1][2] = ostei_deriv1_p_d_p_d;
    simint_osteifunc_array[1][1][2][2][0] = ostei_deriv1_p_d_d_s;
    simint_osteifunc_array[1][1][2][2][1] = ostei_deriv1_p_d_d_p;
    simint_osteifunc_array[1][1][2][2][2] = ostei_deriv1_p_d_d_d;
    simint_osteifunc_array[1][2][0][0][0] = ostei_deriv1_d_s_s_s;
    simint_osteifunc_array[1][2][0][0][1] = ostei_deriv1_d_s_s_p;
    simint_osteifunc_array[1][2][0][0][2] = ostei_deriv1_d_s_s_d;
    simint_osteifunc_array[1][2][0][1][0] = ostei_deriv1_d_s_p_s;
    simint_osteifunc_array[1][2][0][1][1] = ostei_deriv1_d_s_p_p;
    simint_osteifunc_array[1][2][0][1][2] = ostei_deriv1_d_s_p_d;
    simint_osteifunc_array[1][2][0][2][0] = ostei_deriv1_d_s_d_s;
    simint_osteifunc_array[1][2][0][2][1] = ostei_deriv1_d_s_d_p;
    simint_osteifunc_array[1][2][0][2][2] = ostei_deriv1_d_s_d_d;
    simint_osteifunc_array[1][2][1][0][0] = ostei_deriv1_d_p_s_s;
    simint_osteifunc_array[1][2][1][0][1] = ostei_deriv1_d_p_s_p;
    simint_osteifunc_array[1][2][1][0][2] = ostei_deriv1_d_p_s_d;
    simint_osteifunc_array[1][2][1][1][0] = ostei_deriv1_d_p_p_s;
    simint_osteifunc_array[1][2][1][1][1] = ostei_deriv1_d_p_p_p;
    simint_osteifunc_array[1][2][1][1][2] = ostei_deriv1_d_p_p_d;
    simint_osteifunc_array[1][2][1][2][0] = ostei_deriv1_d_p_d_s;
    simint_osteifunc_array[1][2][1][2][1] = ostei_deriv1_d_p_d_p;
    simint_osteifunc_array[1][2][1][2][2] = ostei_deriv1_d_p_d_d;
    simint_osteifunc_array[1][2][2][0][0] = ostei_deriv1_d_d_s_s;
    simint_osteifunc_array[1][2][2][0][1] = ostei_deriv1_d_d_s_p;
    simint_osteifunc_array[1][2][2][0][2] = ostei_deriv1_d_d_s_d;
    simint_osteifunc_array[1][2][2][1][0] = ostei_deriv1_d_d_p_s;
    simint_osteifunc_array[1][2][2][1][1] = ostei_deriv1_d_d_p_p;
    simint_osteifunc_array[1][2][2][1][2] = ostei_deriv1_d_d_p_d;
    simint_osteifunc_array[1][2][2][2][0] = ostei_deriv1_d_d_d_s;
    simint_osteifunc_array[1][2][2][2][1] = ostei_deriv1_d_d_d_p;
    simint_osteifunc_array[1][2][2][2][2] = ostei_deriv1_d_d_d_d;
    #endif

}

