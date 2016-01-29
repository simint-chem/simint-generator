#include "test/simint.hpp"
#include "test/common.hpp"
#include "eri/eri.h"

// TODO - static initialization
static siminterifunc simintfuncs[MAXAM+1][MAXAM+1][MAXAM+1][MAXAM+1];

void Simint_Init(void)
{
    for(int i = 0; i <= MAXAM; i++)
    for(int j = 0; j <= MAXAM; j++)
    for(int k = 0; k <= MAXAM; k++)
    for(int l = 0; l <= MAXAM; l++)
    {
        simintfuncs[i][j][k][l] = siminteri_notyetimplemented;
    }


    simintfuncs[0][0][0][0] = eri_sharedwork_s_s_s_s;
    simintfuncs[0][0][0][1] = eri_sharedwork_s_s_s_p;
    simintfuncs[0][0][0][2] = eri_sharedwork_s_s_s_d;
    simintfuncs[0][0][0][3] = eri_sharedwork_s_s_s_f;
    simintfuncs[0][0][1][0] = eri_sharedwork_s_s_p_s;
    simintfuncs[0][0][1][1] = eri_sharedwork_s_s_p_p;
    simintfuncs[0][0][1][2] = eri_sharedwork_s_s_p_d;
    simintfuncs[0][0][1][3] = eri_sharedwork_s_s_p_f;
    simintfuncs[0][0][2][0] = eri_sharedwork_s_s_d_s;
    simintfuncs[0][0][2][1] = eri_sharedwork_s_s_d_p;
    simintfuncs[0][0][2][2] = eri_sharedwork_s_s_d_d;
    simintfuncs[0][0][2][3] = eri_sharedwork_s_s_d_f;
    simintfuncs[0][0][3][0] = eri_sharedwork_s_s_f_s;
    simintfuncs[0][0][3][1] = eri_sharedwork_s_s_f_p;
    simintfuncs[0][0][3][2] = eri_sharedwork_s_s_f_d;
    simintfuncs[0][0][3][3] = eri_sharedwork_s_s_f_f;
    simintfuncs[0][1][0][0] = eri_sharedwork_s_p_s_s;
    simintfuncs[0][1][0][1] = eri_sharedwork_s_p_s_p;
    simintfuncs[0][1][0][2] = eri_sharedwork_s_p_s_d;
    simintfuncs[0][1][0][3] = eri_sharedwork_s_p_s_f;
    simintfuncs[0][1][1][0] = eri_sharedwork_s_p_p_s;
    simintfuncs[0][1][1][1] = eri_sharedwork_s_p_p_p;
    simintfuncs[0][1][1][2] = eri_sharedwork_s_p_p_d;
    simintfuncs[0][1][1][3] = eri_sharedwork_s_p_p_f;
    simintfuncs[0][1][2][0] = eri_sharedwork_s_p_d_s;
    simintfuncs[0][1][2][1] = eri_sharedwork_s_p_d_p;
    simintfuncs[0][1][2][2] = eri_sharedwork_s_p_d_d;
    simintfuncs[0][1][2][3] = eri_sharedwork_s_p_d_f;
    simintfuncs[0][1][3][0] = eri_sharedwork_s_p_f_s;
    simintfuncs[0][1][3][1] = eri_sharedwork_s_p_f_p;
    simintfuncs[0][1][3][2] = eri_sharedwork_s_p_f_d;
    simintfuncs[0][1][3][3] = eri_sharedwork_s_p_f_f;
    simintfuncs[0][2][0][0] = eri_sharedwork_s_d_s_s;
    simintfuncs[0][2][0][1] = eri_sharedwork_s_d_s_p;
    simintfuncs[0][2][0][2] = eri_sharedwork_s_d_s_d;
    simintfuncs[0][2][0][3] = eri_sharedwork_s_d_s_f;
    simintfuncs[0][2][1][0] = eri_sharedwork_s_d_p_s;
    simintfuncs[0][2][1][1] = eri_sharedwork_s_d_p_p;
    simintfuncs[0][2][1][2] = eri_sharedwork_s_d_p_d;
    simintfuncs[0][2][1][3] = eri_sharedwork_s_d_p_f;
    simintfuncs[0][2][2][0] = eri_sharedwork_s_d_d_s;
    simintfuncs[0][2][2][1] = eri_sharedwork_s_d_d_p;
    simintfuncs[0][2][2][2] = eri_sharedwork_s_d_d_d;
    simintfuncs[0][2][2][3] = eri_sharedwork_s_d_d_f;
    simintfuncs[0][2][3][0] = eri_sharedwork_s_d_f_s;
    simintfuncs[0][2][3][1] = eri_sharedwork_s_d_f_p;
    simintfuncs[0][2][3][2] = eri_sharedwork_s_d_f_d;
    simintfuncs[0][2][3][3] = eri_sharedwork_s_d_f_f;
    simintfuncs[0][3][0][0] = eri_sharedwork_s_f_s_s;
    simintfuncs[0][3][0][1] = eri_sharedwork_s_f_s_p;
    simintfuncs[0][3][0][2] = eri_sharedwork_s_f_s_d;
    simintfuncs[0][3][0][3] = eri_sharedwork_s_f_s_f;
    simintfuncs[0][3][1][0] = eri_sharedwork_s_f_p_s;
    simintfuncs[0][3][1][1] = eri_sharedwork_s_f_p_p;
    simintfuncs[0][3][1][2] = eri_sharedwork_s_f_p_d;
    simintfuncs[0][3][1][3] = eri_sharedwork_s_f_p_f;
    simintfuncs[0][3][2][0] = eri_sharedwork_s_f_d_s;
    simintfuncs[0][3][2][1] = eri_sharedwork_s_f_d_p;
    simintfuncs[0][3][2][2] = eri_sharedwork_s_f_d_d;
    simintfuncs[0][3][2][3] = eri_sharedwork_s_f_d_f;
    simintfuncs[0][3][3][0] = eri_sharedwork_s_f_f_s;
    simintfuncs[0][3][3][1] = eri_sharedwork_s_f_f_p;
    simintfuncs[0][3][3][2] = eri_sharedwork_s_f_f_d;
    simintfuncs[0][3][3][3] = eri_sharedwork_s_f_f_f;
    simintfuncs[1][0][0][0] = eri_sharedwork_p_s_s_s;
    simintfuncs[1][0][0][1] = eri_sharedwork_p_s_s_p;
    simintfuncs[1][0][0][2] = eri_sharedwork_p_s_s_d;
    simintfuncs[1][0][0][3] = eri_sharedwork_p_s_s_f;
    simintfuncs[1][0][1][0] = eri_sharedwork_p_s_p_s;
    simintfuncs[1][0][1][1] = eri_sharedwork_p_s_p_p;
    simintfuncs[1][0][1][2] = eri_sharedwork_p_s_p_d;
    simintfuncs[1][0][1][3] = eri_sharedwork_p_s_p_f;
    simintfuncs[1][0][2][0] = eri_sharedwork_p_s_d_s;
    simintfuncs[1][0][2][1] = eri_sharedwork_p_s_d_p;
    simintfuncs[1][0][2][2] = eri_sharedwork_p_s_d_d;
    simintfuncs[1][0][2][3] = eri_sharedwork_p_s_d_f;
    simintfuncs[1][0][3][0] = eri_sharedwork_p_s_f_s;
    simintfuncs[1][0][3][1] = eri_sharedwork_p_s_f_p;
    simintfuncs[1][0][3][2] = eri_sharedwork_p_s_f_d;
    simintfuncs[1][0][3][3] = eri_sharedwork_p_s_f_f;
    simintfuncs[1][1][0][0] = eri_sharedwork_p_p_s_s;
    simintfuncs[1][1][0][1] = eri_sharedwork_p_p_s_p;
    simintfuncs[1][1][0][2] = eri_sharedwork_p_p_s_d;
    simintfuncs[1][1][0][3] = eri_sharedwork_p_p_s_f;
    simintfuncs[1][1][1][0] = eri_sharedwork_p_p_p_s;
    simintfuncs[1][1][1][1] = eri_sharedwork_p_p_p_p;
    simintfuncs[1][1][1][2] = eri_sharedwork_p_p_p_d;
    simintfuncs[1][1][1][3] = eri_sharedwork_p_p_p_f;
    simintfuncs[1][1][2][0] = eri_sharedwork_p_p_d_s;
    simintfuncs[1][1][2][1] = eri_sharedwork_p_p_d_p;
    simintfuncs[1][1][2][2] = eri_sharedwork_p_p_d_d;
    simintfuncs[1][1][2][3] = eri_sharedwork_p_p_d_f;
    simintfuncs[1][1][3][0] = eri_sharedwork_p_p_f_s;
    simintfuncs[1][1][3][1] = eri_sharedwork_p_p_f_p;
    simintfuncs[1][1][3][2] = eri_sharedwork_p_p_f_d;
    simintfuncs[1][1][3][3] = eri_sharedwork_p_p_f_f;
    simintfuncs[1][2][0][0] = eri_sharedwork_p_d_s_s;
    simintfuncs[1][2][0][1] = eri_sharedwork_p_d_s_p;
    simintfuncs[1][2][0][2] = eri_sharedwork_p_d_s_d;
    simintfuncs[1][2][0][3] = eri_sharedwork_p_d_s_f;
    simintfuncs[1][2][1][0] = eri_sharedwork_p_d_p_s;
    simintfuncs[1][2][1][1] = eri_sharedwork_p_d_p_p;
    simintfuncs[1][2][1][2] = eri_sharedwork_p_d_p_d;
    simintfuncs[1][2][1][3] = eri_sharedwork_p_d_p_f;
    simintfuncs[1][2][2][0] = eri_sharedwork_p_d_d_s;
    simintfuncs[1][2][2][1] = eri_sharedwork_p_d_d_p;
    simintfuncs[1][2][2][2] = eri_sharedwork_p_d_d_d;
    simintfuncs[1][2][2][3] = eri_sharedwork_p_d_d_f;
    simintfuncs[1][2][3][0] = eri_sharedwork_p_d_f_s;
    simintfuncs[1][2][3][1] = eri_sharedwork_p_d_f_p;
    simintfuncs[1][2][3][2] = eri_sharedwork_p_d_f_d;
    simintfuncs[1][2][3][3] = eri_sharedwork_p_d_f_f;
    simintfuncs[1][3][0][0] = eri_sharedwork_p_f_s_s;
    simintfuncs[1][3][0][1] = eri_sharedwork_p_f_s_p;
    simintfuncs[1][3][0][2] = eri_sharedwork_p_f_s_d;
    simintfuncs[1][3][0][3] = eri_sharedwork_p_f_s_f;
    simintfuncs[1][3][1][0] = eri_sharedwork_p_f_p_s;
    simintfuncs[1][3][1][1] = eri_sharedwork_p_f_p_p;
    simintfuncs[1][3][1][2] = eri_sharedwork_p_f_p_d;
    simintfuncs[1][3][1][3] = eri_sharedwork_p_f_p_f;
    simintfuncs[1][3][2][0] = eri_sharedwork_p_f_d_s;
    simintfuncs[1][3][2][1] = eri_sharedwork_p_f_d_p;
    simintfuncs[1][3][2][2] = eri_sharedwork_p_f_d_d;
    simintfuncs[1][3][2][3] = eri_sharedwork_p_f_d_f;
    simintfuncs[1][3][3][0] = eri_sharedwork_p_f_f_s;
    simintfuncs[1][3][3][1] = eri_sharedwork_p_f_f_p;
    simintfuncs[1][3][3][2] = eri_sharedwork_p_f_f_d;
    simintfuncs[1][3][3][3] = eri_sharedwork_p_f_f_f;
    simintfuncs[2][0][0][0] = eri_sharedwork_d_s_s_s;
    simintfuncs[2][0][0][1] = eri_sharedwork_d_s_s_p;
    simintfuncs[2][0][0][2] = eri_sharedwork_d_s_s_d;
    simintfuncs[2][0][0][3] = eri_sharedwork_d_s_s_f;
    simintfuncs[2][0][1][0] = eri_sharedwork_d_s_p_s;
    simintfuncs[2][0][1][1] = eri_sharedwork_d_s_p_p;
    simintfuncs[2][0][1][2] = eri_sharedwork_d_s_p_d;
    simintfuncs[2][0][1][3] = eri_sharedwork_d_s_p_f;
    simintfuncs[2][0][2][0] = eri_sharedwork_d_s_d_s;
    simintfuncs[2][0][2][1] = eri_sharedwork_d_s_d_p;
    simintfuncs[2][0][2][2] = eri_sharedwork_d_s_d_d;
    simintfuncs[2][0][2][3] = eri_sharedwork_d_s_d_f;
    simintfuncs[2][0][3][0] = eri_sharedwork_d_s_f_s;
    simintfuncs[2][0][3][1] = eri_sharedwork_d_s_f_p;
    simintfuncs[2][0][3][2] = eri_sharedwork_d_s_f_d;
    simintfuncs[2][0][3][3] = eri_sharedwork_d_s_f_f;
    simintfuncs[2][1][0][0] = eri_sharedwork_d_p_s_s;
    simintfuncs[2][1][0][1] = eri_sharedwork_d_p_s_p;
    simintfuncs[2][1][0][2] = eri_sharedwork_d_p_s_d;
    simintfuncs[2][1][0][3] = eri_sharedwork_d_p_s_f;
    simintfuncs[2][1][1][0] = eri_sharedwork_d_p_p_s;
    simintfuncs[2][1][1][1] = eri_sharedwork_d_p_p_p;
    simintfuncs[2][1][1][2] = eri_sharedwork_d_p_p_d;
    simintfuncs[2][1][1][3] = eri_sharedwork_d_p_p_f;
    simintfuncs[2][1][2][0] = eri_sharedwork_d_p_d_s;
    simintfuncs[2][1][2][1] = eri_sharedwork_d_p_d_p;
    simintfuncs[2][1][2][2] = eri_sharedwork_d_p_d_d;
    simintfuncs[2][1][2][3] = eri_sharedwork_d_p_d_f;
    simintfuncs[2][1][3][0] = eri_sharedwork_d_p_f_s;
    simintfuncs[2][1][3][1] = eri_sharedwork_d_p_f_p;
    simintfuncs[2][1][3][2] = eri_sharedwork_d_p_f_d;
    simintfuncs[2][1][3][3] = eri_sharedwork_d_p_f_f;
    simintfuncs[2][2][0][0] = eri_sharedwork_d_d_s_s;
    simintfuncs[2][2][0][1] = eri_sharedwork_d_d_s_p;
    simintfuncs[2][2][0][2] = eri_sharedwork_d_d_s_d;
    simintfuncs[2][2][0][3] = eri_sharedwork_d_d_s_f;
    simintfuncs[2][2][1][0] = eri_sharedwork_d_d_p_s;
    simintfuncs[2][2][1][1] = eri_sharedwork_d_d_p_p;
    simintfuncs[2][2][1][2] = eri_sharedwork_d_d_p_d;
    simintfuncs[2][2][1][3] = eri_sharedwork_d_d_p_f;
    simintfuncs[2][2][2][0] = eri_sharedwork_d_d_d_s;
    simintfuncs[2][2][2][1] = eri_sharedwork_d_d_d_p;
    simintfuncs[2][2][2][2] = eri_sharedwork_d_d_d_d;
    simintfuncs[2][2][2][3] = eri_sharedwork_d_d_d_f;
    simintfuncs[2][2][3][0] = eri_sharedwork_d_d_f_s;
    simintfuncs[2][2][3][1] = eri_sharedwork_d_d_f_p;
    simintfuncs[2][2][3][2] = eri_sharedwork_d_d_f_d;
    simintfuncs[2][2][3][3] = eri_sharedwork_d_d_f_f;
    simintfuncs[2][3][0][0] = eri_sharedwork_d_f_s_s;
    simintfuncs[2][3][0][1] = eri_sharedwork_d_f_s_p;
    simintfuncs[2][3][0][2] = eri_sharedwork_d_f_s_d;
    simintfuncs[2][3][0][3] = eri_sharedwork_d_f_s_f;
    simintfuncs[2][3][1][0] = eri_sharedwork_d_f_p_s;
    simintfuncs[2][3][1][1] = eri_sharedwork_d_f_p_p;
    simintfuncs[2][3][1][2] = eri_sharedwork_d_f_p_d;
    simintfuncs[2][3][1][3] = eri_sharedwork_d_f_p_f;
    simintfuncs[2][3][2][0] = eri_sharedwork_d_f_d_s;
    simintfuncs[2][3][2][1] = eri_sharedwork_d_f_d_p;
    simintfuncs[2][3][2][2] = eri_sharedwork_d_f_d_d;
    simintfuncs[2][3][2][3] = eri_sharedwork_d_f_d_f;
    simintfuncs[2][3][3][0] = eri_sharedwork_d_f_f_s;
    simintfuncs[2][3][3][1] = eri_sharedwork_d_f_f_p;
    simintfuncs[2][3][3][2] = eri_sharedwork_d_f_f_d;
    simintfuncs[2][3][3][3] = eri_sharedwork_d_f_f_f;
    simintfuncs[3][0][0][0] = eri_sharedwork_f_s_s_s;
    simintfuncs[3][0][0][1] = eri_sharedwork_f_s_s_p;
    simintfuncs[3][0][0][2] = eri_sharedwork_f_s_s_d;
    simintfuncs[3][0][0][3] = eri_sharedwork_f_s_s_f;
    simintfuncs[3][0][1][0] = eri_sharedwork_f_s_p_s;
    simintfuncs[3][0][1][1] = eri_sharedwork_f_s_p_p;
    simintfuncs[3][0][1][2] = eri_sharedwork_f_s_p_d;
    simintfuncs[3][0][1][3] = eri_sharedwork_f_s_p_f;
    simintfuncs[3][0][2][0] = eri_sharedwork_f_s_d_s;
    simintfuncs[3][0][2][1] = eri_sharedwork_f_s_d_p;
    simintfuncs[3][0][2][2] = eri_sharedwork_f_s_d_d;
    simintfuncs[3][0][2][3] = eri_sharedwork_f_s_d_f;
    simintfuncs[3][0][3][0] = eri_sharedwork_f_s_f_s;
    simintfuncs[3][0][3][1] = eri_sharedwork_f_s_f_p;
    simintfuncs[3][0][3][2] = eri_sharedwork_f_s_f_d;
    simintfuncs[3][0][3][3] = eri_sharedwork_f_s_f_f;
    simintfuncs[3][1][0][0] = eri_sharedwork_f_p_s_s;
    simintfuncs[3][1][0][1] = eri_sharedwork_f_p_s_p;
    simintfuncs[3][1][0][2] = eri_sharedwork_f_p_s_d;
    simintfuncs[3][1][0][3] = eri_sharedwork_f_p_s_f;
    simintfuncs[3][1][1][0] = eri_sharedwork_f_p_p_s;
    simintfuncs[3][1][1][1] = eri_sharedwork_f_p_p_p;
    simintfuncs[3][1][1][2] = eri_sharedwork_f_p_p_d;
    simintfuncs[3][1][1][3] = eri_sharedwork_f_p_p_f;
    simintfuncs[3][1][2][0] = eri_sharedwork_f_p_d_s;
    simintfuncs[3][1][2][1] = eri_sharedwork_f_p_d_p;
    simintfuncs[3][1][2][2] = eri_sharedwork_f_p_d_d;
    simintfuncs[3][1][2][3] = eri_sharedwork_f_p_d_f;
    simintfuncs[3][1][3][0] = eri_sharedwork_f_p_f_s;
    simintfuncs[3][1][3][1] = eri_sharedwork_f_p_f_p;
    simintfuncs[3][1][3][2] = eri_sharedwork_f_p_f_d;
    simintfuncs[3][1][3][3] = eri_sharedwork_f_p_f_f;
    simintfuncs[3][2][0][0] = eri_sharedwork_f_d_s_s;
    simintfuncs[3][2][0][1] = eri_sharedwork_f_d_s_p;
    simintfuncs[3][2][0][2] = eri_sharedwork_f_d_s_d;
    simintfuncs[3][2][0][3] = eri_sharedwork_f_d_s_f;
    simintfuncs[3][2][1][0] = eri_sharedwork_f_d_p_s;
    simintfuncs[3][2][1][1] = eri_sharedwork_f_d_p_p;
    simintfuncs[3][2][1][2] = eri_sharedwork_f_d_p_d;
    simintfuncs[3][2][1][3] = eri_sharedwork_f_d_p_f;
    simintfuncs[3][2][2][0] = eri_sharedwork_f_d_d_s;
    simintfuncs[3][2][2][1] = eri_sharedwork_f_d_d_p;
    simintfuncs[3][2][2][2] = eri_sharedwork_f_d_d_d;
    simintfuncs[3][2][2][3] = eri_sharedwork_f_d_d_f;
    simintfuncs[3][2][3][0] = eri_sharedwork_f_d_f_s;
    simintfuncs[3][2][3][1] = eri_sharedwork_f_d_f_p;
    simintfuncs[3][2][3][2] = eri_sharedwork_f_d_f_d;
    simintfuncs[3][2][3][3] = eri_sharedwork_f_d_f_f;
    simintfuncs[3][3][0][0] = eri_sharedwork_f_f_s_s;
    simintfuncs[3][3][0][1] = eri_sharedwork_f_f_s_p;
    simintfuncs[3][3][0][2] = eri_sharedwork_f_f_s_d;
    simintfuncs[3][3][0][3] = eri_sharedwork_f_f_s_f;
    simintfuncs[3][3][1][0] = eri_sharedwork_f_f_p_s;
    simintfuncs[3][3][1][1] = eri_sharedwork_f_f_p_p;
    simintfuncs[3][3][1][2] = eri_sharedwork_f_f_p_d;
    simintfuncs[3][3][1][3] = eri_sharedwork_f_f_p_f;
    simintfuncs[3][3][2][0] = eri_sharedwork_f_f_d_s;
    simintfuncs[3][3][2][1] = eri_sharedwork_f_f_d_p;
    simintfuncs[3][3][2][2] = eri_sharedwork_f_f_d_d;
    simintfuncs[3][3][2][3] = eri_sharedwork_f_f_d_f;
    simintfuncs[3][3][3][0] = eri_sharedwork_f_f_f_s;
    simintfuncs[3][3][3][1] = eri_sharedwork_f_f_f_p;
    simintfuncs[3][3][3][2] = eri_sharedwork_f_f_f_d;
    simintfuncs[3][3][3][3] = eri_sharedwork_f_f_f_f;
}


int siminteri_notyetimplemented(struct multishell_pair const P,
                                struct multishell_pair const Q,
                                double * const restrict /*dummy*/,
                                double * const restrict /*dummy*/)
{
    printf("****************************\n");
    printf("*** NOT YET IMPLEMENTED! ***\n");
    printf("***  ( %2d %2d | %2d %2d )   ***\n", P.am1, P.am2, Q.am1, Q.am2);
    printf("****************************\n");
    exit(1);
    return 0;
}



TimerType Simint_Integral(struct multishell_pair const P, struct multishell_pair const Q,
                          double * const restrict contwork,
                          double * const restrict integrals)
{
    TimerType ticks0, ticks1;

    CLOCK(ticks0);
    simintfuncs[P.am1][P.am2][Q.am1][Q.am2](P, Q, contwork, integrals);
    CLOCK(ticks1);
    return ticks1 - ticks0;
}

