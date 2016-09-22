#include "test/Simint.hpp"
#include "test/Common.hpp"
#include "simint/eri/eri.h"
#include "simint/simint_init.h"

void Simint_Init(void)
{
    simint_init();
}

void Simint_Finalize(void)
{
    simint_finalize();
}

TimerType Simint_Integral(struct simint_multi_shellpair const * P,
                          struct simint_multi_shellpair const * Q,
                          double * const restrict contwork,
                          double * const restrict integrals)
{
    TimerType ticks0, ticks1;

    CLOCK(ticks0);
    simint_compute_eri_sharedwork(P, Q, contwork, integrals);
    CLOCK(ticks1);
    return ticks1 - ticks0;
}

