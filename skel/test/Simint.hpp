#pragma once

#include <utility>

#include "simint/eri/eri.h"
#include "test/Timer.h"

void Simint_Init(void);

void Simint_Finalize(void);

TimerType Simint_Integral(struct simint_multi_shellpair const * P,
                          struct simint_multi_shellpair const * Q,
                          double * const restrict contwork, double * const restrict integrals);

