#ifndef SIMINT_TEST_SIMINT_HPP
#define SIMINT_TEST_SIMINT_HPP

#include <utility>

#include "simint/eri/eri.h"
#include "test/Timer.h"

void Simint_Init(void);

void Simint_Finalize(void);

TimerType Simint_Integral(struct multishell_pair const P, struct multishell_pair const Q,
                          double * const restrict contwork, double * const restrict integrals);


#endif

