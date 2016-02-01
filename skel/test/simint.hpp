#ifndef SIMINT_TEST_SIMINT_HPP
#define SIMINT_TEST_SIMINT_HPP

#include <utility>

#include "eri/eri.h"
#include "test/timer.h"


// Setting up function pointers and calculating integrals
// using my code
void Simint_Init(void);

TimerType Simint_Integral(struct multishell_pair const P, struct multishell_pair const Q,
                          double * const restrict contwork, double * const restrict integrals);


#endif

