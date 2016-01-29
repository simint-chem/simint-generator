#ifndef SIMINT_TEST_SIMINT_HPP
#define SIMINT_TEST_SIMINT_HPP

#include <utility>

#include "eri/eri.h"
#include "test/timer.h"


// Function pointer typedef
typedef int (*siminterifunc)(struct multishell_pair const, struct multishell_pair const, double * const restrict, double * const restrict);


// Setting up function pointers and calculating integrals
// using my code
void Simint_Init(void);

int siminteri_notyetimplemented(struct multishell_pair const P,
                          struct multishell_pair const Q,
                          double * const restrict /*dummy*/,
                          double * const restrict /*dummy*/);


TimerType Simint_Integral(struct multishell_pair const P, struct multishell_pair const Q,
                          double * const restrict contwork, double * const restrict integrals);


#endif

