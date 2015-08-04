#ifndef SIMINT_TIMER_HPP
#define SIMINT_TIMER_HPP

#include <utility>

typedef unsigned long long TimerType;


#define CLOCK(ticks) do {                                 \
    volatile unsigned int a__, d__;                              \
    __asm__ __volatile__("rdtsc" : "=a" (a__), "=d" (d__) : );   \
    (ticks) = ((TimerType) a__)|(((TimerType)d__)<<32); \
  } while(0)


#endif
