#ifndef SIMINT_TIMER_H
#define SIMINT_TIMER_H

#ifdef __cplusplus
extern "C" {
#endif


typedef unsigned long long TimerType;


#define CLOCK(ticks) do {                                 \
    volatile unsigned int a__, d__;                              \
    __asm__ __volatile__("rdtsc" : "=a" (a__), "=d" (d__) : );   \
    (ticks) = ((TimerType) a__)|(((TimerType)d__)<<32); \
  } while(0)



#ifdef __cplusplus
}
#endif


#endif
