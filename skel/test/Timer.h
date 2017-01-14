#pragma once

#include <chrono>

#ifdef __cplusplus
extern "C" {
#endif


typedef unsigned long long TimerType;


#define CLOCK(ticks, ns) do {                                 \
    volatile unsigned int a__, d__;                              \
    std::chrono::high_resolution_clock::duration dur__ = std::chrono::high_resolution_clock::now().time_since_epoch(); \
    ns = std::chrono::duration_cast<std::chrono::nanoseconds>(dur__).count(); \
    __asm__ __volatile__("rdtsc" : "=a" (a__), "=d" (d__) : );   \
    (ticks) = ((TimerType) a__)|(((TimerType)d__)<<32); \
  } while(0)


#ifdef __cplusplus
}
#endif


