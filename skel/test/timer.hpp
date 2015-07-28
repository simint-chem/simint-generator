#ifndef SIMINT_TIMER_HPP
#define SIMINT_TIMER_HPP

#include <utility>

#define CLOCK(ticks) do {                                 \
    volatile unsigned int a__, d__;                              \
    __asm__ __volatile__("rdtsc" : "=a" (a__), "=d" (d__) : );   \
    (ticks) = ((unsigned long long) a__)|(((unsigned long long)d__)<<32); \
  } while(0)


typedef std::pair<unsigned long long, double> TimerInfo;

inline TimerInfo operator+(const TimerInfo & lhs, const TimerInfo & rhs)
{
    return {lhs.first + rhs.first, lhs.second + rhs.second};
}

inline TimerInfo & operator+=(TimerInfo & lhs, const TimerInfo & rhs)
{
    lhs = lhs + rhs;
    return lhs;
}

#endif
