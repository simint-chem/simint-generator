#!/usr/bin/env python3

#######################################
# Generates the chebyschev polynomials
# if the boys function for a given interval
# to arbitrary precision
#######################################

import mpmath as mp  # arbitrary-precision math

def BoysValue(n, x):
    F = []

    if x == mp.mpf("0"):
        for i in range(0, n+1):
            F.append(mp.mpf(1.0)/(mp.mpf(2.0*i+1)))
    else:
        for i in range(0, n+1):
            N = i+mp.mpf("0.5")
            F.append(mp.gammainc(N, 0, x) * 1.0/(2.0 * mp.power(x, N)))

    return F


def BoysValue_Asymptotic(n, x):
    F = []

    if x == mp.mpf("0"):
        for i in range(0, n+1):
            F.append(mp.mpf(1.0)/(mp.mpf(2.0*i+1)))
    else:
        for i in range(0, n+1):
            val = mp.sqrt( mp.pi() / mp.power(x, 2*i+1) )
            val *= ( mp.fac2(2*i-1) / mp.power("2", i+1) )
            F.append(val)

    return F
