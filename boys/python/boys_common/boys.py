#!/usr/bin/env python3

#######################################
# Generates the chebyschev polynomials
# if the boys function for a given interval
# to arbitrary precision
#######################################

import mpmath as mp  # arbitrary-precision math

# precision of the BoysValue
eps = 1e-16

# we need sqrt(pi)/2
constK = mp.sqrt(mp.pi)/2

# This is the old way of calculating it
def BoysValue_old(n, x):
    F = [None] * (n+1)

    x2 = 2*x
    ex = mp.exp(-x)

    # For small x
    if x < 32.0:
      n2 = 2*n                # n2 = 2*n
      num = mp.fac2(n2-1)        # num = (2*n - 1)!!
      s = 1/(n2+1)
      term1 = 1.0 # just start at some number
      i = 0

      while mp.fabs(term1) > eps:
        i += 1
        num *= x2
        term1 = num / mp.fac2(n2 + 2*i + 1)
        s += term1

      # Downward recursion
      F[n] = s * ex
      m = n - 1
      while m >= 0:
        F[m] = (x2 * F[m+1] + ex) / (2 * m + 1)
        m -= 1

    else: # For large x
      sqx = mp.sqrt(x)
      F[0] = constK * mp.erf(sqx) / sqx
      m = 1
      while m <= n:
        F[m] = ((2*(m-1) + 1) * F[m-1] - ex) / x2
        m += 1

    return F


# Probably a better way
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

