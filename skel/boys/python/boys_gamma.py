#!/usr/bin/env python3

from boys_common.boys import *
import mpmath as mp

mp.mp.dps = 40
tiny = mp.mpf("1e-30")


def LentzIncGamma(a, x, maxiter):
    f = [tiny]
    C = [f[0]]
    D = [mp.mpf(0)]

    # j = 1
    D.append(x + 1.0*D[0])
    C.append(x + 1.0/C[0])
    D[1] = 1.0/D[1]
    delta = C[1]*D[1]
    f.append(f[0]*delta)

    i = 2
    for j in range(1, maxiter):
      jp = mp.mpf(j)
      ja = mp.mpf(jp-a)

      D.append(mp.mpf(0))
      C.append(mp.mpf(0))
      f.append(mp.mpf(0))
      D[i]       = mp.mpf(1) + (ja)*D[i-1]
      C[i]       = 1.0 + (ja)/C[i-1]
      D[i]       = 1.0/D[i]
      delta      = C[i]*D[i]
      f[i]       = f[i-1]*delta
      i += 1

      D.append(mp.mpf(0))
      C.append(mp.mpf(0))
      f.append(mp.mpf(0))
      D[i]     = x + (jp)*D[i-1]
      C[i]     = x + (jp)/C[i-1]
      D[i]     = 1.0/D[i]
      delta    = C[i]*D[i]
      f[i]     = f[i-1]*delta
      i += 1

    print("Delta: {}".format(mp.nstr(delta, 16)))
    return f[-1]


def BoysValue_Gamma(n, x, maxiter):
    F = []

    ex = mp.exp(-x)

    for m in range(0, n+1):
      xn12 = mp.mpf(1.0)/mp.power(x, m+0.5)
      cgam = mp.gammainc(m+0.5)  # "Complete" gamma?

      # Incomplete gamma via continued fractions and Lentz method
      lgam = LentzIncGamma(m+0.5, x, maxiter)

      # put it together
      F.append(mp.mpf(0.5) * (xn12 * cgam - ex*lgam))

    return F



X = 0.01
N = 1
maxiter = 500

oldF = BoysValue(N, X)
newF = BoysValue_Gamma(N, X, maxiter)

sx = mp.sqrt(X)
fac = mp.sqrt(mp.pi)/2.0
print("ERF: {}".format(mp.nstr(fac/sx * mp.erf(sx), 16)))
print("OLD: {}".format(mp.nstr(oldF[0], 16)))
print("NEW: {}".format(mp.nstr(newF[0], 16)))



