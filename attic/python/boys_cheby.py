#!/usr/bin/env python3

#######################################
# Generates the chebyschev polynomials
# if the boys function for a given interval
# to arbitrary precision
#######################################

import mpmath as mp
import numpy as np

def GenCheby(n, bounds):

    bscale = (bounds[1]-bounds[0])/2.0

    # T = coefficients of polynomial
    # xi = roots
    T = mp.matrix(n+1, n+1)  # Add one for T0
    xi = [None]*(n+1)

    T[0,0] = mp.mpf(1)
    xi[0] = None

    T[1,0] = mp.mpf(0)
    T[1,1] = mp.mpf(1)
    xi[1] = [mp.mpf(0)]

    for i in range(2, n+1):
        # Add the 2x part
        for j in range(1, i+1):
            T[i,j] = mp.mpf(2.0) * T[i-1,j-1]

        # Subtract T(i-2)
        for j in range(0, n+1):
            T[i,j] -= T[i-2, j]

        # calculate roots
        # with scaling
        xi[i] = [ bscale * mp.cos(mp.pi * (k - 0.5) / i) for k in range(1, i+1) ]

    # scale Ti for the range
    for i in range(0, n+1):
      for j in range(0, n+1):
          T[i,j] *= mp.power(bscale, -j)

    return (T,xi)


# Polynomial is a + bx + cx**2 + dx**3 + ....
def EvalPolynomial(p, x):
    s = mp.mpf(0.0)
    for i,v in enumerate(p):
        s += v * mp.power(x, i)
    return s    


# n = n value of Boys function
# order = order of polynomial
def Interpolate(n, order, bounds, npoints, func):
    # return value is a dictionary
    ret = {}

    # for shifting the bounds
    bshift = (bounds[0] + bounds[1])/2.0

    # Generate the chebyschev functions for this interval
    (T, xi) = GenCheby(order, bounds)

    # Evaluate all polynomials at roots of highest polynomial
    rts = xi[order]

    #print("Coefs of highest polynomial")
    #print(T[order,:])
    #print()
    #print("Roots of highest polynomial")
    #print(mp.nstr(rts, 8))

    Txi = mp.matrix(order+1, order)
    for k in range(0, order+1):
      Tn = T[k,:]

      for i,x in enumerate(rts):
        for j,a in enumerate(Tn):
          Txi[k, i] += a * mp.power(x, j)

    #print()
    #print("Matrix of root evaluations") 
    #print(mp.nstr(Txi, 8))   


    # evaluate the boys function at roots
    # keeping in mind the shift
    Fxi = [ func(n, x+bshift)[n] for x in rts ]

    #print()
    #print("Value of boys function at (shifted) roots")
    #print(mp.nstr(Fxi, 8))

    # Generate coefficients
    coef = [None] * (order)
    for i in range(0, order):  # note that we are cutting off the last row - it should be zero!
      coef[i] = mp.mpf(0.0)
      
      for j in range(0, order):
        coef[i] += Txi[i,j] * Fxi[j]
      coef[i] *= mp.mpf(2.0)/mp.mpf(order)


    #print()
    #print("Chebyschev Intermediate Coefficients")
    #print(mp.nstr(coef, 8))

    # Apply coefficients to the proper orders and obtain the coefficients of the Chebyschev polynomial
    ChebyMat = mp.matrix(order, order)
    c = [mp.mpf(0)] * order

    for i in range(0, order):
      for j in range(0, order):
        c[j] += T[i,j] * coef[i]

    # subtract 1/2 coef[0] from the x^0 term
    c[0] -= mp.mpf(0.5) * coef[0]

    #print()
    #print("Final polynomial coefficents") 
    #print(mp.nstr(c, 8))

    # return the polynomial
    ret["polynomial"] = c;

    # Generate test points for this interval
    xval = np.linspace(bounds[0], bounds[1], npoints)
    ret["test_x"]   = xval
    ret["test_fx"]  = [ func(n, x)[n] for x in xval ]
    ret["test_ix"]  = [ EvalPolynomial(c, x-bshift) for x in xval ]
    ret["test_abserr"] = [ a - b for a,b in zip(ret["test_ix"], ret["test_fx"]) ]
    ret["test_relerr"] = [ a/b for a,b in zip(ret["test_abserr"], ret["test_fx"]) ]
    ret["max_abserr"]  = max([mp.fabs(a) for a in ret["test_abserr"]])
    ret["max_relerr"]  = max([mp.fabs(a) for a in ret["test_relerr"]])

    return ret;

