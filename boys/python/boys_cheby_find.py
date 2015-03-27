#!/usr/bin/env python3

###############################
# Plots Fn(x), the interpolated
# function, and the error
###############################

import argparse
from boys_common.boys import *
from boys_common.boys_cheby import *
import mpmath as mp
import matplotlib.pyplot as plt


# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("-n",         type=int, required=True,                    help="Boys N value")
parser.add_argument("-o",         type=int, required=True,                    help="Order of polynomial")
parser.add_argument("-p",         action="store_true",                        help="Make the plot")
parser.add_argument("-t",         type=int, required=False, default=250,      help="Number of test points (per bin)")
parser.add_argument("--dps",      type=int, required=False, default=256,      help="Decimal precision/sig figs to use/calculate")
parser.add_argument("--eps",      type=str, required=False, default="5e-14",  help="Max error")
parser.add_argument("x",          type=str,                                   help="Starting x value")
parser.add_argument("sp",         type=str,                                   help="Starting guess at spacing")
args = parser.parse_args()

# Set the dps option
mp.dps = args.dps

# Convert stuff to mpmath
eps = mp.mpf(args.eps)

# order of the polynomial
# We need to calculate args.o + 1 roots,
# etc, so we need one more
order = args.o + 1


# Current spacing
sp = mp.mpf(args.sp)

test_fx = []    # holds value of boys function at test points
test_ix = []    # holds value of interpolated function at test points
test_x = []     # holds x value of test points
test_abserr = []   # Interpolation absolute error (test_ix - test_fx)
test_relerr = []   # Interpolation relative error (test_ix - test_fx)/test_fx


# current x value
x = mp.mpf(args.x)

# and error
err1 = mp.mpf(1.0)

# delta for calculating the derivative
delta = mp.mpf(sp/10.0)

# damping factor
damp = mp.mpf(0.5)

while err1 > eps:
  # Calculate new guess for spacing
  data1 = Interpolate(args.n, order, (x, x+sp), args.t, BoysValue) 
  data2 = Interpolate(args.n, order, (x, x+sp+delta), args.t, BoysValue) 

  err1 = data1["max_abserr"]
  err2 = data2["max_abserr"]
  der = (err2 - err1)/delta

  print("sp = {} err1 = {} der = {}".format(mp.nstr(sp, 6),
                                            mp.nstr(err1, 6),
                                            mp.nstr(der, 6)))

  if err1 < eps:
      print("Found: x = {} spacing = {}".format(mp.nstr(x, 6), mp.nstr(sp, 6)))
      break

  # new guess
  oldsp = sp
  sp = damp * ((eps - err1)/der) + oldsp
  print("step: {}".format(mp.nstr(sp - oldsp, 6)))
  delta = mp.mpf(sp/10.0)
