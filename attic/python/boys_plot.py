#!/usr/bin/env python3

###############################
# Plots Fn(x), the interpolated
# function, and the error
###############################

import argparse
from boys_common.boys_cheby import *
import matplotlib.pyplot as plt


# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("-n",         type=int, required=True,                help="Boys N value")
parser.add_argument("--dps",      type=int, required=False, default=256,  help="Decimal precision/sig figs to use/calculate")
parser.add_argument("-t",         type=int, required=False, default=200,  help="Number of test points")
parser.add_argument("minx",       type=str,                               help="Starting x value")
parser.add_argument("maxx",       type=str,                               help="Ending x value")
args = parser.parse_args()

# Set the dps option
mp.dps = args.dps

# Convert stuff to mpmath
minx = mp.mpf(args.minx)
maxx = mp.mpf(args.maxx)

pts = [ mp.mpf(x) for x in np.linspace(minx, maxx, args.t+1) ]

F = []

# binpoints is a list of tuples (start, x)
for p in pts:
  # Calculate the old way
  #F1.append(BoysValue_old(args.n, p))
  
  # New Way
  F.append(BoysValue(args.n, p)[args.n])

logF = [ 1.0/f for f in F ]

plt.plot(pts, logF, 'r-', pts, F, 'g-')
plt.show()
