#!/usr/bin/env python3

#######################################
# Generates grid for the
# boys function to arbitrary precision
#######################################

import argparse
import sys
from boys_common.boys import *
from mpmath import mp # arbitrary-precision math

# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("--spacing",  type=str, required=True,  default="0.1", help="Space between pre-computed points")
parser.add_argument("--max-n",    type=int, required=True,                 help="Maximum n value to go to")
parser.add_argument("--eps",      type=int, required=False, default=17,    help="Overall precision / accuracy")
args = parser.parse_args()

# Set the dps option
mp.dps = 256

# Convert stuff to mpmath
inc = mp.mpf(args.spacing)
eps = mp.power(10, -args.eps)

cutoffs = [None] * (args.max_n+1)
cutoff_diffs = [None] * (args.max_n+1)

x = mp.mpf(0)
maxdiff = eps * 2.0 # force loop to run

while maxdiff > eps:
  x = x + inc

  F_real = BoysValue(args.max_n, x) 
  F_long = BoysValue_Asymptotic(args.max_n, x)

  diff = [ mp.fabs(F_real[i] - F_long[i]) for i in range(0, len(F_real)) ]
  maxdiff = max(diff)

  for idx,val in enumerate(diff):
      if not cutoffs[idx] and diff[idx] < eps:
          cutoffs[idx] = x
          cutoff_diffs[idx] = val

for idx,val in enumerate(cutoffs):
    print("{:4}  {}".format(idx, mp.nstr(val)))
