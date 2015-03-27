#!/usr/bin/env python3

###############################
###############################

import sys
import argparse
from boys_common.boys import *
from boys_common.boys_cheby import *
import mpmath as mp
import matplotlib.pyplot as plt


# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True,                help="Output file name base (no extension)")
parser.add_argument("-n",         type=int, required=True,                help="Boys N value")
parser.add_argument("-o",         type=int, required=True,                help="Order of polynomial")
parser.add_argument("-s",         type=str, required=False,               help="Starting order list")
parser.add_argument("-b",         type=int, required=False, default=1,    help="Number of bins")
parser.add_argument("-t",         type=int, required=False, default=200,  help="Number of test points (per bin)")
parser.add_argument("-f",         type=str, required=False, default="1",  help="Base factor (ie 1 = base 2, etc)")
parser.add_argument("--dps",      type=int, required=False, default=256,  help="Decimal precision/sig figs to use/calculate")
args = parser.parse_args()

# Set the dps option
mp.dps = args.dps


test_fx = []    # holds value of boys function at test points
test_ix = []    # holds value of interpolated function at test points
test_x = []     # holds x value of test points
test_abserr = []   # Interpolation absolute error (test_ix - test_fx)
test_relerr = []   # Interpolation relative error (test_ix - test_fx)/test_fx

max_abserr = []
max_relerr = []

# split the interval into bins
fact = mp.mpf(args.f)
binpoints = [mp.mpf(0.0), mp.mpf(1.0)]

for i in range(0, args.b):
  binpoints.append(binpoints[-1]+fact*binpoints[-1]+1)
binpoints = list(zip(binpoints[:-1], binpoints[1:]))

# Make orders into lists
# and find max
olst = [args.o] * (args.b+1)
if args.s:
  olst2 =  [int(x) for x in args.s.split(",")]
  for i,v in enumerate(olst2):
    olst[i] = v

maxorder = max(olst)
binpoints = list(zip(binpoints, olst))

# Whole grid
nbin = len(binpoints)
grid = []

# binpoints is a list of tuples (start, x)
for bounds,order in binpoints:
    print("==Interval (order {}): [{} , {}]".format(order, bounds[0], bounds[1]))

    # Re: order of the polynomial
    # We need to calculate args.o + 1 roots,
    # etc, so we need one more
    data = Interpolate(args.n, order+1, bounds, args.t, BoysValue) 
    test_x.extend(data["test_x"])
    test_fx.extend(data["test_fx"])
    test_ix.extend(data["test_ix"])
    test_abserr.extend(data["test_abserr"])
    test_relerr.extend(data["test_relerr"])
    test_relerr.extend(data["test_relerr"])
    max_abserr.append(data["max_abserr"])
    max_relerr.append(data["max_relerr"])

    # pad polynomial with zero
    # and add to grid
    poly = data["polynomial"]
    poly += [mp.mpf(0.0)] * (maxorder + 1 - len(poly))
    grid.append((bounds, order, poly))

# output to file
with open(args.filename + ".c", 'w') as f:
  f.write("/*\n")
  f.write(" Generated with:\n")
  f.write("   " + " ".join(sys.argv[:]))
  f.write("\n")
  f.write("------------------------------------\n")
  f.write("Options for Boys function Fn(x):\n")
  f.write("      DPS: {}\n".format(args.dps))
  f.write("   Factor: {}\n".format(args.f))
  f.write("    Order: {}\n".format(args.o))
  f.write("    Start: {}\n".format(args.s))
  f.write("     nbin: {}\n".format(nbin))
  f.write(" Last bin: {}\n".format(mp.nstr(binpoints[-1][0])))
  f.write("------------------------------------\n")

  f.write("--------------------------------------------------------------------------------------------------------------\n")
  f.write("Error analysis")
  f.write("--------------------------------------------------------------------------------------------------------------\n")
  f.write("{:25} {:10} {:25} {:25}\n".format("Range", "Order", "Max Abserr", "Max Relerr"))
  for i,v in enumerate(binpoints):
    f.write("{:25} {:<10} {:25} {:25}\n".format(mp.nstr(v[0], 4), v[1], mp.nstr(max_abserr[i], 10), mp.nstr(max_relerr[i], 10)))
  f.write("--------------------------------------------------------------------------------------------------------------\n")
  f.write("*/\n\n")

  f.write("const double boys_chebygrid_F{}[{}][{}] = \n".format(args.n, nbin, maxorder+1))
  f.write("{\n")

  for v in grid:
    f.write("/* x = {:15}  order={:<4} */  {{".format(mp.nstr(v[0], 4), v[1]))
    for n in v[2]:
      f.write("{:32}, ".format(mp.nstr(n, 18)))
    f.write("},\n")
  f.write("};\n")

with open(args.filename + ".h", 'w') as f: 
  f.write("#ifndef {}_H\n".format(args.filename.upper()))
  f.write("#define {}_H\n".format(args.filename.upper()))
  f.write("\n")
  f.write("#define BOYS_CHEBY_NBIN {}\n".format(nbin))
  f.write("#define BOYS_CHEBY_MAXX {}\n".format(mp.nstr(grid[-1][0][1], 18)))
  f.write("#define BOYS_CHEBY_ORDER {}\n".format(maxorder))
  f.write("\n")
  f.write("#endif\n")







