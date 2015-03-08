#!/usr/bin/env python3

#######################################
# Generates grid for the
# boys function to arbitrary precision
#######################################

import argparse
import sys
from mpmath import mp # arbitrary-precision math

# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True,               help="Output file name base (no extension)")
parser.add_argument("--max-n",    type=int, required=True,               help="Maximum n value to go to")
parser.add_argument("--max-x",    type=str, required=True,               help="Cutoff for the x value")
parser.add_argument("--spacing",  type=str, required=True,               help="Space between pre-computed points")
parser.add_argument("--dps",      type=int, required=False, default=256, help="Decimal precision/sig figs to use/calculate")
parser.add_argument("--eps",      type=int, required=False, default=32,  help="Overall precision / accuracy")
args = parser.parse_args()

# Set the dps option
mp.dps = args.dps

# Convert stuff to mpmath
inc = mp.mpf(args.spacing)
maxx = mp.mpf(args.max_x)
eps = mp.power(10, -args.eps)
maxn = args.max_n
npoints = int(maxx / inc) + 1

print("------------------------------------")
print("Options for Boys function Fn(x):")
print("    Max n: {}".format(maxn))
print("    Max x: {}".format(maxx))
print("  Spacing: {}".format(inc))
print("  npoints: {}".format(npoints))
print("      DPS: {}".format(args.dps))
print("      EPS: 10^({})".format(-args.eps))
print("         = {}".format(eps))
print("------------------------------------")

# we need sqrt(pi)/2
constK = mp.sqrt(mp.pi)/2


# Start at x=0 and increment up
x = mp.mpf(0)
F = []
pts = []

while x < maxx or (x-maxx < inc):
  F2 = [None]*(maxn+1)  # Creates a list of maxn+1 elements

  x2 = 2*x                   # x2 = 2*x
  ex = mp.exp(-x)
  
  # For small x
  if x < 32.0:
    n2 = 2*maxn                # n2 = 2*n
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
    F2[maxn] = s * ex
    m = maxn - 1
    while m >= 0:
      F2[m] = (x2 * F2[m+1] + ex) / (2 * m + 1)
      m -= 1

  else: # For large x
    sqx = mp.sqrt(x)
    F2[0] = constK * mp.erf(sqx) / sqx
    m = 1
    while m <= maxn:
      F2[m] = ((2*(m-1) + 1) * F2[m-1] - ex) / x2
      m += 1

  F.append(F2)
  pts.append(x)
  x += inc

# Output to file
with open(args.filename + ".c", 'w') as f:
  f.write("/*\n")
  f.write(" Generated with:\n")
  f.write("   " + " ".join(sys.argv[:]))
  f.write("\n")
  f.write("------------------------------------\n")
  f.write("Options for Boys function Fn(x):\n")
  f.write("    Max n: {}\n".format(maxn))
  f.write("    Max x: {}\n".format(maxx))
  f.write("  Spacing: {}\n".format(inc))
  f.write("  npoints: {}\n".format(npoints))
  f.write("      DPS: {}\n".format(args.dps))
  f.write("      EPS: 10^({})\n".format(-args.eps))
  f.write("         = {}\n".format(eps))
  f.write("------------------------------------\n")
  f.write("*/\n\n")

  f.write("double boys_grid[{}][{}] = \n".format(npoints, maxn+1))
  f.write("{\n")

  for p,x in zip(F,pts):
    f.write("/* x = {:12}*/  {{".format(mp.nstr(x, 4)))
    for n in p:
      f.write("{:32}, ".format(mp.nstr(n, 18)))
    f.write("},\n")
  f.write("};\n")

with open(args.filename + ".h", 'w') as f: 
  f.write("#ifndef {}_H\n".format(args.filename.upper()))
  f.write("#define {}_H\n".format(args.filename.upper()))
  f.write("\n")
  f.write("#define BOYS_GRID_MAXN {}\n".format(maxn))
  f.write("#define BOYS_GRID_MAXX {}\n".format(maxx))
  f.write("#define BOYS_GRID_NPOINT {}\n".format(npoints))
  f.write("#define BOYS_GRID_LOOKUPFAC {}\n".format(1.0/inc))
  f.write("#define BOYS_GRID_LOOKUPFAC2 {}\n".format(0.5*inc))
  f.write("\n")
  f.write("#endif\n")
