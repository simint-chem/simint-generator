#!/usr/bin/env python3

#######################################
# Generates an array for the constant
# used in the long-range formula
#
# Eqn. 9.8.9 from Helgaker et al:
# 
# Fn(x) = (2n-1)!!/(2**(n+1)) sqrt(pi/(x**(2n+1)))
#
# So we can pre-compute (2n-1)!!/(2**(n+1)) sqrt(pi) for
# all n that we might need
# After that, Fn(x) = fac * 1/sqrt(x**(2n+1))
#######################################

import argparse
import sys
from mpmath import mp # arbitrary-precision math

# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True,               help="Output file name base (no extension)")
parser.add_argument("--max-n",    type=int, required=True,               help="Maximum n value to go to")
parser.add_argument("--dps",      type=int, required=False, default=256, help="Decimal precision/sig figs to use/calculate")
args = parser.parse_args()

# Set the dps option
mp.dps = args.dps

# Convert stuff to mpmath
maxn = args.max_n

print("------------------------------------")
print("Options for Boys function Fn(x):")
print("    Max n: {}".format(maxn))
print("      DPS: {}".format(args.dps))
print("------------------------------------")


# Generate fac = sqrt(pi)*(2n-1)!!/(2**(n+1))

longfac = [ mp.sqrt(mp.pi) * mp.fac2(2.0*n-1) / mp.power(2.0, n+1) for n in range(0, maxn+1) ]
  


# Output to file
with open(args.filename + ".c", 'w') as f:
  f.write("/*\n")
  f.write(" Generated with:\n")
  f.write("   " + " ".join(sys.argv[:]))
  f.write("\n")
  f.write("------------------------------------\n")
  f.write("Options for long Boys function factor sqrt(pi)*(2n-1)!!/(2**(n+1)):\n")
  f.write("    Max n: {}\n".format(maxn))
  f.write("      DPS: {}\n".format(args.dps))
  f.write("------------------------------------\n")
  f.write("*/\n\n")

  f.write("const double boys_longfac[{}] = \n".format(maxn+1))
  f.write("{\n")

  for i,n in enumerate(longfac):
      f.write("/* n = {:4} */  {:32},\n".format(i, mp.nstr(n, 18)))
  f.write("};\n")

with open(args.filename + ".h", 'w') as f: 
  f.write("#ifndef {}_H\n".format(args.filename.upper()))
  f.write("#define {}_H\n".format(args.filename.upper()))
  f.write("\n")
  f.write("#define BOYS_LONGFAC_MAXN {}\n".format(maxn))
  f.write("\n")
  f.write("#endif\n")
