#!/usr/bin/env python3

#######################################
# Generates some arrays of constants for
# shell related stuff
#######################################

import argparse
import sys
from mpmath import mp # arbitrary-precision math

# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True,               help="Output file name base (no extension)")
parser.add_argument("--max-l",    type=int, required=True,               help="Maximum angular momentum value to go to")
parser.add_argument("--dps",      type=int, required=False, default=256, help="Decimal precision/sig figs to use/calculate")
args = parser.parse_args()

# Set the dps option
mp.dps = args.dps

# Convert stuff to mpmath
maxl = args.max_l

print("------------------------------------")
print("Options for ShellGen constants:")
print("    Max l: {}".format(maxl))
print("      DPS: {}".format(args.dps))
print("------------------------------------")


# For primitive normalization
# c = sqrt(2**l / (pi*sqrt(pi) * (2l-1)!!) )
# afterwards, primitive normalization = c * sqrt(alpha)

pi34 = mp.power(mp.pi, -0.75)
normfac = [None] * (maxl+1)

normfac[0] = pi34;

for l in range(1, maxl+1):
  normfac[l] = normfac[l-1] * mp.sqrt(2.0/(2*l-1))
  

# Output to file
with open(args.filename + ".c", 'w') as f:
  f.write("/*\n")
  f.write(" Generated with:\n")
  f.write("   " + " ".join(sys.argv[:]))
  f.write("\n")
  f.write("------------------------------------\n")
  f.write("Options for ShellGen constants:\n")
  f.write("    Max l: {}\n".format(maxl))
  f.write("      DPS: {}\n".format(args.dps))
  f.write("------------------------------------\n")
  f.write("*/\n\n")

  f.write("const double prim_norm_fac[{}] = \n".format(maxl+1))
  f.write("{\n")

  for i,n in enumerate(normfac):
      f.write("/* l = {:4} */  {:32},\n".format(i, mp.nstr(n, 18)))
  f.write("};\n")

with open(args.filename + ".h", 'w') as f: 
  f.write("#ifndef {}_H\n".format(args.filename.upper()))
  f.write("#define {}_H\n".format(args.filename.upper()))
  f.write("\n")
  f.write("#define SHELL_PRIM_NORMFAC_MAXL {}\n".format(maxl))
  f.write("\n")
  f.write("#endif\n")
