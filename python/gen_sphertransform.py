#!/usr/bin/env python3

import sys
import argparse
from sympy import *

# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True,               help="Output file name base (no extension)")
parser.add_argument("--max-am",    type=int, required=True,              help="Maximum angular momentum value to go to")
args = parser.parse_args()


def IterateCartesian(g, am):
  if g[2] == am:
    newg = (-1, -1, -1)
  elif g[2] < (am - g[0]):
    newg = (g[0], g[1]-1, g[2]+1)
  else:
    newg = (g[0]-1, am-g[0]+1, 0)

  return newg

def IsValid(g, am):
  return (g[0] >= 0 and g[0] <= am and
          g[1] >= 0 and g[1] <= am and
          g[2] >= 0 and g[2] <= am and
          g[0] + g[1] + g[2] == am)



maxam = int(args.max_am)

# For finding indices
cartmap = {}
sphmap = {}

for l in range(0, maxam+1):

    # 0, +-1, +-2, ... 
    sphmap[l] = [ 0 ]

    for x in range(1, l+1):
      sphmap[l].append(x)
      sphmap[l].append(-x)

    g = (l, 0, 0)

    cartmap[l] = []

    while(IsValid(g, l)):
      cartmap[l].append(g)
      g = IterateCartesian(g, l)


########################################################
# Calculate Cartesian -> Spherical
# This algorithm is from Helgaker, Jorgensen, and Olsen
# "Molecular Electronic-Structure Theory
########################################################
tosph = {}
for l in range(0, maxam + 1):
  tosph[l] = {}

  for m in sphmap[l]:
    tosph[l][m] = {}
    absm = abs(m)
  
    d = (1 if m == 0 else 0)
    NmS = 1 / (2**absm * factorial(l))
    NmS *= sqrt( (2*factorial(l+absm)*factorial(l-absm)) / (2**d) )
  
  
    vm2 = (0 if m >= 0 else 1)
    vm  = sympify(vm2)/2
  
    onefourth = sympify("1/4")
  
    for t in range(0, (l-absm)//2 + 1):
      for u in range(0, t + 1):
        for v2 in range(vm2, int(2*(int(absm/2 - vm) + vm)) + 1, 2):
          v = sympify(v2)/2
          lxyz = ( int(2*t+absm-2*(u+v)), int(2*(u+v)), int(l-2*t-absm))
          c = ( (-1)**(t+v-vm) * onefourth**t * binomial(l, t) * binomial(l-t, absm + t) 
                * binomial(t, u) * binomial(absm, v2))
          tosph[l][m][lxyz] = c*NmS

    
#print()
#print("Transformation to Spherical")
#for l in range(-4, maxam+1):
#    cartlist = cartmap[l]
#    sphlist = sphmap[l]
#
#    print()
#    print("l = {}".format(l))
#    for midx,lm in enumerate(sphlist):
#        print("    m = {:3}   (index {:3})".format(str(lm), midx))
#        for k,v in tosph[l][lm].items():
#            cartidx = cartlist.index(k)
#            print("        ({:2}, {:2}, {:2})    {:2} -> {:25}  =  {}".format(k[0], k[1], k[2], cartidx, v.evalf(20), v))
#        print()


srcfile = args.filename + ".h"
with open(srcfile, "w") as f:
    f.write("#pragma once\n\n")

    f.write("/*\n")
    f.write(" Generated with:\n")
    f.write("   " + " ".join(sys.argv[:]))
    f.write("\n")
    f.write("------------------------------------\n")
    f.write("Options for gen_sphtransform:\n")
    f.write("   Max am: {}\n".format(maxam))
    f.write("------------------------------------\n")
    f.write("*/\n\n")

    f.write("struct spherical_transform_info\n")
    f.write("{\n")
    f.write("    int sphidx;  /* Index of the spherical function */\n")
    f.write("    int cartidx; /* Index of the cartesian function */\n")
    f.write("    double coef; /* Coefficient to use with the cartesian function */\n")
    f.write("};\n\n")

    for l in range(0, maxam+1):
        count = 0
        f.write("static const struct spherical_transform_info sinfo_{}[] = {{\n".format(l))
        cartlist = cartmap[l]
        sphlist = sphmap[l]

        f.write("    /* {}  {}    {}*/\n".format("midx", "cartidx", "coef"))
        for midx,lm in enumerate(sphlist):
            for k,v in tosph[l][lm].items():
                cartidx = cartlist.index(k)
                f.write("    {{ {:3}, {:3}, {:>25} /*{}*/ }},  // {}\n".format(midx, cartidx, str(N(v, 18)), v, k))
                count += 1

        f.write("};\n\n")
        f.write("static const int spherical_transform_count_{} = {};\n\n".format(l, count))


