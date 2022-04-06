#!/usr/bin/env python3


import argparse
import sys
import os
import subprocess
import re
import shutil


def UniqueQuartet(q):
  if q[0] < q[1]:
    return False

  if q[2] < q[3]:
    return False

  #if  ( q[0] + q[1] ) < ( q[2] + q[3] ):
  #  return False

  #if (q[0] + q[1]) == (q[2]+q[3]) and (q[0] < q[2]):
  #  return False

  return True


# some helpers
amchar = "spdfghiklmnoqrtuvwxyzabceSPDFGHIKLMNOQRTUVWXYZABCE0123456789"

parser = argparse.ArgumentParser()
parser.add_argument("-l", type=int, required=True, help="Maximum AM")
parser.add_argument("-d", type=int, default=0, help="Derivative level")
parser.add_argument("outfile", type=str, help="Output file")

args = parser.parse_args()


maxam = args.l
der = args.d

print("-------------------------------")
print("Generating Array Filling")
print("Maximum AM: {}".format(maxam))
print("-------------------------------")

valid = {}
for L in range(0, maxam+1):
  valid[L] = set()


for i in range(0, maxam + 1):
  for j in range(0, maxam + 1):
    for k in range(0, maxam + 1):
      for l in range(0, maxam + 1):
        q = (i,j,k,l)
        maxL = max(q)
        valid[maxL].add(q)

# Start the source file
sourcefile = args.outfile
with open(sourcefile, 'w') as f:
  f.write("/*\n")
  f.write(" Generated with:\n")
  f.write("   " + " ".join(sys.argv[:]))
  f.write("\n")
  f.write("*/\n\n")
  f.write("\n\n")
  f.write("#include \"simint/ostei/ostei.h\"\n")
  f.write("#include \"simint/ostei/ostei_init.h\"\n")
  f.write("#include \"simint/ostei/ostei_config.h\"\n")
  if der > 0:
    f.write("#include \"simint/ostei/gen/ostei_deriv{}_generated.h\"\n".format(der))
  else:
    f.write("#include \"simint/ostei/gen/ostei_generated.h\"\n")
  f.write("\n\n")

  f.write("// Stores pointers to the ostei functions\n")
  f.write("#define AMSIZE   SIMINT_OSTEI_MAXAM+1\n")
  f.write("#define DERSIZE  SIMINT_OSTEI_MAXDER+1\n")
  f.write("extern simint_osteifunc simint_osteifunc_array[DERSIZE][AMSIZE][AMSIZE][AMSIZE][AMSIZE];\n")
  f.write("\n\n\n")

  # write the finalize functions
  if der > 0:
    f.write("void simint_ostei_deriv{}_finalize(void)\n".format(der))
  else:
    f.write("void simint_ostei_finalize(void)\n")
  f.write("{\n")
  f.write("    // nothing to do\n")
  f.write("}\n\n\n")

  # Now write the filling function
  if der > 0:
    f.write("void simint_ostei_deriv{}_init(void)\n".format(der))
  else:
    f.write("void simint_ostei_init(void)\n")
  f.write("{\n")

  for L,qset in valid.items():
    if der > 0:
      f.write("    #if SIMINT_OSTEI_DERIV{}_MAXAM >= {}\n".format(der, L))
    else:
      f.write("    #if SIMINT_OSTEI_MAXAM >= {}\n".format(L))
  
    for q in sorted(list(qset)):
      if der > 0:
        fname = "ostei_deriv{}_{}_{}_{}_{}".format(der, amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
      else:
        fname = "ostei_{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
  
      f.write("    simint_osteifunc_array[{}][{}][{}][{}][{}] = {};\n".format(der, q[0], q[1], q[2], q[3], fname))
  
    f.write("    #endif\n\n")

  f.write("}\n\n")
