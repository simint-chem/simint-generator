#!/usr/bin/env python3

#######################################
# Generates the code for the FO method
# from dat files in this directory
#######################################

import argparse
import os
import sys
import re

# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("filename", type=str, help="Output file name base (no extension)")
args = parser.parse_args()

# Load all the *.dat files in this directory
datdir = os.path.dirname(os.path.realpath(__file__))

files = [ f for f in os.listdir(datdir) if f.endswith(".dat") ]

fdat = [None] * len(files)
for i,f in enumerate(files):
  # Get info from filename and read data
  vnm = re.split(r'[_\.]', f)[1:4]
  fdat[int(vnm[0])] = { 
                        'file' : f,
                        'v' : int(vnm[0]),
                        'n' : int(vnm[1]),
                        'm' : int(vnm[2]),
                        'dat' : [ l.strip() for l in open(os.path.join(datdir, f)).readlines()[3:] ] 
                      }

# Output to file
with open(args.filename + ".h", 'w') as f:
  f.write("/*\n")
  f.write("------------------------------------\n")
  f.write(" Generated with:\n")
  f.write("   " + " ".join(sys.argv[:]))
  f.write("\n\n")
  f.write("Dat files:\n")
  for dat in fdat:
    f.write("    {}\n".format(dat['file']))
  f.write("------------------------------------\n")
  f.write("*/\n\n")
  f.write("#ifndef {}_H\n".format(args.filename.upper()))
  f.write("#define {}_H\n".format(args.filename.upper()))
  f.write("\n")
  f.write("#ifdef __cplusplus\n")
  f.write("extern \"C\" {\n");
  f.write("#endif\n")
  f.write("\n")

  f.write("#define BOYS_FO_MAXN {}\n".format(len(fdat)-1))
  f.write("\n")

  f.write("inline void Boys_F_FO(double * const restrict F, int n, double x)\n")
  f.write("{\n")
  f.write("    int idx = n;\n")
  f.write("    switch(n)\n")
  f.write("    {\n")

  for dat in reversed(fdat):
    f.write("        case {}:\n".format(dat['v']))
    f.write("        {\n")
    f.write("            const double num = {}\n".format(dat['dat'][0]))
    for aval in dat['dat'][1:dat['n']+1]:
        f.write("                       + x * ( {}\n".format(aval))
    f.write("                       {};\n".format(")"*dat['n']))
    f.write("\n")
    f.write("            const double den = {}\n".format(1.0))
    for bval in dat['dat'][dat['n']+1:]:
        f.write("                       + x * ( {}\n".format(bval))
    f.write("                       {};\n".format(")"*(dat['m']+1)))

    f.write("\n")
    f.write("            const double frac = num / den;\n")
    f.write("            F[idx--] = pow(frac, {}.5);\n".format(dat['v']))
    f.write("        }\n")
  
  f.write("    }\n")
  f.write("}\n")
  f.write("\n\n")

  # separate F0 function
  f.write("inline double Boys_F0_FO(double x)\n")
  f.write("{\n")

  dat = fdat[0]
  f.write("    const double num = {}\n".format(dat['dat'][0]))
  for aval in dat['dat'][1:dat['n']+1]:
      f.write("               + x * ( {}\n".format(aval))
  f.write("               {};\n".format(")"*dat['n']))
  f.write("\n")
  f.write("    const double den = {}\n".format(1.0))
  for bval in dat['dat'][dat['n']+1:]:
      f.write("               + x * ( {}\n".format(bval))
  f.write("               {};\n".format(")"*(dat['m']+1)))

  f.write("\n")
  f.write("    const double frac = num / den;\n")
  f.write("    return sqrt(frac);\n".format(dat['v']))
  f.write("}\n")
  f.write("\n")
  f.write("#ifdef __cplusplus\n")
  f.write("}\n");
  f.write("#endif\n")
  f.write("\n")
  f.write("#endif\n")
