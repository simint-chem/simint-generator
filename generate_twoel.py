#!/usr/bin/env python3


import argparse
import sys
import os
import subprocess
import re

parser = argparse.ArgumentParser()
parser.add_argument("-l", type=int, required=True, help="Maximum AM")
parser.add_argument("-g", type=str, required=True, help="Path to generator binary")
parser.add_argument("-d", type=str, required=True, help="Directory in which to create files")
parser.add_argument("-p", type=str, required=True, help="Type prefix")
parser.add_argument("-b", type=str, required=True, help="Type of boys function")

parser.add_argument("-f",  required=False, action='store_true', help="Flattened primitive loop")
parser.add_argument("-ve", required=False, type=int, default=1000, help="External VRR for this L value and above")
parser.add_argument("-he", required=False, type=int, default=1000, help="External HRR for this L value and above")
args = parser.parse_args()


amchar = "spdfghijklmnoqrtuvwxyzabceSPDFGHIJKLMNOQRTUVWXYZABCE0123456789"
validboys = [ "FO", "split", "vref" ]

def ValidQuartet(q):
  if q[0] < q[1]:
    return False

  if q[2] < q[3]:
    return False

  if  ( q[0] + q[1] ) < ( q[2] + q[3] ):
    return False

  if q[0] < q[2]:
    return False 

  return True





def QStr(q):
  return "( {} {} | {} {} )".format(q[0], q[1], q[2], q[3])













###################################
# Actual code starts here
###################################

if not os.path.isdir(args.d):
  print("The path \"{}\" does not exist or is not a directory".format(args.d))
  quit(1)

if not os.path.isfile(args.g):
  print("The file \"{}\" does not exist or is not a (binary) file".format(args.g))
  quit(1)

if not args.b in validboys:
  print("Invalid boys function type \"{}\"".format(args.b))
  quit(1)


print("-------------------------------")
print("Maximum AM: {}".format(args.l))
print("Naive combinations: {}".format((args.l+1) ** 4))
print("-------------------------------")

valid = []
invalid = []

print()

for i in range(0, args.l + 1):
  for j in range(0, args.l + 1):
    for k in range(0, args.l + 1):
      for l in range(0, args.l + 1):
        q = [i,j,k,l]
        v = ValidQuartet(q)
        if v:
          valid.append(q)
        else:
          invalid.append(q)

        #print("{}  ->  {}".format(QStr(q), v))


print()
print("Valid: {}".format(len(valid)))
for q in valid:
  print("  {}".format(QStr(q)))

print()
print("Invalid: {}".format(len(invalid)))
#for q in invalid:
#  print("  {}".format(QStr(q)))
print()

print("==========================================================================================")
print()

headerbase = "eri_{}.h".format(args.p)
headerfile = os.path.join(args.d, headerbase)

print()
print("Header file: {}".format(headerfile))
print()


for q in valid:
  filebase = "eri_{}_{}_{}_{}_{}".format(args.p, amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
  outfile = os.path.join(args.d, filebase + ".c")
  logfile = os.path.join(args.d, filebase + ".log")
  print("Creating: {}".format(filebase))
  print("      Output: {}".format(outfile))
  print("     Logfile: {}".format(logfile))

  if sum(q) >= args.ve:
    vrrtype = "External"
  else:
    vrrtype = "Inline";
  print("         VRR: {}".format(vrrtype))

  if sum(q) >= args.he:
    hrrtype = "External"
  else:
    hrrtype = "Inline";
  print("         HRR: {}".format(hrrtype))



  with open(logfile, 'w') as lf:
    cmdline = [args.g];
    cmdline.extend(["-q", str(q[0]), str(q[1]), str(q[2]), str(q[3])])
    cmdline.extend(["-p", args.p])
    cmdline.extend(["-b", args.b])
    cmdline.extend(["-o", outfile])

    if args.f:
        cmdline.append("-f")
    if vrrtype == "External":
        cmdline.append("-ve")
    if hrrtype == "External":
        cmdline.append("-he")

    print()
    print("Command line:")
    print(' '.join(cmdline))
    print()
    ret = subprocess.call(cmdline, stdout=lf, stderr=lf)

    if ret != 0:
      print("\n")
      print("*********************************")
      print("Subprocess returned {} - aborting".format(ret))
      print("*********************************")
      print("\n")
      quit(5)
  print()


with open(headerfile, 'w') as hfile:
  defineline = re.sub('[\W]', '_', headerbase.upper()) 
  hfile.write("#ifndef {}\n".format(defineline))
  hfile.write("#define {}\n".format(defineline))
  hfile.write("\n")

  hfile.write("#ifdef __cplusplus\n")
  hfile.write("extern \"C\" {\n")
  hfile.write("#endif\n")
  hfile.write("\n\n")

  hfile.write("#include \"shell/shell.h\"\n")
  hfile.write("\n\n")

  hfile.write("#define ERI_{}_MAXAM {}\n".format(args.p.upper(), args.l))
  hfile.write("\n\n")

  for q in valid:
    funcname = "eri_{}_{}_{}_{}_{}".format(args.p, amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
    intname = "INT__{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
    if args.f:
        hfile.write("int {}(struct multishell_pair_flat const P, struct multishell_pair_flat const Q, double * const restrict {});\n".format(funcname, intname))
    else:
        hfile.write("int {}(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict {});\n".format(funcname, intname))

  hfile.write("\n\n")

  hfile.write("#ifdef __cplusplus\n")
  hfile.write("}\n")
  hfile.write("#endif\n")
  hfile.write("\n")

  hfile.write("#endif\n")

print("\n")
print("Array filling lines:")
print("\n")
for q in valid:
  funcname = "eri_{}_{}_{}_{}_{}".format(args.p, amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
  print("funcs[{}][{}][{}][{}] = {};".format(q[0], q[1], q[2], q[3], funcname))
print("\n\n")
