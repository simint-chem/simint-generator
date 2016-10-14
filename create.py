#!/usr/bin/env python3


import argparse
import sys
import os
import subprocess
import re
import shutil


def ValidQuartet(q):
  if q[0] < q[1]:
    return False

  if q[2] < q[3]:
    return False

  if  ( q[0] + q[1] ) < ( q[2] + q[3] ):
    return False

  if (q[0] + q[1]) == (q[2]+q[3]) and (q[0] < q[2]):
    return False

  return True





def QStr(q):
  return "( {} {} | {} {} )".format(q[0], q[1], q[2], q[3])



# path to this file
thisfile = os.path.realpath(__file__)
topdir = os.path.dirname(thisfile)

# some helpers
amchar = "spdfghijklmnoqrtuvwxyzabceSPDFGHIJKLMNOQRTUVWXYZABCE0123456789"

parser = argparse.ArgumentParser()
parser.add_argument("-l", type=int, required=True, help="Maximum AM")
parser.add_argument("-g", type=str, required=True, help="Path to directory with generator programs")

parser.add_argument("-c", type=str, required=False, default="", help="CPU Flags (comma separated)")
parser.add_argument("-p", required=False, action='store_true', help="Generate code for permuted angular momentum quartets")

parser.add_argument("-ve", required=False, type=int, default=1000, help="External VRR for this L value and above")
parser.add_argument("-he", required=False, type=int, default=1000, help="External HRR for this L value and above")
parser.add_argument("-s",  required=False, type=int, default=0,    help="Max contracted integral stack size in bytes (per shell quartet)")
parser.add_argument("-i",  required=False, action='store_true', help="Use intrinsics")
parser.add_argument("-S",  required=False, action='store_true', help="Generate scalar code")

parser.add_argument("outdir", type=str, help="Output directory")

args = parser.parse_args()








###################################
# Actual code starts here
###################################
# directory with the files to copy
skeldir = os.path.join(topdir, "skel")


# paths to generator programs
ostei_gen = os.path.join(args.g, "ostei_generator")
hrr_gen = os.path.join(args.g, "ostei_hrr_generator")
vrr_gen = os.path.join(args.g, "ostei_vrr_generator")

if not os.path.isfile(ostei_gen):
  print("The file \"{}\" does not exist or is not a (binary) file".format(ostei_gen))
  quit(1)

if not os.path.isfile(hrr_gen):
  print("The file \"{}\" does not exist or is not a (binary) file".format(hrr_gen))
  quit(1)

if not os.path.isfile(vrr_gen):
  print("The file \"{}\" does not exist or is not a (binary) file".format(vrr_gen))
  quit(1)


if args.S and args.i:
  print("Intrinsics don't make sense with scalar code")
  quit(1)



####################################################
# Create output directory
####################################################
outdir = os.path.join(args.outdir, "simint")
outdir_test = os.path.join(args.outdir, "test")
outdir_examples = os.path.join(args.outdir, "examples")
outdir_cmake = os.path.join(args.outdir, "cmake")

outdir_ostei = os.path.join(outdir, "ostei")
outdir_osteigen = os.path.join(outdir_ostei, "gen")
outdir_vec = os.path.join(outdir, "vectorization")

if os.path.isdir(args.outdir):
  print("WARNING - output directory exists. Overwriting...")
  shutil.rmtree(outdir, ignore_errors=True)
  shutil.rmtree(outdir_test, ignore_errors=True)
  shutil.rmtree(outdir_examples, ignore_errors=True)
  shutil.rmtree(outdir_cmake, ignore_errors=True)


shutil.copytree(os.path.join(skeldir, "simint"),        outdir)
shutil.copytree(os.path.join(skeldir, "test"),          outdir_test)
shutil.copytree(os.path.join(skeldir, "examples"),      outdir_examples)
shutil.copytree(os.path.join(skeldir, "cmake"),         outdir_cmake)
shutil.copy(os.path.join(skeldir, "CMakeLists.txt"),    args.outdir)
shutil.copy(os.path.join(skeldir, "README"),            args.outdir)
shutil.copy(os.path.join(skeldir, "LICENSE"),           args.outdir)
shutil.copy(os.path.join(skeldir, "CHANGELOG"),           args.outdir)



####################################################
# Generate the external HRR source
####################################################
print("-------------------------------")
print("Generating HRR")
print("-------------------------------")

valid = set()
invalid = set()

print()

# TODO - generate too much here?
for i in range(1, args.l+1):
  for j in range(1, i+1):
    valid.add((i,j));
    if args.p:
        valid.add((j, i))

for i in range(args.l+1, 2*args.l):
  for j in range(1, 2*args.l-i+1):
    valid.add((i,j));
    if args.p:
        valid.add((j, i))


print()

print("Valid: {}".format(len(valid)))
for q in valid:
  print("  {}".format(QStr([q[0], q[1], q[0], q[1]])))

print()
print("==========================================================================================")
print()

headerbase = "hrr_generated.h"
headerfile = os.path.join(outdir_osteigen, headerbase)

print()
print("Header file: {}".format(headerfile))
print()


# Start the header file
with open(headerfile, 'w') as hfile:
  hfile.write("#pragma once\n\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("extern \"C\" {\n")
  hfile.write("#endif\n")
  hfile.write("\n\n")
  hfile.write("#include \"simint/vectorization/vectorization.h\"\n")
  hfile.write("\n")


for q in valid:
    filebase = "hrr_{}_{}".format(amchar[q[0]], amchar[q[1]])

    outfile = os.path.join(outdir_osteigen, filebase + ".c")
    logfile = os.path.join(outdir_osteigen, filebase + ".log")
    print("Creating HRR: {}".format(filebase))
    print("      Output: {}".format(outfile))
    print("     Logfile: {}".format(logfile))
    cmdline = [hrr_gen]

    cmdline.extend(["-q", str(q[0]), str(q[1])])
    cmdline.extend(["-o", outfile])
    cmdline.extend(["-oh", headerfile])

    if args.c:
        cmdline.extend(["-c", str(args.c)])
    if args.i:
        cmdline.append("-i")
    if args.S:
        cmdline.append("-S")

    print()
    print("Command line:")
    print(' '.join(cmdline))
    print()

    with open(logfile, 'w') as lf:
      ret = subprocess.call(cmdline, stdout=lf)

    if ret != 0:
      print("\n")
      print("*********************************")
      print("When generating hrr sources")
      print("Subprocess returned {} - aborting".format(ret))
      print("*********************************")
      print("\n")
      quit(1)

# Close out the header file
with open(headerfile, 'a') as hfile:
  hfile.write("\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("}\n")
  hfile.write("#endif\n")
  hfile.write("\n")



####################################################
# Generate the external VRR source
####################################################
print("-------------------------------")
print("Generating VRR")
print("-------------------------------")

valid = set()
invalid = set()

print()

# also need (X s | X s), etc
for i in range(0, args.l*2+1):
  for j in range(0, args.l*2+1):
    if i == 0 and j == 0: # skip ssss
      continue
    valid.add((i, 0, j, 0))

    if args.p:
      valid.add((0, i, 0, j))
      valid.add((i, 0, 0, j))
      valid.add((0, i, j, 0))

print()

print("Valid: {}".format(len(valid)))
for q in valid:
  print("  {}".format(QStr(q)))

print()
print("==========================================================================================")
print()

headerbase = "vrr_generated.h"
headerfile = os.path.join(outdir_osteigen, headerbase)

print()
print("Header file: {}".format(headerfile))
print()

# Start the header file
with open(headerfile, 'w') as hfile:
  hfile.write("#pragma once\n\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("extern \"C\" {\n")
  hfile.write("#endif\n")
  hfile.write("\n\n")
  hfile.write("#include \"simint/vectorization/vectorization.h\"\n")
  hfile.write("\n")


for q in valid:
    filebase = "vrr_{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
    outfile = os.path.join(outdir_osteigen, filebase + ".c")
    logfile = os.path.join(outdir_osteigen, filebase + ".log")
    print("Creating VRR: {}".format(filebase))
    print("      Output: {}".format(outfile))
    print("     Logfile: {}".format(logfile))


    with open(logfile, 'w') as lf:
      cmdline = [vrr_gen]
      cmdline.extend(["-q", str(q[0]), str(q[1]), str(q[2]), str(q[3])])
      cmdline.extend(["-o", outfile])
      cmdline.extend(["-oh", headerfile])

      if args.c:
          cmdline.extend(["-c", str(args.c)])
      if args.i:
          cmdline.append("-i")
      if args.S:
          cmdline.append("-S")

      print()
      print("Command line:")
      print(' '.join(cmdline))
      print()

      ret = subprocess.call(cmdline, stdout=lf)

      if ret != 0:
        print("\n")
        print("*********************************")
        print("When generating vrr sources")
        print("Subprocess returned {} - aborting".format(ret))
        print("*********************************")
        print("\n")
        quit(1)

    print()

# Close out the header file
with open(headerfile, 'a') as hfile:
  hfile.write("\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("}\n")
  hfile.write("#endif\n")
  hfile.write("\n")


####################################################
# Generate the ERI sources and headers
####################################################
print("-------------------------------")
print("Generating ERI")
print("Maximum AM: {}".format(args.l))
print("-------------------------------")

valid = set()
invalid = set()

print()

for i in range(0, args.l + 1):
  for j in range(0, args.l + 1):
    for k in range(0, args.l + 1):
      for l in range(0, args.l + 1):
        q = (i,j,k,l)
        if args.p or ValidQuartet(q):
          valid.add(q)
        else:
          invalid.add(q)

print()
print("Valid: {}".format(len(valid)))
for q in valid:
  print("  {}".format(QStr(q)))

print()
print("==========================================================================================")
print()




# Generate the ostei
headerbase = "ostei_generated.h"
headerfile = os.path.join(outdir_osteigen, headerbase)

# Maximum required contwork
maxworksize = 0  # number of elements

print()
print("Header file: {}".format(headerfile))
print()

# Start the header file
with open(headerfile, 'w') as hfile:
  hfile.write("#pragma once\n\n")
  hfile.write("\n\n")
  hfile.write("#include \"simint/ostei/ostei.h\"\n")
  hfile.write("#include \"simint/ostei/ostei_general.h\"\n")
  hfile.write("#include \"simint/ostei/gen/vrr_generated.h\"\n")
  hfile.write("#include \"simint/ostei/gen/hrr_generated.h\"\n")

  hfile.write("\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("extern \"C\" {\n")
  hfile.write("#endif\n")
  hfile.write("\n")


for q in valid:
  filebase = "ostei_{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
  outfile = os.path.join(outdir_osteigen, filebase + ".c")
  logfile = os.path.join(outdir_osteigen, filebase + ".log")
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
    cmdline = [ostei_gen];
    cmdline.extend(["-q", str(q[0]), str(q[1]), str(q[2]), str(q[3])])
    cmdline.extend(["-o", outfile])
    cmdline.extend(["-oh", headerfile])
    cmdline.extend(["-s", str(args.s)])

    if args.c:
        cmdline.extend(["-c", str(args.c)])
    if vrrtype == "External":
        cmdline.append("-ve")
    if hrrtype == "External":
        cmdline.append("-he")

    if args.i:
        cmdline.append("-i")
    if args.S:
        cmdline.append("-S")

    print()
    print("Command line:")
    print(' '.join(cmdline))
    print()

    ret = subprocess.call(cmdline, stdout=lf, stderr=lf)

    if ret != 0:
      print("\n")
      print("*********************************")
      print("While generating ostei")
      print("Subprocess returned {} - aborting".format(ret))
      print("*********************************")
      print("\n")
      quit(5)


  # reopen the logfile, find contwork
  for line in open(logfile, 'r').readlines():
    if line.startswith("CONTWORK SIZE"):
      maxworksize = max(maxworksize, int(line.split()[2]))

  print()


# Close out the header file

with open(headerfile, 'a') as hfile:
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("}\n")
  hfile.write("#endif\n")
  hfile.write("\n")


######################
# OSTEI config file
######################
headerbase = "ostei_config.h"
headerfile = os.path.join(outdir_ostei, headerbase)
with open(headerfile, 'w') as hfile:
  hfile.write("#pragma once\n\n")
  hfile.write("#include \"simint/vectorization/vectorization.h\"\n")
  hfile.write("\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("extern \"C\" {\n")
  hfile.write("#endif\n")
  hfile.write("\n\n")
  hfile.write("#define SIMINT_OSTEI_MAXAM {}\n".format(args.l))
  hfile.write("#define SIMINT_OSTEI_MAX_WORKSIZE ((SIMINT_SIMD_ROUND(SIMINT_NSHELL_SIMD * {})))\n".format(maxworksize))
  hfile.write("#define SIMINT_OSTEI_MAX_WORKMEM (SIMINT_OSTEI_MAX_WORKSIZE * sizeof(double))\n")
  hfile.write("\n\n")
  hfile.write("\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("}\n")
  hfile.write("#endif\n")
  hfile.write("\n")


####################################################
# Overall config header
####################################################
cpuflags = args.c.lower().split(',')
sinfofile = os.path.join(outdir, "simint_config.h")
with open(sinfofile, 'w') as sf:

  # Parse the cpu flags


  sf.write("#pragma once\n\n")

  if "avx" in cpuflags:
    sf.write("#include \"simint/vectorization/intrinsics_avx.h\"\n")
    sf.write("#define SIMINT_SIMD_LEN 4\n")
    sf.write("#define SIMINT_AVX\n")
  elif "sse2" in cpuflags:
    sf.write("#include \"simint/vectorization/intrinsics_sse.h\"\n")
    sf.write("#define SIMINT_SIMD_LEN 2\n")
    sf.write("#define SIMINT_SSE\n")
  else:
    sf.write("#define SIMINT_SIMD_LEN 1\n")
    sf.write("#define SIMINT_SCALAR\n")

  if args.p:
    sf.write("#define SIMINT_OSTEI_PERMUTATIONS\n")
  else:
    sf.write("#define SIMINT_OSTEI_NO_PERMUTATIONS\n")

  sf.write("\n\n")


####################################################
# Cmake flags file
####################################################
cmflagfile = os.path.join(outdir_cmake, "DefaultFlags_generated.cmake")
with open(cmflagfile, 'w') as sf:

  if "avx" in cpuflags:
    sf.write("include(cmake/DefaultFlags_avx.cmake)")
  elif "sse2" in cpuflags:
    sf.write("include(cmake/DefaultFlags_sse.cmake)")
  else:
    sf.write("include(cmake/DefaultFlags_scalar.cmake)")
  sf.write("\n")



