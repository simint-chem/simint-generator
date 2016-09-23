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
validboys = [ "FO", "split", "vref" ]

parser = argparse.ArgumentParser()
parser.add_argument("-l", type=int, required=True, help="Maximum AM")
parser.add_argument("-g", type=str, required=True, help="Path to directory with generator programs")
parser.add_argument("-d", type=str, required=True, help="Path to dat directory")
parser.add_argument("-b", type=str, required=True, help="Type of boys function")

parser.add_argument("-c", type=str, required=False, default="", help="CPU Flags (comma separated)")
parser.add_argument("-p", required=False, action='store_true', help="Generate code for permuted angular momentum quartets")

parser.add_argument("-et", required=False, action='store_true', help="Use electron transfer rather than second VRR")
parser.add_argument("-ve", required=False, type=int, default=1000, help="External VRR for this L value and above")
parser.add_argument("-ee", required=False, type=int, default=1000, help="External ET for this L value and above")
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
eri_gen = os.path.join(args.g, "eri_generator")
hrr_gen = os.path.join(args.g, "eri_hrr_generator")
et_gen = os.path.join(args.g, "eri_et_generator")
vrr_gen = os.path.join(args.g, "eri_vrr_generator")

if not os.path.isdir(args.d):
  print("The directory \"{}\" does not exist or is not accessible".format(args.d))
  quit(1)

if not os.path.isfile(eri_gen):
  print("The file \"{}\" does not exist or is not a (binary) file".format(eri_gen))
  quit(1)

if not os.path.isfile(hrr_gen):
  print("The file \"{}\" does not exist or is not a (binary) file".format(hrr_gen))
  quit(1)

if not os.path.isfile(et_gen):
  print("The file \"{}\" does not exist or is not a (binary) file".format(et_gen))
  quit(1)

if not os.path.isfile(vrr_gen):
  print("The file \"{}\" does not exist or is not a (binary) file".format(vrr_gen))
  quit(1)

if not args.b in validboys:
  print("Invalid boys function type \"{}\"".format(args.b))
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

outdir_eri = os.path.join(outdir, "eri")
outdir_erigen = os.path.join(outdir_eri, "gen")
outdir_vec = os.path.join(outdir, "vectorization")

if os.path.isdir(args.outdir):
  print("WARNING - output directory exists. Overwriting...")
  shutil.rmtree(outdir, ignore_errors=True)
  shutil.rmtree(outdir_test, ignore_errors=True)
  shutil.rmtree(outdir_examples, ignore_errors=True)


shutil.copytree(os.path.join(skeldir, "simint"),        outdir)
shutil.copytree(os.path.join(skeldir, "test"),          outdir_test)
shutil.copytree(os.path.join(skeldir, "examples"),      outdir_examples)
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
headerfile = os.path.join(outdir_erigen, headerbase)

print()
print("Header file: {}".format(headerfile))
print()


# Start the header file
defineline = re.sub('[\W]', '_', headerbase.upper())
with open(headerfile, 'w') as hfile:
  hfile.write("#ifndef {}\n".format(defineline))
  hfile.write("#define {}\n".format(defineline))
  hfile.write("\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("extern \"C\" {\n")
  hfile.write("#endif\n")
  hfile.write("\n\n")
  hfile.write("#include \"simint/vectorization/vectorization.h\"\n")
  hfile.write("\n")


for q in valid:
    filebase = "hrr_{}_{}".format(amchar[q[0]], amchar[q[1]])

    outfile = os.path.join(outdir_erigen, filebase + ".c")
    logfile = os.path.join(outdir_erigen, filebase + ".log")
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
  hfile.write("#endif // {}\n".format(defineline))
  hfile.write("\n")



####################################################
# Generate the external ET source
####################################################
print("-------------------------------")
print("Generating ET")
print("-------------------------------")

valid = set()
invalid = set()

print()

if args.et:
  for i in range(0, 2*args.l+4):
    for j in range(1, 2*args.l+4):
      valid.add((i,0,j,0))
else:
  print("Skipping ET...")
print()

print("Valid: {}".format(len(valid)))
for q in valid:
  print("  {}".format(QStr(q)))

print()
print("==========================================================================================")
print()

headerbase = "et_generated.h"
headerfile = os.path.join(outdir_erigen, headerbase)

print()
print("Header file: {}".format(headerfile))
print()

# Start the header file
defineline = re.sub('[\W]', '_', headerbase.upper())
with open(headerfile, 'w') as hfile:
  hfile.write("#ifndef {}\n".format(defineline))
  hfile.write("#define {}\n".format(defineline))
  hfile.write("\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("extern \"C\" {\n")
  hfile.write("#endif\n")
  hfile.write("\n\n")
  hfile.write("#include \"simint/vectorization/vectorization.h\"\n")
  hfile.write("\n")


for q in valid:
    filebase = "et_{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
    outfile = os.path.join(outdir_erigen, filebase + ".c")
    logfile = os.path.join(outdir_erigen, filebase + ".log")
    print("Creating ET: {}".format(filebase))
    print("      Output: {}".format(outfile))
    print("     Logfile: {}".format(logfile))

    with open(logfile, 'w') as lf:
      cmdline = [et_gen]
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
  hfile.write("#endif // {}\n".format(defineline))
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

if args.et:
  for i in range(1, args.l*4+1):
    valid.add((i, 0, 0, 0))

    if args.p:
      valid.add((0, 0, i, 0))

    # TODO - more are probably needed here

else:
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
headerfile = os.path.join(outdir_erigen, headerbase)

print()
print("Header file: {}".format(headerfile))
print()

# Start the header file
defineline = re.sub('[\W]', '_', headerbase.upper())
with open(headerfile, 'w') as hfile:
  hfile.write("#ifndef {}\n".format(defineline))
  hfile.write("#define {}\n".format(defineline))
  hfile.write("\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("extern \"C\" {\n")
  hfile.write("#endif\n")
  hfile.write("\n\n")
  hfile.write("#include \"simint/vectorization/vectorization.h\"\n")
  hfile.write("\n")


for q in valid:
    filebase = "vrr_{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
    outfile = os.path.join(outdir_erigen, filebase + ".c")
    logfile = os.path.join(outdir_erigen, filebase + ".log")
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
  hfile.write("#endif // {}\n".format(defineline))
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




# Generate the eri
headerbase = "eri_generated.h"
headerfile = os.path.join(outdir_erigen, headerbase)

# Maximum required contwork
maxworksize = 0  # number of elements

print()
print("Header file: {}".format(headerfile))
print()

# Start the header file
defineline = re.sub('[\W]', '_', headerbase.upper())
with open(headerfile, 'w') as hfile:
  hfile.write("#ifndef {}\n".format(defineline))
  hfile.write("#define {}\n".format(defineline))
  hfile.write("\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("extern \"C\" {\n")
  hfile.write("#endif\n")
  hfile.write("\n\n")
  hfile.write("#include \"simint/vectorization/vectorization.h\"\n")
  hfile.write("#include \"simint/shell/shell.h\"\n")
  hfile.write("\n")


for q in valid:
  filebase = "eri_{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
  outfile = os.path.join(outdir_erigen, filebase + ".c")
  logfile = os.path.join(outdir_erigen, filebase + ".log")
  print("Creating: {}".format(filebase))
  print("      Output: {}".format(outfile))
  print("     Logfile: {}".format(logfile))


  if sum(q) >= args.ve:
    vrrtype = "External"
  else:
    vrrtype = "Inline";
  print("         VRR: {}".format(vrrtype))

  if not args.et:
    ettype = "Skipping"
  elif sum(q) >= args.ee:
    ettype = "External"
  else:
    ettype = "Inline";
  print("          ET: {}".format(ettype))


  if sum(q) >= args.he:
    hrrtype = "External"
  else:
    hrrtype = "Inline";
  print("         HRR: {}".format(hrrtype))


  with open(logfile, 'w') as lf:
    cmdline = [eri_gen];
    cmdline.extend(["-q", str(q[0]), str(q[1]), str(q[2]), str(q[3])])
    cmdline.extend(["-b", args.b])
    cmdline.extend(["-o", outfile])
    cmdline.extend(["-oh", headerfile])
    cmdline.extend(["-s", str(args.s)])
    cmdline.extend(["-d", args.d])

    if args.c:
        cmdline.extend(["-c", str(args.c)])
    if vrrtype == "External":
        cmdline.append("-ve")
    if ettype == "External":
        cmdline.append("-ee")
    if hrrtype == "External":
        cmdline.append("-he")

    if args.i:
        cmdline.append("-i")
    if args.S:
        cmdline.append("-S")
    if args.et:
        cmdline.append("-et")

    print()
    print("Command line:")
    print(' '.join(cmdline))
    print()

    ret = subprocess.call(cmdline, stdout=lf, stderr=lf)

    if ret != 0:
      print("\n")
      print("*********************************")
      print("While generating eri")
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
  hfile.write("\n\n")
  hfile.write("#define SIMINT_ERI_MAXAM {}\n".format(args.l))
  hfile.write("#define SIMINT_ERI_MAX_WORKSIZE ((SIMINT_SIMD_ROUND(SIMINT_NSHELL_SIMD * {})))\n".format(maxworksize))
  hfile.write("#define SIMINT_ERI_MAX_WORKMEM (SIMINT_ERI_MAX_WORKSIZE * sizeof(double))\n")
  hfile.write("\n\n")
  hfile.write("\n")
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("}\n")
  hfile.write("#endif\n")
  hfile.write("\n")
  hfile.write("#endif // {}\n".format(defineline))
  hfile.write("\n")



####################################################
# Overall config header
####################################################
sinfofile = os.path.join(outdir, "simint_config.h")
with open(sinfofile, 'w') as sf:

  # Parse the cpu flags
  cpuflags = args.c.lower().split(',')


  sf.write("#pragma once\n")

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

  if args.et:
    sf.write("#define SIMINT_ERI_USE_ET\n")
  else:
    sf.write("#define SIMINT_ERI_NO_ET\n")

  if args.p:
    sf.write("#define SIMINT_ERI_PERMUTATIONS\n")
  else:
    sf.write("#define SIMINT_ERI_NO_PERMUTATIONS\n")

  sf.write("\n\n")



