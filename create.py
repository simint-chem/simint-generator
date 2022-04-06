#!/usr/bin/env python3


import argparse
import sys
import os
import subprocess
import re
import shutil
from multiprocessing import Pool, set_start_method


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


def QStr(q):
  return "( {} {} | {} {} )".format(q[0], q[1], q[2], q[3])


# The following function allows the use of multiprocess package
def runostei(arg):
  i,q = arg
  filebase = base+"_{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
  outfile = os.path.join(outdir_osteigen, filebase + ".c")
  logfile = os.path.join(outdir_osteigen, filebase + ".log")
  print("Creating: {}\n      Output: {}\n     Logfile: {}\n     {}/{} ".format(filebase,outfile,logfile,i+1,tot))

  with open(logfile, 'w') as lf:
    cmdline = [cmd];
    cmdline.extend(["-q", str(q[0]), str(q[1]), str(q[2]), str(q[3])])
    cmdline.extend(["-o", outfile])
    cmdline.extend(["-oh", headerfile])
    cmdline.extend(["-ve", str(args.ve)])
    cmdline.extend(["-vg", str(args.vg)])
    cmdline.extend(["-he", str(args.he)])
    cmdline.extend(["-hg", str(args.hg)])

    if max(q) >= args.p:
        cmdline.append("-p")

    print()
    print("Command line:")
    print(' '.join(cmdline))
    print()

    ret = subprocess.call(cmdline, stdout=lf, stderr=lf)

    if ret != 0:
      print("\n")
      print("*********************************")
      print("While generating "+base)
      print("Subprocess returned {} - aborting".format(ret))
      print("*********************************")
      print("\n")
      quit(5)

    return



# path to this file
thisfile = os.path.realpath(__file__)
topdir = os.path.dirname(thisfile)

# some helpers
amchar = "spdfghiklmnoqrtuvwxyzabceSPDFGHIKLMNOQRTUVWXYZABCE0123456789"

parser = argparse.ArgumentParser()
parser.add_argument("-n", type=int, required=False, default=1, help="Number of parallel processes to use with generator")
parser.add_argument("-g", type=str, required=True, help="Path to directory with generator programs")
parser.add_argument("-l", type=int, required=True, help="Maximum AM")
parser.add_argument("-p", type=int, required=True, help="Start permuting the slow way when AM of a center reaches this value")
parser.add_argument("-ve", required=False, type=int, default=0, help="External VRR for this L value and above")
parser.add_argument("-vg", required=False, type=int, default=0, help="General VRR for this L value and above")
parser.add_argument("-he", required=False, type=int, default=0, help="External HRR for this L value and above")
parser.add_argument("-hg", required=False, type=int, default=0, help="General HRR for this L value and above")
parser.add_argument("-d", required=False, type=int, default=0, help="Maximum derivative to generate code for")
parser.add_argument("outdir", type=str, help="Output directory")

args = parser.parse_args()

maxam = args.l
derorder = args.d


###################################
# Actual code starts here
###################################
# directory with the files to copy
skeldir = os.path.join(topdir, "skel")


# paths to generator programs
ostei_gen = os.path.join(args.g, "ostei_generator")
ostei_deriv1_gen = os.path.join(args.g, "ostei_deriv1_generator")
hrr_gen = os.path.join(args.g, "ostei_hrr_generator")
vrr_gen = os.path.join(args.g, "ostei_vrr_generator")

if not os.path.isfile(ostei_gen):
  print("The file \"{}\" does not exist or is not a (binary) file".format(ostei_gen))
  quit(1)

if not os.path.isfile(ostei_deriv1_gen):
  print("The file \"{}\" does not exist or is not a (binary) file".format(ostei_deriv1_gen))
  quit(1)

if not os.path.isfile(hrr_gen):
  print("The file \"{}\" does not exist or is not a (binary) file".format(hrr_gen))
  quit(1)

if not os.path.isfile(vrr_gen):
  print("The file \"{}\" does not exist or is not a (binary) file".format(vrr_gen))
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


shutil.copytree(os.path.join(skeldir, "simint"),            outdir)
shutil.copytree(os.path.join(skeldir, "test"),              outdir_test)
shutil.copytree(os.path.join(skeldir, "examples"),          outdir_examples)
shutil.copytree(os.path.join(skeldir, "cmake"),             outdir_cmake)
shutil.copy(os.path.join(skeldir, "CMakeLists.txt"),        args.outdir)
shutil.copy(os.path.join(skeldir, "README"),                args.outdir)
shutil.copy(os.path.join(skeldir, "LICENSE"),               args.outdir)
shutil.copy(os.path.join(skeldir, "CHANGELOG"),             args.outdir)

# We need to create the ostei/gen directory, which
# is empty and therefore not included in the git
# repo (under skel)
if not os.path.exists(outdir_osteigen):
    os.makedirs(outdir_osteigen)


####################################################
# Generate the ERI sources and headers
####################################################
print("-------------------------------")
print("Generating ERI")
print("Maximum AM: {}".format(maxam))
print("-------------------------------")

valid = set()
invalid = set()

print()

for i in range(0, maxam + 1):
  for j in range(0, maxam + 1):
    for k in range(0, maxam + 1):
      for l in range(0, maxam + 1):
        q = (i,j,k,l)
        if UniqueQuartet(q) or max(q) < args.p:
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

# Maximum required work
worksize_bcont = [[0]*(maxam+1) for _ in range(derorder+1)]
worksize_cont  = [[0]*(maxam+1) for _ in range(derorder+1)]
worksize_prim  = [[0]*(maxam+1) for _ in range(derorder+1)]

# Required external HRR and VRR
reqext_hrr = []
reqext_vrr = []

# This keeps track of which L values require what
# files (ERI, VRR, HRR)
filelists = [{ k: [] for k in range(0, maxam+1) }]*(derorder+1)

# This keeps track of the buffer sizes (also by derivative)
memreq = [{}]*(derorder+1)

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

# Define some global variables for all processes
global cmd
global base
global tot
cmd = ostei_gen
base = "ostei"
tot = len(valid)

# Run generator in args.n processes
set_start_method('fork')
with Pool(processes=args.n) as pool:
  pool.map(runostei, [(i,q) for i,q in enumerate(valid)], chunksize=1)

# Update arrays after multiprocess for simplicity
for q in valid:
  filebase = "ostei_{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
  logfile = os.path.join(outdir_osteigen, filebase + ".log")
  filelists[0][max(q)].append(filebase + ".c")

  # reopen the logfile, find work and requirements
  for line in open(logfile, 'r').readlines():
    mq = max(q)
    if line.startswith("WORK SIZE"):
      worksize_bcont[0][mq] = max(worksize_bcont[0][mq], int(line.split()[2]))
      worksize_prim[0][mq] = max(worksize_prim[0][mq], int(line.split()[3]))
      worksize_cont[0][mq] = max(worksize_cont[0][mq], int(line.split()[4]))
    elif line.startswith("SIMINT EXTERNAL HRR"):
      reqam = tuple(line.split()[3:])
      reqext_hrr.append(reqam)
      filelists[0][mq].append( "hrr_{}_{}_{}.c".format(reqam[0],
                                                       amchar[int(reqam[1])],
                                                       amchar[int(reqam[2])]) )
    elif line.startswith("SIMINT EXTERNAL VRR"):
      reqam = tuple(line.split()[3:])
      reqext_vrr.append(reqam)
      filelists[0][mq].append( "vrr_{}_{}_{}_{}_{}.c".format(reqam[0],
                                                             amchar[int(reqam[1])],
                                                             amchar[int(reqam[2])],
                                                             amchar[int(reqam[3])],
                                                             amchar[int(reqam[4])]) )

  print()


# Close out the header file

with open(headerfile, 'a') as hfile:
  hfile.write("#ifdef __cplusplus\n")
  hfile.write("}\n")
  hfile.write("#endif\n")
  hfile.write("\n")

# Overwrite the manually optimized file
shutil.copy(os.path.join(skeldir, "ostei_s_s_s_s.c"), outdir_osteigen)
print("***** Overwrite <OUTDIR>\simint\ostei\gen\ostei_s_s_s_s.c with skel\ostei_s_s_s_s.c *****")

####################################################
# Generate the ERI 1st derivative sources and headers
####################################################
maxam1 = -1
if derorder > 0:
  maxam1 = maxam

print("-------------------------------")
print("Generating ERI 1st Derivatives")
print("Maximum AM: {}".format(maxam1))
print("-------------------------------")

valid = set()
invalid = set()

print()

for i in range(0, maxam1 + 1):
  for j in range(0, maxam1 + 1):
    for k in range(0, maxam1 + 1):
      for l in range(0, maxam1 + 1):
        q = (i,j,k,l)
        if max(q) > maxam1:
          continue
        if UniqueQuartet(q) or max(q) < args.p:
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


# Generate the ostei 1st derivatives
headerbase = "ostei_deriv1_generated.h"
headerfile = os.path.join(outdir_osteigen, headerbase)

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


# Define some global variables for all processes
cmd = ostei_deriv1_gen
base = "ostei_deriv1"
tot = len(valid)

# Run generator in args.n processes
with Pool(processes=args.n) as pool:
    pool.map(runostei, [(i,q) for i,q in enumerate(valid)], chunksize=1)

# Update lists after multiprocess for simplicity
for q in valid:
  filebase = "ostei_deriv1_{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
  logfile = os.path.join(outdir_osteigen, filebase + ".log")
  filelists[1][max(q)].append(filebase + ".c")

  # reopen the logfile, find work and requirements
  for line in open(logfile, 'r').readlines():
    mq = max(q)
    if line.startswith("WORK SIZE"):
      worksize_bcont[1][mq] = max(worksize_bcont[1][mq], int(line.split()[2]))
      worksize_prim[1][mq] = max(worksize_prim[1][mq], int(line.split()[3]))
      worksize_cont[1][mq] = max(worksize_cont[1][mq], int(line.split()[4]))
    elif line.startswith("SIMINT EXTERNAL HRR"):
      reqam = tuple(line.split()[3:])
      reqext_hrr.append(reqam)
      filelists[1][mq].append( "hrr_{}_{}_{}.c".format(reqam[0],
                                                       amchar[int(reqam[1])],
                                                       amchar[int(reqam[2])]) )
    elif line.startswith("SIMINT EXTERNAL VRR"):
      reqam = tuple(line.split()[3:])
      reqext_vrr.append(reqam)
      filelists[1][mq].append( "vrr_{}_{}_{}_{}_{}.c".format(reqam[0],
                                                             amchar[int(reqam[1])],
                                                             amchar[int(reqam[2])],
                                                             amchar[int(reqam[3])],
                                                             amchar[int(reqam[4])]) )

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
headerbase = "ostei_config.h.in"
headerfile = os.path.join(outdir_ostei, headerbase)
with open(headerfile, 'w') as hfile:
  hfile.write("#pragma once\n\n")
  hfile.write("#include \"simint/vectorization/vectorization.h\"\n")
  hfile.write("\n")
  # These are filled in by cmake
  hfile.write("#define SIMINT_OSTEI_MAXAM @SIMINT_MAXAM@\n")
  hfile.write("#define SIMINT_OSTEI_MAXDER @SIMINT_MAXDER@\n")
  hfile.write("#define SIMINT_OSTEI_DERIV1_MAXAM (SIMINT_OSTEI_MAXDER > 0 ? @SIMINT_MAXAM@ : -1)\n")
  hfile.write("\n\n")

  # A function for determining the max work size
  hfile.write("static inline size_t simint_ostei_worksize(int derorder, int maxam)\n")
  hfile.write("{\n")
  hfile.write("    static const size_t nelements[{}][{}] = {{\n".format(derorder+1, maxam+1))

  for d in range(0, derorder+1):
    hfile.write("      {\n")
    for l in range(0, maxam+1):
      max_bcont = max(worksize_bcont[d][:l+1])
      max_cont = max(worksize_cont[d][:l+1])
      max_prim = max(worksize_prim[d][:l+1])
      hfile.write("        (SIMINT_SIMD_ROUND(SIMINT_NSHELL_SIMD*{}) + SIMINT_SIMD_ROUND({}) + SIMINT_SIMD_LEN*{}),\n".format(max_bcont, max_cont, max_prim))
    hfile.write("      },\n")  

  hfile.write("    };\n") # close static array 

  hfile.write("    return nelements[derorder][maxam];\n") 
  hfile.write("}\n\n") # close simint_max_worksize function

  hfile.write("static inline size_t simint_ostei_workmem(int derorder, int maxam)\n")
  hfile.write("{\n")
  hfile.write("    return simint_ostei_worksize(derorder, maxam)*sizeof(double);\n")
  hfile.write("}\n\n")
  hfile.write("\n")



####################################################
# Generate the external HRR source
####################################################
print("-------------------------------")
print("Generating HRR")
print("-------------------------------")


# What do we need
reqext_hrr = set(reqext_hrr)


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


for h in reqext_hrr:
    q = ( int(h[1]), int(h[2]) )
    filebase = "hrr_{}_{}_{}".format(h[0], amchar[q[0]], amchar[q[1]])

    outfile = os.path.join(outdir_osteigen, filebase + ".c")
    logfile = os.path.join(outdir_osteigen, filebase + ".log")
    print("Creating HRR: {}".format(filebase))
    print("      Output: {}".format(outfile))
    print("     Logfile: {}".format(logfile))
    cmdline = [hrr_gen]

    cmdline.extend(["-q", str(q[0]), str(q[1])])
    cmdline.extend(["-o", outfile])
    cmdline.extend(["-oh", headerfile])

    if h[0] == "I" or h[0] == "J":
      cmdline.append("-bra")

    if h[0] == "I":
      cmdline.append("-bra_i")

    if h[0] == "K":
      cmdline.append("-ket_k")

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

# What do we need
reqext_vrr_tmp = set(reqext_vrr)

# Add the permutations that will actually compute stuff
# These are of the form (X s | X s)
reqext_vrr = set()

for v in reqext_vrr_tmp:
  reqext_vrr.add(v)
  if v[1] == "0" and v[3] == "0":
    if v[0] == "J":
      reqext_vrr.add(("I", v[2], v[1], v[4], v[3]))
    else:
      reqext_vrr.add(("K", v[2], v[1], v[4], v[3]))
  elif v[1] == "0":
    if v[0] == "J":
      reqext_vrr.add(("I", v[2], v[1], v[3], v[4]))
    else:
      reqext_vrr.add(("K", v[2], v[1], v[3], v[4]))
  elif v[3] == "0":
    if v[0] == "L":
      reqext_vrr.add(("K", v[1], v[2], v[4], v[3]))
    else:
      reqext_vrr.add(("I", v[1], v[2], v[4], v[3]))


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


for v in reqext_vrr:
    q = ( int(v[1]), int(v[2]), int(v[3]), int(v[4]))
    
    filebase = "vrr_{}_{}_{}_{}_{}".format(v[0], amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
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

      if v[0] == "J":
        cmdline.append("-center_j");
      if v[0] == "K":
        cmdline.append("-center_k");
      if v[0] == "L":
        cmdline.append("-center_l");

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
# Create the separate cmake files for all AM up
# to the max
####################################################

for d in range(0, derorder+1):
    for l,r in filelists[d].items():
      # Make a list of all requirements up to this one
      allr = []
      for l2 in range(0, l+1):
        allr.extend(filelists[d][l2])

      # Write out to file
      fpath = os.path.join(outdir_osteigen, "filelist_d{}_{}".format(d, l))
      with open(fpath, "w") as f:
        for line in allr:
          f.write("ostei/gen/" + line + "\n")


# Make the default AM, highest derivative the same as the max
with open(os.path.join(outdir_cmake, "MaxAvailable.cmake"), 'w') as f:
    f.write("set(SIMINT_MAXDER_AVAILABLE {})\n".format(derorder))
    f.write("set(SIMINT_MAXAM_AVAILABLE {})\n".format(maxam))


