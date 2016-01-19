#!/usr/bin/env python3


import argparse
import sys
import os
import subprocess
import re
import shutil


def ValidQuartet(q):
  return True
  if q[0] < q[1]:
    return False

  if q[2] < q[3]:
    return False

  #if  ( q[0] + q[1] ) < ( q[2] + q[3] ):
  #  return False

  #if q[0] < q[2]:
  #  return False 

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
parser.add_argument("-c", type=str, required=True, default="",   help="CPU Flags (comma separated)")

parser.add_argument("-et", required=False, action='store_true', help="Use electron transfer rather than second VRR")
parser.add_argument("-ve", required=False, type=int, default=1000, help="External VRR for this L value and above")
parser.add_argument("-ee", required=False, type=int, default=1000, help="External ET for this L value and above")
parser.add_argument("-he", required=False, type=int, default=1000, help="External HRR for this L value and above")
parser.add_argument("-etvrr", required=False, action='store_true', help="Use VRR for ( s s | X s ) rather than ET")
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
hrr_gen = os.path.join(args.g, "hrr_generator")
et_gen = os.path.join(args.g, "et_generator")
vrr_gen = os.path.join(args.g, "vrr_generator")

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
outdir_eri = os.path.join(args.outdir, "eri")
outdir_erigen = os.path.join(outdir_eri, "gen")
outdir_vec = os.path.join(args.outdir, "vectorization")
outdir_boys = os.path.join(args.outdir, "boys")
outdir_shell = os.path.join(args.outdir, "shell")
outdir_test = os.path.join(args.outdir, "test")

if os.path.isdir(args.outdir):
  print("WARNING - output directory exists. Overwriting...")
  shutil.rmtree(outdir_eri, ignore_errors=True)
  shutil.rmtree(outdir_vec, ignore_errors=True)
  shutil.rmtree(outdir_boys, ignore_errors=True)
  shutil.rmtree(outdir_shell, ignore_errors=True)
  shutil.rmtree(outdir_test, ignore_errors=True)


shutil.copytree(os.path.join(skeldir, "eri"),           outdir_eri) 
shutil.copytree(os.path.join(skeldir, "vectorization"), outdir_vec)
shutil.copytree(os.path.join(skeldir, "boys"),          outdir_boys)
shutil.copytree(os.path.join(skeldir, "shell"),         outdir_shell)
shutil.copytree(os.path.join(skeldir, "test"),          outdir_test)
shutil.copy(os.path.join(skeldir, "CMakeLists.txt"),    args.outdir)
shutil.copy(os.path.join(skeldir, "constants.h"),       args.outdir)



####################################################
# Generate the external HRR source
####################################################
logfile = os.path.join(outdir_erigen, "hrr.log")

cmdline = [hrr_gen]
cmdline.extend(["-c", str(args.c)])
cmdline.extend(["-o", outdir_erigen])
cmdline.extend(["-L", str(args.l)])
if args.i:
    cmdline.append("-i")
if args.S:
    cmdline.append("-S")
if args.et:
    cmdline.append("-et")
if args.etvrr:
    cmdline.append("-etvrr")

print("Creating HRR sources in {}".format(outdir_erigen))
print("     Logfile: {}".format(logfile))
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



####################################################
# Generate the external ET source
####################################################
if args.et:
    logfile = os.path.join(outdir_erigen, "et.log")

    cmdline = [et_gen]
    cmdline.extend(["-c", str(args.c)])
    cmdline.extend(["-o", outdir_erigen])
    cmdline.extend(["-L", str(args.l*4)])
    if args.i:
        cmdline.append("-i")
    if args.S:
        cmdline.append("-S")
    if args.et:
        cmdline.append("-et")
    if args.etvrr:
        cmdline.append("-etvrr")

    print("Creating ET sources in {}".format(outdir_erigen))
    print("     Logfile: {}".format(logfile))
    print()
    print("Command line:")
    print(' '.join(cmdline))
    print()

    with open(logfile, 'w') as lf:
      ret = subprocess.call(cmdline, stdout=lf)

    if ret != 0:
      print("\n")
      print("*********************************")
      print("When generating et sources")
      print("Subprocess returned {} - aborting".format(ret))
      print("*********************************")
      print("\n")
      quit(1)




####################################################
# Generate the external VRR source
####################################################
logfile = os.path.join(outdir_erigen, "vrr.log")

if args.et:
  if args.etvrr:
    maxL = args.l*2
  else:
    maxL = args.l*4
else:
  maxL = args.l*2


cmdline = [vrr_gen]
cmdline.extend(["-c", str(args.c)])
cmdline.extend(["-o", outdir_erigen])
cmdline.extend(["-L", str(maxL)])
if args.i:
    cmdline.append("-i")
if args.S:
    cmdline.append("-S")
if args.et:
    cmdline.append("-et")
if args.etvrr:
    cmdline.append("-etvrr")

print("Creating VRR sources in {}".format(outdir_erigen))
print("     Logfile: {}".format(logfile))
print()
print("Command line:")
print(' '.join(cmdline))
print()

with open(logfile, 'w') as lf:
  ret = subprocess.call(cmdline, stdout=lf)

if ret != 0:
  print("\n")
  print("*********************************")
  print("When generating vrr sources")
  print("Subprocess returned {} - aborting".format(ret))
  print("*********************************")
  print("\n")
  quit(1)






####################################################
# Generate the ERI sources and headers
####################################################
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




# Generate the eri
headerbase = "eri_generated.h"
headerfile = os.path.join(outdir_erigen, headerbase)

print()
print("Header file: {}".format(headerfile))
print()
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


  if sum(q) >= args.ee:
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
    cmdline.extend(["-s", str(args.s)])
    cmdline.extend(["-c", str(args.c)])
    cmdline.extend(["-d", args.d])

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
    if args.etvrr:
        cmdline.append("-etvrr")

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
  print()


# Create the header file for all eri
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

  hfile.write("#define ERI_MAXAM {}\n".format(args.l))
  hfile.write("\n\n")

  for q in valid:
    funcname = "eri_{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
    intname = "INT__{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
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
  funcname = "eri_{}_{}_{}_{}".format(amchar[q[0]], amchar[q[1]], amchar[q[2]], amchar[q[3]])
  print("funcs[{}][{}][{}][{}] = {};".format(q[0], q[1], q[2], q[3], funcname))
print("\n\n")
