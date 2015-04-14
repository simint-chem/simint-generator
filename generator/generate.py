#!/usr/bin/env python3

import twogen as tg
import argparse

###################################
# Main generator for ERI source
###################################

parser = argparse.ArgumentParser()
parser.add_argument("q",          type=str,                               help="AM quartet string (comma separated)")
#parser.add_argument("--output",   type=str, required=True,                help="Base output directory")
parser.add_argument("-c",         action='store_true', default=False,     help="Print code lines")
args = parser.parse_args()


ams = [ int(x) for x in args.q.split(',') ]
inittargets = tg.MakeInitialTargetQuartets(ams) 
print("Created {} initial target quartets".format(len(inittargets)))
for q in inittargets:
  print(q)


###################################################
# Working backwards
###################################################
###################################################
# 1.) Horizaontal recurrance
#     ( a+b 0 | c+d  ) -> ( a b | c d )
###################################################
hrrlst = tg.Makowski_HRR(inittargets)
tg.HRRSanityCheck(inittargets, hrrlst)

# Now print out the main requirements
# from the HRR step
solvedtargets = []
for h in hrrlst:
  solvedtargets.append(h.target)

topreq = [ x for x in inittargets if (x.bra.right.am == 0 and x.ket.right.am == 0) ]
for h in hrrlst:
  if not h.src1 in topreq and not h.src1 in solvedtargets:
    topreq.append(h.src1)
  if not h.src2 in topreq and not h.src2 in solvedtargets:
    topreq.append(h.src2)

# sort the requirements
topreq = sorted(topreq)
print()
print("---------------------------------------------")
print("Top-level requirements: {}".format(len(topreq)))
print("---------------------------------------------")
for t in topreq:
  print(t)

print()
print("---------------------------------------------")
print("Steps: {}".format(len(hrrlst)))
print("FLOPS: {}".format(2*len(hrrlst)))
print("---------------------------------------------")
for h in hrrlst:
  print(h)

print()

if args.c:
  print("---------------------------------------------")
  print("HRR Code")
  print("---------------------------------------------")
  for h in hrrlst:
    print(h.Code_line())



###################################################
# 2.) Electron Transfer
#     ( a+b+c+d 0 | 0 0 ) -> ( a+b 0 | c+d 0 )
###################################################



###################################################
# 3.) Initial generation
#     ( 0 0 | 0 0 ) > ( a+b+c+d 0 | 0 0 )
###################################################




print()
