#!/usr/bin/env python3

import twogen as tg
import argparse

###################################
# Main generator for ERI source
###################################

parser = argparse.ArgumentParser()
parser.add_argument("q",         type=str,                help="AM quartet string")
#parser.add_argument("--filename", type=str, required=True,                help="Output file name base (no extension)")
#parser.add_argument("-o",         type=int, required=True,                help="Order of polynomial")
#parser.add_argument("-s",         type=str, required=False,               help="Starting order list")
#parser.add_argument("-b",         type=int, required=False, default=1,    help="Number of bins")
#parser.add_argument("-t",         type=int, required=False, default=200,  help="Number of test points (per bin)")
#parser.add_argument("-f",         type=str, required=False, default="1",  help="Base factor (ie 1 = base 2, etc)")
#parser.add_argument("--dps",      type=int, required=False, default=256,  help="Decimal precision/sig figs to use/calculate")
args = parser.parse_args()


ams = [ int(x) for x in args.q.split(',') ]
inittargets = tg.MakeInitialTargetQuartets(ams) 
print("Created {} initial target quartets".format(len(inittargets)))
for q in inittargets:
  print(q)


# list of steps
hrrlst = tg.Makowski_HRR(inittargets)

# sanity
tg.SanityCheck(inittargets, hrrlst)

# Now print out the main requirements
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
print("Steps: {} steps".format(len(hrrlst)))
print("---------------------------------------------")
for h in hrrlst:
  print(h)

print()

