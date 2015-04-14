#!/usr/bin/env python3

from .Gaussian import *
from .Doublet import *
from .Quartet import *
from .HRRStep import *

def MakeInitialTargetQuartets(amlst):
    targets = []
    am1,am2,am3,am4 = amlst
    nam1 = ((am1 + 1) * (am1 + 2)) // 2
    nam2 = ((am2 + 1) * (am2 + 2)) // 2
    nam3 = ((am3 + 1) * (am3 + 2)) // 2
    nam4 = ((am4 + 1) * (am4 + 2)) // 2

    cur1 = Gaussian(am1, 0, 0)
    for i in range(0, nam1):
      cur2 = Gaussian(am2, 0, 0)

      for j in range(0, nam2):
        bra = Doublet("bra", cur1, cur2) 

        cur3 = Gaussian(am3, 0, 0)

        for k in range(0, nam3):
          cur4 = Gaussian(am4, 0, 0)

          for l in range(0, nam4):
            ket = Doublet("ket", cur3, cur4)
            targets.append(Quartet(bra, ket))

            cur4 = cur4.Iterate()

          cur3 = cur3.Iterate()

        cur2 = cur2.Iterate()        

      cur1 = cur1.Iterate()

    return targets


##################################################################
# Sanity checks
# 1.) All requirements are in the list or are in the form (X s|
# 2.) Everything that is solved for is actually used
##################################################################

def SanityCheck(inittargets, hrrlst):
    allreq = [ x for x in inittargets if (x.bra.right.am > 0 or x.ket.right.am > 0) ]
    solvedtargets = []
    for h in hrrlst:
      solvedtargets.append(h.target)
      if not h.src1 in allreq and (h.src1.bra.right.am > 0 or h.src1.ket.right.am > 0):
        allreq.append(h.src1)
      if not h.src2 in allreq and (h.src2.bra.right.am > 0 or h.src2.ket.right.am > 0):
        allreq.append(h.src2)

    notfound = []
    for r in allreq:
      if not r in solvedtargets:
        notfound.append(r)

    notneeded = []
    for s in solvedtargets:
      if not s in allreq and not s in inittargets:
        notneeded.append(s)

    if len(notfound) > 0 or len(notneeded) > 0:
      print("############################################################")
      print("Allreq:")
      for r in allreq:
        print(r)
      print("Solved:")
      for s in solvedtargets:
        print(s)
      print("Not found:")
      for n in notfound:
        print(n)
      print("Not needed:")
      for n in notneeded:
        print(n)

      if len(notfound) > 0:
        raise RuntimeError("Error: Some unresolved requirements!")
      if len(notneeded) > 0:
        raise RuntimeError("Error: Some stuff calculated that is not needed!")
