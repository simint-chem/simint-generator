#!/usr/bin/env python3

from .Gaussian import *
from .Doublet import *
from .Quartet import *
from .HRRStep import *

xyzcomp = ( "x", "y", "z" )

def Makowski_HRR_doublet(t, quartet):
   if t == 'bra':
     doublet = quartet.bra
   elif t == 'ket':
     doublet = quartet.ket
   else:
     raise RuntimeError("Error - unknown type \"{}\"".format(t))

   if doublet.right.am == 0:
     raise RuntimeError("Error - can't recurse to a {} s-function...".format(t))

   # Makowski: Recurse in the lowest non-zero angular
   # component of the second function

   # r will be the value
   # i will be the xyz-index (0=x, 1=y, 2=z)
   s = [ doublet.right.ix, doublet.right.iy, doublet.right.iz ]
   ss = sorted(s)
   r = [ x for x in ss if x > 0 ][0]
   i = 2 - s[::-1].index(r)  # find last occurance of r

   # Now we have the direction
   # figure out the rest
   g1_lst = [doublet.left.ix, doublet.left.iy, doublet.left.iz]
   g2_lst = [doublet.right.ix, doublet.right.iy, doublet.right.iz]
   g1_lst[i] += 1
   g2_lst[i] -= 1
   g1 = Gaussian(g1_lst[0], g1_lst[1], g1_lst[2])
   g2 = Gaussian(g2_lst[0], g2_lst[1], g2_lst[2])
   src1d = Doublet(t, g1, g2)

   g1_lst = [doublet.left.ix, doublet.left.iy, doublet.left.iz]
   g2_lst = [doublet.right.ix, doublet.right.iy, doublet.right.iz]
   g2_lst[i] -= 1
   g1 = Gaussian(g1_lst[0], g1_lst[1], g1_lst[2])
   g2 = Gaussian(g2_lst[0], g2_lst[1], g2_lst[2])
   src2d = Doublet(t, g1, g2)

   if t == 'bra':
     src1q = Quartet(src1d, quartet.ket, quartet.m)
     src2q = Quartet(src2d, quartet.ket, quartet.m)
   else:
     src1q = Quartet(quartet.bra, src1d, quartet.m)
     src2q = Quartet(quartet.bra, src2d, quartet.m)

   xyz = xyzcomp[i]

   return HRRStep(t, quartet, src1q, xyz, src2q)


def NewTargets(t, hrrlst):
  solvedtargets = [ h.target for h in hrrlst ]
  newtargets = []
  for h in hrrlst:
    newtargets.append(h.src1)
    newtargets.append(h.src2)

  # clean out duplicate targets, targets that have been done, and s targets
  # (s on the right-hand of the bra or ket)
  targets = []

  for nt in newtargets:
    if not nt in targets and not nt in solvedtargets:
      if t == 'bra' and nt.bra.right.am > 0:
        targets.append(nt)
      if t == 'ket' and nt.ket.right.am > 0:
        targets.append(nt)

  # probably not necessary, but leads to nice ordering
  targets = sorted(targets)
  targets.reverse()

  return targets


def InitTargets(t, inittargets):
  # clean out duplicate targets, targets that have been done, and s targets
  # (s on the right-hand of the bra or ket)
  targets = []

  for nt in inittargets:
    if not nt in targets:
      if t == 'bra' and nt.bra.right.am > 0:
        targets.append(nt)
      if t == 'ket' and nt.ket.right.am > 0:
        targets.append(nt)

  # probably not necessary, but leads to nice ordering
  targets = sorted(targets)
  targets.reverse()

  return targets


def Makowski_HRR(inittargets):

  hrrlst = []
  targets = InitTargets('bra', inittargets)

  ##############
  # bras first
  ##############
  while len(targets) > 0:
    # find next steps for each target
    for t in targets:
      #prepend
      hrrlst.insert(0, Makowski_HRR_doublet('bra', t))

    targets = NewTargets('bra', hrrlst)


  ##############
  # now do kets
  ##############
  targets = NewTargets('ket', hrrlst)
  if len(hrrlst) == 0:
    # no HRR done yet
    # form of (ss|XY)
    targets = InitTargets('ket', inittargets)

  while len(targets) > 0:
    # find next steps for each target
    for t in targets:
      #prepend
      hrrlst.insert(0, Makowski_HRR_doublet("ket", t))

    targets = NewTargets('ket', hrrlst)

  return hrrlst
