#!/usr/bin/env python3

from .Gaussian import *
from .Doublet import *
from .Quartet import *
from .HRRStep import *

xyzcomp = ( "x", "y", "z" )

def Makowski_HRR_doublet(t, quartet, solvedtargets):
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

    # calculate new targets
    newtargets = set()
    if src1d.right.am > 0:
       newtargets.add(src1q)
    if src2d.right.am > 0:
       newtargets.add(src2q)

    newtargets -= solvedtargets
    solvedtargets.add(quartet) 

    xyz = xyzcomp[i]
 
    return (newtargets, solvedtargets, HRRStep(t, quartet, src1q, xyz, src2q))


# clean out duplicate targets, targets that have been done, and s targets
# (s on the right-hand of the bra or ket)
def NewTargets_bra(hrrlst):
  solvedtargets = { h.target for h in hrrlst }
  newtargets = set()

  for h in hrrlst:
    newtargets.add(h.src1)
    newtargets.add(h.src2)

  targets = { nt for nt in newtargets if nt.bra.right.am > 0 }
  targets -= solvedtargets

  # probably not necessary, but leads to nice ordering
  targets = sorted(list(targets))
  targets.reverse()

  return targets


def NewTargets_ket(hrrlst):
  solvedtargets = { h.target for h in hrrlst }
  newtargets = set()

  for h in hrrlst:
    newtargets.add(h.src1)
    newtargets.add(h.src2)

  targets = { nt for nt in newtargets if nt.ket.right.am > 0 }
  targets -= solvedtargets

  targets = sorted(list(targets))
  targets.reverse()

  return targets



def InitTargets_bra(inittargets):
  targets = { it for it in inittargets if it.bra.right.am > 0 }
  targets = sorted(list(targets))
  targets.reverse()
  return targets


def InitTargets_ket(inittargets):
  targets = { it for it in inittargets if it.ket.right.am > 0 }
  targets = sorted(list(targets))
  targets.reverse()
  return targets



def Makowski_HRR(inittargets):

  hrrlst = []
  solvedtargets = set()

  ##############
  # bras first
  ##############
  targets = InitTargets_bra(inittargets)
  while len(targets) > 0:
    newtargets = set()
    for t in targets:
      (new,sol,hrr) = Makowski_HRR_doublet('bra', t, solvedtargets)
      hrrlst.insert(0, hrr)
      newtargets |= new 
      solvedtargets |= sol
      
    targets = newtargets


  ##############
  # now do kets
  ##############
  if len(hrrlst) == 0:
    # no HRR done yet
    # form of (ss|XY)
    targets = InitTargets_ket(inittargets)
  else:
    targets = NewTargets_ket(hrrlst)

  while len(targets) > 0:
    newtargets = set()
    for t in targets:
      (new,sol,hrr) = Makowski_HRR_doublet('ket', t, solvedtargets)
      hrrlst.insert(0, hrr)
      newtargets |= new 
      solvedtargets |= sol
      
    targets = newtargets

  return hrrlst
