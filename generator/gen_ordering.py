#!/usr/bin/env python3

maxam = 24

class Gaussian:
  def __init__(self, am):
    self.am = am
    self.ijk = [am, 0, 0]

  def __str__(self):
    return "{{ {}, {}, {}}}".format(self.ijk[0], self.ijk[1], self.ijk[2])


  def IsValid(self):
    return (self.ijk[0] >= 0 and self.ijk[0] <= self.am and
            self.ijk[1] >= 0 and self.ijk[1] <= self.am and
            self.ijk[2] >= 0 and self.ijk[2] <= self.am and
            self.ijk[0] + self.ijk[1] + self.ijk[2] == self.am)


def IteratePsi4(g):
  newg = Gaussian(g.am)

  if g.ijk[2] == g.am:
    newg.ijk = [-1, -1, -1]
  elif g.ijk[2] < (g.am - g.ijk[0]):
    newg.ijk = [g.ijk[0], g.ijk[1]-1, g.ijk[2]+1]
  else:
    newg.ijk = [g.ijk[0]-1, g.am-g.ijk[0]+1, 0]

  return newg



def PrintMap(gmap, maxam): 

  print("std::map<const int, const std::vector<ExpList>> gorder_map = {")

  for i in range(0, maxam+1):
    print("{{ {:2}, {{".format(i))

    for g in gmap[i]:
      print("        " + str(g) + ",")
    print("      }")
    print("},")
  print("};")




gmap = {}
for i in range(0, maxam+1):
  g = Gaussian(i)
  gmap[i] = [g]

  g = IteratePsi4(g)
  while g.IsValid(): 
    gmap[i].append(g)
    g = IteratePsi4(g)


PrintMap(gmap, maxam)
