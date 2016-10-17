from .gaussian import Gaussian

def IteratePsi4(g):
  newg = Gaussian(g.am)

  if g.ijk[2] == g.am:
    newg.ijk = [-1, -1, -1]
  elif g.ijk[2] < (g.am - g.ijk[0]):
    newg.ijk = [g.ijk[0], g.ijk[1]-1, g.ijk[2]+1]
  else:
    newg.ijk = [g.ijk[0]-1, g.am-g.ijk[0]+1, 0]

  return newg
