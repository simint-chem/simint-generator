from .gaussian import Gaussian

def GenerateMap(maxam, iterfunc):
  gmap = dict()

  for i in range(0, maxam+1):
    g = Gaussian(i)
    gmap[i] = [g]

    g = iterfunc(g)
    while g.IsValid(): 
      gmap[i].append(g)
      g = iterfunc(g)

  return gmap
