#!/usr/bin/env python3


amchar = "SPDFGHIJKLMNOQRTUVWXYZABCE"

class Gaussian:
  """ A single cartesian gaussian """



  def __init__(self, ix, iy, iz):
    self.am = ix+iy+iz
    self.ix = ix
    self.iy = iy
    self.iz = iz



  def __str__(self):
    return "{}_{}{}{}".format(amchar[self.am], self.ix, self.iy, self.iz)


    
  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return (self.am == other.am and 
              self.ix == other.ix and 
              self.iy == other.iy and
              self.iz == other.iz)
    else:
      return False



  def __ne__(self, other):
    return not self.__eq__(other)



  def __lt__(self, other):
    # note the "reverse" ordering - I like x first of course
    if self.am < other.am:
      return True
    elif self.am == other.am:
      if self.ix > other.ix:
        return True
      elif self.ix == other.ix:
        if self.iy > other.iy:
          return True
        elif self.iy == other.iy:
          if self.iz > other.iz:
            return True

    return False
    
  def __gt__(self, other):
    if self.am > other.am:
      return True
    elif self.am == other.am:
      if self.ix < other.ix:
        return True
      elif self.ix == other.ix:
        if self.iy < other.iy:
          return True
        elif self.iy == other.iy:
          if self.iz < other.iz:
            return True

    return False



  def Iterate(self):
    if(self.iz < (self.am - self.ix)):
      return Gaussian(self.ix, self.iy - 1, self.iz + 1)
    else:
      return Gaussian(self.ix - 1, self.am - self.ix + 1, 0)
