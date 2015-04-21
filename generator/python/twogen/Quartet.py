#!/usr/bin/env python3

class Quartet:
  """ Two doublets """



  def __init__(self, bra, ket, m):
    self.am = bra.am + ket.am
    self.bra = bra
    self.ket = ket
    self.m = m



  def __str__(self):
    return "( {} {} | {} {} )".format(str(self.bra.left), str(self.bra.right),
                                      str(self.ket.left), str(self.ket.right))


  def __hash__(self):
    return hash((self.bra, self.ket, self.m)) # am is not really unique


  def __eq__(self, other):
    #if isinstance(other, self.__class__):
      return (self.bra == other.bra and 
              self.ket == other.ket and 
              self.am == other.am and
              self.m == other.m)
    #else:
    #  return False


  def __lt__(self, other):
    if self.am < other.am:
      return True
    elif self.am == other.am:
      if self.bra < other.bra:
        return True
      elif self.bra == other.bra:
        if self.ket < other.ket:
          return True
        elif self.ket == other.ket:
          if self.m < other.m:
            return True

    return False


  def __gt__(self, other):
    if self.am > other.am:
      return True
    elif self.am == other.am:
      if self.bra > other.bra:
        return True
      elif self.bra == other.bra:
        if self.ket > other.ket:
          return True
        elif self.ket == other.ket:
          if self.m > other.m:
            return True

    return False



  def __ne__(self, other):
    return not self.__eq__(other)




  def Code_variable(self):
    return "{}_{}{}{}_{}{}{}_{}{}{}_{}{}{}_{}".format("Q",
                                                      self.bra.left.ix,  self.bra.left.iy,  self.bra.left.iz,
                                                      self.bra.right.ix, self.bra.right.iy, self.bra.right.iz,
                                                      self.ket.left.ix,  self.ket.left.iy,  self.ket.left.iz,
                                                      self.ket.right.ix, self.ket.right.iy, self.ket.right.iz,
                                                      self.m)
