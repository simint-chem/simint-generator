#!/usr/bin/env python3

class Quartet:
  """ Two doublets """



  def __init__(self, bra, ket):
    self.am = bra.am + ket.am
    self.bra = bra
    self.ket = ket



  def __str__(self):
    return "( {} {} | {} {} )".format(str(self.bra.left), str(self.bra.right),
                                      str(self.ket.left), str(self.ket.right))


  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return (self.bra == other.bra and 
              self.ket == other.ket and 
              self.am == other.am)
    else:
      return False


  def __lt__(self, other):
    if self.am < other.am:
      return True
    elif self.am == other.am:
      if self.bra < other.bra:
        return True
      elif self.bra == other.bra:
        if self.ket < other.ket:
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

    return False



  def __ne__(self, other):
    return not self.__eq__(other)
