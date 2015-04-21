#!/usr/bin/env python3

class Doublet:
  """ A bra or a ket part (2 gaussians ) """



  def __init__(self, t, g_left, g_right):
    self.left = g_left
    self.right = g_right
    self.t = t
    self.am = g_left.am + g_right.am



  de __str__(self):
    if self.t == "bra":
      return "( {} {} |".format(str(self.left), str(self.right)) 
    elif self.t == "ket":
      return "| {} {} )".format(str(self.left), str(self.right)) 
    else:
      raise RuntimeError("Unknown doublet type \"{}\"".format(self.t))


  def __hash__(self):
    return hash((self.left, self.right, self.t))  # am isn't really unique


  def __eq__(self, other):
    #if isinstance(other, self.__class__):
      return (self.left == other.left and 
              self.right == other.right and 
              self.t == other.t and
              self.am == other.am)
    #else:
    #  return False



  def __lt__(self, other):
    if self.t != other.t:
      raise RuntimeError("Error - can't compare a bra with a ket!")
    if self.am < other.am:
      return True
    elif self.am == other.am:
      if self.left < other.left:
        return True
      elif self.left == other.left:
        if self.right < other.right:
          return True

    return False



  def __gt__(self, other):
    if self.t != other.t:
      raise RuntimeError("Error - can't compare a bra with a ket!")
    if self.am > other.am:
      return True
    elif self.am == other.am:
      if self.left > other.left:
        return True
      elif self.left == other.left:
        if self.right > other.right:
          return True

    return False


  def __ne__(self, other):
    return not self.__eq__(other)
