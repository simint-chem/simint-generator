#!/usr/bin/env python3

class HRRStep:
  """A class holding a single step of the horizontal recurrance relation"""


  def __init__(self, t, target, src1, xyz, src2):
    self.target = target
    self.src1 = src1
    self.src2 = src2
    self.xyz = xyz
    self.t = t


  def __str__(self):
    return "{} = {} + {}{}".format(str(self.target), str(self.src1), str(self.xyz), str(self.src2))


  def __eq__(self, other):
    #if isinstance(other, self.__class__):
      return (self.target == other.target and 
              self.src1 == other.src1 and 
              self.src2 == other.src2 and
              self.xyz == other.xyz and
              self.t == other.t)
    #else:
    #  return False



  def __ne__(self, other):
    return not self.__eq__(other)


  def Code_line(self):

    comp = self.xyz.upper()

    if self.t == 'bra':
      sub = "ab"
    elif self.t == 'ket':
      sub = "cd"
    else:
      raise RuntimeError("Error - unknown type {}".format(t))

    s = "const double {} = {} + {}{} * {};".format(self.target.Code_variable(),
                                                   self.src1.Code_variable(),
                                                   comp, sub,
                                                   self.src2.Code_variable())

    return s
    
