
class Gaussian:
  def __init__(self, am):
      self.am = am
      self.ijk = [am, 0, 0]

  def __str__(self):
      return "{{ {:3}, {:3}, {:3}}}".format(self.ijk[0], self.ijk[1], self.ijk[2])


  def hashstr(self):
      i = 0
      i = i | (self.ijk[0] << 0x00)
      i = i | (self.ijk[1] << 0x08)
      i = i | (self.ijk[2] << 0x10)
      return str(i)

  def IsValid(self):
    return (self.ijk[0] >= 0 and self.ijk[0] <= self.am and
            self.ijk[1] >= 0 and self.ijk[1] <= self.am and
            self.ijk[2] >= 0 and self.ijk[2] <= self.am and
            self.ijk[0] + self.ijk[1] + self.ijk[2] == self.am)

