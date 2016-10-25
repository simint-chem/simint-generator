#!/usr/bin/env python3

#######################################
# Generates shell ordering
#######################################

import argparse
import sys
from ordering_common import Gaussian, IteratePsi4, GenerateMap


def FindIndex(glist, ijk):
  for i,g in enumerate(glist):
    if ijk == g.ijk:
      return i
  return -1


parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True,               help="Output file name base (no extension)")
parser.add_argument("--max-am",    type=int, required=True,               help="Maximum angular momentum value to go to")
args = parser.parse_args()
maxam = args.max_am

# we need one more for incrementing the very last AM
gmap = GenerateMap(maxam+1, IteratePsi4)
print(gmap.keys())

print("------------------------------------")
print("Options for gen_recurinfo:")
print("   Max am: {}".format(maxam))
print("------------------------------------")


# number of elements in the array
numel = 0
for k in range(0, maxam+1):  # Slice off last value
  numel += len(gmap[k])

# Output to file
with open(args.filename + ".c", 'w') as f:
  f.write("/*\n")
  f.write(" Generated with:\n")
  f.write("   " + " ".join(sys.argv[:]))
  f.write("\n")
  f.write("------------------------------------\n")
  f.write("Options for gen_recurinfo:\n")
  f.write("   Max am: {}\n".format(maxam))
  f.write("------------------------------------\n")
  f.write("*/\n\n")

  f.write("#include \"simint/ostei/recur_lookup.h\"\n\n")

  f.write("extern struct RecurInfo const recurinfo_array[{}] = {{\n".format(numel))

  starts = []
  count = 0

  for k in range(0, maxam+1):  # Slice off last value
    v = gmap[k]

    starts.append(count)
    count += len(v)

    for vi in v:
      # exponents
      f.write("    { ")
      f.write("{{ {:>2}, {:>2}, {:>2} }}, ".format(vi.ijk[0], vi.ijk[1], vi.ijk[2]))

      # What is the optimal direction to recurse down
      if k == 0:
        d = -1
      else:
        # smallest nonzero value, blah blah (see Makowski)
        srt = sorted(vi.ijk)
        m = next((x for i, x in enumerate(srt) if x > 0), None) # Smallest non-zero
        d = 2 - list(reversed(vi.ijk)).index(m)  # Favoring the right-hand side

      f.write("{:>2}, ".format(str(d)))

      if k == 0:
        idx_ijk_x_m1 = -1
        idx_ijk_y_m1 = -1
        idx_ijk_z_m1 = -1
        idx_ijk_x_m2 = -1
        idx_ijk_y_m2 = -1
        idx_ijk_z_m2 = -1
        idx_ijk_x_p1 =  0
        idx_ijk_y_p1 =  1
        idx_ijk_z_p1 =  2
      else:
        # decrement i, j, and k by 1, 2, and increment by 1
        ijk_x_m1 = [ vi.ijk[0]-1, vi.ijk[1]  , vi.ijk[2]   ]
        ijk_y_m1 = [ vi.ijk[0]  , vi.ijk[1]-1, vi.ijk[2]   ]
        ijk_z_m1 = [ vi.ijk[0]  , vi.ijk[1]  , vi.ijk[2]-1 ]
        ijk_x_m2 = [ vi.ijk[0]-2, vi.ijk[1]  , vi.ijk[2]   ]
        ijk_y_m2 = [ vi.ijk[0]  , vi.ijk[1]-2, vi.ijk[2]   ]
        ijk_z_m2 = [ vi.ijk[0]  , vi.ijk[1]  , vi.ijk[2]-2 ]
        ijk_x_p1 = [ vi.ijk[0]+1, vi.ijk[1]  , vi.ijk[2]   ]
        ijk_y_p1 = [ vi.ijk[0]  , vi.ijk[1]+1, vi.ijk[2]   ]
        ijk_z_p1 = [ vi.ijk[0]  , vi.ijk[1]  , vi.ijk[2]+1 ]

        # indices of the above
        idx_ijk_x_m1 = FindIndex(gmap[k-1], ijk_x_m1)
        idx_ijk_y_m1 = FindIndex(gmap[k-1], ijk_y_m1)
        idx_ijk_z_m1 = FindIndex(gmap[k-1], ijk_z_m1)
        idx_ijk_x_p1 = FindIndex(gmap[k+1], ijk_x_p1)
        idx_ijk_y_p1 = FindIndex(gmap[k+1], ijk_y_p1)
        idx_ijk_z_p1 = FindIndex(gmap[k+1], ijk_z_p1)

        if k > 1:
          idx_ijk_x_m2 = FindIndex(gmap[k-2], ijk_x_m2)
          idx_ijk_y_m2 = FindIndex(gmap[k-2], ijk_y_m2)
          idx_ijk_z_m2 = FindIndex(gmap[k-2], ijk_z_m2)
        else:
          idx_ijk_x_m2 = -1
          idx_ijk_y_m2 = -1
          idx_ijk_z_m2 = -1

      f.write("{ ")
      f.write("{{ {:>2}, {:>2}, {:>2} }}, ".format(idx_ijk_x_m1, idx_ijk_x_m2, idx_ijk_x_p1))
      f.write("{{ {:>2}, {:>2}, {:>2} }}, ".format(idx_ijk_y_m1, idx_ijk_y_m2, idx_ijk_y_p1))
      f.write("{{ {:>2}, {:>2}, {:>2} }}  ".format(idx_ijk_z_m1, idx_ijk_z_m2, idx_ijk_z_p1))
      f.write("} ")

      f.write("},\n")

  f.write("};\n\n\n")

  # now the indexing map
  f.write("extern int const am_recur_map[{}] = {{\n".format(maxam+1))
  for v in starts:
    f.write("    {:>7},\n".format(v))
  f.write("};\n\n")


print()
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
print("!! REMEMBER TO CHANGE THE HEADER FILE!")
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
