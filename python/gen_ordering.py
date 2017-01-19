#!/usr/bin/env python3

#######################################
# Generates shell ordering
#######################################

import argparse
import sys
from ordering_common import Gaussian, IteratePsi4, GenerateMap

parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True,               help="Output file name base (no extension)")
parser.add_argument("--max-am",    type=int, required=True,               help="Maximum angular momentum value to go to")
args = parser.parse_args()
maxam = args.max_am

gmap = GenerateMap(maxam, IteratePsi4)

print("------------------------------------")
print("Options for gen_ordering:")
print("   Max am: {}".format(maxam))
print("------------------------------------")


# Output to file
with open(args.filename, 'w') as f:
  f.write("/*\n")
  f.write(" Generated with:\n")
  f.write("   " + " ".join(sys.argv[:]))
  f.write("\n")
  f.write("------------------------------------\n")
  f.write("Options for gen_ordering:\n")
  f.write("   Max am: {}\n".format(maxam))
  f.write("------------------------------------\n")
  f.write("*/\n\n")

  f.write("#include <array>\n")
  f.write("#include <map>\n")
  f.write("#include <vector>\n")
  f.write("\n")
  f.write("typedef std::array<int, 3> ExpList;\n")
  f.write("\n")
  f.write("std::map<const int, const std::vector<ExpList>> gorder_map = {\n")

  for i in range(0, maxam+1):
    f.write("{{ {:2}, {{\n".format(i))

    for g in gmap[i]:
      f.write("        " + str(g) + ",\n")
    f.write("      }\n")
    f.write("},\n")
  f.write("};\n")

  f.write("\n\n")

  f.write("std::map<const int, const std::vector<int>> gindex_map = {\n")

  for i in range(0, maxam+1):
    f.write("{{ {:2}, {{\n".format(i))

    for g in gmap[i]:
      f.write("        /*" + str(g) + "*/ " + g.hashstr() + ",\n")
    f.write("      }\n")
    f.write("},\n")
  f.write("};\n")
