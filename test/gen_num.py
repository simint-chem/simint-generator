#!/usr/bin/env python3

import argparse
import os
import random
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-r", action='store_true', help="Use random numbers")
parser.add_argument("-a", type=float, default=50.0, help="Maximum alpha value")
parser.add_argument("-c", type=float, default=2.0, help="Maximum coefficient value")
parser.add_argument("-x", type=float, default=0.5, help="Maximum xyz magnitude")
parser.add_argument("-o", type=str, required=True, help="Output file")
parser.add_argument("-ns", type=int, default=10, help="Number of shells")
parser.add_argument("-np", type=int, default=10, help="Number of primitives per shell")
args = parser.parse_args()

print("-------------------------------")
print("  Max alpha: {}".format(args.a))
print("  Max coeff: {}".format(args.c))
print("  Max coord: {}".format(args.x))
print("# of Shells: {}".format(args.ns))
print(" # of Prims: {}".format(args.np))
print("     Output: {}".format(args.o))
print("-------------------------------")

# Generate the data
ndat = args.ns * args.np
if args.r:
    alpha_dat = list(np.random.uniform(0, args.a, ndat))
    coeff_dat = list(np.random.uniform(0, args.c, ndat))
    x_dat = list(np.random.uniform(-args.x, args.x, args.ns))
    y_dat = list(np.random.uniform(-args.x, args.x, args.ns))
    z_dat = list(np.random.uniform(-args.x, args.x, args.ns))
else:
    alpha_dat = list(np.linspace(0, args.a, ndat))
    coeff_dat = list(np.linspace(0, args.c, ndat))
    x_dat = list(np.linspace(-args.x, args.x, args.ns))
    y_dat = list(np.linspace(-args.x, args.x, args.ns))
    z_dat = list(np.linspace(-args.x, args.x, args.ns))
    random.shuffle(alpha_dat)
    random.shuffle(coeff_dat)
    random.shuffle(x_dat)
    random.shuffle(y_dat)
    random.shuffle(z_dat)

# Write the file!
with open(args.o, 'w') as f:
    f.write("{} {}\n".format(args.ns, args.np))
    for i in range(0, args.ns):
        f.write("{} {} {}\n".format(x_dat[i], y_dat[i], z_dat[i]))
        
        for i in range(0, args.np):
            f.write("{} {}\n".format(alpha_dat[i], coeff_dat[i]))

