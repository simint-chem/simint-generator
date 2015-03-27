#!/usr/bin/env python3

###############################
# Plots Fn(x), the interpolated
# function, and the error
###############################

import argparse
from boys_common.boys import *
from boys_common.boys_cheby import *
import mpmath as mp
import matplotlib.pyplot as plt


# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("-n",         type=int, required=True,                help="Boys N value")
parser.add_argument("-o",         type=int, required=True,                help="Order of polynomial")
parser.add_argument("-s",         type=str, required=False,               help="Starting order list")
parser.add_argument("-b",         type=int, required=False, default=1,    help="Number of bins")
parser.add_argument("-p",         action="store_true",                    help="Make the plot")
parser.add_argument("-t",         type=int, required=False, default=200,  help="Number of test points (per bin)")
parser.add_argument("-f",         type=str, required=False, default="1",  help="Base factor (ie 1 = base 2, etc)")
parser.add_argument("--dps",      type=int, required=False, default=256,  help="Decimal precision/sig figs to use/calculate")
args = parser.parse_args()

# Set the dps option
mp.dps = args.dps


test_fx = []    # holds value of boys function at test points
test_ix = []    # holds value of interpolated function at test points
test_x = []     # holds x value of test points
test_abserr = []   # Interpolation absolute error (test_ix - test_fx)
test_relerr = []   # Interpolation relative error (test_ix - test_fx)/test_fx

max_abserr = []
max_relerr = []

# split the interval into bins
fact = mp.mpf(args.f)
binpoints = [mp.mpf(0.0), mp.mpf(1.0)]

for i in range(0, args.b):
  binpoints.append(binpoints[-1]+fact*binpoints[-1]+1)
binpoints = list(zip(binpoints[:-1], binpoints[1:]))

olst = [args.o] * (args.b+1)
if args.s:
  olst2 =  [int(x) for x in args.s.split(",")]
  for i,v in enumerate(olst2):
    olst[i] = v

binpoints = list(zip(binpoints, olst))

# binpoints is a list of tuples (start, x)
for bounds,order in binpoints:
    print("==Interval (order {}): [{} , {}]".format(order, bounds[0], bounds[1]))

    # Re: order of the polynomial
    # We need to calculate args.o + 1 roots,
    # etc, so we need one more
    data = Interpolate(args.n, order+1, bounds, args.t, BoysValue) 
    test_x.extend(data["test_x"])
    test_fx.extend(data["test_fx"])
    test_ix.extend(data["test_ix"])
    test_abserr.extend(data["test_abserr"])
    test_relerr.extend(data["test_relerr"])
    max_abserr.append(data["max_abserr"])
    max_relerr.append(data["max_relerr"])

print()
print("------------------------------------------------------------------------------------------------------------------------")
print("Error analysis")
print("------------------------------------------------------------------------------------------------------------------------")
print("{:25} {:10} {:25} {:25}".format("Range", "Order", "Max Abserr", "Max Relerr"))
for i,v in enumerate(binpoints):
  print("{:25} {:<10} {:25} {:25}".format(mp.nstr(v[0], 4), v[1], mp.nstr(max_abserr[i], 10), mp.nstr(max_relerr[i], 10)))
print("------------------------------------------------------------------------------------------------------------------------")


if args.p:
    
    f, axarr = plt.subplots(3, sharex=True)
    axarr[0].set_title("Boys Function with Interpolation: n = {}".format(args.n))
    axarr[0].get_yaxis().set_major_formatter(plt.FormatStrFormatter("%4.2e"))
    axarr[0].plot(test_x, test_fx, 'g.', test_x, test_ix, 'b-') 

    axarr[1].set_title("Absolute error")
    axarr[1].axhline(0, color='gray')
    axarr[1].get_yaxis().set_major_formatter(plt.FormatStrFormatter("%4.2e"))
    axarr[1].plot(test_x, test_abserr, 'r-')

    axarr[2].set_title("Relative error")
    axarr[2].axhline(0, color='gray')
    axarr[2].get_yaxis().set_major_formatter(plt.FormatStrFormatter("%4.2e"))
    axarr[2].plot(test_x, test_relerr, 'k-')

    f.subplots_adjust(hspace=0.4)

    plt.show()
