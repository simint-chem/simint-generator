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
parser.add_argument("-b",         type=int, required=False, default=1,    help="Number of bins")
parser.add_argument("-p",         action="store_true",                    help="Make the plot")
parser.add_argument("-t",         type=int, required=False, default=200,  help="Number of test points (per bin)")
parser.add_argument("--dps",      type=int, required=False, default=256,  help="Decimal precision/sig figs to use/calculate")
parser.add_argument("minx",       type=str,                               help="Starting x value")
parser.add_argument("maxx",       type=str,                               help="Ending x value")
args = parser.parse_args()

# Set the dps option
mp.dps = args.dps

# Convert stuff to mpmath
minx = mp.mpf(args.minx)
maxx = mp.mpf(args.maxx)

# order of the polynomial
# We need to calculate args.o + 1 roots,
# etc, so we need one more
order = args.o + 1


test_fx = []    # holds value of boys function at test points
test_ix = []    # holds value of interpolated function at test points
test_x = []     # holds x value of test points
test_abserr = []   # Interpolation absolute error (test_ix - test_fx)
test_relerr = []   # Interpolation relative error (test_ix - test_fx)/test_fx


# split the interval into bins
binpoints = np.linspace(minx, maxx, args.b+1)
binpoints = list(zip(binpoints[:-1], binpoints[1:]))

# binpoints is a list of tuples (start, x)
for bounds in binpoints:
    print("==Interval [{} , {}]".format(bounds[0], bounds[1]))

    data = Interpolate(args.n, order, bounds, args.t, BoysValue) 
    test_x.extend(data["test_x"])
    test_fx.extend(data["test_fx"])
    test_ix.extend(data["test_ix"])
    test_abserr.extend(data["test_abserr"])
    test_relerr.extend(data["test_relerr"])
    print("Max absolute error in this interval: {}".format(mp.nstr(data["max_abserr"])))
    print("Max relative error in this interval: {}".format(mp.nstr(data["max_relerr"])))
    print()


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
