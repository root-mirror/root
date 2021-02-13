## \file
## \ingroup tutorial_dataframe
## \notebook -nodraw
## Use just-in-time-compiled Filters and Defines for quick prototyping.
##
## This tutorial illustrates how to use jit-compiling features of RDataFrame
## to define data using C++ code in a Python script
##
## \macro_code
## \macro_output
##
## \date October 2017
## \author Guilherme Amadio (CERN)

import ROOT

## We will inefficiently calculate an approximation of pi by generating
## some data and doing very simple filtering and analysis on it.

## We start by creating an empty dataframe where we will insert 10 million
## random points in a square of side 2.0 (that is, with an inscribed unit
## circle).

npoints = 10000000
df = ROOT.RDataFrame(npoints)

## Define what data we want inside the dataframe. We do not need to define p
## as an array, but we do it here to demonstrate how to use jitting with RDataFrame

pidf = df.Define("x", "gRandom->Uniform(-1.0, 1.0)") \
         .Define("y", "gRandom->Uniform(-1.0, 1.0)") \
         .Define("p", "std::array<double, 2> v{x, y}; return v;") \
         .Define("r", "double r2 = 0.0; for (auto&& w : p) r2 += w*w; return sqrt(r2);")

## Now we have a dataframe with columns x, y, p (which is a point based on x
## and y), and the radius r = sqrt(x*x + y*y). In order to approximate pi, we
## need to know how many of our data points fall inside the circle of radius
## one compared with the total number of points. The ratio of the areas is
##
##     A_circle / A_square = pi r*r / l * l, where r = 1.0, and l = 2.0
##
## Therefore, we can approximate pi with 4 times the number of points inside
## the unit circle over the total number of points:

incircle = pidf.Filter("r <= 1.0").Count().GetValue()

pi_approx = 4.0 * incircle / npoints

print("pi is approximately equal to %g" % (pi_approx))
