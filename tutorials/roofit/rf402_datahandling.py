## \file
## \ingroup tutorial_roofit
## \notebook
## Data and categories: tools for manipulation of (un)binned datasets
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

from __future__ import print_function
import ROOT
import math

# WVE Add reduction by range

# Binned (RooDataHist) and unbinned datasets (RooDataSet) share
# many properties and inherit from a common abstract base class
# (RooAbsData), provides an interface for all operations
# that can be performed regardless of the data format

x = ROOT.RooRealVar("x", "x", -10, 10)
y = ROOT.RooRealVar("y", "y", 0, 40)
c = ROOT.RooCategory("c", "c")
c.defineType("Plus", +1)
c.defineType("Minus", -1)

# Basic operations on unbinned datasetss
# --------------------------------------------------------------

# ROOT.RooDataSet is an unbinned dataset (a collection of points in
# N-dimensional space)
d = ROOT.RooDataSet("d", "d", ROOT.RooArgSet(x, y, c))

# Unlike ROOT.RooAbsArgs (ROOT.RooAbsPdf, ROOT.RooFormulaVar,....) datasets are not attached to
# the variables they are constructed from. Instead they are attached to an internal
# clone of the supplied set of arguments

# Fill d with dummy values
for i in range(1000):
    x.setVal(i / 50 - 10)
    y.setVal(math.sqrt(1.0 * i))
    if (i % 2):
        c.setLabel("Plus")
    else:
        c.setLabel("Minus")

    # We must explicitly refer to x,y, here to pass the values because
    # d is not linked to them (as explained above)
    if i < 3:
        print(x, y, c)
        print(type(x))
    d.add(ROOT.RooArgSet(x, y, c))

d.Print("v")
print("")

# The get() function returns a pointer to the internal copy of the RooArgSet(x,y,c)
# supplied in the constructor
row = d.get()
row.Print("v")
print("")

# Get with an argument loads a specific data point in row and returns
# a pointer to row argset. get() always returns the same pointer, unless
# an invalid row number is specified. In that case a null ptr is returned
d.get(900).Print("v")
print("")

# Reducing, appending and merging
# -------------------------------------------------------------

# The reduce() function returns a dataset which is a subset of the
# original
print("\n >> d1 has only columns x,c")
d1 = d.reduce(ROOT.RooArgSet(x, c))
d1.Print("v")

print("\n >> d2 has only column y")
d2 = d.reduce(ROOT.RooArgSet(y))
d2.Print("v")

print("\n >> d3 has only the points with y>5.17")
d3 = d.reduce("y>5.17")
d3.Print("v")

print("\n >> d4 has only columns x, for data points with y>5.17")
d4 = d.reduce(ROOT.RooArgSet(x, c), "y>5.17")
d4.Print("v")

# The merge() function adds two data set column-wise
print("\n >> merge d2(y) with d1(x,c) to form d1(x,c,y)")
d1.merge(d2)
d1.Print("v")

# The append() function addes two datasets row-wise
print("\n >> append data points of d3 to d1")
d1.append(d3)
d1.Print("v")

# Operations on binned datasets
# ---------------------------------------------------------

# A binned dataset can be constructed empty, an unbinned dataset, or
# from a ROOT native histogram (TH1,2,3)

print(">> construct dh (binned) from d(unbinned) but only take the x and y dimensions, ")
print(">> the category 'c' will be projected in the filling process")

# The binning of real variables (like x,y) is done using their fit range
# 'get/setRange()' and number of specified fit bins 'get/setBins()'.
# Category dimensions of binned datasets get one bin per defined category
# state
x.setBins(10)
y.setBins(10)
dh = ROOT.RooDataHist("dh", "binned version of d", ROOT.RooArgSet(x, y), d)
dh.Print("v")

yframe = y.frame(ROOT.RooFit.Bins(10), ROOT.RooFit.Title(
    "Operations on binned datasets"))
dh.plotOn(yframe)  # plot projection of 2D binned data on y

# Examine the statistics of a binned dataset
print(">> number of bins in dh   : ", dh.numEntries())
print(">> sum of weights in dh   : ", dh.sum(ROOT.kFALSE))
# accounts for bin volume
print(">> integral over histogram: ", dh.sum(ROOT.kTRUE))

# Locate a bin from a set of coordinates and retrieve its properties
x.setVal(0.3)
y.setVal(20.5)
print(">> retrieving the properties of the bin enclosing coordinate (x,y) = (0.3,20.5) bin center:")
# load bin center coordinates in internal buffer
dh.get(ROOT.RooArgSet(x, y)).Print("v")
print(" weight = ", dh.weight())  # return weight of last loaded coordinates

# Reduce the 2-dimensional binned dataset to a 1-dimensional binned dataset
#
# All reduce() methods are interfaced in RooAbsData. All reduction techniques
# demonstrated on unbinned datasets can be applied to binned datasets as
# well.
print(">> Creating 1-dimensional projection on y of dh for bins with x>0")
dh2 = dh.reduce(ROOT.RooArgSet(y), "x>0")
dh2.Print("v")

# Add dh2 to yframe and redraw
dh2.plotOn(yframe, ROOT.RooFit.LineColor(ROOT.kRed),
           ROOT.RooFit.MarkerColor(ROOT.kRed))

# Saving and loading from file
# -------------------------------------------------------

# Datasets can be persisted with ROOT I/O
print("\n >> Persisting d via ROOT I/O")
f = ROOT.TFile("rf402_datahandling.root", "RECREATE")
d.Write()
f.ls()

# To read back in future session:
# > ROOT.TFile f("rf402_datahandling.root")
# > d = (ROOT.RooDataSet*) f.FindObject("d")

c = ROOT.TCanvas("rf402_datahandling", "rf402_datahandling", 600, 600)
ROOT.gPad.SetLeftMargin(0.15)
yframe.GetYaxis().SetTitleOffset(1.4)
yframe.Draw()

c.SaveAs("rf402_datahandling.png")
