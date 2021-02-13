## \file
## \ingroup tutorial_roofit
## \notebook
## Basic functionality: examples on normalization and integration of pdfs, construction
## of cumulative distribution functions from monodimensional pdfs
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

from __future__ import print_function
import ROOT

# Set up model
# ---------------------

# Create observables x,y
x = ROOT.RooRealVar("x", "x", -10, 10)

# Create pdf gaussx(x,-2,3)
gx = ROOT.RooGaussian(
    "gx", "gx", x, ROOT.RooFit.RooConst(-2), ROOT.RooFit.RooConst(3))

# Retrieve raw & normalized values of RooFit pdfs
# --------------------------------------------------------------------------------------------------

# Return 'raw' unnormalized value of gx
print("gx = ", gx.getVal())

# Return value of gx normalized over x in range [-10,10]
nset = ROOT.RooArgSet(x)
print("gx_Norm[x] = ", gx.getVal(nset))

# Create object representing integral over gx
# which is used to calculate  gx_Norm[x] == gx / gx_Int[x]
igx = gx.createIntegral(ROOT.RooArgSet(x))
print("gx_Int[x] = ", igx.getVal())

# Integrate normalized pdf over subrange
# ----------------------------------------------------------------------------

# Define a range named "signal" in x from -5,5
x.setRange("signal", -5, 5)

# Create an integral of gx_Norm[x] over x in range "signal"
# ROOT.This is the fraction of of pdf gx_Norm[x] which is in the
# range named "signal"
xset = ROOT.RooArgSet(x)
igx_sig = gx.createIntegral(xset, ROOT.RooFit.NormSet(xset), ROOT.RooFit.Range("signal"))
print("gx_Int[x|signal]_Norm[x] = ", igx_sig.getVal())

# Construct cumulative distribution function from pdf
# -----------------------------------------------------------------------------------------------------

# Create the cumulative distribution function of gx
# i.e. calculate Int[-10,x] gx(x') dx'
gx_cdf = gx.createCdf(ROOT.RooArgSet(x))

# Plot cdf of gx versus x
frame = x.frame(ROOT.RooFit.Title("cdf of Gaussian pdf"))
gx_cdf.plotOn(frame)

# Draw plot on canvas
c = ROOT.TCanvas("rf110_normintegration",
                 "rf110_normintegration", 600, 600)
ROOT.gPad.SetLeftMargin(0.15)
frame.GetYaxis().SetTitleOffset(1.6)
frame.Draw()

c.SaveAs("rf110_normintegration.png")
