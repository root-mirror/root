## \file
## \ingroup tutorial_roofit
## \notebook
## Multidimensional models: utility functions classes available for use in tailoring of
## composite (multidimensional) pdfs
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT

# Create observables, parameters
# -----------------------------------------------------------

# Create observables
x = ROOT.RooRealVar("x", "x", -5, 5)
y = ROOT.RooRealVar("y", "y", -5, 5)

# Create parameters
a0 = ROOT.RooRealVar("a0", "a0", -1.5, -5, 5)
a1 = ROOT.RooRealVar("a1", "a1", -0.5, -1, 1)
sigma = ROOT.RooRealVar("sigma", "width of gaussian", 0.5)

# Using RooFormulaVar to tailor pdf
# -----------------------------------------------------------------------

# Create interpreted function f(y) = a0 - a1*sqrt(10*abs(y))
fy_1 = ROOT.RooFormulaVar(
    "fy_1", "a0-a1*sqrt(10*abs(y))", ROOT.RooArgList(y, a0, a1))

# Create gauss(x,f(y),s)
model_1 = ROOT.RooGaussian(
    "model_1", "Gaussian with shifting mean", x, fy_1, sigma)

# Using RooPolyVar to tailor pdf
# -----------------------------------------------------------------------

# Create polynomial function f(y) = a0 + a1*y
fy_2 = ROOT.RooPolyVar("fy_2", "fy_2", y, ROOT.RooArgList(a0, a1))

# Create gauss(x,f(y),s)
model_2 = ROOT.RooGaussian(
    "model_2", "Gaussian with shifting mean", x, fy_2, sigma)

# Using RooAddition to tailor pdf
# -----------------------------------------------------------------------

# Create sum function f(y) = a0 + y
fy_3 = ROOT.RooAddition("fy_3", "a0+y", ROOT.RooArgList(a0, y))

# Create gauss(x,f(y),s)
model_3 = ROOT.RooGaussian(
    "model_3", "Gaussian with shifting mean", x, fy_3, sigma)

# Using RooProduct to tailor pdf
# -----------------------------------------------------------------------

# Create product function f(y) = a1*y
fy_4 = ROOT.RooProduct("fy_4", "a1*y", ROOT.RooArgList(a1, y))

# Create gauss(x,f(y),s)
model_4 = ROOT.RooGaussian(
    "model_4", "Gaussian with shifting mean", x, fy_4, sigma)

# Plot all pdfs
# ----------------------------

# Make two-dimensional plots in x vs y
hh_model_1 = model_1.createHistogram("hh_model_1", x, ROOT.RooFit.Binning(
    50), ROOT.RooFit.YVar(y, ROOT.RooFit.Binning(50)))
hh_model_2 = model_2.createHistogram("hh_model_2", x, ROOT.RooFit.Binning(
    50), ROOT.RooFit.YVar(y, ROOT.RooFit.Binning(50)))
hh_model_3 = model_3.createHistogram("hh_model_3", x, ROOT.RooFit.Binning(
    50), ROOT.RooFit.YVar(y, ROOT.RooFit.Binning(50)))
hh_model_4 = model_4.createHistogram("hh_model_4", x, ROOT.RooFit.Binning(
    50), ROOT.RooFit.YVar(y, ROOT.RooFit.Binning(50)))
hh_model_1.SetLineColor(ROOT.kBlue)
hh_model_2.SetLineColor(ROOT.kBlue)
hh_model_3.SetLineColor(ROOT.kBlue)
hh_model_4.SetLineColor(ROOT.kBlue)

# Make canvas and draw ROOT.RooPlots
c = ROOT.TCanvas("rf302_utilfuncs", "rf302_utilfuncs", 800, 800)
c.Divide(2, 2)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.20)
hh_model_1.GetZaxis().SetTitleOffset(2.5)
hh_model_1.Draw("surf")
c.cd(2)
ROOT.gPad.SetLeftMargin(0.20)
hh_model_2.GetZaxis().SetTitleOffset(2.5)
hh_model_2.Draw("surf")
c.cd(3)
ROOT.gPad.SetLeftMargin(0.20)
hh_model_3.GetZaxis().SetTitleOffset(2.5)
hh_model_3.Draw("surf")
c.cd(4)
ROOT.gPad.SetLeftMargin(0.20)
hh_model_4.GetZaxis().SetTitleOffset(2.5)
hh_model_4.Draw("surf")

c.SaveAs("rf302_utilfuncs.png")
