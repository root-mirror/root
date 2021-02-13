## \file
## \ingroup tutorial_roofit
## \notebook -nodraw
## Data and categories: latex printing of lists and sets of RooArgSets
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT


# Setup composite pdf
# --------------------------------------

# Declare observable x
x = ROOT.RooRealVar("x", "x", 0, 10)

# Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and
# their parameters
mean = ROOT.RooRealVar("mean", "mean of gaussians", 5)
sigma1 = ROOT.RooRealVar("sigma1", "width of gaussians", 0.5)
sigma2 = ROOT.RooRealVar("sigma2", "width of gaussians", 1)
sig1 = ROOT.RooGaussian("sig1", "Signal component 1", x, mean, sigma1)
sig2 = ROOT.RooGaussian("sig2", "Signal component 2", x, mean, sigma2)

# Sum the signal components into a composite signal pdf
sig1frac = ROOT.RooRealVar(
    "sig1frac", "fraction of component 1 in signal", 0.8, 0., 1.)
sig = ROOT.RooAddPdf(
    "sig", "Signal", ROOT.RooArgList(sig1, sig2), ROOT.RooArgList(sig1frac))

# Build Chebychev polynomial pdf
a0 = ROOT.RooRealVar("a0", "a0", 0.5, 0., 1.)
a1 = ROOT.RooRealVar("a1", "a1", -0.2, 0., 1.)
bkg1 = ROOT.RooChebychev("bkg1", "Background 1",
                         x, ROOT.RooArgList(a0, a1))

# Build expontential pdf
alpha = ROOT.RooRealVar("alpha", "alpha", -1)
bkg2 = ROOT.RooExponential("bkg2", "Background 2", x, alpha)

# Sum the background components into a composite background pdf
bkg1frac = ROOT.RooRealVar(
    "sig1frac", "fraction of component 1 in background", 0.2, 0., 1.)
bkg = ROOT.RooAddPdf(
    "bkg", "Signal", ROOT.RooArgList(bkg1, bkg2), ROOT.RooArgList(sig1frac))

# Sum the composite signal and background
bkgfrac = ROOT.RooRealVar("bkgfrac", "fraction of background", 0.5, 0., 1.)
model = ROOT.RooAddPdf(
    "model", "g1+g2+a", ROOT.RooArgList(bkg, sig), ROOT.RooArgList(bkgfrac))

# Make list of parameters before and after fit
# ----------------------------------------------------------------------------------------

# Make list of model parameters
params = model.getParameters(ROOT.RooArgSet(x))

# Save snapshot of prefit parameters
initParams = params.snapshot()

# Do fit to data, obtain error estimates on parameters
data = model.generate(ROOT.RooArgSet(x), 1000)
model.fitTo(data)

# Print LateX table of parameters of pdf
# --------------------------------------------------------------------------

# Print parameter list in LaTeX for (one column with names, column with
# values)
params.printLatex()

# Print parameter list in LaTeX for (names values|names values)
params.printLatex(ROOT.RooFit.Columns(2))

# Print two parameter lists side by side (name values initvalues)
params.printLatex(ROOT.RooFit.Sibling(initParams))

# Print two parameter lists side by side (name values initvalues|name
# values initvalues)
params.printLatex(ROOT.RooFit.Sibling(initParams), ROOT.RooFit.Columns(2))

# Write LaTex table to file
params.printLatex(ROOT.RooFit.Sibling(initParams),
                  ROOT.RooFit.OutputFile("rf407_latextables.tex"))
