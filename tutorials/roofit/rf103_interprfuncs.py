## \file
## \ingroup tutorial_roofit
## \notebook
## Basic functionality: interpreted functions and pdfs
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT

# Generic interpreted pdf
# ------------------------------

# Declare observable x
x = ROOT.RooRealVar("x", "x", -20, 20)

# Construct generic pdf from interpreted expression
# ------------------------------------------------------

# ROOT.To construct a proper pdf, the formula expression is explicitly normalized internally by dividing
# it by a numeric integral of the expresssion over x in the range [-20,20]
#
alpha = ROOT.RooRealVar("alpha", "alpha", 5, 0.1, 10)
genpdf = ROOT.RooGenericPdf(
    "genpdf",
    "genpdf",
    "(1+0.1*abs(x)+sin(sqrt(abs(x*alpha+0.1))))",
    ROOT.RooArgList(
        x,
        alpha))

# Sample, fit and plot generic pdf
# ---------------------------------------------------------------

# Generate a toy dataset from the interpreted pdf
data = genpdf.generate(ROOT.RooArgSet(x), 10000)

# Fit the interpreted pdf to the generated data
genpdf.fitTo(data)

# Make a plot of the data and the pdf overlaid
xframe = x.frame(ROOT.RooFit.Title("Interpreted expression pdf"))
data.plotOn(xframe)
genpdf.plotOn(xframe)

# Standard pdf adjust with interpreted helper function
# ------------------------------------------------------------------------------------------------------------
# Make a gauss(x,sqrt(mean2),sigma) from a standard ROOT.RooGaussian                                               #
#
# Construct standard pdf with formula replacing parameter
# ------------------------------------------------------------------------------------------------------------

# Construct parameter mean2 and sigma
mean2 = ROOT.RooRealVar("mean2", "mean^2", 10, 0, 200)
sigma = ROOT.RooRealVar("sigma", "sigma", 3, 0.1, 10)

# Construct interpreted function mean = sqrt(mean^2)
mean = ROOT.RooFormulaVar(
    "mean", "mean", "sqrt(mean2)", ROOT.RooArgList(mean2))

# Construct a gaussian g2(x,sqrt(mean2),sigma)
g2 = ROOT.RooGaussian("g2", "h2", x, mean, sigma)

# Generate toy data
# ---------------------------------

# Construct a separate gaussian g1(x,10,3) to generate a toy Gaussian
# dataset with mean 10 and width 3
g1 = ROOT.RooGaussian("g1", "g1", x, ROOT.RooFit.RooConst(
    10), ROOT.RooFit.RooConst(3))
data2 = g1.generate(ROOT.RooArgSet(x), 1000)

# Fit and plot tailored standard pdf
# -------------------------------------------------------------------

# Fit g2 to data from g1
r = g2.fitTo(data2, ROOT.RooFit.Save())  # ROOT.RooFitResult
r.Print()

# Plot data on frame and overlay projection of g2
xframe2 = x.frame(ROOT.RooFit.Title("Tailored Gaussian pdf"))
data2.plotOn(xframe2)
g2.plotOn(xframe2)

# Draw all frames on a canvas
c = ROOT.TCanvas("rf103_interprfuncs", "rf103_interprfuncs", 800, 400)
c.Divide(2)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
xframe.GetYaxis().SetTitleOffset(1.4)
xframe.Draw()
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
xframe2.GetYaxis().SetTitleOffset(1.4)
xframe2.Draw()

c.SaveAs("rf103_interprfuncs.png")
