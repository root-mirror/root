## \file
## \ingroup tutorial_roofit
## \notebook
## Multidimensional models: simple uncorrelated multi-dimensional pdfs
##
## `pdf = gauss(x,mx,sx) * gauss(y,my,sy)`
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT


# Create component pdfs in x and y
# ----------------------------------------------------------------

# Create two pdfs gaussx(x,meanx,sigmax) gaussy(y,meany,sigmay) and its
# variables
x = ROOT.RooRealVar("x", "x", -5, 5)
y = ROOT.RooRealVar("y", "y", -5, 5)

meanx = ROOT.RooRealVar("mean1", "mean of gaussian x", 2)
meany = ROOT.RooRealVar("mean2", "mean of gaussian y", -2)
sigmax = ROOT.RooRealVar("sigmax", "width of gaussian x", 1)
sigmay = ROOT.RooRealVar("sigmay", "width of gaussian y", 5)

gaussx = ROOT.RooGaussian("gaussx", "gaussian PDF", x, meanx, sigmax)
gaussy = ROOT.RooGaussian("gaussy", "gaussian PDF", y, meany, sigmay)

# Construct uncorrelated product pdf
# -------------------------------------------------------------------

# Multiply gaussx and gaussy into a two-dimensional pdf gaussxy
gaussxy = ROOT.RooProdPdf(
    "gaussxy", "gaussx*gaussy", ROOT.RooArgList(gaussx, gaussy))

# Sample pdf, plot projection on x and y
# ---------------------------------------------------------------------------

# Generate 10000 events in x and y from gaussxy
data = gaussxy.generate(ROOT.RooArgSet(x, y), 10000)

# Plot x distribution of data and projection of gaussxy x = Int(dy)
# gaussxy(x,y)
xframe = x.frame(ROOT.RooFit.Title("X projection of gauss(x)*gauss(y)"))
data.plotOn(xframe)
gaussxy.plotOn(xframe)

# Plot x distribution of data and projection of gaussxy y = Int(dx)
# gaussxy(x,y)
yframe = y.frame(ROOT.RooFit.Title("Y projection of gauss(x)*gauss(y)"))
data.plotOn(yframe)
gaussxy.plotOn(yframe)

# Make canvas and draw ROOT.RooPlots
c = ROOT.TCanvas("rf304_uncorrprod", "rf304_uncorrprod", 800, 400)
c.Divide(2)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
xframe.GetYaxis().SetTitleOffset(1.4)
xframe.Draw()
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
yframe.GetYaxis().SetTitleOffset(1.4)
yframe.Draw()

c.SaveAs("rf304_uncorrprod.png")
