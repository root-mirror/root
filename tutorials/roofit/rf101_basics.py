## \file
## \ingroup tutorial_roofit
## \notebook
## This tutorial illustrates the basic features of RooFit.
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT

# Set up model
# ---------------------
# Declare variables x,mean,sigma with associated name, title, initial
# value and allowed range
x = ROOT.RooRealVar("x", "x", -10, 10)
mean = ROOT.RooRealVar("mean", "mean of gaussian", 1, -10, 10)
sigma = ROOT.RooRealVar("sigma", "width of gaussian", 1, 0.1, 10)

# Build gaussian pdf in terms of x,mean and sigma
gauss = ROOT.RooGaussian("gauss", "gaussian PDF", x, mean, sigma)

# Construct plot frame in 'x'
xframe = x.frame(ROOT.RooFit.Title("Gaussian pdf"))  # RooPlot

# Plot model and change parameter values
# ---------------------------------------------------------------------------
# Plot gauss in frame (i.e. in x)
gauss.plotOn(xframe)

# Change the value of sigma to 3
sigma.setVal(3)

# Plot gauss in frame (i.e. in x) and draw frame on canvas
gauss.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kRed))

# Generate events
# -----------------------------
# Generate a dataset of 1000 events in x from gauss
data = gauss.generate(ROOT.RooArgSet(x), 10000)  # ROOT.RooDataSet

# Make a second plot frame in x and draw both the
# data and the pdf in the frame
xframe2 = x.frame(ROOT.RooFit.Title(
    "Gaussian pdf with data"))  # RooPlot
data.plotOn(xframe2)
gauss.plotOn(xframe2)

# Fit model to data
# -----------------------------
# Fit pdf to data
gauss.fitTo(data)

# Print values of mean and sigma (that now reflect fitted values and
# errors)
mean.Print()
sigma.Print()

# Draw all frames on a canvas
c = ROOT.TCanvas("rf101_basics", "rf101_basics", 800, 400)
c.Divide(2)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
xframe.GetYaxis().SetTitleOffset(1.6)
xframe.Draw()
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
xframe2.GetYaxis().SetTitleOffset(1.6)
xframe2.Draw()

c.SaveAs("rf101_basics.png")
