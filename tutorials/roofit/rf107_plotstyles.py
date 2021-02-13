## \file
## \ingroup tutorial_roofit
## \notebook
## Basic functionality: demonstration of various plotting styles of data, functions in a RooPlot
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT


# Set up model
# ---------------------

# Create observables
x = ROOT.RooRealVar("x", "x", -10, 10)

# Create Gaussian
sigma = ROOT.RooRealVar("sigma", "sigma", 3, 0.1, 10)
mean = ROOT.RooRealVar("mean", "mean", -3, -10, 10)
gauss = ROOT.RooGaussian("gauss", "gauss", x, mean, sigma)

# Generate a sample of 100 events with sigma=3
data = gauss.generate(ROOT.RooArgSet(x), 100)

# Fit pdf to data
gauss.fitTo(data)

# Make plot frames
# -------------------------------

# Make four plot frames to demonstrate various plotting features
frame1 = x.frame(ROOT.RooFit.Name("xframe"), ROOT.RooFit.Title(
    "Red Curve / SumW2 Histo errors"), ROOT.RooFit.Bins(20))
frame2 = x.frame(ROOT.RooFit.Name("xframe"), ROOT.RooFit.Title(
    "Dashed Curve / No XError bars"), ROOT.RooFit.Bins(20))
frame3 = x.frame(ROOT.RooFit.Name("xframe"), ROOT.RooFit.Title(
    "Filled Curve / Blue Histo"), ROOT.RooFit.Bins(20))
frame4 = x.frame(ROOT.RooFit.Name("xframe"), ROOT.RooFit.Title(
    "Partial Range / Filled Bar chart"), ROOT.RooFit.Bins(20))

# Data plotting styles
# ---------------------------------------

# Use sqrt(sum(weights^2)) error instead of Poisson errors
data.plotOn(frame1, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

# Remove horizontal error bars
data.plotOn(frame2, ROOT.RooFit.XErrorSize(0))

# Blue markers and error bors
data.plotOn(frame3, ROOT.RooFit.MarkerColor(
    ROOT.kBlue), ROOT.RooFit.LineColor(ROOT.kBlue))

# Filled bar chart
data.plotOn(
    frame4,
    ROOT.RooFit.DrawOption("B"),
    ROOT.RooFit.DataError(
        ROOT.RooAbsData.ErrorType(2)),
    ROOT.RooFit.XErrorSize(0),
    ROOT.RooFit.FillColor(
        ROOT.kGray))

# Function plotting styles
# -----------------------------------------------

# Change line color to red
gauss.plotOn(frame1, ROOT.RooFit.LineColor(ROOT.kRed))

# Change line style to dashed
gauss.plotOn(frame2, ROOT.RooFit.LineStyle(ROOT.kDashed))

# Filled shapes in green color
gauss.plotOn(frame3, ROOT.RooFit.DrawOption("F"),
             ROOT.RooFit.FillColor(ROOT.kOrange), ROOT.RooFit.MoveToBack())

#
gauss.plotOn(frame4, ROOT.RooFit.Range(-8, 3),
             ROOT.RooFit.LineColor(ROOT.kMagenta))

c = ROOT.TCanvas("rf107_plotstyles", "rf107_plotstyles", 800, 800)
c.Divide(2, 2)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
frame1.GetYaxis().SetTitleOffset(1.6)
frame1.Draw()
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
frame2.GetYaxis().SetTitleOffset(1.6)
frame2.Draw()
c.cd(3)
ROOT.gPad.SetLeftMargin(0.15)
frame3.GetYaxis().SetTitleOffset(1.6)
frame3.Draw()
c.cd(4)
ROOT.gPad.SetLeftMargin(0.15)
frame4.GetYaxis().SetTitleOffset(1.6)
frame4.Draw()

c.SaveAs("rf107_plotstyles.png")
