## \file
## \ingroup tutorial_roofit
## \notebook
## Basic functionality: adding boxes with parameters to RooPlots and decorating with arrows, etc...
##
## \macro_code
##
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT

# Set up model
# ---------------------

# Create observables
x = ROOT.RooRealVar("x", "x", -10, 10)

# Create Gaussian
sigma = ROOT.RooRealVar("sigma", "sigma", 1, 0.1, 10)
mean = ROOT.RooRealVar("mean", "mean", -3, -10, 10)
gauss = ROOT.RooGaussian("gauss", "gauss", x, mean, sigma)

# Generate a sample of 1000 events with sigma=3
data = gauss.generate(ROOT.RooArgSet(x), 1000)

# Fit pdf to data
gauss.fitTo(data)

# Plot pdf and data
# -------------------------------------

# Overlay projection of gauss on data
frame = x.frame(ROOT.RooFit.Name("xframe"), ROOT.RooFit.Title(
    "RooPlot with decorations"), ROOT.RooFit.Bins(40))
data.plotOn(frame)
gauss.plotOn(frame)

# Add box with pdf parameters
# -----------------------------------------------------

# Left edge of box starts at 55% of Xaxis)
gauss.paramOn(frame, ROOT.RooFit.Layout(0.55))

# Add box with data statistics
# -------------------------------------------------------

# X size of box is from 55% to 99% of Xaxis range, of box is at 80% of
# Yaxis range)
data.statOn(frame, ROOT.RooFit.Layout(0.55, 0.99, 0.8))

# Add text and arrow
# -----------------------------------

# Add text to frame
txt = ROOT.TText(2, 100, "Signal")
txt.SetTextSize(0.04)
txt.SetTextColor(ROOT.kRed)
frame.addObject(txt)

# Add arrow to frame
arrow = ROOT.TArrow(2, 100, -1, 50, 0.01, "|>")
arrow.SetLineColor(ROOT.kRed)
arrow.SetFillColor(ROOT.kRed)
arrow.SetLineWidth(3)
frame.addObject(arrow)

# Persist frame with all decorations in ROOT file
# ---------------------------------------------------------------------------------------------

f = ROOT.TFile("rf106_plotdecoration.root", "RECREATE")
frame.Write()
f.Close()

# To read back and plot frame with all decorations in clean root session do
# root> ROOT.TFile f("rf106_plotdecoration.root")
# root>  xframe.Draw()

c = ROOT.TCanvas("rf106_plotdecoration", "rf106_plotdecoration", 600, 600)
ROOT.gPad.SetLeftMargin(0.15)
frame.GetYaxis().SetTitleOffset(1.6)
frame.Draw()

c.SaveAs("rf106_plotdecoration.png")
