## \file
## \ingroup tutorial_roofit
## \notebook
## Addition and convolution: fitting and plotting in sub ranges
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

from __future__ import print_function
import ROOT

# Set up model
# ---------------------

# Construct observables x
x = ROOT.RooRealVar("x", "x", -10, 10)

# Construct gaussx(x,mx,1)
mx = ROOT.RooRealVar("mx", "mx", 0, -10, 10)
gx = ROOT.RooGaussian("gx", "gx", x, mx, ROOT.RooFit.RooConst(1))

# px = 1 (flat in x)
px = ROOT.RooPolynomial("px", "px", x)

# model = f*gx + (1-f)px
f = ROOT.RooRealVar("f", "f", 0., 1.)
model = ROOT.RooAddPdf(
    "model", "model", ROOT.RooArgList(gx, px), ROOT.RooArgList(f))

# Generated 10000 events in (x,y) from pdf model
modelData = model.generate(ROOT.RooArgSet(x), 10000)

# Fit full range
# ---------------------------

# Fit pdf to all data
r_full = model.fitTo(modelData, ROOT.RooFit.Save(ROOT.kTRUE))

# Fit partial range
# ----------------------------------

# Define "signal" range in x as [-3,3]
x.setRange("signal", -3, 3)

# Fit pdf only to data in "signal" range
r_sig = model.fitTo(modelData, ROOT.RooFit.Save(
    ROOT.kTRUE), ROOT.RooFit.Range("signal"))

# Plot/print results
# ---------------------------------------

# Make plot frame in x and add data and fitted model
frame = x.frame(ROOT.RooFit.Title("Fitting a sub range"))
modelData.plotOn(frame)
model.plotOn(
    frame, ROOT.RooFit.Range("Full"), ROOT.RooFit.LineStyle(
        ROOT.kDashed), ROOT.RooFit.LineColor(
            ROOT.kRed)) # Add shape in full ranged dashed
model.plotOn(frame)  # By default only fitted range is shown

# Print fit results
print("result of fit on all data ")
r_full.Print()
print("result of fit in in signal region (note increased error on signal fraction)")
r_sig.Print()

# Draw frame on canvas
c = ROOT.TCanvas("rf203_ranges", "rf203_ranges", 600, 600)
ROOT.gPad.SetLeftMargin(0.15)
frame.GetYaxis().SetTitleOffset(1.4)
frame.Draw()

c.SaveAs("rf203_ranges.png")
