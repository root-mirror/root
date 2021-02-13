## \file
## \ingroup tutorial_roofit
## \notebook
## Multidimensional models: working with parameterized ranges in a fit.
## This an example of a fit with an acceptance that changes per-event
##
## `pdf = exp(-t/tau)` with `t[tmin,5]`
##
## where `t` and `tmin` are both observables in the dataset
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT


# Define observables and decay pdf
# ---------------------------------------------------------------

# Declare observables
t = ROOT.RooRealVar("t", "t", 0, 5)
tmin = ROOT.RooRealVar("tmin", "tmin", 0, 0, 5)

# Make parameterized range in t : [tmin,5]
t.setRange(tmin, ROOT.RooFit.RooConst(t.getMax()))

# Make pdf
tau = ROOT.RooRealVar("tau", "tau", -1.54, -10, -0.1)
model = ROOT.RooExponential("model", "model", t, tau)

# Create input data
# ------------------------------------

# Generate complete dataset without acceptance cuts (for reference)
dall = model.generate(ROOT.RooArgSet(t), 10000)

# Generate a (fake) prototype dataset for acceptance limit values
tmp = ROOT.RooGaussian("gmin", "gmin", tmin, ROOT.RooFit.RooConst(
    0), ROOT.RooFit.RooConst(0.5)).generate(ROOT.RooArgSet(tmin), 5000)

# Generate dataset with t values that observe (t>tmin)
dacc = model.generate(ROOT.RooArgSet(t), ROOT.RooFit.ProtoData(tmp))

# Fit pdf to data in acceptance region
# -----------------------------------------------------------------------

r = model.fitTo(dacc, ROOT.RooFit.Save())

# Plot fitted pdf on full and accepted data
# ---------------------------------------------------------------------------------

# Make plot frame, datasets and overlay model
frame = t.frame(ROOT.RooFit.Title("Fit to data with per-event acceptance"))
dall.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kRed),
            ROOT.RooFit.LineColor(ROOT.kRed))
model.plotOn(frame)
dacc.plotOn(frame)

# Print fit results to demonstrate absence of bias
r.Print("v")

c = ROOT.TCanvas("rf314_paramranges", "rf314_paramranges", 600, 600)
ROOT.gPad.SetLeftMargin(0.15)
frame.GetYaxis().SetTitleOffset(1.6)
frame.Draw()

c.SaveAs("rf314_paramranges.png")
