## \file
## \ingroup tutorial_roofit
## \notebook
## Multidimensional models: projecting pdf and data ranges in continuous observables
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT

# Create 3D pdf and data
# -------------------------------------------

# Create observables
x = ROOT.RooRealVar("x", "x", -5, 5)
y = ROOT.RooRealVar("y", "y", -5, 5)
z = ROOT.RooRealVar("z", "z", -5, 5)

# Create signal pdf gauss(x)*gauss(y)*gauss(z)
gx = ROOT.RooGaussian(
    "gx", "gx", x, ROOT.RooFit.RooConst(0), ROOT.RooFit.RooConst(1))
gy = ROOT.RooGaussian(
    "gy", "gy", y, ROOT.RooFit.RooConst(0), ROOT.RooFit.RooConst(1))
gz = ROOT.RooGaussian(
    "gz", "gz", z, ROOT.RooFit.RooConst(0), ROOT.RooFit.RooConst(1))
sig = ROOT.RooProdPdf("sig", "sig", ROOT.RooArgList(gx, gy, gz))

# Create background pdf poly(x)*poly(y)*poly(z)
px = ROOT.RooPolynomial("px", "px", x, ROOT.RooArgList(
    ROOT.RooFit.RooConst(-0.1), ROOT.RooFit.RooConst(0.004)))
py = ROOT.RooPolynomial("py", "py", y, ROOT.RooArgList(
    ROOT.RooFit.RooConst(0.1), ROOT.RooFit.RooConst(-0.004)))
pz = ROOT.RooPolynomial("pz", "pz", z)
bkg = ROOT.RooProdPdf("bkg", "bkg", ROOT.RooArgList(px, py, pz))

# Create composite pdf sig+bkg
fsig = ROOT.RooRealVar("fsig", "signal fraction", 0.1, 0., 1.)
model = ROOT.RooAddPdf(
    "model", "model", ROOT.RooArgList(
        sig, bkg), ROOT.RooArgList(fsig))

data = model.generate(ROOT.RooArgSet(x, y, z), 20000)

# Project pdf and data on x
# -------------------------------------------------

# Make plain projection of data and pdf on x observable
frame = x.frame(ROOT.RooFit.Title(
    "Projection of 3D data and pdf on X"), ROOT.RooFit.Bins(40))
data.plotOn(frame)
model.plotOn(frame)

# Project pdf and data on x in signal range
# ----------------------------------------------------------------------------------

# Define signal region in y and z observables
y.setRange("sigRegion", -1, 1)
z.setRange("sigRegion", -1, 1)

# Make plot frame
frame2 = x.frame(ROOT.RooFit.Title(
    "Same projection on X in signal range of (Y,Z)"), ROOT.RooFit.Bins(40))

# Plot subset of data in which all observables are inside "sigRegion"
# For observables that do not have an explicit "sigRegion" range defined (e.g. observable)
# an implicit definition is used that is identical to the full range (i.e.
# [-5,5] for x)
data.plotOn(frame2, ROOT.RooFit.CutRange("sigRegion"))

# Project model on x, projected observables (y,z) only in "sigRegion"
model.plotOn(frame2, ROOT.RooFit.ProjectionRange("sigRegion"))

c = ROOT.TCanvas("rf311_rangeplot", "rf310_rangeplot", 800, 400)
c.Divide(2)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
frame.GetYaxis().SetTitleOffset(1.4)
frame.Draw()
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
frame2.GetYaxis().SetTitleOffset(1.4)
frame2.Draw()

c.SaveAs("rf311_rangeplot.png")
