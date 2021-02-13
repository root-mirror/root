## \file
## \ingroup tutorial_roofit
## \notebook
## Likelihood and minimization: setting up a multi-core parallelized unbinned maximum likelihood fit
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
model = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(
    sig, bkg), ROOT.RooArgList(fsig))

# Generate large dataset
data = model.generate(ROOT.RooArgSet(x, y, z), 200000)

# Parallel fitting
# -------------------------------

# In parallel mode the likelihood calculation is split in N pieces,
# that are calculated in parallel and added a posteriori before passing
# it back to MINUIT.

# Use four processes and time results both in wall time and CPU time
model.fitTo(data, ROOT.RooFit.NumCPU(4), ROOT.RooFit.Timer(ROOT.kTRUE))

# Parallel MC projections
# ----------------------------------------------

# Construct signal, likelihood projection on (y,z) observables and
# likelihood ratio
sigyz = sig.createProjection(ROOT.RooArgSet(x))
totyz = model.createProjection(ROOT.RooArgSet(x))
llratio_func = ROOT.RooFormulaVar(
    "llratio", "log10(@0)-log10(@1)", ROOT.RooArgList(sigyz, totyz))

# Calculate likelihood ratio for each event, subset of events with high
# signal likelihood
data.addColumn(llratio_func)
dataSel = data.reduce(ROOT.RooFit.Cut("llratio>0.7"))

# Make plot frame and plot data
frame = x.frame(ROOT.RooFit.Title(
    "Projection on X with LLratio(y,z)>0.7"), ROOT.RooFit.Bins(40))
dataSel.plotOn(frame)

# Perform parallel projection using MC integration of pdf using given input dataSet.
# In self mode the data-weighted average of the pdf is calculated by splitting the
# input dataset in N equal pieces and calculating in parallel the weighted average
# one each subset. The N results of those calculations are then weighted into the
# final result

# Use four processes
model.plotOn(frame, ROOT.RooFit.ProjWData(dataSel), ROOT.RooFit.NumCPU(4))

c = ROOT.TCanvas("rf603_multicpu", "rf603_multicpu", 600, 600)
ROOT.gPad.SetLeftMargin(0.15)
frame.GetYaxis().SetTitleOffset(1.6)
frame.Draw()

c.SaveAs("rf603_multicpu.png")
