## \file
## \ingroup tutorial_roofit
## \notebook
## Multidimensional models: usage of full pdf with per-event errors
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT

# B-physics pdf with per-event Gaussian resolution
# ----------------------------------------------------------------------------------------------

# Observables
dt = ROOT.RooRealVar("dt", "dt", -10, 10)
dterr = ROOT.RooRealVar("dterr", "per-event error on dt", 0.01, 10)

# Build a gaussian resolution model scaled by the per-error =
# gauss(dt,bias,sigma*dterr)
bias = ROOT.RooRealVar("bias", "bias", 0, -10, 10)
sigma = ROOT.RooRealVar(
    "sigma", "per-event error scale factor", 1, 0.1, 10)
gm = ROOT.RooGaussModel(
    "gm1", "gauss model scaled bt per-event error", dt, bias, sigma, dterr)

# Construct decay(dt) (x) gauss1(dt|dterr)
tau = ROOT.RooRealVar("tau", "tau", 1.548)
decay_gm = ROOT.RooDecay("decay_gm", "decay", dt,
                         tau, gm, ROOT.RooDecay.DoubleSided)

# Construct empirical pdf for per-event error
# -----------------------------------------------------------------

# Use landau pdf to get empirical distribution with long tail
pdfDtErr = ROOT.RooLandau("pdfDtErr", "pdfDtErr", dterr, ROOT.RooFit.RooConst(
    1), ROOT.RooFit.RooConst(0.25))
expDataDterr = pdfDtErr.generate(ROOT.RooArgSet(dterr), 10000)

# Construct a histogram pdf to describe the shape of the dtErr distribution
expHistDterr = expDataDterr.binnedClone()
pdfErr = ROOT.RooHistPdf(
    "pdfErr", "pdfErr", ROOT.RooArgSet(dterr), expHistDterr)

# Construct conditional product decay_dm(dt|dterr)*pdf(dterr)
# ----------------------------------------------------------------------------------------------------------------------

# Construct production of conditional decay_dm(dt|dterr) with empirical
# pdfErr(dterr)
model = ROOT.RooProdPdf(
    "model",
    "model",
    ROOT.RooArgSet(pdfErr),
    ROOT.RooFit.Conditional(
        ROOT.RooArgSet(decay_gm),
        ROOT.RooArgSet(dt)))

# (Alternatively you could also use the landau shape pdfDtErr)
# ROOT.RooProdPdf model("model", "model",pdfDtErr,
# ROOT.RooFit.Conditional(decay_gm,dt))

# Sample, fit and plot product model
# ------------------------------------------------------------------

# Specify external dataset with dterr values to use model_dm as
# conditional pdf
data = model.generate(ROOT.RooArgSet(dt, dterr), 10000)

# Fit conditional decay_dm(dt|dterr)
# ---------------------------------------------------------------------

# Specify dterr as conditional observable
model.fitTo(data)

# Plot conditional decay_dm(dt|dterr)
# ---------------------------------------------------------------------

# Make two-dimensional plot of conditional pdf in (dt,dterr)
hh_model = model.createHistogram("hh_model", dt, ROOT.RooFit.Binning(
    50), ROOT.RooFit.YVar(dterr, ROOT.RooFit.Binning(50)))
hh_model.SetLineColor(ROOT.kBlue)

# Make projection of data an dt
frame = dt.frame(ROOT.RooFit.Title("Projection of model(dt|dterr) on dt"))
data.plotOn(frame)
model.plotOn(frame)

# Draw all frames on canvas
c = ROOT.TCanvas("rf307_fullpereventerrors",
                 "rf307_fullpereventerrors", 800, 400)
c.Divide(2)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.20)
hh_model.GetZaxis().SetTitleOffset(2.5)
hh_model.Draw("surf")
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
frame.GetYaxis().SetTitleOffset(1.6)
frame.Draw()

c.SaveAs("rf307_fullpereventerrors.png")
