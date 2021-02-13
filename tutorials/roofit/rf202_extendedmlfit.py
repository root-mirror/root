## \file
## \ingroup tutorial_roofit
## \notebook
## Addition and convolution: setting up an extended maximum likelihood fit
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT

# Set up component pdfs
# ---------------------------------------

# Declare observable x
x = ROOT.RooRealVar("x", "x", 0, 10)

# Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and
# their parameters
mean = ROOT.RooRealVar("mean", "mean of gaussians", 5)
sigma1 = ROOT.RooRealVar("sigma1", "width of gaussians", 0.5)
sigma2 = ROOT.RooRealVar("sigma2", "width of gaussians", 1)

sig1 = ROOT.RooGaussian("sig1", "Signal component 1", x, mean, sigma1)
sig2 = ROOT.RooGaussian("sig2", "Signal component 2", x, mean, sigma2)

# Build Chebychev polynomial pdf
a0 = ROOT.RooRealVar("a0", "a0", 0.5, 0., 1.)
a1 = ROOT.RooRealVar("a1", "a1", -0.2, 0., 1.)
bkg = ROOT.RooChebychev("bkg", "Background", x, ROOT.RooArgList(a0, a1))

# Sum the signal components into a composite signal pdf
sig1frac = ROOT.RooRealVar(
    "sig1frac", "fraction of component 1 in signal", 0.8, 0., 1.)
sig = ROOT.RooAddPdf(
    "sig", "Signal", ROOT.RooArgList(sig1, sig2), ROOT.RooArgList(sig1frac))

# Method 1 - Construct extended composite model
# -------------------------------------------------------------------

# Sum the composite signal and background into an extended pdf
# nsig*sig+nbkg*bkg
nsig = ROOT.RooRealVar("nsig", "number of signal events", 500, 0., 10000)
nbkg = ROOT.RooRealVar(
    "nbkg", "number of background events", 500, 0, 10000)
model = ROOT.RooAddPdf(
    "model",
    "(g1+g2)+a",
    ROOT.RooArgList(
        bkg,
        sig),
    ROOT.RooArgList(
        nbkg,
        nsig))

# Sample, fit and plot extended model
# ---------------------------------------------------------------------

# Generate a data sample of expected number events in x from model
# = model.expectedEvents() = nsig+nbkg
data = model.generate(ROOT.RooArgSet(x))

# Fit model to data, ML term automatically included
model.fitTo(data)

# Plot data and PDF overlaid, expected number of events for pdf projection normalization
# rather than observed number of events (==data.numEntries())
xframe = x.frame(ROOT.RooFit.Title("extended ML fit example"))
data.plotOn(xframe)
model.plotOn(xframe, ROOT.RooFit.Normalization(
    1.0, ROOT.RooAbsReal.RelativeExpected))

# Overlay the background component of model with a dashed line
ras_bkg = ROOT.RooArgSet(bkg)
model.plotOn(
    xframe, ROOT.RooFit.Components(ras_bkg), ROOT.RooFit.LineStyle(
        ROOT.kDashed), ROOT.RooFit.Normalization(
            1.0, ROOT.RooAbsReal.RelativeExpected))

# Overlay the background+sig2 components of model with a dotted line
ras_bkg_sig2 = ROOT.RooArgSet(bkg, sig2)
model.plotOn(
    xframe, ROOT.RooFit.Components(ras_bkg_sig2), ROOT.RooFit.LineStyle(
        ROOT.kDotted), ROOT.RooFit.Normalization(
            1.0, ROOT.RooAbsReal.RelativeExpected))

# Print structure of composite pdf
model.Print("t")


# Method 2 - Construct extended components first
# ---------------------------------------------------------------------

# Associated nsig/nbkg as expected number of events with sig/bkg
esig = ROOT.RooExtendPdf("esig", "extended signal pdf", sig, nsig)
ebkg = ROOT.RooExtendPdf("ebkg", "extended background pdf", bkg, nbkg)

# Sum extended components without coefs
# -------------------------------------------------------------------------

# Construct sum of two extended pdf (no coefficients required)
model2 = ROOT.RooAddPdf("model2", "(g1+g2)+a", ROOT.RooArgList(ebkg, esig))

# Draw the frame on the canvas
c = ROOT.TCanvas("rf202_extendedmlfit", "rf202_extendedmlfit", 600, 600)
ROOT.gPad.SetLeftMargin(0.15)
xframe.GetYaxis().SetTitleOffset(1.4)
xframe.Draw()

c.SaveAs("rf202_extendedmlfit.png")
