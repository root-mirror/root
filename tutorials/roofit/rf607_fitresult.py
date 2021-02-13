## \file
## \ingroup tutorial_roofit
## \notebook
## Likelihood and minimization: demonstration of options of the RooFitResult class
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

from __future__ import print_function
import ROOT


# Create pdf, data
# --------------------------------

# Declare observable x
x = ROOT.RooRealVar("x", "x", 0, 10)

# Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and
# their parameters
mean = ROOT.RooRealVar("mean", "mean of gaussians", 5, -10, 10)
sigma1 = ROOT.RooRealVar("sigma1", "width of gaussians", 0.5, 0.1, 10)
sigma2 = ROOT.RooRealVar("sigma2", "width of gaussians", 1, 0.1, 10)

sig1 = ROOT.RooGaussian("sig1", "Signal component 1", x, mean, sigma1)
sig2 = ROOT.RooGaussian("sig2", "Signal component 2", x, mean, sigma2)

# Build Chebychev polynomial pdf
a0 = ROOT.RooRealVar("a0", "a0", 0.5, 0., 1.)
a1 = ROOT.RooRealVar("a1", "a1", -0.2)
bkg = ROOT.RooChebychev("bkg", "Background", x, ROOT.RooArgList(a0, a1))

# Sum the signal components into a composite signal pdf
sig1frac = ROOT.RooRealVar(
    "sig1frac", "fraction of component 1 in signal", 0.8, 0., 1.)
sig = ROOT.RooAddPdf(
    "sig", "Signal", ROOT.RooArgList(sig1, sig2), ROOT.RooArgList(sig1frac))

# Sum the composite signal and background
bkgfrac = ROOT.RooRealVar("bkgfrac", "fraction of background", 0.5, 0., 1.)
model = ROOT.RooAddPdf(
    "model", "g1+g2+a", ROOT.RooArgList(bkg, sig), ROOT.RooArgList(bkgfrac))

# Generate 1000 events
data = model.generate(ROOT.RooArgSet(x), 1000)

# Fit pdf to data, save fit result
# -------------------------------------------------------------

# Perform fit and save result
r = model.fitTo(data, ROOT.RooFit.Save())

# Print fit results
# ---------------------------------

# Summary printing: Basic info plus final values of floating fit parameters
r.Print()

# Verbose printing: Basic info, of constant parameters, and
# final values of floating parameters, correlations
r.Print("v")

# Visualize correlation matrix
# -------------------------------------------------------

# Construct 2D color plot of correlation matrix
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
hcorr = r.correlationHist()

# Visualize ellipse corresponding to single correlation matrix element
frame = ROOT.RooPlot(sigma1, sig1frac, 0.45, 0.60, 0.65, 0.90)
frame.SetTitle("Covariance between sigma1 and sig1frac")
r.plotOn(frame, sigma1, sig1frac, "ME12ABHV")

# Access fit result information
# ---------------------------------------------------------

# Access basic information
print("EDM = ", r.edm())
print("-log(L) minimum = ", r.minNll())

# Access list of final fit parameter values
print("final value of floating parameters")
r.floatParsFinal().Print("s")

# Access correlation matrix elements
print("correlation between sig1frac and a0 is  ", r.correlation(
    sig1frac, a0))
print("correlation between bkgfrac and mean is ", r.correlation(
    "bkgfrac", "mean"))

# Extract covariance and correlation matrix as ROOT.TMatrixDSym
cor = r.correlationMatrix()
cov = r.covarianceMatrix()

# Print correlation, matrix
print("correlation matrix")
cor.Print()
print("covariance matrix")
cov.Print()

# Persist fit result in root file
# -------------------------------------------------------------

# Open ROOT file save save result
f = ROOT.TFile("rf607_fitresult.root", "RECREATE")
r.Write("rf607")
f.Close()

# In a clean ROOT session retrieve the persisted fit result as follows:
# r = gDirectory.Get("rf607")

c = ROOT.TCanvas("rf607_fitresult", "rf607_fitresult", 800, 400)
c.Divide(2)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
hcorr.GetYaxis().SetTitleOffset(1.4)
hcorr.Draw("colz")
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
frame.GetYaxis().SetTitleOffset(1.6)
frame.Draw()

c.SaveAs("rf607_fitresult.png")
