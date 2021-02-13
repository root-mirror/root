## \file
## \ingroup tutorial_roofit
## \notebook
## Addition and convolution: decay function pdfs with optional B physics effects (mixing
## and CP violation) that can be analytically convolved with e.g. Gaussian resolution functions
##
## ```
## pdf1 = decay(t,tau) (x) delta(t)
## pdf2 = decay(t,tau) (x) gauss(t,m,s)
## pdf3 = decay(t,tau) (x) (f*gauss1(t,m1,s1) + (1-f)*gauss2(t,m1,s1))
## ```
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT

# B-physics pdf with truth resolution
# ---------------------------------------------------------------------

# Variables of decay pdf
dt = ROOT.RooRealVar("dt", "dt", -10, 10)
tau = ROOT.RooRealVar("tau", "tau", 1.548)

# Build a truth resolution model (delta function)
tm = ROOT.RooTruthModel("tm", "truth model", dt)

# Construct decay(t) (x) delta(t)
decay_tm = ROOT.RooDecay("decay_tm", "decay", dt,
                         tau, tm, ROOT.RooDecay.DoubleSided)

# Plot pdf (dashed)
frame = dt.frame(ROOT.RooFit.Title("Bdecay (x) resolution"))
decay_tm.plotOn(frame, ROOT.RooFit.LineStyle(ROOT.kDashed))

# B-physics pdf with Gaussian resolution
# ----------------------------------------------------------------------------

# Build a gaussian resolution model
bias1 = ROOT.RooRealVar("bias1", "bias1", 0)
sigma1 = ROOT.RooRealVar("sigma1", "sigma1", 1)
gm1 = ROOT.RooGaussModel("gm1", "gauss model 1", dt, bias1, sigma1)

# Construct decay(t) (x) gauss1(t)
decay_gm1 = ROOT.RooDecay("decay_gm1", "decay",
                          dt, tau, gm1, ROOT.RooDecay.DoubleSided)

# Plot pdf
decay_gm1.plotOn(frame)

# B-physics pdf with double Gaussian resolution
# ------------------------------------------------------------------------------------------

# Build another gaussian resolution model
bias2 = ROOT.RooRealVar("bias2", "bias2", 0)
sigma2 = ROOT.RooRealVar("sigma2", "sigma2", 5)
gm2 = ROOT.RooGaussModel("gm2", "gauss model 2", dt, bias2, sigma2)

# Build a composite resolution model f*gm1+(1-f)*gm2
gm1frac = ROOT.RooRealVar("gm1frac", "fraction of gm1", 0.5)
gmsum = ROOT.RooAddModel(
    "gmsum",
    "sum of gm1 and gm2",
    ROOT.RooArgList(
        gm1,
        gm2),
    ROOT.RooArgList(gm1frac))

# Construct decay(t) (x) (f*gm1 + (1-f)*gm2)
decay_gmsum = ROOT.RooDecay(
    "decay_gmsum", "decay", dt, tau, gmsum, ROOT.RooDecay.DoubleSided)

# Plot pdf (red)
decay_gmsum.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))

# Draw all frames on canvas
c = ROOT.TCanvas("rf209_anaconv", "rf209_anaconv", 600, 600)
ROOT.gPad.SetLeftMargin(0.15)
frame.GetYaxis().SetTitleOffset(1.6)
frame.Draw()

c.SaveAs("rf209_anaconv.png")
