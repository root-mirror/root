## \file
## \ingroup tutorial_roofit
## \notebook
## Multidimensional models: multi-dimensional pdfs with conditional pdfs in product
##
## `pdf = gauss(x,f(y),sx | y ) * gauss(y,ms,sx)`    with `f(y) = a0 + a1*y`
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT

# Create conditional pdf gx(x|y)
# -----------------------------------------------------------

# Create observables
x = ROOT.RooRealVar("x", "x", -5, 5)
y = ROOT.RooRealVar("y", "y", -5, 5)

# Create function f(y) = a0 + a1*y
a0 = ROOT.RooRealVar("a0", "a0", -0.5, -5, 5)
a1 = ROOT.RooRealVar("a1", "a1", -0.5, -1, 1)
fy = ROOT.RooPolyVar("fy", "fy", y, ROOT.RooArgList(a0, a1))

# Create gaussx(x,f(y),sx)
sigmax = ROOT.RooRealVar("sigma", "width of gaussian", 0.5)
gaussx = ROOT.RooGaussian(
    "gaussx", "Gaussian in x with shifting mean in y", x, fy, sigmax)

# Create pdf gy(y)
# -----------------------------------------------------------

# Create gaussy(y,0,5)
gaussy = ROOT.RooGaussian(
    "gaussy",
    "Gaussian in y",
    y,
    ROOT.RooFit.RooConst(0),
    ROOT.RooFit.RooConst(3))

# Create product gx(x|y)*gy(y)
# -------------------------------------------------------

# Create gaussx(x,sx|y) * gaussy(y)
model = ROOT.RooProdPdf(
    "model",
    "gaussx(x|y)*gaussy(y)",
    ROOT.RooArgSet(gaussy),
    ROOT.RooFit.Conditional(
        ROOT.RooArgSet(gaussx),
        ROOT.RooArgSet(x)))

# Sample, fit and plot product pdf
# ---------------------------------------------------------------

# Generate 1000 events in x and y from model
data = model.generate(ROOT.RooArgSet(x, y), 10000)

# Plot x distribution of data and projection of model x = Int(dy)
# model(x,y)
xframe = x.frame()
data.plotOn(xframe)
model.plotOn(xframe)

# Plot x distribution of data and projection of model y = Int(dx)
# model(x,y)
yframe = y.frame()
data.plotOn(yframe)
model.plotOn(yframe)

# Make two-dimensional plot in x vs y
hh_model = model.createHistogram("hh_model", x, ROOT.RooFit.Binning(
    50), ROOT.RooFit.YVar(y, ROOT.RooFit.Binning(50)))
hh_model.SetLineColor(ROOT.kBlue)

# Make canvas and draw ROOT.RooPlots
c = ROOT.TCanvas("rf305_condcorrprod", "rf05_condcorrprod", 1200, 400)
c.Divide(3)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
xframe.GetYaxis().SetTitleOffset(1.6)
xframe.Draw()
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
yframe.GetYaxis().SetTitleOffset(1.6)
yframe.Draw()
c.cd(3)
ROOT.gPad.SetLeftMargin(0.20)
hh_model.GetZaxis().SetTitleOffset(2.5)
hh_model.Draw("surf")

c.SaveAs("rf305_condcorrprod.png")
