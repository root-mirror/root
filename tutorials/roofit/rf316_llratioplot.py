## \file
## \ingroup tutorial_roofit
## \notebook
## Multidimensional models: using the likelihood ratio techique to construct a signal
## enhanced one-dimensional projection of a multi-dimensional pdf
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

# Define projected signal likelihood ratio
# ----------------------------------------------------------------------------------

# Calculate projection of signal and total likelihood on (y,z) observables
# i.e. integrate signal and composite model over x
sigyz = sig.createProjection(ROOT.RooArgSet(x))
totyz = model.createProjection(ROOT.RooArgSet(x))

# Construct the log of the signal / signal+background probability
llratio_func = ROOT.RooFormulaVar(
    "llratio", "log10(@0)-log10(@1)", ROOT.RooArgList(sigyz, totyz))

# Plot data with a LL ratio cut
# -------------------------------------------------------

# Calculate the llratio value for each event in the dataset
data.addColumn(llratio_func)

# Extract the subset of data with large signal likelihood
dataSel = data.reduce(ROOT.RooFit.Cut("llratio>0.7"))

# Make plot frame
frame2 = x.frame(ROOT.RooFit.Title(
    "Same projection on X with LLratio(y,z)>0.7"), ROOT.RooFit.Bins(40))

# Plot select data on frame
dataSel.plotOn(frame2)

# Make MC projection of pdf with same LL ratio cut
# ---------------------------------------------------------------------------------------------

# Generate large number of events for MC integration of pdf projection
mcprojData = model.generate(ROOT.RooArgSet(x, y, z), 10000)

# Calculate LL ratio for each generated event and select MC events with
# llratio)0.7
mcprojData.addColumn(llratio_func)
mcprojDataSel = mcprojData.reduce(ROOT.RooFit.Cut("llratio>0.7"))

# Project model on x, projected observables (y,z) with Monte Carlo technique
# on set of events with the same llratio cut as was applied to data
model.plotOn(frame2, ROOT.RooFit.ProjWData(mcprojDataSel))

c = ROOT.TCanvas("rf316_llratioplot", "rf316_llratioplot", 800, 400)
c.Divide(2)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
frame.GetYaxis().SetTitleOffset(1.4)
frame.Draw()
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
frame2.GetYaxis().SetTitleOffset(1.4)
frame2.Draw()
c.SaveAs("rf316_llratioplot.png")
