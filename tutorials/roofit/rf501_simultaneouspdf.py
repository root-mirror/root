## \file
## \ingroup tutorial_roofit
## \notebook
## Organization and simultaneous fits: using simultaneous pdfs to describe simultaneous
## fits to multiple datasets
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT


# Create model for physics sample
# -------------------------------------------------------------

# Create observables
x = ROOT.RooRealVar("x", "x", -8, 8)

# Construct signal pdf
mean = ROOT.RooRealVar("mean", "mean", 0, -8, 8)
sigma = ROOT.RooRealVar("sigma", "sigma", 0.3, 0.1, 10)
gx = ROOT.RooGaussian("gx", "gx", x, mean, sigma)

# Construct background pdf
a0 = ROOT.RooRealVar("a0", "a0", -0.1, -1, 1)
a1 = ROOT.RooRealVar("a1", "a1", 0.004, -1, 1)
px = ROOT.RooChebychev("px", "px", x, ROOT.RooArgList(a0, a1))

# Construct composite pdf
f = ROOT.RooRealVar("f", "f", 0.2, 0., 1.)
model = ROOT.RooAddPdf(
    "model", "model", ROOT.RooArgList(gx, px), ROOT.RooArgList(f))

# Create model for control sample
# --------------------------------------------------------------

# Construct signal pdf.
# NOTE that sigma is shared with the signal sample model
mean_ctl = ROOT.RooRealVar("mean_ctl", "mean_ctl", -3, -8, 8)
gx_ctl = ROOT.RooGaussian("gx_ctl", "gx_ctl", x, mean_ctl, sigma)

# Construct the background pdf
a0_ctl = ROOT.RooRealVar("a0_ctl", "a0_ctl", -0.1, -1, 1)
a1_ctl = ROOT.RooRealVar("a1_ctl", "a1_ctl", 0.5, -0.1, 1)
px_ctl = ROOT.RooChebychev(
    "px_ctl", "px_ctl", x, ROOT.RooArgList(a0_ctl, a1_ctl))

# Construct the composite model
f_ctl = ROOT.RooRealVar("f_ctl", "f_ctl", 0.5, 0., 1.)
model_ctl = ROOT.RooAddPdf(
    "model_ctl",
    "model_ctl",
    ROOT.RooArgList(
        gx_ctl,
        px_ctl),
    ROOT.RooArgList(f_ctl))

# Generate events for both samples
# ---------------------------------------------------------------

# Generate 1000 events in x and y from model
data = model.generate(ROOT.RooArgSet(x), 100)
data_ctl = model_ctl.generate(ROOT.RooArgSet(x), 2000)

# Create index category and join samples
# ---------------------------------------------------------------------------

# Define category to distinguish physics and control samples events
sample = ROOT.RooCategory("sample", "sample")
sample.defineType("physics")
sample.defineType("control")

# Construct combined dataset in (x,sample)
combData = ROOT.RooDataSet(
    "combData",
    "combined data",
    ROOT.RooArgSet(x),
    ROOT.RooFit.Index(sample),
    ROOT.RooFit.Import(
        "physics",
        data),
    ROOT.RooFit.Import(
        "control",
        data_ctl))

# Construct a simultaneous pdf in (x, sample)
# -----------------------------------------------------------------------------------

# Construct a simultaneous pdf using category sample as index
simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", sample)

# Associate model with the physics state and model_ctl with the control
# state
simPdf.addPdf(model, "physics")
simPdf.addPdf(model_ctl, "control")

# Perform a simultaneous fit
# ---------------------------------------------------

# Perform simultaneous fit of model to data and model_ctl to data_ctl
simPdf.fitTo(combData)

# Plot model slices on data slices
# ----------------------------------------------------------------

# Make a frame for the physics sample
frame1 = x.frame(ROOT.RooFit.Bins(30), ROOT.RooFit.Title("Physics sample"))

# Plot all data tagged as physics sample
combData.plotOn(frame1, ROOT.RooFit.Cut("sample==sample::physics"))

# Plot "physics" slice of simultaneous pdf.
# NB: You *must* project the sample index category with data using ProjWData
# as a RooSimultaneous makes no prediction on the shape in the index category
# and can thus not be integrated
# NB2: The sampleSet *must* be named. It will not work to pass this as a temporary
# because python will delete it. The same holds for fitTo() and plotOn() below.
sampleSet = ROOT.RooArgSet(sample)
simPdf.plotOn(frame1, ROOT.RooFit.Slice(sample, "physics"), ROOT.RooFit.Components(
    "px"), ROOT.RooFit.ProjWData(sampleSet, combData), ROOT.RooFit.LineStyle(ROOT.kDashed))

# The same plot for the control sample slice
frame2 = x.frame(ROOT.RooFit.Bins(30), ROOT.RooFit.Title("Control sample"))
combData.plotOn(frame2, ROOT.RooFit.Cut("sample==sample::control"))
simPdf.plotOn(frame2, ROOT.RooFit.Slice(sample, "control"),
                ROOT.RooFit.ProjWData(sampleSet, combData))
simPdf.plotOn(frame2, ROOT.RooFit.Slice(sample, "control"), ROOT.RooFit.Components(
    "px_ctl"), ROOT.RooFit.ProjWData(sampleSet, combData), ROOT.RooFit.LineStyle(ROOT.kDashed))

c = ROOT.TCanvas("rf501_simultaneouspdf",
                 "rf501_simultaneouspdf", 800, 400)
c.Divide(2)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
frame1.GetYaxis().SetTitleOffset(1.4)
frame1.Draw()
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
frame2.GetYaxis().SetTitleOffset(1.4)
frame2.Draw()

c.SaveAs("rf501_simultaneouspdf.png")
