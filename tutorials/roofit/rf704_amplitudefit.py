## \file
## \ingroup tutorial_roofit
## \notebook
## Special pdf's: using a pdf defined by a sum of real-valued amplitude components
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT


# Setup 2D amplitude functions
# -------------------------------------------------------

# Observables
t = ROOT.RooRealVar("t", "time", -1., 15.)
cosa = ROOT.RooRealVar("cosa", "cos(alpha)", -1., 1.)

# Use ROOT.RooTruthModel to obtain compiled implementation of sinh/cosh
# modulated decay functions
tau = ROOT.RooRealVar("tau", "#tau", 1.5)
deltaGamma = ROOT.RooRealVar("deltaGamma", "deltaGamma", 0.3)
tm = ROOT.RooTruthModel("tm", "tm", t)
coshGBasis = ROOT.RooFormulaVar(
    "coshGBasis",
    "exp(-@0/ @1)*cosh(@0*@2/2)",
    ROOT.RooArgList(
        t,
        tau,
        deltaGamma))
sinhGBasis = ROOT.RooFormulaVar(
    "sinhGBasis",
    "exp(-@0/ @1)*sinh(@0*@2/2)",
    ROOT.RooArgList(
        t,
        tau,
        deltaGamma))
coshGConv = tm.convolution(coshGBasis, t)
sinhGConv = tm.convolution(sinhGBasis, t)

# Construct polynomial amplitudes in cos(a)
poly1 = ROOT.RooPolyVar(
    "poly1",
    "poly1",
    cosa,
    ROOT.RooArgList(
        ROOT.RooFit.RooConst(0.5),
        ROOT.RooFit.RooConst(0.2),
        ROOT.RooFit.RooConst(0.2)),
    0)
poly2 = ROOT.RooPolyVar("poly2", "poly2", cosa, ROOT.RooArgList(
    ROOT.RooFit.RooConst(1), ROOT.RooFit.RooConst(-0.2), ROOT.RooFit.RooConst(3)), 0)

# Construct 2D amplitude as uncorrelated product of amp(t)*amp(cosa)
ampl1 = ROOT.RooProduct("ampl1", "amplitude 1",
                        ROOT.RooArgList(poly1, coshGConv))
ampl2 = ROOT.RooProduct("ampl2", "amplitude 2",
                        ROOT.RooArgList(poly2, sinhGConv))

# Construct amplitude sum pdf
# -----------------------------------------------------

# Amplitude strengths
f1 = ROOT.RooRealVar("f1", "f1", 1, 0, 2)
f2 = ROOT.RooRealVar("f2", "f2", 0.5, 0, 2)

# Construct pdf
pdf = ROOT.RooRealSumPdf("pdf", "pdf", ROOT.RooArgList(
    ampl1, ampl2), ROOT.RooArgList(f1, f2))

# Generate some toy data from pdf
data = pdf.generate(ROOT.RooArgSet(t, cosa), 10000)

# Fit pdf to toy data with only amplitude strength floating
pdf.fitTo(data)

# Plot amplitude sum pdf
# -------------------------------------------

# Make 2D plots of amplitudes
hh_cos = ampl1.createHistogram("hh_cos", t, ROOT.RooFit.Binning(
    50), ROOT.RooFit.YVar(cosa, ROOT.RooFit.Binning(50)))
hh_sin = ampl2.createHistogram("hh_sin", t, ROOT.RooFit.Binning(
    50), ROOT.RooFit.YVar(cosa, ROOT.RooFit.Binning(50)))
hh_cos.SetLineColor(ROOT.kBlue)
hh_sin.SetLineColor(ROOT.kRed)

# Make projection on t, data, and its components
# Note component projections may be larger than sum because amplitudes can
# be negative
frame1 = t.frame()
data.plotOn(frame1)
pdf.plotOn(frame1)
# workaround, see https://root.cern.ch/phpBB3/viewtopic.php?t=7764
ras_ampl1 = ROOT.RooArgSet(ampl1)
pdf.plotOn(frame1, ROOT.RooFit.Components(ras_ampl1),
           ROOT.RooFit.LineStyle(ROOT.kDashed))
ras_ampl2 = ROOT.RooArgSet(ampl2)
pdf.plotOn(frame1, ROOT.RooFit.Components(ras_ampl2), ROOT.RooFit.LineStyle(
    ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

# Make projection on cosa, data, and its components
# Note that components projection may be larger than sum because
# amplitudes can be negative
frame2 = cosa.frame()
data.plotOn(frame2)
pdf.plotOn(frame2)
pdf.plotOn(frame2, ROOT.RooFit.Components(ras_ampl1),
           ROOT.RooFit.LineStyle(ROOT.kDashed))
pdf.plotOn(frame2, ROOT.RooFit.Components(ras_ampl2), ROOT.RooFit.LineStyle(
    ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))

c = ROOT.TCanvas("rf704_amplitudefit", "rf704_amplitudefit", 800, 800)
c.Divide(2, 2)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
frame1.GetYaxis().SetTitleOffset(1.4)
frame1.Draw()
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
frame2.GetYaxis().SetTitleOffset(1.4)
frame2.Draw()
c.cd(3)
ROOT.gPad.SetLeftMargin(0.20)
hh_cos.GetZaxis().SetTitleOffset(2.3)
hh_cos.Draw("surf")
c.cd(4)
ROOT.gPad.SetLeftMargin(0.20)
hh_sin.GetZaxis().SetTitleOffset(2.3)
hh_sin.Draw("surf")

c.SaveAs("rf704_amplitudefit.png")
