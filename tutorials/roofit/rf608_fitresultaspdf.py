## \file
## \ingroup tutorial_roofit
## \notebook
## Likelihood and minimization: representing the parabolic approximation of the fit as a
## multi-variate Gaussian on the parameters of the fitted pdf
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT


# Create model and dataset
# -----------------------------------------------

# Observable
x = ROOT.RooRealVar("x", "x", -20, 20)

# Model (intentional strong correlations)
mean = ROOT.RooRealVar("mean", "mean of g1 and g2", 0, -1, 1)
sigma_g1 = ROOT.RooRealVar("sigma_g1", "width of g1", 2)
g1 = ROOT.RooGaussian("g1", "g1", x, mean, sigma_g1)

sigma_g2 = ROOT.RooRealVar("sigma_g2", "width of g2", 4, 3.0, 5.0)
g2 = ROOT.RooGaussian("g2", "g2", x, mean, sigma_g2)

frac = ROOT.RooRealVar("frac", "frac", 0.5, 0.0, 1.0)
model = ROOT.RooAddPdf(
    "model", "model", ROOT.RooArgList(
        g1, g2), ROOT.RooArgList(frac))

# Generate 1000 events
data = model.generate(ROOT.RooArgSet(x), 1000)

# Fit model to data
# ----------------------------------

r = model.fitTo(data, ROOT.RooFit.Save())

# Create MV Gaussian pdf of fitted parameters
# ------------------------------------------------------------------------------------

parabPdf = r.createHessePdf(ROOT.RooArgSet(frac, mean, sigma_g2))

# Some exercises with the parameter pdf
# -----------------------------------------------------------------------------

# Generate 100K points in the parameter space, from the MVGaussian pdf
d = parabPdf.generate(ROOT.RooArgSet(mean, sigma_g2, frac), 100000)

# Sample a 3-D histogram of the pdf to be visualized as an error
# ellipsoid using the GLISO draw option
hh_3d = parabPdf.createHistogram("mean,sigma_g2,frac", 25, 25, 25)
hh_3d.SetFillColor(ROOT.kBlue)

# Project 3D parameter pdf down to 3 permutations of two-dimensional pdfs
# The integrations corresponding to these projections are performed analytically
# by the MV Gaussian pdf
pdf_sigmag2_frac = parabPdf.createProjection(ROOT.RooArgSet(mean))
pdf_mean_frac = parabPdf.createProjection(ROOT.RooArgSet(sigma_g2))
pdf_mean_sigmag2 = parabPdf.createProjection(ROOT.RooArgSet(frac))

# Make 2D plots of the 3 two-dimensional pdf projections
hh_sigmag2_frac = pdf_sigmag2_frac.createHistogram("sigma_g2,frac", 50, 50)
hh_mean_frac = pdf_mean_frac.createHistogram("mean,frac", 50, 50)
hh_mean_sigmag2 = pdf_mean_sigmag2.createHistogram("mean,sigma_g2", 50, 50)
hh_mean_frac.SetLineColor(ROOT.kBlue)
hh_sigmag2_frac.SetLineColor(ROOT.kBlue)
hh_mean_sigmag2.SetLineColor(ROOT.kBlue)

# Draw the 'sigar'
ROOT.gStyle.SetCanvasPreferGL(True)
ROOT.gStyle.SetPalette(1)
c1 = ROOT.TCanvas("rf608_fitresultaspdf_1", "rf608_fitresultaspdf_1", 600, 600)
hh_3d.Draw("gliso")

c1.SaveAs("rf608_fitresultaspdf_1.png")

# Draw the 2D projections of the 3D pdf
c2 = ROOT.TCanvas("rf608_fitresultaspdf_2",
                  "rf608_fitresultaspdf_2", 900, 600)
c2.Divide(3, 2)
c2.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
hh_mean_sigmag2.GetZaxis().SetTitleOffset(1.4)
hh_mean_sigmag2.Draw("surf3")
c2.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
hh_sigmag2_frac.GetZaxis().SetTitleOffset(1.4)
hh_sigmag2_frac.Draw("surf3")
c2.cd(3)
ROOT.gPad.SetLeftMargin(0.15)
hh_mean_frac.GetZaxis().SetTitleOffset(1.4)
hh_mean_frac.Draw("surf3")

# Draw the distributions of parameter points sampled from the pdf
tmp1 = d.createHistogram(mean, sigma_g2, 50, 50)
tmp2 = d.createHistogram(sigma_g2, frac, 50, 50)
tmp3 = d.createHistogram(mean, frac, 50, 50)

c2.cd(4)
ROOT.gPad.SetLeftMargin(0.15)
tmp1.GetZaxis().SetTitleOffset(1.4)
tmp1.Draw("lego3")
c2.cd(5)
ROOT.gPad.SetLeftMargin(0.15)
tmp2.GetZaxis().SetTitleOffset(1.4)
tmp2.Draw("lego3")
c2.cd(6)
ROOT.gPad.SetLeftMargin(0.15)
tmp3.GetZaxis().SetTitleOffset(1.4)
tmp3.Draw("lego3")

c2.SaveAs("rf608_fitresultaspdf_2.png")
