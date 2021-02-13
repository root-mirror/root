## \file
## \ingroup tutorial_roofit
## \notebook
## Multidimensional models: making 2/3 dimensional plots of pdfs and datasets
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT

# Create 2D model and dataset
# -----------------------------------------------------

# Create observables
x = ROOT.RooRealVar("x", "x", -5, 5)
y = ROOT.RooRealVar("y", "y", -5, 5)

# Create parameters
a0 = ROOT.RooRealVar("a0", "a0", -3.5, -5, 5)
a1 = ROOT.RooRealVar("a1", "a1", -1.5, -1, 1)
sigma = ROOT.RooRealVar("sigma", "width of gaussian", 1.5)

# Create interpreted function f(y) = a0 - a1*sqrt(10*abs(y))
fy = ROOT.RooFormulaVar("fy", "a0-a1*sqrt(10*abs(y))",
                        ROOT.RooArgList(y, a0, a1))

# Create gauss(x,f(y),s)
model = ROOT.RooGaussian(
    "model", "Gaussian with shifting mean", x, fy, sigma)

# Sample dataset from gauss(x,y)
data = model.generate(ROOT.RooArgSet(x, y), 10000)

# Make 2D plots of data and model
# -------------------------------------------------------------

# Create and fill ROOT 2D histogram (20x20 bins) with contents of dataset
# hh_data = data.createHistogram("hh_data",x, ROOT.RooFit.Binning(20), ROOT.RooFit.YVar(y, ROOT.RooFit.Binning(20)))
# hh_data = data.createHistogram("x,y", 20, 20) # does not work, see
# https://root.cern.ch/phpBB3/viewtopic.php?t=16648
hh_data = ROOT.RooAbsData.createHistogram(data, "x,y", x, ROOT.RooFit.Binning(
    20), ROOT.RooFit.YVar(y, ROOT.RooFit.Binning(20)))

# Create and fill ROOT 2D histogram (50x50 bins) with sampling of pdf
# hh_pdf = model.createHistogram("hh_model",x, ROOT.RooFit.Binning(50), ROOT.RooFit.YVar(y, ROOT.RooFit.Binning(50)))
hh_pdf = model.createHistogram("x,y", 50, 50)
hh_pdf.SetLineColor(ROOT.kBlue)

# Create 3D model and dataset
# -----------------------------------------------------

# Create observables
z = ROOT.RooRealVar("z", "z", -5, 5)

gz = ROOT.RooGaussian(
    "gz", "gz", z, ROOT.RooFit.RooConst(0), ROOT.RooFit.RooConst(2))
model3 = ROOT.RooProdPdf("model3", "model3", ROOT.RooArgList(model, gz))

data3 = model3.generate(ROOT.RooArgSet(x, y, z), 10000)

# Make 3D plots of data and model
# -------------------------------------------------------------

# Create and fill ROOT 2D histogram (8x8x8 bins) with contents of dataset
# hh_data3 = data3.createHistogram("hh_data3", x, ROOT.RooFit.Binning(8), ROOT.RooFit.YVar(y, ROOT.RooFit.Binning(8)), ROOT.RooFit.ZVar(z, ROOT.RooFit.Binning(8)))
hh_data3 = ROOT.RooAbsData.createHistogram(
    data3, "hh_data3", x, ROOT.RooFit.Binning(8), ROOT.RooFit.YVar(
        y, ROOT.RooFit.Binning(8)), ROOT.RooFit.ZVar(
            z, ROOT.RooFit.Binning(8)))

# Create and fill ROOT 2D histogram (20x20x20 bins) with sampling of pdf
hh_pdf3 = model3.createHistogram(
    "hh_model3", x, ROOT.RooFit.Binning(20), ROOT.RooFit.YVar(
        y, ROOT.RooFit.Binning(20)), ROOT.RooFit.ZVar(
            z, ROOT.RooFit.Binning(20)))
hh_pdf3.SetFillColor(ROOT.kBlue)

c1 = ROOT.TCanvas("rf309_2dimplot", "rf309_2dimplot", 800, 800)
c1.Divide(2, 2)
c1.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
hh_data.GetZaxis().SetTitleOffset(1.4)
hh_data.Draw("lego")
c1.cd(2)
ROOT.gPad.SetLeftMargin(0.20)
hh_pdf.GetZaxis().SetTitleOffset(2.5)
hh_pdf.Draw("surf")
c1.cd(3)
ROOT.gPad.SetLeftMargin(0.15)
hh_data.GetZaxis().SetTitleOffset(1.4)
hh_data.Draw("box")
c1.cd(4)
ROOT.gPad.SetLeftMargin(0.15)
hh_pdf.GetZaxis().SetTitleOffset(2.5)
hh_pdf.Draw("cont3")
c1.SaveAs("rf309_2dimplot.png")

c2 = ROOT.TCanvas("rf309_3dimplot", "rf309_3dimplot", 800, 400)
c2.Divide(2)
c2.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
hh_data3.GetZaxis().SetTitleOffset(1.4)
hh_data3.Draw("lego")
c2.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
hh_pdf3.GetZaxis().SetTitleOffset(1.4)
hh_pdf3.Draw("iso")
c2.SaveAs("rf309_3dimplot.png")
