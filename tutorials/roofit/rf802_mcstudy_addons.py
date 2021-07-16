## \ingroup tutorial_roofit
## \notebook
##
## 'VALIDATION AND MC STUDIES' RooFit tutorial macro #802
##
## RooMCStudy: using separate fit and generator models, the chi^2 calculator model
##
## \macro_code
##
## \date February 2018
## \author Clemens Lange


import ROOT


# Create model
# -----------------------

# Observables, parameters
x = ROOT.RooRealVar("x", "x", -10, 10)
x.setBins(10)
mean = ROOT.RooRealVar("mean", "mean of gaussian", 0, -2.0, 1.8)
sigma = ROOT.RooRealVar("sigma", "width of gaussian", 5, 1, 10)

# Create Gaussian pdf
gauss = ROOT.RooGaussian("gauss", "gaussian PDF", x, mean, sigma)

# Create manager with chi^2 add-on module
# ----------------------------------------------------------------------------

# Create study manager for binned likelihood fits of a Gaussian pdf in 10
# bins
mcs = ROOT.RooMCStudy(gauss, ROOT.RooArgSet(x), ROOT.RooFit.Silence(), ROOT.RooFit.Binned())

# Add chi^2 calculator module to mcs
chi2mod = ROOT.RooChi2MCSModule()
mcs.addModule(chi2mod)

# Generate 1000 samples of 1000 events
mcs.generateAndFit(2000, 1000)

# Fill histograms with distributions chi2 and prob(chi2,ndf) that
# are calculated by ROOT.RooChiMCSModule
hist_chi2 = ROOT.RooAbsData.createHistogram(mcs.fitParDataSet(), "chi2")
hist_prob = ROOT.RooAbsData.createHistogram(mcs.fitParDataSet(), "prob")

# Create manager with separate fit model
# ----------------------------------------------------------------------------

# Create alternate pdf with shifted mean
mean2 = ROOT.RooRealVar("mean2", "mean of gaussian 2", 2.0)
gauss2 = ROOT.RooGaussian("gauss2", "gaussian PDF2", x, mean2, sigma)

# Create study manager with separate generation and fit model. ROOT.This configuration
# is set up to generate bad fits as the fit and generator model have different means
# and the mean parameter is not floating in the fit
mcs2 = ROOT.RooMCStudy(
    gauss2, ROOT.RooArgSet(x), ROOT.RooFit.FitModel(gauss), ROOT.RooFit.Silence(), ROOT.RooFit.Binned()
)

# Add chi^2 calculator module to mcs
chi2mod2 = ROOT.RooChi2MCSModule()
mcs2.addModule(chi2mod2)

# Generate 1000 samples of 1000 events
mcs2.generateAndFit(2000, 1000)

# Request a the pull plot of mean. The pulls will be one-sided because
# `mean` is limited to 1.8.
# Note that RooFit will have trouble to compute the pulls because the parameters
# are called `mean` in the fit, but `mean2` in the generator model. It is not obvious
# that these are related. RooFit will nevertheless compute pulls, but complain that
# this is risky.
pullMeanFrame = mcs2.plotPull(mean)

# Fill histograms with distributions chi2 and prob(chi2,ndf) that
# are calculated by ROOT.RooChiMCSModule
hist2_chi2 = ROOT.RooAbsData.createHistogram(mcs2.fitParDataSet(), "chi2")
hist2_prob = ROOT.RooAbsData.createHistogram(mcs2.fitParDataSet(), "prob")
hist2_chi2.SetLineColor(ROOT.kRed)
hist2_prob.SetLineColor(ROOT.kRed)

c = ROOT.TCanvas("rf802_mcstudy_addons", "rf802_mcstudy_addons", 800, 400)
c.Divide(3)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
hist_chi2.GetYaxis().SetTitleOffset(1.4)
hist_chi2.Draw()
hist2_chi2.Draw("esame")
c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
hist_prob.GetYaxis().SetTitleOffset(1.4)
hist_prob.Draw()
hist2_prob.Draw("esame")
c.cd(3)
pullMeanFrame.Draw()

c.SaveAs("rf802_mcstudy_addons.png")

# Make RooMCStudy object available on command line after
# macro finishes
ROOT.gDirectory.Add(mcs)
