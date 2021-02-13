## \file
## \ingroup tutorial_roofit
## \notebook
## Data and categories: demonstration of real-discrete mapping functions
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT


# Define pdf in x, sample dataset in x
# ------------------------------------------------------------------------

# Define a dummy PDF in x
x = ROOT.RooRealVar("x", "x", 0, 10)
a = ROOT.RooArgusBG("a", "argus(x)", x, ROOT.RooFit.RooConst(
    10), ROOT.RooFit.RooConst(-1))

# Generate a dummy dataset
data = a.generate(ROOT.RooArgSet(x), 10000)

# Create a threshold real -> cat function
# --------------------------------------------------------------------------

# A RooThresholdCategory is a category function that maps regions in a real-valued
# input observable observables to state names. At construction time a 'default'
# state name must be specified to which all values of x are mapped that are not
# otherwise assigned
xRegion = ROOT.RooThresholdCategory(
    "xRegion", "region of x", x, "Background")

# Specify thresholds and state assignments one-by-one.
# Each statement specifies that all values _below_ the given value
# (and above any lower specified threshold) are mapped to the
# category state with the given name
#
# Background | SideBand | Signal | SideBand | Background
#           4.23       5.23     8.23       9.23
xRegion.addThreshold(4.23, "Background")
xRegion.addThreshold(5.23, "SideBand")
xRegion.addThreshold(8.23, "Signal")
xRegion.addThreshold(9.23, "SideBand")

# Use threshold function to plot data regions
# ----------------------------------------------

# Add values of threshold function to dataset so that it can be used as
# observable
data.addColumn(xRegion)

# Make plot of data in x
xframe = x.frame(ROOT.RooFit.Title(
    "Demo of threshold and binning mapping functions"))
data.plotOn(xframe)

# Use calculated category to select sideband data
data.plotOn(
    xframe,
    ROOT.RooFit.Cut("xRegion==xRegion::SideBand"),
    ROOT.RooFit.MarkerColor(
        ROOT.kRed),
    ROOT.RooFit.LineColor(
        ROOT.kRed))

# Create a binning real -> cat function
# ----------------------------------------------------------------------

# A RooBinningCategory is a category function that maps bins of a (named) binning definition
# in a real-valued input observable observables to state names. The state names are automatically
# constructed from the variable name, binning name and the bin number. If no binning name
# is specified the default binning is mapped

x.setBins(10, "coarse")
xBins = ROOT.RooBinningCategory("xBins", "coarse bins in x", x, "coarse")

# Use binning function for tabulation and plotting
# -----------------------------------------------------------------------------------------------

# Print table of xBins state multiplicity. Note that xBins does not need to be an observable in data
# it can be a function of observables in data as well
xbtable = data.table(xBins)
xbtable.Print("v")

# Add values of xBins function to dataset so that it can be used as
# observable
xb = data.addColumn(xBins)

# Define range "alt" as including bins 1,3,5,7,9
xb.setRange(
    "alt",
    "x_coarse_bin1,x_coarse_bin3,x_coarse_bin5,x_coarse_bin7,x_coarse_bin9")

# Construct subset of data matching range "alt" but only for the first
# 5000 events and plot it on the frame
dataSel = data.reduce(ROOT.RooFit.CutRange(
    "alt"), ROOT.RooFit.EventRange(0, 5000))
dataSel.plotOn(xframe, ROOT.RooFit.MarkerColor(ROOT.kGreen),
               ROOT.RooFit.LineColor(ROOT.kGreen))

c = ROOT.TCanvas("rf405_realtocatfuncs", "rf405_realtocatfuncs", 600, 600)
xframe.SetMinimum(0.01)
ROOT.gPad.SetLeftMargin(0.15)
xframe.GetYaxis().SetTitleOffset(1.4)
xframe.Draw()

c.SaveAs("rf405_realtocatfuncs.png")
