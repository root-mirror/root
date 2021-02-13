## \file
## \ingroup tutorial_roofit
## \notebook
## Numeric algorithm tuning: caching of slow numeric integrals and parameterizations of slow numeric integrals
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import sys
import ROOT


def getWorkspace(mode):
    # Create, save or load workspace with pdf
    # -----------------------------------------------------------------------------------
    #
    # Mode = 0 : Create workspace for plain running (no integral caching)
    # Mode = 1 : Generate workspace with precalculated integral and store it on file
    # Mode = 2 : Load previously stored workspace from file

    w = ROOT.RooWorkspace()

    if mode != 2:
        # Create empty workspace workspace
        w = ROOT.RooWorkspace("w", 1)

        # Make a difficult to normalize  pdf in 3 dimensions that is
        # integrated numerically.
        w.factory(
            "EXPR::model('1/((x-a)*(x-a)+0.01)+1/((y-a)*(y-a)+0.01)+1/((z-a)*(z-a)+0.01)',x[-1,1],y[-1,1],z[-1,1],a[-5,5])")

    if mode == 1:
        # Instruct model to precalculate normalization integral that integrate at least
        # two dimensions numerically. In self specific case the integral value for
        # all values of parameter 'a' are stored in a histogram and available for use
        # in subsequent fitting and plotting operations (interpolation is
        # applied)

        # w.pdf("model").setNormValueCaching(3)
        w.pdf("model").setStringAttribute("CACHEPARMINT", "x:y:z")

        # Evaluate pdf once to trigger filling of cache
        normSet = ROOT.RooArgSet(w.var("x"), w.var("y"), w.var("z"))
        w.pdf("model").getVal(normSet)
        w.writeToFile("rf903_numintcache.root")

    if (mode == 2):
        # Load preexisting workspace from file in mode==2
        f = ROOT.TFile("rf903_numintcache.root")
        w = f.Get("w")

    # Return created or loaded workspace
    return w


mode = 0
# Mode = 0 : Run plain fit (slow)
# Mode = 1 : Generate workspace with precalculated integral and store it on file (prepare for accelerated running)
# Mode = 2 : Run fit from previously stored workspace including cached
# integrals (fast, run in mode=1 first)

# Create, save or load workspace with pdf
# -----------------------------------------------------------------------------------

# Make/load workspace, here in mode 1
w = getWorkspace(mode)
if mode == 1:
    # Show workspace that was created
    w.Print()

    # Show plot of cached integral values
    hhcache = w.expensiveObjectCache().getObj(1)
    if (hhcache):
        ROOT.TCanvas("rf903_numintcache", "rf903_numintcache", 600, 600)
        hhcache.createHistogram("a").Draw()
    else:
        ROOT.RooFit.Error("rf903_numintcache",
                          "Cached histogram is not existing in workspace")
        sys.exit()

# Use pdf from workspace for generation and fitting
# -----------------------------------------------------------------------------------

# ROOT.This is always slow (need to find maximum function value
# empirically in 3D space)
d = w.pdf("model").generate(
    ROOT.RooArgSet(
        w.var("x"),
        w.var("y"),
        w.var("z")),
    1000)

# ROOT.This is slow in mode 0, fast in mode 1
w.pdf("model").fitTo(
    d, ROOT.RooFit.Verbose(
        ROOT.kTRUE), ROOT.RooFit.Timer(
            ROOT.kTRUE))

# Projection on x (always slow as 2D integral over Y, at fitted value of a
# is not cached)
framex = w.var("x").frame(ROOT.RooFit.Title("Projection of 3D model on X"))
d.plotOn(framex)
w.pdf("model").plotOn(framex)

# Draw x projection on canvas
c = ROOT.TCanvas("rf903_numintcache", "rf903_numintcache", 600, 600)
framex.Draw()

c.SaveAs("rf903_numintcache.png")

# Make workspace available on command line after macro finishes
ROOT.gDirectory.Add(w)
