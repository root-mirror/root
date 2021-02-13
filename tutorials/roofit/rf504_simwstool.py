## \file
## \ingroup tutorial_roofit
## \notebook -nodraw
## Organization and simultaneous fits: using RooSimWSTool to construct a simultaneous pdf
## that is built of variations of an input pdf
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT


# Create master pdf
# ---------------------------------

# Construct gauss(x,m,s)
x = ROOT.RooRealVar("x", "x", -10, 10)
m = ROOT.RooRealVar("m", "m", 0, -10, 10)
s = ROOT.RooRealVar("s", "s", 1, -10, 10)
gauss = ROOT.RooGaussian("g", "g", x, m, s)

# Construct poly(x,p0)
p0 = ROOT.RooRealVar("p0", "p0", 0.01, 0., 1.)
poly = ROOT.RooPolynomial("p", "p", x, ROOT.RooArgList(p0))

# model = f*gauss(x) + (1-f)*poly(x)
f = ROOT.RooRealVar("f", "f", 0.5, 0., 1.)
model = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(
    gauss, poly), ROOT.RooArgList(f))

# Create category observables for splitting
# ----------------------------------------------------------------------------------

# Define two categories that can be used for splitting
c = ROOT.RooCategory("c", "c")
c.defineType("run1")
c.defineType("run2")

d = ROOT.RooCategory("d", "d")
d.defineType("foo")
d.defineType("bar")

# Set up SimWSTool
# -----------------------------

# Import ingredients in a workspace
w = ROOT.RooWorkspace("w", "w")
w.Import(ROOT.RooArgSet(model, c, d))

# Make Sim builder tool
sct = ROOT.RooSimWSTool(w)

# Build a simultaneous model with one split
# ---------------------------------------------------------------------------------

# Construct a simultaneous pdf with the following form
#
# model_run1(x) = f*gauss_run1(x,m_run1,s) + (1-f)*poly
# model_run2(x) = f*gauss_run2(x,m_run2,s) + (1-f)*poly
# simpdf(x,c) = model_run1(x) if c=="run1"
#             = model_run2(x) if c=="run2"
#
# Returned pdf is owned by the workspace
model_sim = sct.build("model_sim", "model",
                      ROOT.RooFit.SplitParam("m", "c"))

# Print tree structure of model
model_sim.Print("t")

# Adjust model_sim parameters in workspace
w.var("m_run1").setVal(-3)
w.var("m_run2").setVal(+3)

# Print contents of workspace
w.Print("v")

# Build a simultaneous model with product split
# -----------------------------------------------------------------------------------------

# Build another simultaneous pdf using a composite split in states c X d
model_sim2 = sct.build("model_sim2", "model",
                       ROOT.RooFit.SplitParam("p0", "c,d"))

# Print tree structure of self model
model_sim2.Print("t")
