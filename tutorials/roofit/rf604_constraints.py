## \file
## \ingroup tutorial_roofit
## \notebook -nodraw
## Likelihood and minimization: fitting with constraints
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

from __future__ import print_function
import ROOT


# Create model and dataset
# ----------------------------------------------

# Construct a Gaussian pdf
x = ROOT.RooRealVar("x", "x", -10, 10)

m = ROOT.RooRealVar("m", "m", 0, -10, 10)
s = ROOT.RooRealVar("s", "s", 2, 0.1, 10)
gauss = ROOT.RooGaussian("gauss", "gauss(x,m,s)", x, m, s)

# Construct a flat pdf (polynomial of 0th order)
poly = ROOT.RooPolynomial("poly", "poly(x)", x)

# model = f*gauss + (1-f)*poly
f = ROOT.RooRealVar("f", "f", 0.5, 0., 1.)
model = ROOT.RooAddPdf(
    "model",
    "model",
    ROOT.RooArgList(
        gauss,
        poly),
    ROOT.RooArgList(f))

# Generate small dataset for use in fitting below
d = model.generate(ROOT.RooArgSet(x), 50)

# Create constraint pdf
# -----------------------------------------

# Construct Gaussian constraint pdf on parameter f at 0.8 with
# resolution of 0.1
fconstraint = ROOT.RooGaussian(
    "fconstraint",
    "fconstraint",
    f,
    ROOT.RooFit.RooConst(0.8),
    ROOT.RooFit.RooConst(0.1))

# Method 1 - add internal constraint to model
# -------------------------------------------------------------------------------------

# Multiply constraint term with regular pdf using ROOT.RooProdPdf
# Specify in fitTo() that internal constraints on parameter f should be
# used

# Multiply constraint with pdf
modelc = ROOT.RooProdPdf(
    "modelc", "model with constraint", ROOT.RooArgList(model, fconstraint))

# Fit model (without use of constraint term)
r1 = model.fitTo(d, ROOT.RooFit.Save())

# Fit modelc with constraint term on parameter f
r2 = modelc.fitTo(
    d,
    ROOT.RooFit.Constrain(
        ROOT.RooArgSet(f)),
    ROOT.RooFit.Save())

# Method 2 - specify external constraint when fitting
# ------------------------------------------------------------------------------------------

# Construct another Gaussian constraint pdf on parameter f at 0.8 with
# resolution of 0.1
fconstext = ROOT.RooGaussian("fconstext", "fconstext", f, ROOT.RooFit.RooConst(
    0.2), ROOT.RooFit.RooConst(0.1))

# Fit with external constraint
r3 = model.fitTo(d, ROOT.RooFit.ExternalConstraints(
    ROOT.RooArgSet(fconstext)), ROOT.RooFit.Save())

# Print the fit results
print("fit result without constraint (data generated at f=0.5)")
r1.Print("v")
print("fit result with internal constraint (data generated at f=0.5, is f=0.8+/-0.2)")
r2.Print("v")
print("fit result with (another) external constraint (data generated at f=0.5, is f=0.2+/-0.1)")
r3.Print("v")
