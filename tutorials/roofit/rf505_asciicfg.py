## \file rf505_asciicfg.py
## \ingroup tutorial_roofit
## \notebook -nodraw
## Organization and simultaneous fits: reading and writing ASCII configuration files
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

from __future__ import print_function
import ROOT


# Create pdf
# ------------------

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

# Fit model to toy data
# -----------------------------------------

d = model.generate(ROOT.RooArgSet(x), 1000)
model.fitTo(d)

# Write parameters to ASCII file
# -----------------------------------------------------------

# Obtain set of parameters
params = model.getParameters(ROOT.RooArgSet(x))

# Write parameters to file
params.writeToFile("rf505_asciicfg_example.txt")

# Read parameters from ASCII file
# ----------------------------------------------------------------

# Read parameters from file
params.readFromFile("rf505_asciicfg_example.txt")
params.Print("v")

configFile = ROOT.gROOT.GetTutorialDir().Data() + "/roofit/rf505_asciicfg.txt"

# Read parameters from section 'Section2' of file
params.readFromFile(configFile, "", "Section2")
params.Print("v")

# Read parameters from section 'Section3' of file. Mark all
# variables that were processed with the "READ" attribute
params.readFromFile(configFile, "READ", "Section3")

# Print the list of parameters that were not read from Section3
print("The following parameters of the were _not_ read from Section3: ",
      params.selectByAttrib("READ", ROOT.kFALSE))

# Read parameters from section 'Section4' of file, contains
# 'include file' statement of rf505_asciicfg_example.txt
# so that we effective read the same
params.readFromFile(configFile, "", "Section4")
params.Print("v")
