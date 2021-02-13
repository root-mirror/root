## \file
## \ingroup tutorial_roofit
## \notebook -nodraw
## Numeric algorithm tuning: configuration and customization of how numeric (partial) integrals are executed
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

from __future__ import print_function
import ROOT


# Adjust global 1D integration precision
# ----------------------------------------------------------------------------

# Print current global default configuration for numeric integration
# strategies
ROOT.RooAbsReal.defaultIntegratorConfig().Print("v")

# Example: Change global precision for 1D integrals from 1e-7 to 1e-6
#
# The relative epsilon (change as fraction of current best integral estimate) and
# absolute epsilon (absolute change w.r.t last best integral estimate) can be specified
# separately. For most pdf integrals the relative change criterium is the most important,
# however for certain non-pdf functions that integrate out to zero a separate absolute
# change criterium is necessary to declare convergence of the integral
#
# NB: ROOT.This change is for illustration only. In general the precision should be at least 1e-7
# for normalization integrals for MINUIT to succeed.
#
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-6)
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-6)

# N u m e r i c   i n t e g r a t i o n   o f   l a n d a u   p d f
# ------------------------------------------------------------------

# Construct pdf without support for analytical integrator for
# demonstration purposes
x = ROOT.RooRealVar("x", "x", -10, 10)
landau = ROOT.RooLandau("landau", "landau", x,
                        ROOT.RooFit.RooConst(0), ROOT.RooFit.RooConst(0.1))

# Activate debug-level messages for topic integration to be able to follow
# actions below
ROOT.RooMsgService.instance().addStream(
    ROOT.RooFit.DEBUG, ROOT.RooFit.Topic(ROOT.RooFit.Integration))

# Calculate integral over landau with default choice of numeric integrator
intLandau = landau.createIntegral(ROOT.RooArgSet(x))
val = intLandau.getVal()
print(" [1] int_dx landau(x) = ", val)  # setprecision(15)

# Same with custom configuration
# -----------------------------------------------------------

# Construct a custom configuration which uses the adaptive Gauss-Kronrod technique
# for closed 1D integrals
customConfig = ROOT.RooNumIntConfig(
    ROOT.RooAbsReal.defaultIntegratorConfig())
integratorGKNotExisting = customConfig.method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D")
if (integratorGKNotExisting) :
   print("WARNING: RooAdaptiveGaussKronrodIntegrator is not existing because ROOT is built without Mathmore support")

# Calculate integral over landau with custom integral specification
intLandau2 = landau.createIntegral(
    ROOT.RooArgSet(x), ROOT.RooFit.NumIntConfig(customConfig))
val2 = intLandau2.getVal()
print(" [2] int_dx landau(x) = ", val2)

# Adjusting default config for a specific pdf
# -------------------------------------------------------------------------------------

# Another possibility: associate custom numeric integration configuration
# as default for object 'landau'
landau.setIntegratorConfig(customConfig)

# Calculate integral over landau custom numeric integrator specified as
# object default
intLandau3 = landau.createIntegral(ROOT.RooArgSet(x))
val3 = intLandau3.getVal()
print(" [3] int_dx landau(x) = ", val3)

# Another possibility: Change global default for 1D numeric integration
# strategy on finite domains
if (not integratorGKNotExisting) :
   ROOT.RooAbsReal.defaultIntegratorConfig().method1D().setLabel(
       "RooAdaptiveGaussKronrodIntegrator1D")

# Adjusting parameters of a speficic technique
# ---------------------------------------------------------------------------------------

# Adjust maximum number of steps of ROOT.RooIntegrator1D in the global
# default configuration
   ROOT.RooAbsReal.defaultIntegratorConfig().getConfigSection(
       "RooIntegrator1D").setRealValue("maxSteps", 30)

# Example of how to change the parameters of a numeric integrator
# (Each config section is a ROOT.RooArgSet with ROOT.RooRealVars holding real-valued parameters
#  and ROOT.RooCategories holding parameters with a finite set of options)
   customConfig.getConfigSection(
       "RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 50)
   customConfig.getConfigSection(
       "RooAdaptiveGaussKronrodIntegrator1D").setCatLabel("method", "15Points")

# Example of how to print set of possible values for "method" category
   customConfig.getConfigSection(
       "RooAdaptiveGaussKronrodIntegrator1D").find("method").Print("v")
