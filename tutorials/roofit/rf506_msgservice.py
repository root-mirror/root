## \file
## \ingroup tutorial_roofit
## \notebook -nodraw
## Organization and simultaneous fits: tuning and customizing the ROOT.RooFit message logging facility
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

import ROOT

# Create pdf
# --------------------

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

data = model.generate(ROOT.RooArgSet(x), 10)

# Print configuration of message service
# ------------------------------------------

# Print streams configuration
ROOT.RooMsgService.instance().Print()

# Adding integration topic to existing INFO stream
# ---------------------------------------------------

# Print streams configuration
ROOT.RooMsgService.instance().Print()

# Add Integration topic to existing INFO stream
ROOT.RooMsgService.instance().getStream(1).addTopic(ROOT.RooFit.Integration)

# Construct integral over gauss to demonstrate message stream
igauss = gauss.createIntegral(ROOT.RooArgSet(x))
igauss.Print()

# Print streams configuration in verbose, also shows inactive streams
ROOT.RooMsgService.instance().Print()

# Remove stream
ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Integration)

# Examples of pdf value tracing
# -----------------------------------------------------------------------

# Show DEBUG level message on function tracing, ROOT.RooGaussian only
ROOT.RooMsgService.instance().addStream(
    ROOT.RooFit.DEBUG,
    ROOT.RooFit.Topic(
        ROOT.RooFit.Tracing),
    ROOT.RooFit.ClassName("RooGaussian"))

# Perform a fit to generate some tracing messages
model.fitTo(data, ROOT.RooFit.Verbose(ROOT.kTRUE))

# Reset message service to default stream configuration
ROOT.RooMsgService.instance().reset()

# Show DEBUG level message on function tracing on all objects, output to
# file
ROOT.RooMsgService.instance().addStream(
    ROOT.RooFit.DEBUG,
    ROOT.RooFit.Topic(
        ROOT.RooFit.Tracing),
    ROOT.RooFit.OutputFile("rf506_debug.log"))

# Perform a fit to generate some tracing messages
model.fitTo(data, ROOT.RooFit.Verbose(ROOT.kTRUE))

# Reset message service to default stream configuration
ROOT.RooMsgService.instance().reset()

# Example of another debugging stream
# ---------------------------------------------------------------------

# Show DEBUG level messages on client/server link state management
ROOT.RooMsgService.instance().addStream(
    ROOT.RooFit.DEBUG, ROOT.RooFit.Topic(ROOT.RooFit.LinkStateMgmt))
ROOT.RooMsgService.instance().Print("v")

# Clone composite pdf g to trigger some link state management activity
gprime = gauss.cloneTree()
gprime.Print()

# Reset message service to default stream configuration
ROOT.RooMsgService.instance().reset()
