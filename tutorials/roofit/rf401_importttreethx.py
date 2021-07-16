## \file
## \ingroup tutorial_roofit
## \notebook
##
## 'DATA AND CATEGORIES' RooFit tutorial macro #401
##
## Overview of advanced option for importing data from ROOT ROOT.TTree and ROOT.THx histograms
## Basic import options are demonstrated in rf102_dataimport.py
##
## \macro_code
##
## \date February 2018
## \author Clemens Lange
## \author Wouter Verkerke (C version)

import ROOT
from array import array


def makeTH1(name, mean, sigma):
    """Create ROOT TH1 filled with a Gaussian distribution."""

    hh = ROOT.TH1D(name, name, 100, -10, 10)
    for i in range(1000):
        hh.Fill(ROOT.gRandom.Gaus(mean, sigma))

    return hh


def makeTTree():
    """Create ROOT ROOT.TTree filled with a Gaussian distribution in x and a uniform distribution in y."""

    tree = ROOT.TTree("tree", "tree")
    px = array("d", [0])
    py = array("d", [0])
    pz = array("d", [0])
    pi = array("i", [0])
    tree.Branch("x", px, "x/D")
    tree.Branch("y", py, "y/D")
    tree.Branch("z", pz, "z/D")
    tree.Branch("i", pi, "i/I")
    for i in range(100):
        px[0] = ROOT.gRandom.Gaus(0, 3)
        py[0] = ROOT.gRandom.Uniform() * 30 - 15
        pz[0] = ROOT.gRandom.Gaus(0, 5)
        pi[0] = i % 3
        tree.Fill()

    return tree


# Import multiple TH1 into a RooDataHist
# ----------------------------------------------------------

# Create thee ROOT ROOT.TH1 histograms
hh_1 = makeTH1("hh1", 0, 3)
hh_2 = makeTH1("hh2", -3, 1)
hh_3 = makeTH1("hh3", +3, 4)

# Declare observable x
x = ROOT.RooRealVar("x", "x", -10, 10)

# Create category observable c that serves as index for the ROOT histograms
c = ROOT.RooCategory("c", "c")
c.defineType("SampleA")
c.defineType("SampleB")
c.defineType("SampleC")

# Create a binned dataset that imports contents of all ROOT.TH1 mapped by
# index category c
dh = ROOT.RooDataHist(
    "dh", "dh", ROOT.RooArgList(x), Index=c, Import={"SampleA": hh_1, "SampleB": hh_2, "SampleC": hh_3}
)
dh.Print()

dh2 = ROOT.RooDataHist(
    "dh", "dh", ROOT.RooArgList(x), Index=c, Import={"SampleA": hh_1, "SampleB": hh_2, "SampleC": hh_3}
)
dh2.Print()

# Importing a ROOT TTree into a RooDataSet with cuts
# --------------------------------------------------------------------

tree = makeTTree()

# Define observables y,z
y = ROOT.RooRealVar("y", "y", -10, 10)
z = ROOT.RooRealVar("z", "z", -10, 10)

# Import only observables (y,z)
ds = ROOT.RooDataSet("ds", "ds", ROOT.RooArgSet(x, y), Import=tree)
ds.Print()

# Import observables (x,y,z) but only event for which (y+z<0) is ROOT.True
# Import observables (x,y,z) but only event for which (y+z<0) is ROOT.True
ds2 = ROOT.RooDataSet("ds2", "ds2", ROOT.RooArgSet(x, y, z), Import=tree, Cut="y+z<0")
ds2.Print()

# Importing integer ROOT TTree branches
# ---------------------------------------------------------------

# Import integer tree branch as ROOT.RooRealVar
i = ROOT.RooRealVar("i", "i", 0, 5)
ds3 = ROOT.RooDataSet("ds3", "ds3", ROOT.RooArgSet(i, x), Import=tree)
ds3.Print()

# Define category i
icat = ROOT.RooCategory("i", "i", {"State0": 0, "State1": 1})

# Import integer tree branch as ROOT.RooCategory (only events with i==0 and i==1
# will be imported as those are the only defined states)
ds4 = ROOT.RooDataSet("ds4", "ds4", ROOT.RooArgSet(icat, x), Import=tree)
ds4.Print()

# Import multiple RooDataSets into a RooDataSet
# ----------------------------------------------------------------------------------------

# Create three ROOT.RooDataSets in (y,z)
dsA = ds2.reduce(ROOT.RooArgSet(x, y), "z<-5")
dsB = ds2.reduce(ROOT.RooArgSet(x, y), "abs(z)<5")
dsC = ds2.reduce(ROOT.RooArgSet(x, y), "z>5")

# Create a dataset that imports contents of all the above datasets mapped
# by index category c
dsABC = ROOT.RooDataSet(
    "dsABC",
    "dsABC",
    ROOT.RooArgSet(x, y),
    ROOT.RooFit.Import("SampleA", dsA),
    ROOT.RooFit.Import("SampleB", dsB),
    Index=c,
    Import=("SampleC", dsC),
)

dsABC.Print()
