## \file
## \ingroup tutorial_roofit
## \notebook -nodraw
## Data and categories: working with ROOT.RooCategory objects to describe discrete variables
##
## \macro_code
##
## \date February 2018
## \authors Clemens Lange, Wouter Verkerke (C++ version)

from __future__ import print_function
import ROOT


# Construct a category with labels
# --------------------------------------------

# Define a category with labels only
tagCat = ROOT.RooCategory("tagCat", "Tagging category")
tagCat.defineType("Lepton")
tagCat.defineType("Kaon")
tagCat.defineType("NetTagger-1")
tagCat.defineType("NetTagger-2")
tagCat.Print()

# Construct a category with labels and indices
# ------------------------------------------------

# Define a category with explicitly numbered states
b0flav = ROOT.RooCategory("b0flav", "B0 flavour eigenstate", {"B0": -1, "B0bar": 1})
b0flav.Print()

# Generate dummy data for tabulation demo
# ------------------------------------------------

# Generate a dummy dataset
x = ROOT.RooRealVar("x", "x", 0, 10)
data = ROOT.RooPolynomial("p", "p", x).generate(ROOT.RooArgSet(x, b0flav, tagCat), 10000)

# Print tables of category contents of datasets
# --------------------------------------------------

# Tables are equivalent of plots for categories
btable = data.table(b0flav)
btable.Print()
btable.Print("v")

# Create table for subset of events matching cut expression
ttable = data.table(tagCat, "x>8.23")
ttable.Print()
ttable.Print("v")

# Create table for all (tagCat x b0flav) state combinations
bttable = data.table(ROOT.RooArgSet(tagCat, b0flav))
bttable.Print("v")

# Retrieve number of events from table
# Number can be non-integer if source dataset has weighed events
nb0 = btable.get("B0")
print("Number of events with B0 flavor is ", nb0)

# Retrieve fraction of events with "Lepton" tag
fracLep = ttable.getFrac("Lepton")
print("Fraction of events tagged with Lepton tag is ", fracLep)

# Defining ranges for plotting, fitting on categories
# ------------------------------------------------------------------------------------------------------

# Define named range as comma separated list of labels
tagCat.setRange("good", "Lepton,Kaon")

# Or add state names one by one
tagCat.addToRange("soso", "NetTagger-1")
tagCat.addToRange("soso", "NetTagger-2")

# Use category range in dataset reduction specification
goodData = data.reduce(CutRange="good")
goodData.table(tagCat).Print("v")
