## \file
## \ingroup tutorial_roofit
## \notebook -nodraw
## Using the RooCustomizer to create multiple PDFs that share a lot of properties, but have unique parameters for each category.
## As an extra complication, some of the new parameters need to be functions of a mass parameter.
##
## \macro_code
## \macro_output
##
## \date June 2021
## \author Harshal Shende, Stephan Hageboeck (C++ version)

import ROOT

E = ROOT.RooRealVar("Energy", "Energy", 0, 3000)

meanG = ROOT.RooRealVar("meanG", "meanG", 100.0, 0.0, 3000.0)
sigmaG = ROOT.RooRealVar("sigmaG", "sigmaG", 3.0)
gauss = ROOT.RooGaussian("gauss", "gauss", E, meanG, sigmaG)

pol1 = ROOT.RooRealVar("pol1", "Constant of the polynomial", 1, -10, 10)
linear = ROOT.RooPolynomial("linear", "linear", E, pol1)

yieldSig = ROOT.RooRealVar("yieldSig", "yieldSig", 1, 0, 1.0e4)
yieldBkg = ROOT.RooRealVar("yieldBkg", "yieldBkg", 1, 0, 1.0e4)

model = ROOT.RooAddPdf("model", "S + B model", ROOT.RooArgList(gauss, linear), ROOT.RooArgList(yieldSig, yieldBkg))

print("The proto model before customisation:\n")
model.Print("T")  # "T" prints the model as a tree


# Build the categories
sample = ROOT.RooCategory("sample", "sample", {"Sample1": 1, "Sample2": 2, "Sample3": 3})


# Start to customise the proto model that was defined above.
# ---------------------------------------------------------------------------

# We need two sets for bookkeeping of PDF nodes:
newLeafs = ROOT.RooArgSet()
allCustomiserNodes = ROOT.RooArgSet()


# 1. Each sample should have its own mean for the gaussian
# The customiser will make copies of `meanG` for each category.
# These will all appear in the set `newLeafs`, which will own the new nodes.
cust = ROOT.RooCustomizer(model, sample, newLeafs, allCustomiserNodes)
cust.splitArg(meanG, sample)


# 2. Each sample should have its own signal yield, but there is an extra complication:
# We need the yields 1 and 2 to be a function of the variable "mass".
# For this, we pre-define nodes with exacly the names that the customiser would have created automatically,
# that is, "<nodeName>_<categoryName>", and we register them in the set of customiser nodes.
# The customiser will pick them up instead of creating new ones.
# If we don't provide one (e.g. for "yieldSig_Sample3"), it will be created automatically by cloning `yieldSig`.
mass = ROOT.RooRealVar("M", "M", 1, 0, 12000)
yield1 = ROOT.RooFormulaVar("yieldSig_Sample1", "Signal yield in the first sample", "M/3.360779", mass)
yield2 = ROOT.RooFormulaVar("yieldSig_Sample2", "Signal yield in the second sample", "M/2", mass)
allCustomiserNodes.add(yield1)
allCustomiserNodes.add(yield2)

# Instruct the customiser to replace all yieldSig nodes for each sample:
cust.splitArg(yieldSig, sample)


# Now we can start building the PDFs for all categories:
pdf1 = cust.build("Sample1")
pdf2 = cust.build("Sample2")
pdf3 = cust.build("Sample3")

# And we inspect the two PDFs
print("\nPDF 1 with a yield depending on M:\n")
pdf1.Print("T")
print("\nPDF 2 with a yield depending on M:\n")
pdf2.Print("T")
print("\nPDF 3 with a free yield:\n")
pdf3.Print("T")

print("\nThe following leafs have been created automatically while customising:\n")
newLeafs.Print("V")

#  If we needed to set reasonable values for the means of the gaussians, this could be done as follows:
meanG1 = allCustomiserNodes["meanG_Sample1"]
meanG1.setVal(200)
meanG2 = allCustomiserNodes["meanG_Sample2"]
meanG2.setVal(300)

print(
    "\nThe following leafs have been used while customising\n\t(partial overlap with the set of automatically created leaves.\n\ta new customiser for a different PDF could reuse them if necessary.):"
)
allCustomiserNodes.Print("V")

del cust
