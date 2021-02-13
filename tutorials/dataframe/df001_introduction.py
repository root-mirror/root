## \file
## \ingroup tutorial_dataframe
## \notebook -nodraw
## Basic usage of RDataFrame from python.
##
## This tutorial illustrates the basic features of the RDataFrame class,
## a utility which allows to interact with data stored in TTrees following
## a functional-chain like approach.
##
## \macro_code
## \macro_output
##
## \date May 2017
## \author Danilo Piparo (CERN)

import ROOT

# A simple helper function to fill a test tree: this makes the example stand-alone.
def fill_tree(treeName, fileName):
    df = ROOT.RDataFrame(10)
    df.Define("b1", "(double) rdfentry_")\
      .Define("b2", "(int) rdfentry_ * rdfentry_").Snapshot(treeName, fileName)

# We prepare an input tree to run on
fileName = "df001_introduction_py.root"
treeName = "myTree"
fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, fileName)

# Operations on the dataframe
# We now review some *actions* which can be performed on the data frame.
# All actions but ForEach return a TActionResultPtr<T>. The series of
# operations on the data frame is not executed until one of those pointers
# is accessed.
# But first of all, let us we define now our cut-flow with two strings.
# Filters can be expressed as strings. The content must be C++ code. The
# name of the variables must be the name of the branches. The code is
# just in time compiled.
cutb1 = 'b1 < 5.'
cutb1b2 = 'b2 % 2 && b1 < 4.'

# `Count` action
# The `Count` allows to retrieve the number of the entries that passed the
# filters. Here we show how the automatic selection of the column kicks
# in in case the user specifies none.
entries1 = d.Filter(cutb1) \
            .Filter(cutb1b2) \
            .Count();

print("%s entries passed all filters" %entries1.GetValue())

entries2 = d.Filter("b1 < 5.").Count();
print("%s entries passed all filters" %entries2.GetValue())

# `Min`, `Max` and `Mean` actions
# These actions allow to retrieve statistical information about the entries
# passing the cuts, if any.
b1b2_cut = d.Filter(cutb1b2)
minVal = b1b2_cut.Min('b1')
maxVal = b1b2_cut.Max('b1')
meanVal = b1b2_cut.Mean('b1')
nonDefmeanVal = b1b2_cut.Mean("b2")
print("The mean is always included between the min and the max: %s <= %s <= %s" %(minVal.GetValue(), meanVal.GetValue(), maxVal.GetValue()))

# `Histo1D` action
# The `Histo1D` action allows to fill an histogram. It returns a TH1F filled
# with values of the column that passed the filters. For the most common
# types, the type of the values stored in the column is automatically
# guessed.
hist = d.Filter(cutb1).Histo1D('b1')
print("Filled h %s times, mean: %s" %(hist.GetEntries(), hist.GetMean()))

# Express your chain of operations with clarity!
# We are discussing an example here but it is not hard to imagine much more
# complex pipelines of actions acting on data. Those might require code
# which is well organised, for example allowing to conditionally add filters
# or again to clearly separate filters and actions without the need of
# writing the entire pipeline on one line. This can be easily achieved.
# We'll show this re-working the `Count` example:
cutb1_result = d.Filter(cutb1);
cutb1b2_result = d.Filter(cutb1b2);
cutb1_cutb1b2_result = cutb1_result.Filter(cutb1b2)

# Now we want to count:
evts_cutb1_result = cutb1_result.Count()
evts_cutb1b2_result = cutb1b2_result.Count()
evts_cutb1_cutb1b2_result = cutb1_cutb1b2_result.Count()

print("Events passing cutb1: %s" %evts_cutb1_result.GetValue())
print("Events passing cutb1b2: %s" %evts_cutb1b2_result.GetValue())
print("Events passing both: %s" %evts_cutb1_cutb1b2_result.GetValue())

# Calculating quantities starting from existing columns
# Often, operations need to be carried out on quantities calculated starting
# from the ones present in the columns. We'll create in this example a third
# column the values of which are the sum of the *b1* and *b2* ones, entry by
# entry. The way in which the new quantity is defined is via a callable.
# It is important to note two aspects at this point:
# - The value is created on the fly only if the entry passed the existing
# filters.
# - The newly created column behaves as the one present on the file on disk.
# - The operation creates a new value, without modifying anything. De facto,
# this is like having a general container at disposal able to accommodate
# any value of any type.
# Let's dive in an example:
entries_sum = d.Define('sum', 'b2 + b1') \
               .Filter('sum > 4.2') \
               .Count()
print(entries_sum.GetValue())
