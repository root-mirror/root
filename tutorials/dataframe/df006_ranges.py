## \file
## \ingroup tutorial_dataframe
## \notebook -nodraw
## Use Range to limit the amount of data processed.
##
## This tutorial shows how to express the concept of ranges when working with the RDataFrame.
##
## \macro_code
## \macro_output
##
## \date March 2017
## \author Danilo Piparo (CERN)

import ROOT

def fill_tree(treeName, fileName):
    df = ROOT.RDataFrame(100)
    df.Define("b1", "(int) rdfentry_")\
      .Define("b2", "(float) rdfentry_ * rdfentry_").Snapshot(treeName, fileName)


# We prepare an input tree to run on
fileName = "df006_ranges_py.root"
treeName = "myTree"

fill_tree(treeName, fileName)

# We read the tree from the file and create a RDataFrame.
d = ROOT.RDataFrame(treeName, fileName)

# ## Usage of ranges
# Now we'll count some entries using ranges
c_all = d.Count()

# This is how you can express a range of the first 30 entries
d_0_30 = d.Range(0, 30)
c_0_30 = d_0_30.Count()

# This is how you pick all entries from 15 onwards
d_15_end = d.Range(15, 0)
c_15_end = d_15_end.Count()

# We can use a stride too, in this case we pick an event every 3
d_15_end_3 = d.Range(15, 0, 3)
c_15_end_3 = d_15_end_3.Count()

# The Range is a 1st class citizen in the RDataFrame graph:
# not only actions (like Count) but also filters and new columns can be added to it.
d_0_50 = d.Range(0, 50)
c_0_50_odd_b1 = d_0_50.Filter("1 == b1 % 2").Count()

# An important thing to notice is that the counts of a filter are relative to the
# number of entries a filter "sees". Therefore, if a Range depends on a filter,
# the Range will act on the entries passing the filter only.
c_0_3_after_even_b1 = d.Filter("0 == b1 % 2").Range(0, 3).Count()

# Ok, time to wrap up: let's print all counts!
print("Usage of ranges:")
print(" - All entries:", c_all.GetValue())
print(" - Entries from 0 to 30:", c_0_30.GetValue())
print(" - Entries from 15 onwards:", c_15_end.GetValue())
print(" - Entries from 15 onwards in steps of 3:", c_15_end_3.GetValue())
print(" - Entries from 0 to 50, odd only:", c_0_50_odd_b1.GetValue())
print(" - First three entries of all even entries:", c_0_3_after_even_b1.GetValue())
