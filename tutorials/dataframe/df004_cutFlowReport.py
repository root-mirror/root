## \file
## \ingroup tutorial_dataframe
## \notebook
## Display cut/Filter efficiencies with RDataFrame.
##
## This tutorial shows how to get information about the efficiency of the filters
## applied
##
## \macro_code
## \macro_output
##
## \date May 2017
## \author Danilo Piparo (CERN)

import ROOT

def fill_tree(treeName, fileName):
    df = ROOT.RDataFrame(50)
    df.Define("b1", "(double) rdfentry_")\
      .Define("b2", "(int) rdfentry_ * rdfentry_").Snapshot(treeName, fileName)

# We prepare an input tree to run on
fileName = 'df004_cutFlowReport_py.root'
treeName = 'myTree'
fill_tree(treeName, fileName)

# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, fileName)

# ## Define cuts and create the report
# An optional string parameter name can be passed to the Filter method to create a named filter.
# Named filters work as usual, but also keep track of how many entries they accept and reject.
filtered1 = d.Filter('b1 > 25', 'Cut1')
filtered2 = d.Filter('0 == b2 % 2', 'Cut2')

augmented1 = filtered2.Define('b3', 'b1 / b2')
filtered3 = augmented1.Filter('b3 < .5','Cut3')

# Statistics are retrieved through a call to the Report method:
# when Report is called on the main RDataFrame object, it retrieves stats for
# all named filters declared up to that point. When called on a stored chain
# state (i.e. a chain/graph node), it retrieves stats for all named filters in
# the section of the chain between the main RDataFrame and that node (included).
# Stats are printed in the same order as named filters have been added to the
# graph, and refer to the latest event-loop that has been run using the relevant
# RDataFrame.
print('Cut3 stats:')
filtered3.Report()
print('All stats:')
allCutsReport = d.Report()
allCutsReport.Print()
