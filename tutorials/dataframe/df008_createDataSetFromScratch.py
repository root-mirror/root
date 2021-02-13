## \file
## \ingroup tutorial_dataframe
## \notebook -nodraw
## Create data from scratch with RDataFrame.
##
## This tutorial shows how to create a dataset from scratch with RDataFrame
##
## \macro_code
##
## \date June 2017
## \author Danilo Piparo (CERN)

import ROOT

# We create an empty data frame of 100 entries
df = ROOT.RDataFrame(100)

# We now fill it with random numbers
ROOT.gRandom.SetSeed(1)
df_1 = df.Define("rnd", "gRandom->Gaus()")

# And we write out the dataset on disk
df_1.Snapshot("randomNumbers", "df008_createDataSetFromScratch_py.root")
