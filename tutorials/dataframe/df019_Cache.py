## \file
## \ingroup tutorial_dataframe
## \notebook -draw
## This tutorial shows how the content of a data frame can be cached in memory
## in form of a data frame. The content of the columns is stored in memory in
## contiguous slabs of memory and is "ready to use", i.e. no ROOT IO operation
## is performed.
## All steps in the caching are lazy, i.e. the cached data frame is actually filled
## only when the event loop is triggered on it.
##
## \macro_code
##
## \date June 2018
## \author Danilo Piparo

import ROOT
RDataFrame = ROOT.ROOT.RDataFrame
import os

# We create a data frame on top of the hsimple example
hsimplePath = os.path.join(str(ROOT.gROOT.GetTutorialDir().Data()), "hsimple.root")
df = RDataFrame("ntuple", hsimplePath)

#We apply a simple cut and define a new column
df_cut = df.Filter("py > 0.f")\
           .Define("px_plus_py", "px + py")

# We cache the content of the dataset. Nothing has happened yet: the work to accomplish
# has been described.
df_cached = df_cut.Cache()

h = df_cached.Histo1D("px_plus_py")

# Now the event loop on the cached dataset is triggered. This event triggers the loop
# on the `df` data frame lazily.
h.Draw()