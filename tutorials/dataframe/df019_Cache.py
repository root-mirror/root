## \file
## \ingroup tutorial_dataframe
## \notebook -draw
## Cache a processed RDataFrame in memory for further usage.
##
## This tutorial shows how the content of a data frame can be cached in memory
## in form of a data frame. The content of the columns is stored in memory in
## contiguous slabs of memory and is "ready to use", i.e. no ROOT IO operation
## is performed.
##
## Creating a cached data frame storing all of its content deserialised and uncompressed
## in memory is particularly useful when dealing with datasets of a moderate size
## (small enough to fit the RAM) over which several explorative loops need to be
## performed at as fast as possible. In addition, caching can be useful when no file
## on disk needs to be created as a side effect of checkpointing part of the analysis.
##
## All steps in the caching are lazy, i.e. the cached data frame is actually filled
## only when the event loop is triggered on it.
##
## \macro_code
## \macro_image
##
## \date June 2018
## \author Danilo Piparo (CERN)

import ROOT
import os

# We create a data frame on top of the hsimple example
hsimplePath = os.path.join(str(ROOT.gROOT.GetTutorialDir().Data()), "hsimple.root")
df = ROOT.RDataFrame("ntuple", hsimplePath)

# We apply a simple cut and define a new column
df_cut = df.Filter("py > 0.f")\
           .Define("px_plus_py", "px + py")

# We cache the content of the dataset. Nothing has happened yet: the work to accomplish
# has been described.
df_cached = df_cut.Cache()

h = df_cached.Histo1D("px_plus_py")

# Now the event loop on the cached dataset is triggered by accessing the histogram.
# This event triggers the loop on the `df` data frame lazily.
c = ROOT.TCanvas()
h.Draw()
c.SaveAs("df019_Cache.png")

print("Saved figure to df019_Cache.png")
