## \file
## \ingroup tutorial_dataframe
## \notebook -draw
## Fill a TGraph using RDataFrame.
##
## \macro_code
## \macro_image
##
## \date July 2018
## \authors Enrico Guiraud, Danilo Piparo (CERN), Massimo Tumolo (Politecnico di Torino)

import ROOT

ROOT.ROOT.EnableImplicitMT(2)
d = ROOT.RDataFrame(160)

# Create a trivial parabola
dd = d.Alias("x", "rdfentry_").Define("y", "x*x")

graph = dd.Graph("x", "y")

# This tutorial is ran with multithreading enabled. The order in which points are inserted is not known, so to have a meaningful representation points are sorted.
c = ROOT.TCanvas()
graph.Sort()
graph.Draw("APL")
c.SaveAs("df021_createTGraph.png")

print("Saved figure to df021_createTGraph.png")
