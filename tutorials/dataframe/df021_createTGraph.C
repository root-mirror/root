/// \file
/// \ingroup tutorial_dataframe
/// \notebook -draw
/// Fill a TGraph using RDataFrame.
///
/// \macro_code
/// \macro_image
///
/// \date July 2018
/// \authors Enrico Guiraud, Danilo Piparo (CERN), Massimo Tumolo (Politecnico di Torino)



void df021_createTGraph()
{
   ROOT::EnableImplicitMT(2);

   ROOT::RDataFrame d(160);

   // Create a trivial parabola
   auto dd = d.Alias("x", "rdfentry_").Define("y", "x*x");

   auto graph = dd.Graph("x", "y");

   // This tutorial is ran with multithreading enabled. The order in which points are inserted is not known, so to have a meaningful representation points are sorted.
   graph->Sort();
   auto c = new TCanvas();
   graph->DrawClone("APL");
}
