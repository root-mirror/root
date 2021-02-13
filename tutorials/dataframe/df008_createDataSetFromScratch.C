/// \file
/// \ingroup tutorial_dataframe
/// \notebook -nodraw
/// Create data from scratch with RDataFrame.
///
/// This tutorial shows how to create a dataset from scratch with RDataFrame
///
/// \macro_code
///
/// \date June 2017
/// \author Danilo Piparo (CERN)

void df008_createDataSetFromScratch()
{
   // We create an empty data frame of 100 entries
   ROOT::RDataFrame df(100);

   // We now fill it with random numbers
   gRandom->SetSeed(1);
   auto df_1 = df.Define("rnd", []() { return gRandom->Gaus(); });

   // And we write out the dataset on disk
   df_1.Snapshot("randomNumbers", "df008_createDataSetFromScratch.root");
}
