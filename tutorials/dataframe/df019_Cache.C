/// \file
/// \ingroup tutorial_dataframe
/// \notebook -draw
/// Cache a processed RDataFrame in memory for further usage.
///
/// This tutorial shows how the content of a data frame can be cached in memory
/// in form of a data frame. The content of the columns is stored in memory in
/// contiguous slabs of memory and is "ready to use", i.e. no ROOT IO operation
/// is performed.
///
/// Creating a cached data frame storing all of its content deserialised and uncompressed
/// in memory is particularly useful when dealing with datasets of a moderate size
/// (small enough to fit the RAM) over which several explorative loops need to be
/// performed at as fast as possible. In addition, caching can be useful when no file
/// on disk needs to be created as a side effect of checkpointing part of the analysis.
///
/// All steps in the caching are lazy, i.e. the cached data frame is actually filled
/// only when the event loop is triggered on it.
///
/// \macro_code
/// \macro_image
///
/// \date June 2018
/// \author Danilo Piparo (CERN)

void df019_Cache()
{
   // We create a data frame on top of the hsimple example
   auto hsimplePath = gROOT->GetTutorialDir();
   hsimplePath += "/hsimple.root";
   ROOT::RDataFrame df("ntuple", hsimplePath.Data());

   // We apply a simple cut and define a new column
   auto df_cut = df.Filter([](float py) { return py > 0.f; }, {"py"})
                    .Define("px_plus_py", [](float px, float py) { return px + py; }, {"px", "py"});

   // We cache the content of the dataset. Nothing has happened yet: the work to accomplish
   // has been described. As for `Snapshot`, the types and columns can be written out explicitly
   // or left for the jitting to handle (`df_cached` is intentionally unused - it shows how to
   // to create a *cached* data frame specifying column types explicitly):
   auto df_cached = df_cut.Cache<float, float>({"px_plus_py", "py"});
   auto df_cached_implicit = df_cut.Cache();
   auto h = df_cached_implicit.Histo1D<float>("px_plus_py");

   // Now the event loop on the cached dataset is triggered. This event triggers the loop
   // on the `df` data frame lazily.
   h->DrawCopy();
}
