/// \file
/// \ingroup tutorial_multicore
/// \notebook
/// Demonstrate how to activate and use the implicit parallelisation of TTree::GetEntry.
/// Such parallelisation creates one task per top-level branch of the tree being read.
/// In this example, most of the branches are floating point numbers, which are very fast to read.
///  This parallelisation can be used, though, on bigger trees with many (complex) branches, which
///  are more likely to benefit from speedup gains.
///
/// \macro_code
///
/// \date 26/09/2016
/// \author Enric Tejedor

int imt001_parBranchProcessing()
{
   // First enable implicit multi-threading globally, so that the implicit parallelisation is on.
   // The parameter of the call specifies the number of threads to use.
   int nthreads = 4;
   ROOT::EnableImplicitMT(nthreads);

   // Open the file containing the tree
   auto file = TFile::Open("http://root.cern/files/h1/dstarmb.root");

   // Get the tree
   auto tree = file->Get<TTree>("h42");

   const auto nEntries = tree->GetEntries();

   // Read the branches in parallel.
   // Note that the interface does not change, the parallelisation is internal
   for (auto i : ROOT::TSeqUL(nEntries)) {
      tree->GetEntry(i); // parallel read
   }

   // IMT parallelisation can be disabled for a specific tree
   tree->SetImplicitMT(false);

   // If now GetEntry is invoked on the tree, the reading is sequential
   for (auto i : ROOT::TSeqUL(nEntries)) {
      tree->GetEntry(i); // sequential read
   }

   // Parallel reading can be re-enabled
   tree->SetImplicitMT(true);

   // IMT can be also disabled globally.
   // As a result, no tree will run GetEntry in parallel
   ROOT::DisableImplicitMT();

   // This is still sequential: the global flag is disabled, even if the
   // flag for this particular tree is enabled
   for (auto i : ROOT::TSeqUL(nEntries)) {
      tree->GetEntry(i); // sequential read
   }

   return 0;
}
