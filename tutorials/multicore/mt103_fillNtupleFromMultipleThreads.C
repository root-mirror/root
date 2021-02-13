/// \file
/// \ingroup tutorial_multicore
/// \notebook
/// Fill the same TNtuple from different threads.
/// This tutorial illustrates the basics of how it's possible with ROOT
/// to write simultaneously to a single output file using TBufferMerger.
///
/// \macro_code
///
/// \date May 2017
/// \author Guilherme Amadio

void mt103_fillNtupleFromMultipleThreads()
{
   // Avoid unnecessary output
   gROOT->SetBatch();

   // Make ROOT thread-safe
   ROOT::EnableThreadSafety();

   // Total number of events
   const size_t nEntries = 65535;

   // Match number of threads to what the hardware can do
   const size_t nWorkers = 4;

   // Split work in equal parts
   const size_t nEventsPerWorker = nEntries / nWorkers;

   // Create the TBufferMerger: this class orchestrates the parallel writing
   auto fileName = "mt103_fillNtupleFromMultipleThreads.root";
   ROOT::Experimental::TBufferMerger merger(fileName);

   // Define what each worker will do
   // We obtain from a merger a TBufferMergerFile, which is nothing more than
   // a file which is held in memory and that flushes to the TBufferMerger its
   // content.
   auto work_function = [&](int seed) {
      auto f = merger.GetFile();
      TNtuple ntrand("ntrand", "Random Numbers", "r");

      TRandom rnd(seed);
      for (auto i : ROOT::TSeqI(nEntries))
         ntrand.Fill(rnd.Gaus());
      f->Write();
   };

   // Create worker threads
   std::vector<std::thread> workers;

   for (auto i : ROOT::TSeqI(nWorkers))
      workers.emplace_back(work_function, i + 1); // seed==0 means random seed :)

   // Make sure workers are done
   for (auto &&worker : workers)
      worker.join();
}
