/// \file
/// \ingroup tutorial_multicore
/// \notebook
/// Fill histograms in parallel and write them on file.
/// This example expresses the parallelism of the mt001_fillHistos.C tutorial
/// with multiprocessing techniques.
///
/// \macro_code
///
/// \date January 2016
/// \author Danilo Piparo

// Total amount of numbers
const UInt_t nNumbers = 20000000U;

// The number of workers
const UInt_t nThreads = 4U;

Int_t mtbb001_fillHistos()
{
   // We define our work item
   auto workItem = [](UInt_t workerID) {
      // One generator, file and ntuple per worker
      TRandom3 workerRndm(workerID); // Change the seed
      TFile f(Form("myFile_mtbb001_%u.root", workerID), "RECREATE");
      TH1F h(Form("myHisto_%u", workerID), "The Histogram", 64, -4, 4);
      for (UInt_t i = 0; i < nNumbers; ++i) {
         h.Fill(workerRndm.Gaus());
      }
      h.Write();
      return 0;
   };

   // Create the pool of threads
   ROOT::TThreadExecutor pool(nThreads);

   // Fill the pool with work
   pool.Map(workItem, ROOT::TSeqI(nThreads));
   return 0;
}
