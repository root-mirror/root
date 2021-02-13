/// \file
/// \ingroup tutorial_multicore
/// \notebook
/// Fill n-tuples in distinct workers.
/// This tutorial illustrates the basics of how it's possible with ROOT to
/// offload heavy operations on multiple processes and how it's possible to write
/// simultaneously multiple files. The operation performed in this case is the
/// creation of random gaussian numbers.
///
/// \macro_code
///
/// \date January 2016
/// \author Danilo Piparo

// Some useful constants and functions

// Total amount of numbers
const UInt_t nNumbers = 20000000U;

// The number of workers
const UInt_t nThreads = 4U;

// We split the work in equal parts
const auto workSize = nNumbers / nThreads;

// A simple function to fill ntuples randomly
void fillRandom(TNtuple &ntuple, TRandom3 &rndm, UInt_t n)
{
   for (auto i : ROOT::TSeqI(n))
      ntuple.Fill(rndm.Gaus());
}

Int_t mtbb101_fillNtuples()
{
   // No nuisance for batch execution
   gROOT->SetBatch();

   // Perform the operation sequentially ---------------------------------------

   // Create a random generator and and Ntuple to hold the numbers
   TRandom3 rndm(1);
   TFile ofile("mtbb101_singleCore.root", "RECREATE");
   TNtuple randomNumbers("singleCore", "Random Numbers", "r");
   fillRandom(randomNumbers, rndm, nNumbers);
   randomNumbers.Write();
   ofile.Close();

   // We now go MP! ------------------------------------------------------------

   // We define our work item
   auto workItem = [](UInt_t workerID) {
      // One generator, file and ntuple per worker
      TRandom3 workerRndm(workerID); // Change the seed
      TFile ofile(Form("mtbb101_multiCore_%u.root", workerID), "RECREATE");
      TNtuple workerRandomNumbers("multiCore", "Random Numbers", "r");
      fillRandom(workerRandomNumbers, workerRndm, workSize);
      workerRandomNumbers.Write();
      return 0;
   };

   // Create the pool of workers
   ROOT::TThreadExecutor pool(nThreads);

   // Fill the pool with work
   pool.Map(workItem, ROOT::TSeqI(nThreads));

   return 0;
}
