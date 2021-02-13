/// \file
/// \ingroup tutorial_multicore
/// \notebook
/// Fill n-tuples in distinct workers.
/// This tutorial illustrates the basics of how it's possible with ROOT to
/// offload heavy operations on multiple threads and how it's possible to write
/// simultaneously multiple files. The operation performed in this case is the
/// creation of random gaussian numbers.
/// NOTE: this code can be executed in a macro, ACLiC'ed or not, but not yet at
/// the command line prompt.
///
/// \macro_code
///
/// \date January 2016
/// \author Danilo Piparo

// Some useful constants and functions

// Total amount of numbers
const UInt_t nNumbers = 20000000U;

// The number of workers
const UInt_t nWorkers = 4U;

// We split the work in equal parts
const auto workSize = nNumbers / nWorkers;

// A simple function to fill ntuples randomly
void fillRandom(TNtuple &ntuple, TRandom3 &rndm, UInt_t n)
{
   for (auto i : ROOT::TSeqI(n))
      ntuple.Fill(rndm.Gaus());
}

Int_t mt101_fillNtuples()
{

   // No nuisance for batch execution
   gROOT->SetBatch();

   // Perform the operation sequentially ---------------------------------------

   // Create a random generator and and Ntuple to hold the numbers
   TRandom3 rndm(1);
   TFile ofile("mt101_singleCore.root", "RECREATE");
   TNtuple randomNumbers("singleCore", "Random Numbers", "r");
   fillRandom(randomNumbers, rndm, nNumbers);
   randomNumbers.Write();
   ofile.Close();

   // We now go MT! ------------------------------------------------------------

   // The first, fundamental operation to be performed in order to make ROOT
   // thread-aware.
   ROOT::EnableThreadSafety();

   // We define our work item
   auto workItem = [](UInt_t workerID) {
      // One generator, file and ntuple per worker
      TRandom3 workerRndm(workerID); // Change the seed
      TFile ofile(Form("mt101_multiCore_%u.root", workerID), "RECREATE");
      TNtuple workerRandomNumbers("multiCore", "Random Numbers", "r");
      fillRandom(workerRandomNumbers, workerRndm, workSize);
      workerRandomNumbers.Write();
      return 0;
   };

   // Create the collection which will hold the threads, our "pool"
   std::vector<std::thread> workers;

   // Fill the "pool" with workers
   for (auto workerID : ROOT::TSeqI(nWorkers)) {
      workers.emplace_back(workItem, workerID);
   }

   // Now join them
   for (auto &&worker : workers)
      worker.join();

   return 0;
}
