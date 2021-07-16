/// \file
/// \ingroup tutorial_ntuple
/// \notebook
/// Example of multi-threaded writes using multuple REntry objects
///
/// \macro_image
/// \macro_code
///
/// \date July 2021
/// \author The ROOT Team

// NOTE: The RNTuple classes are experimental at this point.
// Functionality, interface, and data format is still subject to changes.
// Do not use for real data!

// Until C++ runtime modules are universally used, we explicitly load the ntuple library.  Otherwise
// triggering autoloading from the use of templated types would require an exhaustive enumeration
// of "all" template instances in the LinkDef file.
R__LOAD_LIBRARY(ROOTNTuple)

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>

#include <TCanvas.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TSystem.h>

#include <atomic>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>
#include <utility>

// Import classes from experimental namespace for the time being
using REntry = ROOT::Experimental::REntry;
using RNTupleModel = ROOT::Experimental::RNTupleModel;
using RNTupleReader = ROOT::Experimental::RNTupleReader;
using RNTupleWriter = ROOT::Experimental::RNTupleWriter;

// Where to store the ntuple of this example
constexpr char const *kNTupleFileName = "ntpl007_mtFill.root";

// Number of parallel threads to fill the ntuple
constexpr int kNWriterThreads = 4;

// Number of events to generate is kNEventsPerThread * kNWriterThreads
constexpr int kNEventsPerThread = 25000;

// Thread function to generate and write events
void FillData(std::unique_ptr<REntry> entry, RNTupleWriter *ntuple) {
   // Protect the ntuple->Fill() call
   static std::mutex gLock;

   static std::atomic<std::uint32_t> gThreadId;
   auto threadId = gThreadId++;

   auto prng = std::make_unique<TRandom3>();
   prng->SetSeed();

   auto id = entry->Get<std::uint32_t>("id");
   auto vpx = entry->Get<std::vector<float>>("vpx");
   auto vpy = entry->Get<std::vector<float>>("vpy");
   auto vpz = entry->Get<std::vector<float>>("vpz");

   for (int i = 0; i < kNEventsPerThread; i++) {
      vpx->clear();
      vpy->clear();
      vpz->clear();

      *id = threadId * kNEventsPerThread + i;

      int npx = static_cast<int>(prng->Rndm(1) * 15);
      // Set the field data for the current event
      for (int j = 0; j < npx; ++j) {
         float px, py, pz;
         prng->Rannor(px, py);
         pz = px*px + py*py;

         vpx->emplace_back(px);
         vpy->emplace_back(py);
         vpz->emplace_back(pz);
      }

      std::lock_guard<std::mutex> guard(gLock);
      ntuple->Fill(*entry);
   }
}

// Generate kNEvents with multiple threads in kNTupleFileName
void Write()
{
   // Create the data model
   auto model = RNTupleModel::Create();
   model->MakeField<std::uint32_t>("id");
   model->MakeField<std::vector<float>>("vpx");
   model->MakeField<std::vector<float>>("vpy");
   model->MakeField<std::vector<float>>("vpz");

   // We hand-over the data model to a newly created ntuple of name "NTuple", stored in kNTupleFileName
   auto ntuple = RNTupleWriter::Recreate(std::move(model), "NTuple", kNTupleFileName);

   std::vector<std::unique_ptr<REntry>> entries;
   std::vector<std::thread> threads;
   for (int i = 0; i < kNWriterThreads; ++i)
      entries.emplace_back(ntuple->CreateEntry());
   for (int i = 0; i < kNWriterThreads; ++i)
      threads.emplace_back(FillData, std::move(entries[i]), ntuple.get());
   for (int i = 0; i < kNWriterThreads; ++i)
      threads[i].join();

   // The ntuple unique pointer goes out of scope here.  On destruction, the ntuple flushes unwritten data to disk
   // and closes the attached ROOT file.
}


// For all of the events, histogram only one of the written vectors
void Read()
{
   auto ntuple = RNTupleReader::Open("NTuple", kNTupleFileName);
   // TODO(jblomer): the "inner name" of the vector should become "vpx._0"
   auto viewVpx = ntuple->GetView<float>("vpx.float");

   TCanvas *c2 = new TCanvas("c2", "Dynamic Filling Example", 200, 10, 700, 500);

   TH1F h("h", "This is the px distribution", 100, -4, 4);
   h.SetFillColor(48);
   // Iterate through all values of vpx in all events
   for (auto i : viewVpx.GetFieldRange())
      h.Fill(viewVpx(i));
   // Prevent the histogram from disappearing
   h.DrawCopy();

   auto nEvents = ntuple->GetNEntries();
   auto viewId = ntuple->GetView<std::uint32_t>("id");
   TH2F hFillSequence("","Entry Id vs Event Id", 100, 0, nEvents, 100, 0, nEvents);
   for (auto i : ntuple->GetEntryRange())
      hFillSequence.Fill(i, viewId(i));
   hFillSequence.DrawCopy();
}


void ntpl007_mtFill()
{
   Write();
   Read();
}
