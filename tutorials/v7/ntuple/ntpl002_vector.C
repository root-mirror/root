/// \file
/// \ingroup tutorial_ntuple
/// \notebook
/// Write and read STL vectors with RNTuple.  Adapted from the hvector tree tutorial.
///
/// \macro_image
/// \macro_code
///
/// \date April 2019
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

#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <utility>

// Import classes from experimental namespace for the time being
using RNTupleModel = ROOT::Experimental::RNTupleModel;
using RNTupleReader = ROOT::Experimental::RNTupleReader;
using RNTupleWriter = ROOT::Experimental::RNTupleWriter;

// Where to store the ntuple of this example
constexpr char const* kNTupleFileName = "ntpl002_vector.root";

// Update the histogram GUI every so many fills
constexpr int kUpdateGuiFreq = 1000;

// Number of events to generate
constexpr int kNEvents = 25000;

// Generate kNEvents with vectors in kNTupleFileName
void Write()
{
   // We create a unique pointer to an empty data model
   auto model = RNTupleModel::Create();

   // Creating fields of std::vector is the same as creating fields of simple types.  As a result, we get
   // shared pointers of the given type
   std::shared_ptr<std::vector<float>> fldVpx = model->MakeField<std::vector<float>>("vpx");
   auto fldVpy   = model->MakeField<std::vector<float>>("vpy");
   auto fldVpz   = model->MakeField<std::vector<float>>("vpz");
   auto fldVrand = model->MakeField<std::vector<float>>("vrand");

   // We hand-over the data model to a newly created ntuple of name "F", stored in kNTupleFileName
   // In return, we get a unique pointer to an ntuple that we can fill
   auto ntuple = RNTupleWriter::Recreate(std::move(model), "F", kNTupleFileName);

   TH1F hpx("hpx", "This is the px distribution", 100, -4, 4);
   hpx.SetFillColor(48);

   auto c1 = new TCanvas("c1", "Dynamic Filling Example", 200, 10, 700, 500);

   gRandom->SetSeed();
   for (int i = 0; i < kNEvents; i++) {
      int npx = static_cast<int>(gRandom->Rndm(1) * 15);

      fldVpx->clear();
      fldVpy->clear();
      fldVpz->clear();
      fldVrand->clear();

      // Set the field data for the current event
      for (int j = 0; j < npx; ++j) {
         float px, py, pz;
         gRandom->Rannor(px, py);
         pz = px*px + py*py;
         auto random = gRandom->Rndm(1);

         hpx.Fill(px);

         fldVpx->emplace_back(px);
         fldVpy->emplace_back(py);
         fldVpz->emplace_back(pz);
         fldVrand->emplace_back(random);
      }

      // Gui updates
      if (i && (i % kUpdateGuiFreq) == 0) {
         if (i == kUpdateGuiFreq) hpx.Draw();
         c1->Modified();
         c1->Update();
         if (gSystem->ProcessEvents())
            break;
      }

      ntuple->Fill();
   }

   hpx.DrawCopy();

   // The ntuple unique pointer goes out of scope here.  On destruction, the ntuple flushes unwritten data to disk
   // and closes the attached ROOT file.
}


// For all of the events, histogram only one of the written vectors
void Read()
{
   // Get a unique pointer to an empty RNTuple model
   auto model = RNTupleModel::Create();

   // We only define the fields that are needed for reading
   auto fldVpx = model->MakeField<std::vector<float>>("vpx");

   // Create an ntuple without imposing a specific data model.  We could generate the data model from the ntuple
   // but here we prefer the view because we only want to access a single field
   auto ntuple = RNTupleReader::Open(std::move(model), "F", kNTupleFileName);

   // Quick overview of the ntuple's key meta-data
   ntuple->PrintInfo();

   std::cout << "Entry number 42 in JSON format:" << std::endl;
   ntuple->Show(41);
   // In a future version of RNTuple, there will be support for ntuple->Scan()

   TCanvas *c2 = new TCanvas("c2", "Dynamic Filling Example", 200, 10, 700, 500);
   TH1F h("h", "This is the px distribution", 100, -4, 4);
   h.SetFillColor(48);

   // Iterate through all the events using i as event number and as an index for accessing the view
   for (auto entryId : *ntuple) {
      ntuple->LoadEntry(entryId);

      for (auto px : *fldVpx) {
         h.Fill(px);
      }

      if (entryId && (entryId % kUpdateGuiFreq) == 0) {
         if (entryId == kUpdateGuiFreq) h.Draw();
         c2->Modified();
         c2->Update();
         if (gSystem->ProcessEvents())
            break;
      }
   }

   // Prevent the histogram from disappearing
   h.DrawCopy();
}


void ntpl002_vector()
{
   Write();
   Read();
}
