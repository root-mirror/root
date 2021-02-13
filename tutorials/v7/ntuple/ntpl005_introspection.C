/// \file
/// \ingroup tutorial_ntuple
/// \notebook
/// Write and read an RNTuple from a user-defined class.  Adapted from tv3.C
/// Illustrates various RNTuple introspection methods.
///
/// \macro_image
/// \macro_code
///
/// \date April 2020
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
#include <ROOT/RNTupleOptions.hxx>

#include <Compression.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TRandom.h>
#include <TSystem.h>

#include <cassert>

// Import classes from experimental namespace for the time being
using ENTupleInfo = ROOT::Experimental::ENTupleInfo;
using ENTupleShowFormat = ROOT::Experimental::ENTupleShowFormat;
using RNTupleModel = ROOT::Experimental::RNTupleModel;
using RNTupleReader = ROOT::Experimental::RNTupleReader;
using RNTupleWriter = ROOT::Experimental::RNTupleWriter;
using RNTupleWriteOptions = ROOT::Experimental::RNTupleWriteOptions;

constexpr char const* kNTupleFileName = "ntpl005_introspection.root";

// Store entries of type Vector3 in the ntuple
class Vector3 {
private:
   double fX = 0;
   double fY = 0;
   double fZ = 0;

public:
   double x() const { return fX; }
   double y() const { return fY; }
   double z() const { return fZ; }

   void SetXYZ(double x, double y, double z) {
      fX = x;
      fY = y;
      fZ = z;
   }
};


void Generate()
{
   auto model = RNTupleModel::Create();
   auto fldVector3 = model->MakeField<Vector3>("v3");

   // Explicitly enforce a certain compression algorithm
   RNTupleWriteOptions options;
   options.SetCompression(ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose);

   auto ntuple = RNTupleWriter::Recreate(std::move(model), "Vector3", kNTupleFileName, options);
   TRandom r;
   for (unsigned int i = 0; i < 500000; ++i) {
      fldVector3->SetXYZ(r.Gaus(0,1), r.Landau(0,1), r.Gaus(100,10));
      ntuple->Fill();
   }
}


void ntpl005_introspection() {
   Generate();

   auto ntuple = RNTupleReader::Open("Vector3", kNTupleFileName);

   // Display the schema of the ntuple
   ntuple->PrintInfo();

   // Display information about the storage layout of the data
   ntuple->PrintInfo(ENTupleInfo::kStorageDetails);

   // Display the first entry; no model was defined, hence we use kCompleteJSON to force showing all fields
   ntuple->Show(0, ENTupleShowFormat::kCompleteJSON);

   // Collect I/O runtime counters when processing the data set.
   // Maintaining the counters comes with a small performance overhead, so it has to be explicitly enabled
   ntuple->EnableMetrics();

   // Plot the y components of vector3
   TCanvas *c1 = new TCanvas("c1","RNTuple Demo", 10, 10, 600, 800);
   c1->Divide(1,2);
   c1->cd(1);
   TH1F h1("x", "x component of Vector3", 100, -3, 3);
   {
      /// We enclose viewX in a scope in order to indicate to the RNTuple when we are not
      /// anymore interested in v3.fX
      auto viewX = ntuple->GetView<double>("v3.fX");
      for (auto i : ntuple->GetEntryRange()) {
         h1.Fill(viewX(i));
      }
   }
   h1.DrawCopy();

   c1->cd(2);
   TH1F h2("y", "y component of Vector3", 100, -5, 20);
   auto viewY = ntuple->GetView<double>("v3.fY");
   for (auto i : ntuple->GetEntryRange()) {
      h2.Fill(viewY(i));
   }
   h2.DrawCopy();

   // Display the I/O operation statistics performed by the RNTuple reader
   ntuple->PrintInfo(ENTupleInfo::kMetrics);

   // We read 2 out of the 3 Vector3 members and thus should have requested approximately 2/3 of the file
   FileStat_t fileStat;
   auto retval = gSystem->GetPathInfo(kNTupleFileName, fileStat);
   assert(retval == 0);
   float fileSize = static_cast<float>(fileStat.fSize);
   float nbytesRead = ntuple->GetMetrics().GetCounter("RNTupleReader.RPageSourceFile.szReadPayload")->GetValueAsInt() +
                      ntuple->GetMetrics().GetCounter("RNTupleReader.RPageSourceFile.szReadOverhead")->GetValueAsInt();

   std::cout << "File size:      " << fileSize / 1024. / 1024. << " MiB" << std::endl;
   std::cout << "Read from file: " << nbytesRead / 1024. / 1024. << " MiB" << std::endl;
   std::cout << "Ratio:          " << nbytesRead / fileSize << std::endl;
}
