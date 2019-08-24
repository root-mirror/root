#include "gtest/gtest.h"

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleDS.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RPageStorageRaw.hxx>

#include <TRandom3.h>

#include <cstdio>
#include <memory>
#include <string>
#include <vector>
#include <utility>

using RNTupleModel = ROOT::Experimental::RNTupleModel;
using RNTupleReader = ROOT::Experimental::RNTupleReader;
using RNTupleWriter = ROOT::Experimental::RNTupleWriter;
using RPageSinkRaw = ROOT::Experimental::Detail::RPageSinkRaw;
using RPageSource = ROOT::Experimental::Detail::RPageSource;
using RPageSourceRaw = ROOT::Experimental::Detail::RPageSourceRaw;

namespace {

/**
 * An RAII wrapper around an open temporary file on disk. It cleans up the guarded file when the wrapper object
 * goes out of scope.
 */
class FileRaii {
private:
   std::string fPath;
public:
   explicit FileRaii(const std::string &path) : fPath(path) { }
   FileRaii(const FileRaii&) = delete;
   FileRaii& operator=(const FileRaii&) = delete;
   ~FileRaii() { std::remove(fPath.c_str()); }
   std::string GetPath() const { return fPath; }
};

} // anonymous namespace

TEST(RNTuple, Basics)
{
   FileRaii fileGuard("test_ntuple_rawfile.ntuple");

   auto model = RNTupleModel::Create();
   auto wrPt = model->MakeField<float>("pt", 42.0);

   {
      RNTupleWriter ntuple(std::move(model), std::make_unique<RPageSinkRaw>("f", fileGuard.GetPath()));
      ntuple.Fill();
      ntuple.CommitCluster();
      *wrPt = 24.0;
      ntuple.Fill();
      *wrPt = 12.0;
      ntuple.Fill();
   }

   RNTupleReader ntuple(std::make_unique<RPageSourceRaw>("f", fileGuard.GetPath()));
   EXPECT_EQ(3U, ntuple.GetNEntries());
   auto rdPt = ntuple.GetModel()->GetDefaultEntry()->Get<float>("pt");

   ntuple.LoadEntry(0);
   EXPECT_EQ(42.0, *rdPt);
   ntuple.LoadEntry(1);
   EXPECT_EQ(24.0, *rdPt);
   ntuple.LoadEntry(2);
   EXPECT_EQ(12.0, *rdPt);
}


TEST(RNTuple, Extended)
{
   FileRaii fileGuard("test_ntuple_rawfile_ext.ntuple");

   auto model = RNTupleModel::Create();
   auto wrVector = model->MakeField<std::vector<double>>("vector");

   TRandom3 rnd(42);
   double chksumWrite = 0.0;
   {
      RNTupleWriter ntuple(std::move(model), std::make_unique<RPageSinkRaw>("f", fileGuard.GetPath()));
      constexpr unsigned int nEvents = 32000;
      for (unsigned int i = 0; i < nEvents; ++i) {
         auto nVec = 1 + floor(rnd.Rndm() * 1000.);
         wrVector->resize(nVec);
         for (unsigned int n = 0; n < nVec; ++n) {
            auto val = 1 + rnd.Rndm()*1000. - 500.;
            (*wrVector)[n] = val;
            chksumWrite += val;
         }
         ntuple.Fill();
         if (i % 1000 == 0)
            ntuple.CommitCluster();
      }
   }

   RNTupleReader ntuple(std::make_unique<RPageSourceRaw>("f", fileGuard.GetPath()));
   auto rdVector = ntuple.GetModel()->GetDefaultEntry()->Get<std::vector<double>>("vector");

   double chksumRead = 0.0;
   for (auto entryId : ntuple) {
      ntuple.LoadEntry(entryId);
      for (auto v : *rdVector)
         chksumRead += v;
   }

   EXPECT_EQ(chksumRead, chksumWrite);

   ROOT::EnableImplicitMT();
   auto rdf = ROOT::Experimental::MakeNTupleDataFrame("f", fileGuard.GetPath());
   auto minLenght = *rdf.Define("vecSize", "vector.size()").Min("vecSize");
   EXPECT_GE(1, minLenght);
   EXPECT_LE(minLenght, 1000);
}
