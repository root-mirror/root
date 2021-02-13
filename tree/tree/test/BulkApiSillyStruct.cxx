#include <stdio.h>

#include "Bytes.h"
#include "Rtypes.h"
#include "SillyStruct.h"
#include "TBranch.h"
#include "TBufferFile.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "ROOT/TTreeReaderFast.hxx"
#include "ROOT/TTreeReaderValueFast.hxx"
#include "ROOT/TIOFeatures.hxx"

#include "gtest/gtest.h"

class BulkApiSillyStructTest : public ::testing::Test {
public:
   static constexpr Long64_t fClusterSize = 1e5;
   static constexpr Long64_t fEventCount = 4e6;
   const std::string fFileName = "BulkApiSillyStruct.root";

protected:
   virtual void SetUp()
   {
      TFile *hfile = new TFile(fFileName.c_str(), "RECREATE", "TTree silly-struct benchmark");
      hfile->SetCompressionLevel(0); // No compression at all.

      TTree *tree = new TTree("T", "A ROOT tree of silly-struct branches.");
      tree->SetBit(TTree::kOnlyFlushAtCluster);
      tree->SetAutoFlush(fClusterSize);
      ROOT::TIOFeatures features;
      features.Set(ROOT::Experimental::EIOFeatures::kGenerateOffsetMap);
      tree->SetIOFeatures(features);

      SillyStruct ss;
      TBranch *branch = tree->Branch("myEvent", &ss, 32000, 99);
      branch->SetAutoDelete(kFALSE);

      Int_t nb = 0;
      for (Long64_t ev = 0; ev < fEventCount; ev++) {
         ss.i = ev;
         ss.f = ev;
         ss.d = ev;
         nb += tree->Fill();
      } 

      hfile = tree->GetCurrentFile();
      hfile->Write();
      tree->Print();
      printf("Successful write of all events, nb = %d.\n", nb);

      delete hfile;
   }
};

constexpr Long64_t BulkApiSillyStructTest::fClusterSize;
constexpr Long64_t BulkApiSillyStructTest::fEventCount;

TEST_F(BulkApiSillyStructTest, stdReadStruct)
{
   auto hfile  = TFile::Open(fFileName.c_str());

   TTreeReader myReader("T", hfile);
   TTreeReaderValue<SillyStruct>  ss(myReader, "myEvent");

   int    evI = 0;
   float  evF = 0.0;
   double evD = 0.0;
   while (myReader.Next()) {
      ASSERT_TRUE(ss->f == evF);
      ASSERT_TRUE(ss->i == evI);
      ASSERT_TRUE(ss->d == evD);
      evI++;
      evF++;
      evD++;
   }
}

TEST_F(BulkApiSillyStructTest, stdReadSplitBranch)
{
   auto hfile  = TFile::Open(fFileName.c_str());

   TTreeReader myReader("T", hfile);
   TTreeReaderValue<float>        myF(myReader, "f");
   TTreeReaderValue<int>          myI(myReader, "i");
   TTreeReaderValue<double>       myD(myReader, "d");

   int    evI = 0;
   float  evF = 0.0;
   double evD = 0.0;
   while (myReader.Next()) {
      ASSERT_TRUE(*myF == evF);
      ASSERT_TRUE(*myI == evI);
      ASSERT_TRUE(*myD == evD);
      evI++;
      evF++;
      evD++;
   }
}

TEST_F(BulkApiSillyStructTest, fastRead)
{
   TBufferFile bufF(TBuffer::kWrite, 10000);
   TBufferFile bufI(TBuffer::kWrite, 10000);
   TBufferFile bufD(TBuffer::kWrite, 10000);
   auto hfile  = TFile::Open(fFileName.c_str());
   auto tree = dynamic_cast<TTree*>(hfile->Get("T"));
   ASSERT_TRUE(tree);
   TBranch *branchF = tree->GetBranch("f");
   TBranch *branchI = tree->GetBranch("i");
   TBranch *branchD = tree->GetBranch("d");

   Long64_t events = fEventCount;
   Int_t count = std::min(fClusterSize, fEventCount);
   Long64_t evt_idx = 0;
   int      evI = 0;
   float    evF = 0.0;
   double   evD = 0.0;

   while (events) {
      auto countF = branchF->GetBulkRead().GetBulkEntries(evt_idx, bufF);
      auto countI = branchI->GetBulkRead().GetBulkEntries(evt_idx, bufI);
      auto countD = branchD->GetBulkRead().GetBulkEntries(evt_idx, bufD);
      ASSERT_EQ(countF, count);
      ASSERT_EQ(countI, count);
      ASSERT_EQ(countD, count);

      if (events > count) {
         events -= count;
      } else {
         events = 0;
      }
      float *float_buf = reinterpret_cast<float*>(bufF.GetCurrent());
      double *double_buf = reinterpret_cast<double*>(bufD.GetCurrent());
      int *int_buf = reinterpret_cast<int*>(bufI.GetCurrent());
      for (Int_t idx = 0; idx < count; idx++) {
         ASSERT_EQ(float_buf[idx], evF);
         ASSERT_EQ(int_buf[idx], evI);
         ASSERT_EQ(double_buf[idx], evD);
         evF++;
         evI++;
         evD++;
         evt_idx++;
      }
   }
}
