#include <stdio.h>

#include "Bytes.h"
#include "TBranch.h"
#include "TBufferFile.h"
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "ROOT/TTreeReaderFast.hxx"
#include "ROOT/TTreeReaderValueFast.hxx"


#include "gtest/gtest.h"

class BulkApiTest : public ::testing::Test {
public:
   static constexpr Long64_t fEventCount = (Long64_t)1e7;
   const std::string fFileName = "BulkApiTest.root";

protected:
   virtual void SetUp()
   {
      auto hfile = TFile::Open(fFileName.c_str(), "recreate", "TTree float micro benchmark ROOT file");
      hfile->SetCompressionLevel(0); // No compression at all.

      // Otherwise, we keep with the current ROOT defaults.
      auto tree = new TTree("T", "A ROOT tree of floats.");
      float f = 2;
      tree->Branch("myFloat", &f, 320000, 1);
      for (Long64_t ev = 0; ev < fEventCount; ev++) {
         tree->Fill();
         ++f;
      }
      hfile = tree->GetCurrentFile();
      hfile->Write();
      tree->Print();
      printf("Successful write of all events.\n");

      // We also want a copy the TTree with a basket stored alongside the TTree
      // rather than on its own.
      for (Long64_t ev = 0; ev < 10; ev++) {
         tree->Fill();
         ++f;
      }
      hfile->WriteTObject(tree, "TwithBasket");
      hfile->Close();

      delete hfile;
   }
};

TEST_F(BulkApiTest, stdRead)
{
   auto hfile = TFile::Open(fFileName.c_str());
   printf("Starting read of file %s.\n", fFileName.c_str());
   TStopwatch sw;

   printf("Using standard read APIs.\n");
   // Read via standard APIs.
   TTreeReader myReader("T", hfile);
   TTreeReaderValue<float> myF(myReader, "myFloat");
   Long64_t idx = 0;
   float idx_f = 1;
   auto events = fEventCount;
   sw.Start();
   while (myReader.Next()) {
      if (R__unlikely(idx == events)) {
         break;
      }
      idx_f++;
      if (R__unlikely((idx < 16000000) && (abs((*myF) - idx_f) > std::numeric_limits<float>::epsilon()))) {
         printf("Incorrect value on myFloat branch: %f, expected %f (event %lld)\n", *myF, idx_f, idx);
         ASSERT_TRUE(false);
      }
      idx++;
   }
   sw.Stop();
   printf("TTreeReader: Successful read of all events.\n");
   printf("TTreeReader: Total elapsed time (seconds) for standard APIs: %.2f\n", sw.RealTime());
   delete hfile;
}

void SimpleReadFunc(const char *filename, const char *treename)
{
   auto hfile = TFile::Open(filename);
   printf("Starting read of file %s.\n", filename);
   TStopwatch sw;

   printf("Using inline bulk read APIs.\n");
   TBufferFile branchbuf(TBuffer::kWrite, 32 * 1024);
   auto tree = hfile->Get<TTree>(treename);
   ASSERT_TRUE(tree);

   TBranch *branchF = tree->GetBranch("myFloat");
   ASSERT_TRUE(branchF);
   branchF->GetListOfBaskets()->ls();

   float idx_f = 1;
   Long64_t evt_idx = 0;
   Long64_t events = tree->GetEntries();
   while (events) {
      auto count = branchF->GetBulkRead().GetEntriesSerialized(evt_idx, branchbuf);
      ASSERT_GE(count, 0);
      events = events > count ? (events - count) : 0;

      auto entry = reinterpret_cast<float *>(branchbuf.GetCurrent());
      for (Int_t idx = 0; idx < count; idx++) {
         idx_f++;
         auto tmp = *reinterpret_cast<Int_t *>(&entry[idx]);
         auto tmp_ptr = reinterpret_cast<char *>(&tmp);
         frombuf(tmp_ptr, entry + idx);

         if (R__unlikely((evt_idx < 16000000) && (entry[idx] != idx_f))) {
            printf("in %s Incorrect value on myFloat branch: %f (event %lld)\n", treename, entry[idx], evt_idx + idx);
            ASSERT_TRUE(false);
         }
      }
      evt_idx += count;
   }
   sw.Stop();
   printf("GetEntriesSerialized: Successful read of all events in %s.\n", treename);
   printf("GetEntriesSerialized: Total elapsed time (seconds) for bulk APIs: %.2f\n", sw.RealTime());
   delete hfile;
}

TEST_F(BulkApiTest, simpleRead)
{
   SimpleReadFunc(fFileName.c_str(), "T");
}

TEST_F(BulkApiTest, simpleReadTrailingBasket)
{
   SimpleReadFunc(fFileName.c_str(), "TwithBasket");
}

void SimpleBulkReadFunc(const char *filename, const char *treename)
{
   auto hfile = TFile::Open(filename);
   printf("Starting read of file %s.\n", filename);
   TStopwatch sw;

   printf("Using outlined bulk read APIs.\n");
   TBufferFile branchbuf(TBuffer::kWrite, 32 * 1024);
   auto tree = hfile->Get<TTree>(treename);
   ASSERT_TRUE(tree);

   TBranch *branchF = tree->GetBranch("myFloat");
   ASSERT_TRUE(branchF);
   branchF->GetListOfBaskets()->ls();

   float idx_f = 1;
   Long64_t evt_idx = 0;
   Long64_t events = tree->GetEntries();
   while (events) {
      auto count = branchF->GetBulkRead().GetBulkEntries(evt_idx, branchbuf);
      ASSERT_GE(count, 0);
      events = events > count ? (events - count) : 0;

      auto entry = reinterpret_cast<float *>(branchbuf.GetCurrent());
      for (Int_t idx = 0; idx < count; idx++) {
         idx_f++;
         if (R__unlikely((evt_idx < 16000000) && (entry[idx] != idx_f))) {
            printf("in %s Incorrect value on myFloat branch: %f (event %lld)\n", treename, entry[idx], evt_idx + idx);
            ASSERT_TRUE(false);
         }
      }
      evt_idx += count;
   }
   sw.Stop();
   printf("GetBulkEntries: Successful read of all events in %s.\n", treename);
   printf("GetBulkEntries: Total elapsed time (seconds) for bulk APIs: %.2f\n", sw.RealTime());
   delete hfile;
}

TEST_F(BulkApiTest, simpleBulkRead)
{
   SimpleBulkReadFunc(fFileName.c_str(), "T");
}

TEST_F(BulkApiTest, simpleBulkReadTrailingBasket)
{
   SimpleBulkReadFunc(fFileName.c_str(), "TwithBasket");
}

TEST_F(BulkApiTest, fastRead)
{
   auto hfile = TFile::Open(fFileName.c_str());
   printf("Starting read of file %s.\n", fFileName.c_str());
   TStopwatch sw;

   printf("Using TTreeReaderFast.\n");
   ROOT::Experimental::TTreeReaderFast myReader("T", hfile);
   ROOT::Experimental::TTreeReaderValueFast<float> myF(myReader, "myFloat");
   myReader.SetEntry(0);
   if (ROOT::Internal::TTreeReaderValueBase::kSetupMatch != myF.GetSetupStatus()) {
      printf("TTreeReaderValueFast<float> failed to initialize.  Status code: %d\n", myF.GetSetupStatus());
      ASSERT_TRUE(false);
   }
   if (myReader.GetEntryStatus() != TTreeReader::kEntryValid) {
      printf("TTreeReaderFast failed to initialize.  Entry status: %d\n", myReader.GetEntryStatus());
      ASSERT_TRUE(false);
   }
   auto events = fEventCount;
   Long64_t idx = 0;
   float idx_f = 1;
   for (auto reader_idx : myReader) {
      ASSERT_LT(reader_idx, events);
      ASSERT_EQ(reader_idx, idx);
      idx_f++;
      if (R__unlikely((idx < 16000000) && (abs((*myF) - idx_f) > std::numeric_limits<float>::epsilon()))) {
         printf("Incorrect value on myFloat branch: %f, expected %f (event %lld)\n", *myF, idx_f, idx);
         ASSERT_TRUE(false);
      }
      idx++;
   }
   sw.Stop();
   printf("TTreeReaderFast: Successful read of all events.\n");
   printf("TTreeReaderFast: Total elapsed time (seconds) for bulk APIs: %.2f\n", sw.RealTime());
   delete hfile;
}

TEST_F(BulkApiTest, BulkInMem)
{
   TTree t("t","t");
   int i = 3;
   t.Branch("fAlpha",&i);
   for(i = 0; i < 100000; ++i)
      t.Fill();

   TBufferFile buf(TBuffer::EMode::kWrite, 32*1024);
   auto &r = t.GetBranch("fAlpha")->GetBulkRead();

   {
      auto s = r.GetBulkEntries(0, buf);
      ASSERT_EQ(7982, s) << "Did not read the expected number of entries.";
      s = r.GetBulkEntries(0, buf);
      ASSERT_EQ(7982, s) << "Did not read the expected number of entries.";
   }

   int iteration = 0;
   Long64_t nentries = t.GetEntries();
   Long64_t event_idx = 0;
   while (nentries) {
      auto s = r.GetBulkEntries(event_idx, buf);
      if (iteration < 12)
         ASSERT_EQ(7982, s) << "Did not read the expected number of entries.";
      else
         ASSERT_EQ(4216, s) << "Did not read the expected number of entries.";
      nentries -= s;
      event_idx += s;
      ++iteration;
      if (s < 0)
         break;
   }
}
