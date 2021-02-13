#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBasket.h"
#include "TBranchElement.h"
#include "TLeafElement.h"
#include "TRandom.h"

#include "gtest/gtest.h"

#include "ElementStruct.h"

class TOffsetGeneration : public ::testing::Test {
protected:
   static constexpr int fEventCount = 10000;

   virtual void SetUp()
   {
      TRandom *random = new TRandom(837);
      auto file = new TFile("TOffsetGeneration1.root", "RECREATE");
      auto tree = new TTree("tree", "A test tree");
      ROOT::TIOFeatures features;
      features.Set(ROOT::Experimental::EIOFeatures::kGenerateOffsetMap);
      tree->SetIOFeatures(features);
      tree->SetBit(TTree::kOnlyFlushAtCluster);
      tree->SetAutoFlush(10);
      Int_t sample[10];
      Int_t elem = 1;
      tree->Branch("elem", &elem, "elem/I");
      tree->Branch("sample", &sample, "sample[elem]/I");

      for (Int_t ev = 0; ev < fEventCount; ev++) {
         sample[0] = random->Gaus(100, 7);
         tree->Fill();
      }
      file->Write();
      delete tree;
      delete file;

      file = new TFile("TOffsetGeneration2.root", "RECREATE");
      tree = new TTree("tree", "A test tree");
      tree->SetBit(TTree::kOnlyFlushAtCluster);
      tree->SetAutoFlush(10);
      tree->Branch("elem", &elem, "elem/I");
      tree->Branch("sample", &sample, "sample[elem]/I");

      for (Int_t ev = 0; ev < fEventCount; ev++) {
         sample[0] = random->Gaus(100, 7);
         tree->Fill();
      }
      file->Write();

      file = new TFile("TOffsetGeneration3.root", "RECREATE");
      tree = new TTree("tree", "A test tree");
      tree->SetBit(TTree::kOnlyFlushAtCluster);
      tree->SetAutoFlush(5000);
      ElementStruct sample2;
      sample2.i = 1;
      double d[10];
      sample2.d = d;
      tree->Branch("sample", &sample2, 32*1024, 99);

      for (Int_t ev = 0; ev < fEventCount; ev++) {
         sample2.d[0] = random->Gaus(100, 7);
         tree->Fill();
      }
      file->Write();
      file->Close();
      delete file;

      file = new TFile("TOffsetGeneration4.root", "RECREATE");
      auto tree2 = new TTree("tree2", "A test tree");
      features.Set(ROOT::Experimental::EIOFeatures::kGenerateOffsetMap);
      tree2->SetIOFeatures(features);
      tree2->SetBit(TTree::kOnlyFlushAtCluster);
      tree2->SetAutoFlush(5000);
      sample2.i = 1;
      sample2.d = d;
      tree2->Branch("sample", &sample2, 32*1024, 99);

      for (Int_t ev = 0; ev < fEventCount; ev++) {
         sample2.d[0] = random->Gaus(100, 7);
         tree2->Fill();
      }
      file->Write();

      delete tree2;
      delete file;

      file = new TFile("TOffsetGeneration5.root", "RECREATE");
      tree = new TTree("tree", "A test tree");
      tree2 = new TTree("tree2", "A test tree");
      features.Set(ROOT::Experimental::EIOFeatures::kGenerateOffsetMap);
      tree->SetIOFeatures(features);
      tree->SetBit(TTree::kOnlyFlushAtCluster);
      tree2->SetBit(TTree::kOnlyFlushAtCluster);
      tree->SetAutoFlush(5000);
      tree2->SetAutoFlush(5000);
      sample2.i = 1;
      sample2.d = d;
      auto br = tree->Branch("sample", &sample2, 32*1024, 99);
      br->SetAutoDelete(kFALSE);
      auto br2 = tree2->Branch("sample", &sample2, 32*1024, 99);
      br2->SetAutoDelete(kFALSE);
      tree->Branch("elem", &elem, "elem/I");
      tree->Branch("sample2", &sample, "sample2[elem]/I");

      for (Int_t ev = 0; ev < 10; ev++) {
         sample2.i = ev;
         elem = ev;
         for (Int_t idx = 0; idx < ev; idx++) {
            sample2.d[idx] = random->Gaus(100, 7);
            sample[idx] = random->Gaus(100, 7);
         }
         tree->Fill();
         tree2->Fill();
      }
      file->Write();
      delete random;
      delete file;
   }
};

TEST_F(TOffsetGeneration, offsetArrayValues)
{
   std::unique_ptr<TFile> file(new TFile("TOffsetGeneration5.root"));
   auto tree = static_cast<TTree *>(file->Get("tree"));
   auto br = tree->GetBranch("d");
   auto basket = br->GetBasket(0);
   auto tree2 = static_cast<TTree *>(file->Get("tree2"));
   auto br2 = tree2->GetBranch("d");
   auto basket2 = br2->GetBasket(0);
   auto br3 = tree->GetBranch("sample2");
   auto basket3 = br3->GetBasket(0);
   Int_t *offsetArray = basket->GetEntryOffset();
   Int_t *offsetArray2 = basket2->GetEntryOffset();
   Int_t *offsetArray3 = basket3->GetEntryOffset();

   Int_t lastOffset = offsetArray[0];
   Int_t lastOffsetPrim = offsetArray3[0];
   for (Int_t idx = 0; idx < 10; idx++) {
      //printf("Event #%d; generated offset: %d, right offset: %d, primitive offset: %d\n", idx, offsetArray[idx], offsetArray2[idx], offsetArray3[idx]);
      Int_t curOffset = offsetArray[idx];
      Int_t curOffsetPrim = offsetArray3[idx];
      if (idx) {
         ASSERT_EQ(curOffset - lastOffset, (idx - 1) * 8 + 1);
         ASSERT_EQ(curOffsetPrim - lastOffsetPrim, (idx - 1) * 4);
      }
      lastOffset = offsetArray[idx];
      lastOffsetPrim = offsetArray3[idx];
      ASSERT_EQ(offsetArray[idx], offsetArray2[idx]);
   }
}

TEST_F(TOffsetGeneration, primitiveTest)
{
   std::unique_ptr<TFile> file(new TFile("TOffsetGeneration1.root"));
   auto tree = static_cast<TTree *>(file->Get("tree"));
   auto br = tree->GetBranch("sample");
   ASSERT_TRUE(br->GetTotalSize() < fEventCount * 14);

   file.reset(new TFile("TOffsetGeneration2.root"));
   tree = static_cast<TTree *>(file->Get("tree"));
   br = tree->GetBranch("sample");
   ASSERT_TRUE(br->GetTotalSize() > fEventCount * 14);
}

TEST_F(TOffsetGeneration, elementsTest)
{
   std::unique_ptr<TFile> file(new TFile("TOffsetGeneration3.root"));
   auto tree = static_cast<TTree *>(file->Get("tree"));
   auto br = tree->GetBranch("d");
   ASSERT_TRUE(br->GetTotalSize() > fEventCount * 10);

   file.reset(new TFile("TOffsetGeneration4.root"));
   tree = static_cast<TTree *>(file->Get("tree2"));
   br = tree->GetBranch("d");
   TClass *expectedClass = nullptr;
   EDataType expectedType;
   ASSERT_FALSE(br->GetExpectedType(expectedClass, expectedType));

   // Verifies splitting is working correctly.
   TLeaf *leaf = static_cast<TBranchElement*>(br)->FindLeaf("d");
   ASSERT_TRUE(leaf);

   ASSERT_TRUE(br->GetTotalSize() < fEventCount * 10);
}
