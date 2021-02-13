#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"

#include "gtest/gtest.h"

class TChainSaveTest : public ::testing::Test {
public:
   static const Int_t nev1 = 100, nev2 = 100;

protected:
   void SetUp() override
   {
      TFile file1("chaintree1.root", "RECREATE");
      TTree tree1("testtree", "A test tree");

      TFile file2("chaintree2.root", "RECREATE");
      TTree tree2("testtree", "A test tree");

      TFile filefr("chaintree-fr.root", "RECREATE");
      TTree treefr("testtree_fr", "A friend test tree");

      Int_t x = 0, y = 0;
      tree1.Branch("x", &x);
      tree2.Branch("x", &x);
      treefr.Branch("y", &y);

      for (Int_t ev = 0; ev < nev1; ev++) {
         x = ev;
         y = x % 2;
         tree1.Fill();
         treefr.Fill();
      }
      file1.Write();

      for (Int_t ev = 0; ev < nev2; ev++) {
         x = ev + nev1;
         y = x % 2;
         tree2.Fill();
         treefr.Fill();
      }
      file2.Write();
      filefr.Write();

      TChain chain("testtree", "Test chain");
      chain.AddFile("chaintree1.root");
      chain.AddFile("chaintree2.root");
      chain.AddFriend("testtree_fr", "chaintree-fr.root");
      chain.GetEntries();
      chain.SetMarkerStyle(7);
      chain.SetMarkerColor(4);
      
      chain.SaveAs("testchain.C", "friend");
   }
   
   void TearDown() override
   {
      gSystem->Unlink("chaintree1.root");
      gSystem->Unlink("chaintree2.root");
      gSystem->Unlink("chaintree-fr.root");
      gSystem->Unlink("testchain.C");
   }
   
};

// Read TChain from a C++ file and check entries
TEST_F(TChainSaveTest, readTest)
{
   gROOT->Macro("testchain.C");

   TChain *pchain = dynamic_cast<TChain *>(gROOT->GetListOfDataSets()->FindObject("testtree"));

   ASSERT_NE(pchain, nullptr);
   EXPECT_EQ(pchain->GetEntries(), TChainSaveTest::nev1 + TChainSaveTest::nev2);

   // Check marker properties
   EXPECT_EQ(pchain->GetMarkerStyle(), 7);
   EXPECT_EQ(pchain->GetMarkerColor(), 4);

   // Check values
   Int_t nev = pchain->Draw("y:x", "", "goff");
   Double_t *x = pchain->GetV2(), *y = pchain->GetV1();
   ASSERT_NE(x, nullptr);
   ASSERT_NE(y, nullptr);

   for (Int_t i = 0; i < nev; i++) {
      EXPECT_EQ((Int_t)x[i], i);
      EXPECT_EQ((Int_t)y[i], i % 2);
   }
}
