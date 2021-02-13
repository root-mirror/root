#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom.h"

#include "gtest/gtest.h"

class TBranchTest : public ::testing::Test {
protected:
   virtual void SetUp()
   {
      TRandom *random = new TRandom(837);
      TFile *file = new TFile("TBranchTestTree.root", "RECREATE");
      TTree *tree = new TTree("tree", "A test tree");
      TTree *tree2 = new TTree("tree2", "A test tree");
      tree->SetAutoSave(10);
      Float_t data = 0;
      tree->Branch("branch", &data);
      tree2->Branch("branch", &data);

      for (Int_t ev = 0; ev < 100; ev++) {
         data = random->Gaus(100, 7);
         tree->Fill();
         tree2->Fill();
         if (ev % 10 == 7) {
            tree->FlushBaskets(false);
         }
         if (ev % 10 == 6) {
            tree2->FlushBaskets();
         }
      }
      file->Write();
      delete random;
      delete tree;
      delete tree2;
      delete file;
   }
};

TEST_F(TBranchTest, clusterIteratorTest)
{
   std::unique_ptr<TFile> file(new TFile("TBranchTestTree.root"));
   TTree *tree = (TTree *)file->Get("tree2");

   auto iter = tree->GetClusterIterator(0);
   Long64_t last_start = 7;
   auto current = iter.Next();
   ASSERT_TRUE(current == 0);
   while (current < tree->GetEntries()) {
      current = iter.Next();
      ASSERT_TRUE(last_start == current);
      last_start += 10;
      last_start = std::min(last_start, tree->GetEntries());
   }
   ASSERT_TRUE(current == tree->GetEntries());
}

TEST_F(TBranchTest, nonePreviousTest)
{
   TFile *file = new TFile("TBranchTestTree.root");
   TTree *tree = (TTree *)file->Get("tree");
   TBranch *branch = tree->GetBranch("branch");

   tree->SetClusterPrefetch(false);
   tree->SetMaxVirtualSize(0);

   // Checks for normal behavior when change is
   // not being used
   branch->GetEntry(0);
   ASSERT_TRUE(branch->GetListOfBaskets()->At(0));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(1));

   branch->GetEntry(10);
   ASSERT_FALSE(branch->GetListOfBaskets()->At(0));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(1));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(2));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(3));
   delete file;
}

TEST_F(TBranchTest, onePreviousTest)
{
   TFile *file = new TFile("TBranchTestTree.root");
   TTree *tree = (TTree *)file->Get("tree");
   TBranch *branch = tree->GetBranch("branch");

   // Checks to make sure only first basket is loaded
   branch->GetEntry(0);
   ASSERT_TRUE(branch->GetListOfBaskets()->At(0));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(1));

   tree->SetClusterPrefetch(true);

   // Checks to make sure the whole cluster is loaded
   branch->GetEntry(10);
   ASSERT_FALSE(branch->GetListOfBaskets()->At(0));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(1));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(2));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(3));

   // Checks to make sure clusters are being removed
   // from memory, with none retained
   branch->GetEntry(20);
   ASSERT_FALSE(branch->GetListOfBaskets()->At(2));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(3));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(4));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(5));

   tree->SetMaxVirtualSize(-1);

   // Checks to make sure previous cluster is retained in memory
   branch->GetEntry(30);
   ASSERT_TRUE(branch->GetListOfBaskets()->At(4));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(5));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(6));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(7));

   // Checks to make sure clusters are being removed
   // from memory
   branch->GetEntry(40);
   ASSERT_FALSE(branch->GetListOfBaskets()->At(4));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(5));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(6));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(7));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(8));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(9));
   delete file;
}

TEST_F(TBranchTest, twoPreviousTest)
{
   TFile *file = new TFile("TBranchTestTree.root");
   TTree *tree = (TTree *)file->Get("tree");
   TBranch *branch = tree->GetBranch("branch");

   tree->SetMaxVirtualSize(-2);

   // Checks to make sure only first basket is loaded
   branch->GetEntry(0);
   ASSERT_TRUE(branch->GetListOfBaskets()->At(0));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(1));

   // Checks to make sure previous cluster is retained in memory
   branch->GetEntry(10);
   ASSERT_TRUE(branch->GetListOfBaskets()->At(0));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(1));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(2));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(3));

   // Checks to make sure previous cluster is retained in memory
   branch->GetEntry(19);
   ASSERT_TRUE(branch->GetListOfBaskets()->At(0));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(1));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(2));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(3));

   tree->SetClusterPrefetch(true);

   // Checks to make sure previous cluster is retained in memory
   branch->GetEntry(20);
   ASSERT_TRUE(branch->GetListOfBaskets()->At(0));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(1));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(2));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(3));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(4));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(5));

   // Checks to make sure clusters are being removed
   // from memory
   branch->GetEntry(30);
   ASSERT_FALSE(branch->GetListOfBaskets()->At(0));
   ASSERT_FALSE(branch->GetListOfBaskets()->At(1));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(2));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(3));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(4));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(5));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(6));
   ASSERT_TRUE(branch->GetListOfBaskets()->At(7));
   delete file;
}

bool nocomp(int mode = 0)
{
   TString filename("f_compression_");
   filename += mode;
   filename += ".root";

   // create a new file
   TFile *f_write = TFile::Open(filename, "recreate", "", 101);

   if (mode == 1)
      f_write->cd();                 // Create regular TTree.
   TTree *t = new TTree("t", "t"); // Create a new DISK RESIDENT tree (when mode == 0)
   if (mode == 2)
      t->SetDirectory(f_write);      // Attach after creation of TTree.

   // fill the tree
   Float_t v;
   t->Branch("v", &v, "v/F");
   if (mode == 3)
      t->SetDirectory(f_write);      // Attach after creation of TBranch.
   for (Int_t i = 0; i < 999999; i++)
   {
      v = gRandom->Gaus(0., 1.);
      t->Fill();
   }

   if (mode == 4)
      t->SetDirectory(f_write);      // Attach after filling of TBranch.

   std::cout << "Branch compression: " << t->GetBranch("v")->GetCompressionLevel() << '\n';

   // save the tree to a file
   f_write->cd();
   t->Write();
   delete f_write; // automatically deletes "t", too

   TFile *f_read = TFile::Open(filename, "READ");
   if (!f_read || f_read->IsZombie())
      return false;
   t = f_read->Get<TTree>("t");
   if (!t)
      return false;

   return (t->GetZipBytes() != t->GetTotBytes());
}

TEST_F(TBranchTest, branchInheritsCompression)
{
   for(int mode = 4; mode >= 0; --mode)
      ASSERT_TRUE(nocomp(mode)) << "Failed for mode: " << mode;
}