#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "gtest/gtest.h"

#include <string>


// ROOT-10702
TEST(TTreeReaderRegressions, CompositeTypeWithNameClash)
{
   struct Int { int x; };
	gInterpreter->Declare("struct Int { int x; };");

   const auto fname = "ttreereader_compositetypewithnameclash.root";

   {
      TFile f(fname, "recreate");
      Int i{-1};
      int x = 1;
      TTree t("t", "t");
      const auto toJit = "((TTree*)" + std::to_string(reinterpret_cast<std::size_t>(&t)) + ")->Branch(\"i\", (Int*)" +
                         std::to_string(reinterpret_cast<std::size_t>(&i)) + ");";
      gInterpreter->ProcessLine(toJit.c_str());
      t.Branch("x", &x);
      t.Fill();
      t.Write();
      f.Close();
   }

   TFile f(fname);
   TTreeReader r("t", &f);
   TTreeReaderValue<int> iv(r, "i.x");
   TTreeReaderValue<int> xv(r, "x");
   r.Next();
   EXPECT_EQ(xv.GetSetupStatus(), 0);
   if (xv.GetSetupStatus() == 0) {
      EXPECT_EQ(*xv, 1);
   }
   EXPECT_EQ(iv.GetSetupStatus(), 0);
   if (iv.GetSetupStatus() == 0) {
      EXPECT_EQ(*iv, -1);
   }

   gSystem->Unlink(fname);
}

// Regression test for https://github.com/root-project/root/issues/6993
TEST(TTreeReaderRegressions, AutoloadedFriends)
{
   const auto fname = "treereaderautoloadedfriends.root";
   {
      // write a TTree and its friend to the same file:
      // when t1 is read back, it automatically also loads its friend
      TFile f(fname, "recreate");
      TTree t1("t1", "t1");
      TTree t2("t2", "t2");
      int x = 42;
      t2.Branch("x", &x);
      t1.Fill();
      t2.Fill();
      t1.AddFriend(&t2);
      t1.Write();
      t2.Write();
   }

   // reading t2.x via TTreeReader segfaults
   TChain c("t1");
   c.Add(fname);
   c.LoadTree(0);
   TTreeReader r(&c);
   TTreeReaderValue<int> rv(r, "t2.x");
   ASSERT_TRUE(r.Next());
   EXPECT_EQ(*rv, 42);
   EXPECT_FALSE(r.Next());

   gSystem->Unlink(fname);
}

// ROOT-10824
TEST(TTreeReaderRegressions, IndexedFriend)
{
   const auto fname = "treereader_fillindexedfriend.root";

   {
      TFile f(fname, "recreate");
      // Create main tree
      TTree mainTree("mainTree", "mainTree");
      int idx;
      mainTree.Branch("idx", &idx);
      float x;
      mainTree.Branch("x", &x);

      idx = 1;
      x = 1.f;
      mainTree.Fill();
      idx = 1;
      x = 2.f;
      mainTree.Fill();
      idx = 2;
      x = 10.f;
      mainTree.Fill();
      idx = 2;
      x = 20.f;
      mainTree.Fill();
      mainTree.Write();

      // Create aux tree
      TTree auxTree("auxTree", "auxTree");
      auxTree.Branch("idx", &idx);
      std::string s;
      auxTree.Branch("s", &s);
      idx = 1;
      s = "small";
      auxTree.Fill();
      idx = 2;
      s = "big";
      auxTree.Fill();
      auxTree.Write();
      f.Close();
   }

   auto checkTreeReader = [](TTreeReader &r, TTreeReaderValue<float> &rx, TTreeReaderValue<std::string> &rs) {
      ASSERT_TRUE(r.Next());
      EXPECT_EQ(*rx, 1.f);
      EXPECT_EQ(*rs, "small");
      ASSERT_TRUE(r.Next());
      EXPECT_EQ(*rx, 2.f);
      EXPECT_EQ(*rs, "small");
      ASSERT_TRUE(r.Next());
      EXPECT_EQ(*rx, 10.f);
      EXPECT_EQ(*rs, "big");
      ASSERT_TRUE(r.Next());
      EXPECT_EQ(*rx, 20.f);
      EXPECT_EQ(*rs, "big");
      ASSERT_FALSE(r.Next());
   };

   // Test reading back with TTreeReader+TTree
   {
      TFile f(fname);
      auto mainTree = f.Get<TTree>("mainTree");
      auto auxTree = f.Get<TTree>("auxTree");

      auxTree->BuildIndex("idx");
      mainTree->AddFriend(auxTree);

      TTreeReader r(mainTree);
      TTreeReaderValue<float> rx(r, "x");
      TTreeReaderValue<std::string> rs(r, "auxTree.s");
      checkTreeReader(r, rx, rs);
   }

   // Test reading back with TTreeReader+TChain
   {
      TChain mainChain("mainTree", "mainTree");
      mainChain.Add(fname);
      TChain auxChain("auxTree", "auxTree");
      auxChain.Add(fname);

      auxChain.BuildIndex("idx");
      mainChain.AddFriend(&auxChain);

      TTreeReader r(&mainChain);
      TTreeReaderValue<float> rx(r, "x");
      TTreeReaderValue<std::string> rs(r, "auxTree.s");
      checkTreeReader(r, rx, rs);
   }

   gSystem->Unlink(fname);
}
