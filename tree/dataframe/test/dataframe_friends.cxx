#include "ROOT/RDataFrame.hxx"
#include "ROOT/TSeq.hxx"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "gtest/gtest.h"

#include <algorithm> // std::equal

// fixture that creates two files with two trees of 10 events each. One has branch `x`, the other branch `y`, both ints.
class RDFAndFriends : public ::testing::Test {
protected:
   constexpr static auto kFile1 = "test_tdfandfriends.root";
   constexpr static auto kFile2 = "test_tdfandfriends2.root";
   constexpr static auto kFile3 = "test_tdfandfriends3.root";
   constexpr static auto kFile4 = "test_tdfandfriends4.root";
   constexpr static auto kFile5 = "test_tdfandfriends5.root";
   constexpr static ULong64_t kSizeSmall = 4;
   constexpr static ULong64_t kSizeBig = 10000;
   static void SetUpTestCase()
   {
      ROOT::RDataFrame d(kSizeSmall);
      d.Define("x", [] { return 1; }).Snapshot<int>("t", kFile1, {"x"});
      d.Define("y", [] { return 2; }).Snapshot<int>("t2", kFile2, {"y"});

      TFile f(kFile3, "RECREATE");
      TTree t("t3", "t3");
      float arr[4];
      t.Branch("arr", arr, "arr[4]/F");
      for (auto i : ROOT::TSeqU(4)) {
         for (auto j : ROOT::TSeqU(4)) {
            arr[j] = i + j;
         }
         t.Fill();
      }
      t.Write();

      ROOT::RDataFrame d2(kSizeBig);
      d2.Define("x", [] { return 4; }).Snapshot<int>("t", kFile4, {"x"});
      d2.Define("y", [] { return 5; }).Snapshot<int>("t2", kFile5, {"y"});
   }

   static void TearDownTestCase()
   {
      for (auto fileName : {kFile1, kFile2, kFile3, kFile4, kFile5})
         gSystem->Unlink(fileName);
   }
};

TEST_F(RDFAndFriends, FriendByFile)
{
   TFile f1(kFile1);
   auto t1 = f1.Get<TTree>("t");
   t1->AddFriend("t2", kFile2);
   ROOT::RDataFrame d(*t1);
   auto x = d.Min<int>("x");
   auto t = d.Take<int>("y");
   EXPECT_EQ(*x, 1);
   for (auto v : t)
      EXPECT_EQ(v, 2);
}

TEST_F(RDFAndFriends, FriendByPointer)
{
   TFile f1(kFile1);
   auto t1 = f1.Get<TTree>("t");
   TFile f2(kFile2);
   auto t2 = f2.Get<TTree>("t2");
   t1->AddFriend(t2);
   ROOT::RDataFrame d(*t1);
   auto x = d.Min<int>("x");
   auto t = d.Take<int>("y");
   EXPECT_EQ(*x, 1);
   for (auto v : t)
      EXPECT_EQ(v, 2);
}

TEST_F(RDFAndFriends, FriendArrayByFile)
{
   TFile f1(kFile1);
   auto t1 = f1.Get<TTree>("t");
   t1->AddFriend("t3", kFile3);
   ROOT::RDataFrame d(*t1);

   int i(0);
   auto checkArr = [&i](ROOT::VecOps::RVec<float> av) {
      auto ifloat = float(i);
      EXPECT_EQ(ifloat, av[0]);
      EXPECT_EQ(ifloat + 1, av[1]);
      EXPECT_EQ(ifloat + 2, av[2]);
      EXPECT_EQ(ifloat + 3, av[3]);
      i++;
   };
   d.Foreach(checkArr, {"arr"});
}

TEST_F(RDFAndFriends, FriendArrayByPointer)
{
   TFile f1(kFile1);
   auto t1 = f1.Get<TTree>("t");
   TFile f3(kFile3);
   auto t3 = f3.Get<TTree>("t3");
   t1->AddFriend(t3);
   ROOT::RDataFrame d(*t1);

   int i(0);
   auto checkArr = [&i](ROOT::VecOps::RVec<float> av) {
      auto ifloat = float(i);
      EXPECT_EQ(ifloat, av[0]);
      EXPECT_EQ(ifloat + 1, av[1]);
      EXPECT_EQ(ifloat + 2, av[2]);
      EXPECT_EQ(ifloat + 3, av[3]);
      i++;
   };
   d.Foreach(checkArr, {"arr"});
}

TEST_F(RDFAndFriends, QualifiedBranchName)
{
   TFile f1(kFile1);
   auto t1 = f1.Get<TTree>("t");
   t1->AddFriend("t2", kFile2);
   ROOT::RDataFrame d(*t1);
   auto x = d.Min<int>("x");
   EXPECT_EQ(*x, 1);
   auto t = d.Take<int>("t2.y");
   for (auto v : t)
      EXPECT_EQ(v, 2);
}

TEST_F(RDFAndFriends, FromDefine)
{
   TFile f1(kFile1);
   auto t1 = f1.Get<TTree>("t");
   t1->AddFriend("t2", kFile2);
   ROOT::RDataFrame d(*t1);

   auto m = d.Define("yy", [](int y) { return y * y; }, {"y"}).Mean("yy");
   EXPECT_DOUBLE_EQ(*m, 4.);
}

TEST_F(RDFAndFriends, FromJittedDefine)
{
   TFile f1(kFile1);
   auto t1 = f1.Get<TTree>("t");
   t1->AddFriend("t2", kFile2);
   ROOT::RDataFrame d(*t1);

   auto m = d.Define("yy", "y * y").Mean("yy");
   EXPECT_DOUBLE_EQ(*m, 4.);
}

// NOW MT!-------------
#ifdef R__USE_IMT

TEST_F(RDFAndFriends, FriendMT)
{
   ROOT::EnableImplicitMT(4u);

   TFile f1(kFile4);
   auto t1 = f1.Get<TTree>("t");
   t1->AddFriend("t2", kFile5);
   ROOT::RDataFrame d(*t1);
   auto x = d.Min<int>("x");
   auto t = d.Take<int>("y");
   EXPECT_EQ(*x, 4);
   for (auto v : t)
      EXPECT_EQ(v, 5);
   ROOT::DisableImplicitMT();
}

TEST_F(RDFAndFriends, FriendAliasMT)
{
   ROOT::EnableImplicitMT(4u);
   TFile f1(kFile1);
   auto t1 = f1.Get<TTree>("t");
   TFile f2(kFile4);
   auto t2 = f2.Get<TTree>("t");
   t1->AddFriend(t2, "myfriend");
   ROOT::RDataFrame d(*t1);
   auto x = d.Min<int>("x");
   auto t = d.Take<int>("myfriend.x");
   EXPECT_EQ(*x, 1);
   for (auto v : t)
      EXPECT_EQ(v, 4);
   ROOT::DisableImplicitMT();
}

TEST_F(RDFAndFriends, FriendChainMT)
{
   ROOT::EnableImplicitMT(4u);
   TChain c1("t");
   c1.AddFile(kFile1);
   c1.AddFile(kFile4);
   c1.AddFile(kFile1);
   c1.AddFile(kFile4);
   TChain c2("t2");
   c2.AddFile(kFile2);
   c2.AddFile(kFile5);
   c2.AddFile(kFile2);
   c2.AddFile(kFile5);
   c1.AddFriend(&c2);

   ROOT::RDataFrame d(c1);
   auto c = d.Count();
   EXPECT_EQ(*c, 2 * (kSizeSmall + kSizeBig));
   auto x = d.Min<int>("x");
   auto y = d.Max<int>("y");
   EXPECT_EQ(*x, 1);
   EXPECT_EQ(*y, 5);
   ROOT::DisableImplicitMT();
}

// ROOT-9559
void FillIndexedFriend(const char *mainfile, const char *auxfile)
{
   // Start by creating main Tree
   TFile f(mainfile, "RECREATE");
   TTree mainTree("mainTree", "mainTree");
   int idx;
   mainTree.Branch("idx", &idx);
   int x;
   mainTree.Branch("x", &x);

   idx = 1;
   x = 1;
   mainTree.Fill();
   idx = 1;
   x = 2;
   mainTree.Fill();
   idx = 1;
   x = 3;
   mainTree.Fill();
   idx = 2;
   x = 4;
   mainTree.Fill();
   idx = 2;
   x = 5;
   mainTree.Fill();
   mainTree.Write();
   f.Close();

   // And aux tree
   TFile f2(auxfile, "RECREATE");
   TTree auxTree("auxTree", "auxTree");
   auxTree.Branch("idx", &idx);
   int y;
   auxTree.Branch("y", &y);
   idx = 2;
   y = 5;
   auxTree.Fill();
   idx = 1;
   y = 7;
   auxTree.Fill();
   auxTree.Write();
   f2.Close();
}

TEST(RDFAndFriendsNoFixture, IndexedFriend)
{
   auto mainFile = "IndexedFriend_main.root";
   auto auxFile = "IndexedFriend_aux.root";
   FillIndexedFriend(mainFile, auxFile);

   TChain mainChain("mainTree", "mainTree");
   mainChain.Add(mainFile);
   TChain auxChain("auxTree", "auxTree");
   auxChain.Add(auxFile);

   auxChain.BuildIndex("idx");
   mainChain.AddFriend(&auxChain);

   auto df = ROOT::RDataFrame(mainChain);
   auto x = df.Take<int>("x");
   auto y = df.Take<int>("auxTree.y");

   std::vector<int> refx{{1,2,3,4,5}};
   EXPECT_TRUE(std::equal(x->begin(), x->end(), refx.begin()));
   std::vector<int> refy{{7,7,7,5,5}};
   EXPECT_TRUE(std::equal(y->begin(), y->end(), refy.begin()));

   gSystem->Unlink(mainFile);
   gSystem->Unlink(auxFile);
}

// Test for https://github.com/root-project/root/issues/6741
TEST(RDFAndFriendsNoFixture, AutomaticFriendsLoad)
{
   const auto fname = "rdf_automaticfriendsloadtest.root";
   {
      // write a TTree and its friend to the same file
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
      f.Close();
   }
   EXPECT_EQ(ROOT::RDataFrame("t1", fname).Max<int>("t2.x").GetValue(), 42);

   gSystem->Unlink(fname);
}

#endif // R__USE_IMT
