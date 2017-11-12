#include "ROOT/TDataFrame.hxx"
#include "ROOT/TSeq.hxx"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "gtest/gtest.h"

using namespace ROOT::Experimental;

// fixture that creates two files with two trees of 10 events each. One has branch `x`, the other branch `y`, both ints.
class TDFAndFriends : public ::testing::Test {
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
      TDataFrame d(kSizeSmall);
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

      TDataFrame d2(kSizeBig);
      d2.Define("x", [] { return 4; }).Snapshot<int>("t", kFile4, {"x"});
      d2.Define("y", [] { return 5; }).Snapshot<int>("t2", kFile5, {"y"});
   }

   static void TearDownTestCase()
   {
      for (auto fileName : {kFile1, kFile2, kFile3})
         gSystem->Unlink(fileName);
   }
};

TEST_F(TDFAndFriends, FriendByFile)
{
   TFile f1(kFile1);
   TTree *t1 = static_cast<TTree *>(f1.Get("t"));
   t1->AddFriend("t2", kFile2);
   TDataFrame d(*t1);
   auto x = d.Min<int>("x");
   auto t = d.Take<int>("y");
   EXPECT_EQ(*x, 1);
   for (auto v : t)
      EXPECT_EQ(v, 2);
}

TEST_F(TDFAndFriends, FriendByPointer)
{
   TFile f1(kFile1);
   TTree *t1 = static_cast<TTree *>(f1.Get("t"));
   TFile f2(kFile2);
   TTree *t2 = static_cast<TTree *>(f2.Get("t2"));
   t1->AddFriend(t2);
   TDataFrame d(*t1);
   auto x = d.Min<int>("x");
   auto t = d.Take<int>("y");
   EXPECT_EQ(*x, 1);
   for (auto v : t)
      EXPECT_EQ(v, 2);
}

TEST_F(TDFAndFriends, FriendArrayByFile)
{
   TFile f1(kFile1);
   TTree *t1 = static_cast<TTree *>(f1.Get("t"));
   t1->AddFriend("t3", kFile3);
   TDataFrame d(*t1);

   int i(0);
   auto checkArr = [&i](TDF::TArrayBranch<float> av) {
      auto ifloat = float(i);
      EXPECT_EQ(ifloat, av[0]);
      EXPECT_EQ(ifloat + 1, av[1]);
      EXPECT_EQ(ifloat + 2, av[2]);
      EXPECT_EQ(ifloat + 3, av[3]);
      i++;
   };
   d.Foreach(checkArr, {"arr"});
}

TEST_F(TDFAndFriends, FriendArrayByPointer)
{
   TFile f1(kFile1);
   TTree *t1 = static_cast<TTree *>(f1.Get("t"));
   TFile f3(kFile3);
   TTree *t3 = static_cast<TTree *>(f3.Get("t3"));
   t1->AddFriend(t3);
   TDataFrame d(*t1);

   int i(0);
   auto checkArr = [&i](TDF::TArrayBranch<float> av) {
      auto ifloat = float(i);
      EXPECT_EQ(ifloat, av[0]);
      EXPECT_EQ(ifloat + 1, av[1]);
      EXPECT_EQ(ifloat + 2, av[2]);
      EXPECT_EQ(ifloat + 3, av[3]);
      i++;
   };
   d.Foreach(checkArr, {"arr"});
}

TEST_F(TDFAndFriends, QualifiedBranchName)
{
   TFile f1(kFile1);
   TTree *t1 = static_cast<TTree *>(f1.Get("t"));
   t1->AddFriend("t2", kFile2);
   TDataFrame d(*t1);
   auto x = d.Min<int>("x");
   EXPECT_EQ(*x, 1);
   auto t = d.Take<int>("t2.y");
   for (auto v : t)
      EXPECT_EQ(v, 2);
}

TEST_F(TDFAndFriends, FromDefine)
{
   TFile f1(kFile1);
   TTree *t1 = static_cast<TTree *>(f1.Get("t"));
   t1->AddFriend("t2", kFile2);
   TDataFrame d(*t1);

   auto m = d.Define("yy", [](int y) { return y * y; }, {"y"}).Mean("yy");
   EXPECT_DOUBLE_EQ(*m, 4.);
}

TEST_F(TDFAndFriends, FromJittedDefine)
{
   TFile f1(kFile1);
   TTree *t1 = static_cast<TTree *>(f1.Get("t"));
   t1->AddFriend("t2", kFile2);
   TDataFrame d(*t1);

   auto m = d.Define("yy", "y * y").Mean("yy");
   EXPECT_DOUBLE_EQ(*m, 4.);
}

// NOW MT!-------------
#ifdef R__USE_IMT

TEST_F(TDFAndFriends, FriendMT)
{
   auto nSlots = 4U;
   ROOT::EnableImplicitMT(nSlots);

   TFile f1(kFile4);
   TTree *t1 = static_cast<TTree *>(f1.Get("t"));
   t1->AddFriend("t2", kFile5);
   TDataFrame d(*t1);
   auto x = d.Min<int>("x");
   auto t = d.Take<int>("y");
   EXPECT_EQ(*x, 4);
   for (auto v : t)
      EXPECT_EQ(v, 5);
}

TEST_F(TDFAndFriends, FriendAliasMT)
{
   TFile f1(kFile4);
   TTree *t1 = static_cast<TTree *>(f1.Get("t"));
   TFile f2(kFile4);
   TTree *t2 = static_cast<TTree *>(f2.Get("t"));
   t1->AddFriend(t2, "myfriend");
   TDataFrame d(*t1);
   auto x = d.Min<int>("x");
   auto t = d.Take<int>("myfriend.x");
   EXPECT_EQ(*x, 4);
   for (auto v : t)
      EXPECT_EQ(v, 4);
}

TEST_F(TDFAndFriends, FriendChainMT)
{
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

   TDataFrame d(c1);
   auto c = d.Count();
   EXPECT_EQ(*c, 2 * (kSizeSmall + kSizeBig));
   auto x = d.Min<int>("x");
   auto y = d.Max<int>("y");
   EXPECT_EQ(*x, 1);
   EXPECT_EQ(*y, 5);
}

#endif // R__USE_IMT
