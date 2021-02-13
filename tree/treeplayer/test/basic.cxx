#include "ROOT/RMakeUnique.hxx"
#include "TEntryListArray.h"

#include "TChain.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TSystem.h"

#include "gtest/gtest.h"
#include <stdlib.h>

#include "RErrorIgnoreRAII.hxx"

std::unique_ptr<TTree> MakeTree() {
   double x[3]{};
   float z = 0.;
   struct {
      unsigned int ny = 100;
      int y[100];
   } yData;
   std::string str;
   Double32_t Double32 = 12.;
   Float16_t Float16 = -12.;

   auto tree = std::make_unique<TTree>("T", "test tree");
   tree->Branch("one", &x, "x[3]/D");
   tree->Branch("two", &yData, "ny/i:y[ny]/I");
   tree->Branch("three", &z, "z");
   tree->Branch("str", &str);
   tree->Branch("d32", &Double32);
   tree->Branch("f16", &Float16);
   tree->Branch("0.2.0.energy", &z);

   x[1] = 42.;
   yData.ny = 42;
   yData.y[0] = 17;
   str = "first";
   tree->Fill();

   x[2] = 43.;
   yData.ny = 5;
   yData.y[4] = 7;
   str = "";
   tree->Fill();

   for (int entry = 2; entry < 20; ++entry) {
      z = entry * (1 - 2 * (entry % 2)); // +entry for even, -entry for odd
      str = std::string(entry, '*');
      tree->Fill();
   }

   tree->ResetBranchAddresses();

   return tree;
}

TEST(TTreeReaderBasic, Interfaces) {
   auto tree = MakeTree();

   TTreeReader tr(tree.get());
   TTreeReaderArray<double> x(tr, "one.x");
   TTreeReaderArray<int> y(tr, "two.y");
   TTreeReaderValue<unsigned int> ny(tr, "two.ny");
   TTreeReaderValue<std::string> nstr(tr, "str");

   // Before reading data:
   EXPECT_NE(tr.begin(), tr.end());
   //MISSING: EXPECT_EQ(2, tr.end() - tr.begin());
   EXPECT_EQ(-1, tr.GetCurrentEntry());
   EXPECT_EQ(20, tr.GetEntries(false));
   EXPECT_EQ(TTreeReader::kEntryNotLoaded, tr.GetEntryStatus());
   EXPECT_EQ(tree.get(), tr.GetTree());
   EXPECT_FALSE(tr.IsChain());

   EXPECT_EQ(nullptr, ny.GetAddress());
   EXPECT_FALSE(ny.IsValid());
   EXPECT_STREQ("two.ny", ny.GetBranchName());
   EXPECT_EQ(ROOT::Internal::TTreeReaderValueBase::kReadNothingYet, ny.GetReadStatus());
   EXPECT_EQ(ROOT::Internal::TTreeReaderValueBase::kSetupNotSetup, ny.GetSetupStatus());
   EXPECT_EQ(ROOT::Internal::TTreeReaderValueBase::kReadNothingYet, ny.ProxyRead());

   // Skip to second entry:
   EXPECT_TRUE(tr.Next());
   EXPECT_TRUE(tr.Next());

   EXPECT_NE(tr.begin(), tr.end());
   //MISSING: EXPECT_EQ(2, tr.end() - tr.begin());
   EXPECT_EQ(1, tr.GetCurrentEntry());
   EXPECT_EQ(20, tr.GetEntries(false));


   EXPECT_EQ(5u, *ny);
   EXPECT_NE(nullptr, ny.GetAddress());
   EXPECT_TRUE(ny.IsValid());
   EXPECT_EQ(ny.GetAddress(), ny.Get());
   EXPECT_STREQ("two", ny.GetBranchName());
   EXPECT_STREQ("ny", ny.GetLeaf()->GetName());
   EXPECT_EQ(ROOT::Internal::TTreeReaderValueBase::kReadSuccess, ny.GetReadStatus());
   EXPECT_EQ(ROOT::Internal::TTreeReaderValueBase::kSetupMatch, ny.GetSetupStatus());
   EXPECT_EQ(ROOT::Internal::TTreeReaderValueBase::kReadSuccess, ny.ProxyRead());

   EXPECT_EQ(3u, x.GetSize());
   EXPECT_EQ(5u, y.GetSize());
   EXPECT_DOUBLE_EQ(43., x[2]);
   EXPECT_EQ(7, y[4]);

   for (int entry = 2; entry < 20; ++entry)
      EXPECT_TRUE(tr.Next());

   EXPECT_FALSE(tr.Next());
   //FAILS: EXPECT_FALSE(ny.IsValid());
}


TEST(TTreeReaderBasic, ErrorProbing) {
   TTreeReader tr("doesNotExist", gROOT);
   EXPECT_EQ(TTreeReader::kEntryNoTree, tr.GetEntryStatus());
   EXPECT_EQ(nullptr, tr.GetTree());

   tr.SetTree((TTree*)nullptr);
   EXPECT_EQ(TTreeReader::kEntryNoTree, tr.GetEntryStatus());
   EXPECT_EQ(nullptr, tr.GetTree());

   TTreeReaderValue<double> val(tr, "branchDoesNotExist");
   EXPECT_FALSE(val.IsValid());
}


TEST(TTreeReaderBasic, Range) {
   auto tree = MakeTree();
   TTreeReader tr(tree.get());

   EXPECT_EQ(TTreeReader::kEntryValid, tr.SetEntriesRange(5, 8));

   EXPECT_TRUE(tr.Next());
   EXPECT_EQ(5, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, tr.GetEntryStatus());

   EXPECT_TRUE(tr.Next());
   EXPECT_EQ(6, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, tr.GetEntryStatus());

   EXPECT_TRUE(tr.Next());
   EXPECT_EQ(7, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, tr.GetEntryStatus());

   // Reached end:
   EXPECT_FALSE(tr.Next());
   EXPECT_EQ(8, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryBeyondEnd, tr.GetEntryStatus());

   // Read beyond end:
   EXPECT_FALSE(tr.Next());
   EXPECT_EQ(TTreeReader::kEntryBeyondEnd, tr.GetEntryStatus());
   EXPECT_EQ(9, tr.GetCurrentEntry());


   // Restart, now with different entries.
   tr.Restart();
   EXPECT_EQ(TTreeReader::kEntryValid, tr.SetEntriesRange(0, 2));

   EXPECT_TRUE(tr.Next());
   EXPECT_EQ(0, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, tr.GetEntryStatus());

   EXPECT_TRUE(tr.Next());
   EXPECT_EQ(1, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, tr.GetEntryStatus());

   // Reached end:
   EXPECT_FALSE(tr.Next());
   EXPECT_EQ(2, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryBeyondEnd, tr.GetEntryStatus());

}

TEST(TTreeReaderBasic, InvalidRange) {
   auto tree = MakeTree();
   TTreeReader tr(tree.get());

   {
      RErrorIgnoreRAII errIgnRAAII;
      EXPECT_EQ(TTreeReader::kEntryNotFound, tr.SetEntriesRange(tree->GetEntries(), 0));
   }

   // Is SetEntriesRange() simply ignored as it should be?
   EXPECT_TRUE(tr.Next());
   EXPECT_EQ(0, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, tr.GetEntryStatus());
}


TEST(TTreeReaderBasic, OneEntryRange) {
   auto tree = MakeTree();
   TTreeReader tr(tree.get());

   EXPECT_EQ(TTreeReader::kEntryValid, tr.SetEntriesRange(1, 2));

   EXPECT_TRUE(tr.Next());
   EXPECT_EQ(1, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, tr.GetEntryStatus());

   // Read beyond end:
   EXPECT_FALSE(tr.Next());
   EXPECT_EQ(TTreeReader::kEntryBeyondEnd, tr.GetEntryStatus());
   EXPECT_EQ(2, tr.GetCurrentEntry());

   for (int entry = 2; entry < tree->GetEntries(); ++entry)
      EXPECT_FALSE(tr.Next());
}


TEST(TTreeReaderBasic, ZeroEntryRange) {
   auto tree = MakeTree();
   TTreeReader tr(tree.get());

   // end is ignored:
   EXPECT_EQ(TTreeReader::kEntryValid, tr.SetEntriesRange(18, 18));

   EXPECT_TRUE(tr.Next());
   EXPECT_EQ(18, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, tr.GetEntryStatus());

   EXPECT_TRUE(tr.Next());
   EXPECT_EQ(19, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, tr.GetEntryStatus());

   // Read beyond end:
   EXPECT_FALSE(tr.Next());
   // As the TTree only has up to entry 19, 20 is kEntryNotFound:
   EXPECT_EQ(TTreeReader::kEntryNotFound, tr.GetEntryStatus());
   EXPECT_EQ(20, tr.GetCurrentEntry());
}


TEST(TTreeReaderBasic, InvertedEntryRange) {
   auto tree = MakeTree();
   TTreeReader tr(tree.get());

   EXPECT_EQ(TTreeReader::kEntryValid, tr.SetEntriesRange(18, 3));

   EXPECT_TRUE(tr.Next());
   EXPECT_EQ(18, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, tr.GetEntryStatus());

   EXPECT_TRUE(tr.Next());
   EXPECT_EQ(19, tr.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, tr.GetEntryStatus());

   // Read beyond end:
   EXPECT_FALSE(tr.Next());
   // As the TTree only has up to entry 19, 20 is kEntryNotFound:
   EXPECT_EQ(TTreeReader::kEntryNotFound, tr.GetEntryStatus());
   EXPECT_EQ(20, tr.GetCurrentEntry());
}


TEST(TTreeReaderBasic, EntryList) {
   // See https://root.cern.ch/phpBB3/viewtopic.php?f=3&t=22850&p=100796
   auto tree = MakeTree();
   EXPECT_EQ(9, tree->Draw(">>negZ","three.z<0", "entrylistarray"));
   TEntryListArray* selected = (TEntryListArray*)gDirectory->Get("negZ");
   TTreeReader aReader(tree.get(), selected);

   EXPECT_EQ(9, aReader.GetEntries(false));

   TTreeReaderValue<float> z(aReader, "three.z");

   int count = 0;
   while (aReader.Next()) {
      // Make sure all Next()-ed entries have z<0 as selected by the entry list.
      EXPECT_LT(*z, 0);
      ++count;
   }
   EXPECT_EQ(9, count);

   aReader.Restart();
   EXPECT_EQ(TTreeReader::kEntryValid, aReader.SetEntriesRange(0, 2));

   EXPECT_TRUE(aReader.Next());
   EXPECT_EQ(0, aReader.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, aReader.GetEntryStatus());

   EXPECT_EQ(9, aReader.GetEntries(false));
}

TEST(TTreeReaderBasic, EntryListBeyondEnd) {
   auto tree = MakeTree();
   TEntryList selected;
   // Add the last valid entry; a subsequent Next() must return false.
   selected.Enter(tree->GetEntries() - 1);

   TTreeReader aReader(tree.get(), &selected);

   EXPECT_EQ(1, aReader.GetEntries(false));

   EXPECT_TRUE(aReader.Next());
   EXPECT_EQ(0, aReader.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryValid, aReader.GetEntryStatus());

   EXPECT_FALSE(aReader.Next());
   EXPECT_EQ(1, aReader.GetCurrentEntry());
   EXPECT_EQ(TTreeReader::kEntryNotFound, aReader.GetEntryStatus());
}

// two files treereader_entrylists{1,2}.root with branch "x" and values {0,1} and {2,3}
struct InputFilesRAII {
   const char *fNames[2] = {"treereader_entrylists1.root", "treereader_entrylists2.root"};

   InputFilesRAII()
   {
      int x = 0;
      for (const auto fName : fNames) {
         TFile f(fName, "recreate");
         TTree t("t", "t");
         t.Branch("x", &x);
         t.Fill();
         ++x;
         t.Fill();
         ++x;
         t.Write();
      }
   }

   ~InputFilesRAII()
   {
      gSystem->Unlink(fNames[0]);
      gSystem->Unlink(fNames[1]);
   }
};

void TestWithEntryList(TEntryList &l, const InputFilesRAII &files)
{
   TChain c("t");
   c.Add(files.fNames[0]);
   c.Add(files.fNames[1]);
   if (l.GetLists()) // must associate with TChain if there are sub-entrylists
      c.SetEntryList(&l);

   TTreeReader r(&c, &l);
   TTreeReaderValue<int> xv(r, "x");
   const std::vector<int> values({0, 2});
   int i = 0;
   while (r.Next())
      EXPECT_EQ(*xv, values.at(i++));
}

TEST(TTreeReaderBasic, TChainWithGlobalEntryList)
{
   InputFilesRAII files;

   // if using a global list with a TChain, no need to specify treename and filename to the TEntryList (but harmless)
   TEntryList l;
   l.Enter(0);
   l.Enter(2);

   TestWithEntryList(l, files);
}

TEST(TTreeReaderBasic, TChainWithSubEntryLists)
{
   InputFilesRAII files;

   TEntryList l1("l1", "l1", "t", files.fNames[0]);
   l1.Enter(0);
   TEntryList l2("l2", "l2", "t", files.fNames[1]);
   l2.Enter(0);

   TEntryList l;
   l.Add(&l1);
   l.Add(&l2);

   TestWithEntryList(l, files);
}

TEST(TTreeReaderBasic, EntryListAndEntryRange)
{
   TTree t("t", "t");
   int x = 0;
   t.Branch("x", &x);
   for (x = 0; x < 10; ++x)
      t.Fill();

   TEntryList l;
   t.SetEntryList(&l);
   l.Enter(1);
   l.Enter(2);
   l.Enter(8);
   l.Enter(9);

   TTreeReader r(&t, &l);
   r.SetEntriesRange(1, 3);
   TTreeReaderValue<int> rv(r, "x");
   EXPECT_TRUE(r.Next());
   EXPECT_EQ(*rv, 2);
   EXPECT_TRUE(r.Next());
   EXPECT_EQ(*rv, 8);
   EXPECT_FALSE(r.Next());
}

TEST(TTreeReaderBasic, TChainWithSubEntryListsAndEntryRange)
{
   InputFilesRAII files;

   TEntryList l1("l1", "l1", "t", files.fNames[0]);
   l1.Enter(1);
   TEntryList l2("l2", "l2", "t", files.fNames[1]);
   l2.Enter(0);
   l2.Enter(1);

   TEntryList l;
   l.Add(&l1);
   l.Add(&l2);

   TChain c("t");
   c.Add(files.fNames[0]);
   c.Add(files.fNames[1]);
   c.SetEntryList(&l); // strictly required!

   TTreeReader r(&c, &l);
   TTreeReaderValue<int> xv(r, "x");
   r.SetEntriesRange(1, 2);
   EXPECT_TRUE(r.Next());
   EXPECT_EQ(*xv, 2);
   EXPECT_FALSE(r.Next());
}

TEST(TTreeReaderBasic, Values) {
   auto tree = MakeTree();
   TTreeReader tr(tree.get());
   TTreeReaderArray<double> x(tr, "one.x");
   TTreeReaderArray<int> y(tr, "two.y");
   TTreeReaderValue<unsigned int> ny(tr, "two.ny");
   TTreeReaderValue<std::string> str(tr, "str");
   TTreeReaderValue<Double32_t> d32(tr, "d32");
   TTreeReaderValue<Float16_t> f16(tr, "f16");

   // Check values for first entry.
   EXPECT_TRUE(tr.Next());

   EXPECT_DOUBLE_EQ(42, x[1]);
   EXPECT_EQ(42u, *ny);
   EXPECT_EQ(17, y[0]);
   EXPECT_STREQ("first", str->c_str());
   EXPECT_FLOAT_EQ(12, *d32);
   EXPECT_FLOAT_EQ(-12, *f16);
}

// #PR 3692
TEST(TTreeReaderBasic, InfLoop)
{
   auto tree = MakeTree();
   TTreeReader tr(tree.get());
   TTreeReaderArray<float> x(tr, "0.2.0.energy");
   tr.Next();
}

// ROOT-10019
TEST(TTreeReaderBasic, DisappearingBranch)
{

   auto createFile = [](const char *fileName, int ncols) {
      // auto r = ROOT::RDataFrame(1).Define("col0",[](){return 0;}).Snapshot<int>("t","f1.root",{"col0"});
      // r->Define("col1",[](){return 0;}).Snapshot<int,int>("t","f0.root",{"col0","col1"});
      TFile f(fileName, "RECREATE");
      TTree t("t", "t");
      int i = 42;
      t.Branch("col0", &i);
      if (ncols == 2)
         t.Branch("col1", &i);
      t.Fill();
      t.Write();
      f.Close();
   };
   createFile("DisappearingBranch0.root", 2);
   createFile("DisappearingBranch1.root", 1);

   TChain c("t");
   c.Add("DisappearingBranch*.root");
   TTreeReader r(&c);
   TTreeReaderValue<int> rv(r, "col1");
   r.Next();
   EXPECT_EQ(*rv,42);
   EXPECT_FALSE(r.Next());
   
   // Make the warnings fatal.
   gErrorAbortLevel = kWarning;
   TChain c2("t");
   c2.Add("DisappearingBranch0.root");
   c2.Add("DisappearingBranch0.root");
   c2.Add("DisappearingBranch0.root");
   TTreeReader r2(&c2);
   EXPECT_EQ(r2.GetEntries(true),3);
   EXPECT_FALSE(r2.SetEntry(0));
   gErrorAbortLevel = kFatal;

   gSystem->Unlink("DisappearingBranch0.root");
   gSystem->Unlink("DisappearingBranch1.root");
}
