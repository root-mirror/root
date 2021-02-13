#include "ROOT/RDataFrame.hxx"
#include "TChain.h"
#include "TEntryList.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "gtest/gtest.h"

#include <algorithm> // std::sort
#include <vector>

// Write to disk one file for each filename, with TTree "t" with nEntries of branch "e" with increasing values
// starting at valueStart
void MakeInputFile(const std::string &filename, int nEntries, int valueStart = 0)
{
   int e = valueStart;
   const auto treename = "t";
   auto d = ROOT::RDataFrame(nEntries).Define("e", [&] { return e++; }).Snapshot<int>(treename, filename, {"e"});
}

class RMTRAII {
   bool fIsMT;

public:
   RMTRAII(bool isMT) : fIsMT(isMT)
   {
      if (fIsMT)
         ROOT::EnableImplicitMT();
   }
   ~RMTRAII()
   {
      if (fIsMT)
         ROOT::DisableImplicitMT();
   }
};

void TestTreeWithEntryList(bool isMT = false)
{
   const auto nEntries = 10;
   const auto treename = "t";
   const auto filename = "rdfentrylist.root";
   MakeInputFile(filename, nEntries);

   RMTRAII gomt(isMT);

   // do NOT pass treename and filename to the TEntryList here, see ROOT-10775
   TEntryList elist("e", "e");
   TFile f(filename);
   auto t = f.Get<TTree>(treename);
   t->SetEntryList(&elist);
   ASSERT_TRUE(elist.GetLists() == nullptr) << "Failure in setting up the TEntryList";

   // entries must be added to the TEntryList AFTER it is associated to the TTree, otherwise subentrylists are created
   elist.Enter(0);
   elist.Enter(nEntries - 1);

   auto entries = ROOT::RDataFrame(*t).Take<int>("e").GetValue();
   std::sort(entries.begin(), entries.end()); // could be out of order in MT runs
   EXPECT_EQ(entries, std::vector<int>({0, nEntries - 1}));

   gSystem->Unlink(filename);
}

void TestChainWithEntryList(bool isMT = false)
{
   // The test might seem contrived, because it is. We want to check:
   // - that we read the correct entries from the correct trees, so we need to read different entries per tree
   // - that TTreeProcessorMT re-builds TEntryLists correctly for each file/cluster, so we need a different number
   //   of entries and we need to select, via TEntryList, entries of tree2 that do not exist in tree1

   const auto nEntries = 10;
   const auto treename = "t";
   const auto file1 = "rdfentrylist1.root";
   const auto file2 = "rdfentrylist2.root";
   MakeInputFile(file1, nEntries);
   MakeInputFile(file2, nEntries * 2, /*valueStart=*/nEntries);

   RMTRAII gomt(isMT);

   TEntryList elist1("e", "e", treename, file1);
   elist1.Enter(0);
   elist1.Enter(nEntries - 1);
   TEntryList elist2("e", "e", treename, file2);
   elist2.Enter(nEntries + 1);
   elist2.Enter(2 * nEntries - 1);

   // make a TEntryList that contains two TEntryLists in its list of TEntryLists,
   // as required by TChain (see TEntryList's doc)
   TEntryList elists;
   elists.Add(&elist1);
   elists.Add(&elist2);

   TChain c(treename);
   c.Add(file1, nEntries);
   c.Add(file2, nEntries * 2);
   c.SetEntryList(&elists);

   auto entries = ROOT::RDataFrame(c).Take<int>("e").GetValue();
   std::sort(entries.begin(), entries.end()); // could be out of order in MT runs
   EXPECT_EQ(entries, std::vector<int>({0, nEntries - 1, 2 * nEntries + 1, 3 * nEntries - 1}));

   gSystem->Unlink(file1);
   gSystem->Unlink(file2);
}

TEST(RDFEntryList, Chain)
{
   TestChainWithEntryList();
}

TEST(RDFEntryList, Tree)
{
   TestTreeWithEntryList();
}

#ifdef R__USE_IMT
TEST(RDFEntryList, ChainMT)
{
   TestChainWithEntryList(true);
}

TEST(RDFEntryList, TreeMT)
{
   TestTreeWithEntryList(true);
}
#endif
