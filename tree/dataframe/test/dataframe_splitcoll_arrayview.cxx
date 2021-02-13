// this test was roottest's root/dataframe/test_splitcoll_arrayview.C
#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "ROOT/TSeq.hxx"
#include "ROOT/RVec.hxx"
#include "TSystem.h"
#include "TROOT.h"
#include "gtest/gtest.h"

#include "TwoFloats.h"

#include <iostream>

void fill_tree(const char *filename, const char *treeName)
{
   TFile f(filename, "RECREATE");
   TTree t(treeName, treeName);
   std::vector<TwoFloats> v;
   t.Branch("v", "vector<TwoFloats>", &v, 32000, 2);
   for (auto i : ROOT::TSeqI(10)) {
      v.emplace_back(i);
      t.Fill();
   }
   t.Write();
   f.Close();
}

void test_splitcoll_arrayview(const std::string &fileName, const std::string &treeName)
{
   ROOT::RDataFrame d1(treeName, fileName, {"v.a"});
   auto c1 = d1.Filter([](ROOT::VecOps::RVec<float> d) {
                  float ex_v = 0.f;
                  for (auto v : d) {
                     EXPECT_DOUBLE_EQ(v, ex_v);
                     ex_v += 1.f;
                  }
                  return d[0] > 5;
               })
                .Count();
   EXPECT_EQ(*c1, 0ull);

   ROOT::RDataFrame d2(treeName, fileName, {"v"});
   auto c2 = d2.Filter([](ROOT::VecOps::RVec<TwoFloats> d) {
                  int q = 0;
                  float ex_a = 0.f;
                  for (auto v : d) {
                     EXPECT_DOUBLE_EQ(v.a, ex_a);
                     q += int(v.a);
                     ex_a += 1.f;
                  }
                  return 0 == q % 3;
               })
                .Count();
   EXPECT_EQ(*c2, 7ull);
}

TEST(RDFSimpleTests, SplitCollectionArrayView)
{
   auto fileName = "myfile_test_splitcoll_arrayview.root";
   auto treeName = "myTree";
   fill_tree(fileName, treeName);
   test_splitcoll_arrayview(fileName, treeName);
   gSystem->Unlink(fileName);
}
