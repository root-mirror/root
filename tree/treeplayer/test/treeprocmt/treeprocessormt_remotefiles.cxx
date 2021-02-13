#include <ROOT/TTreeProcessorMT.hxx>
#include <TTree.h>

#include <atomic>

#include "gtest/gtest.h"

class TTree;
class TTreeReader;

TEST(TreeProcessorMT, PathName)
{
   auto fname = "root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/ZZTo4mu.root";
   auto f = std::unique_ptr<TFile>(TFile::Open(fname));
   ASSERT_TRUE(f != nullptr) << "Could not open remote file\n";
   auto tree = f->Get<TTree>("Events");
   ROOT::TTreeProcessorMT p(*tree);
   std::atomic<unsigned int> n(0U);
   auto func = [&n](TTreeReader &t) {
      while (t.Next())
         n++;
   };
   p.Process(func);
   EXPECT_EQ(n.load(), 1499064U) << "Wrong number of events processed!\n";
}
