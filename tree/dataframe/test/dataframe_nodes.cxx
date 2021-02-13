#include "ROOTUnitTestSupport.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RSlotStack.hxx>
#include <TStatistic.h> // To check reading of columns with types which are mothers of the column type
#include <TSystem.h>

#include <mutex>
#include <thread>
#include <stdexcept> // std::runtime_error

#include "gtest/gtest.h"

#ifndef NDEBUG

TEST(RDataFrameNodes, RSlotStackGetOneTooMuch)
{
   auto theTest = []() {
      unsigned int n(2);
      ROOT::Internal::RDF::RSlotStack s(n);

      std::vector<std::thread> ts;

      for (unsigned int i = 0; i < 3; ++i) {
         ts.emplace_back([&s]() { s.GetSlot(); });
      }

      for (auto &&t : ts)
         t.join();
   };

   EXPECT_DEATH(theTest(), "Trying to pop a slot from an empty stack!");
}

TEST(RDataFrameNodes, RSlotStackPutBackTooMany)
{
   auto theTest = []() {
      ROOT::Internal::RDF::RSlotStack s(1);
      s.ReturnSlot(0);
   };

   EXPECT_DEATH(theTest(), "Trying to put back a slot to a full stack!");
}

#endif

TEST(RDataFrameNodes, RLoopManagerGetLoopManagerUnchecked)
{
   ROOT::Detail::RDF::RLoopManager lm(nullptr, {});
   ASSERT_EQ(&lm, lm.GetLoopManagerUnchecked());
}

TEST(RDataFrameNodes, RLoopManagerJitWrongCode)
{
   ROOT::Detail::RDF::RLoopManager lm(nullptr, {});
   lm.ToJitExec("souble d = 3.14");
   EXPECT_THROW(lm.Run(), std::runtime_error) << "Bogus C++ code was jitted and nothing was detected!";
}

TEST(RDataFrameNodes, DoubleEvtLoop)
{
   ROOT::RDataFrame d1(4);
   auto d = d1.Define("x", []() { return 2; });

   std::vector<std::string> files{"f1.root", "f2.root"};

   for (auto &f : files)
      d.Snapshot<int>("t1", f, {"x"});

   ROOT::RDataFrame tdf("t1", files);
   *tdf.Count();

   // Check that this is not printed
   // Warning in <TTreeReader::SetEntryBase()>: The current tree in the TChain t1 has changed (e.g. by TTree::Process)
   // even though TTreeReader::SetEntry() was called, which switched the tree again. Did you mean to call
   // TTreeReader::SetLocalEntry()?

   ROOT_EXPECT_NODIAG(*tdf.Count());

   for (auto &f : files)
      gSystem->Unlink(f.c_str());
}

// ROOT-9736
TEST(RDataFrameNodes, InheritanceOfDefines)
{
   ROOT::RDataFrame df(1);
   const auto nBinsExpected = 42;
   // Read the TH1F as a TH1
   df.Define("b", [&]() { return TH1F("b", "b", nBinsExpected, 0, 1); })
      .Foreach([&](TH1 &h) { EXPECT_EQ(h.GetNbinsX(), nBinsExpected);}, {"b"});

   const auto ofileName = "InheritanceOfDefines.root";

   const auto val = 42.;
   auto createStat = [&val]() {
      TStatistic t;
      t.Fill(val);
      return t;
   };

   // Read as TObject from disk a TStatistics object
   auto checkStat = [&val](TObject &o) { EXPECT_EQ(val, ((TStatistic *)&o)->GetMean()); };
   ROOT::RDataFrame(1).Define("x", createStat).Snapshot<TStatistic>("t", ofileName, {"x"})->Foreach(checkStat, {"x"});
   gSystem->Unlink(ofileName);
}
