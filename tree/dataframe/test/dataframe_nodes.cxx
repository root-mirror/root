#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RSlotStack.hxx>
#include <TSystem.h>

#include <mutex>
#include <thread>

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

TEST(RDataFrameNodes, RLoopManagerJit)
{
   ROOT::Detail::RDF::RLoopManager lm(nullptr, {});
   lm.ToJitExec("souble d = 3.14");
   int ret(1);
   try {
      testing::internal::CaptureStderr();
      lm.Run();
   } catch (const std::runtime_error &) {
      ret = 0;
   }
   EXPECT_EQ(0, ret) << "Bogus C++ code was jitted and nothing was detected!";
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

   testing::internal::CaptureStdout();
   *tdf.Count();
   auto output = testing::internal::GetCapturedStdout();
   EXPECT_STREQ("", output.c_str()) << "An error was printed: " << output << std::endl;

   for (auto &f : files)
      gSystem->Unlink(f.c_str());
}

// ROOT-9736
TEST(RDataFrameNodes, InheritanceOfCustomColumns) 
{
   ROOT::RDataFrame df(1);
   int nBins = -1;
   const auto nBinsExpected = 42;
   df.Define("b", [](){return TH1F("b","b",nBinsExpected,0,1);})
      .Foreach([&nBins](TH1& h){nBins = h.GetNbinsX();}, {"b"});
   EXPECT_EQ(nBins, nBinsExpected);
}