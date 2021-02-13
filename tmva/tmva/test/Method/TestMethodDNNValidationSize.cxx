#include "gtest/gtest.h"

#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/MethodDNN.h"

#include "Rtypes.h"

#include <exception>
#include <memory>

using namespace std;

namespace TMVA {

struct TestMethodDNNValidationSize {

   TestMethodDNNValidationSize(UInt_t nEvents)
   {
      if (nEvents < 4)
         throw runtime_error("Use min 4 events");

      fFactory = unique_ptr<Factory>(new Factory("", "Silent:AnalysisType=Classification"));
      fDataLoader = shared_ptr<DataLoader>(new DataLoader("dataset"));

      fDataLoader->AddVariable("x", 'F');

      // Add nEvents events
      for (UInt_t i = 0; i < nEvents - 3; ++i) {
         fDataLoader->AddEvent("cls1", Types::kTraining, {0}, 1);
      }
      fDataLoader->AddEvent("cls1", Types::kTraining, {1}, 1);

      fDataLoader->AddEvent("cls2", Types::kTraining, {0}, 1);
      fDataLoader->AddEvent("cls2", Types::kTraining, {1}, 1);

      // To aviod triggering crash if there is no test tree
      fDataLoader->AddEvent("cls1", Types::kTesting, {0}, 1);
      fDataLoader->AddEvent("cls2", Types::kTesting, {0}, 1);

      fDataLoader->PrepareTrainingAndTestTree("", "");
   };

   size_t GetProcessedValidationOption(TString options)
   {
#ifdef DNNCPU
      const TString defaultOptions = "!H:!V:Layout=RELU|50:Architecture=CPU:";
#elif defined DNNCUDA
      const TString defaultOptions = "!H:!V:Layout=RELU|50:Architecture=GPU:";
#else
#error "This should not happen. Can only compile with CPU or CUDA implementations."
#endif

      IMethod *m = fFactory->BookMethod(fDataLoader.get(), Types::kDNN, "DNN", defaultOptions + options);
      MethodDNN *mdnn = dynamic_cast<MethodDNN *>(m);

      size_t numValidationSamples = mdnn->GetNumValidationSamples();
      fFactory->DeleteAllMethods();
      fFactory->fMethodsMap.clear();

      return numValidationSamples;
   };

private:
   std::unique_ptr<TMVA::Factory> fFactory;
   std::shared_ptr<TMVA::DataLoader> fDataLoader;
};

} // End namespace TMVA

TEST(MethodDnnValidationSize, RelativePercent)
{
   const UInt_t nEvents = 100;
   TMVA::TestMethodDNNValidationSize context{nEvents};

   size_t nValidSamples;

   // 20.0%
   nValidSamples = context.GetProcessedValidationOption("ValidationSize=20.0%");
   EXPECT_EQ(nValidSamples, (size_t)(nEvents * 0.2));

   // 20%
   nValidSamples = context.GetProcessedValidationOption("ValidationSize=20%");
   EXPECT_EQ(nValidSamples, (size_t)(nEvents * 0.2));

   // 2e1%
   nValidSamples = context.GetProcessedValidationOption("ValidationSize=2e1%");
   EXPECT_EQ(nValidSamples, (size_t)(nEvents * 0.2));

   // 100% - all events in training set -> error
   EXPECT_THROW(context.GetProcessedValidationOption("ValidationSize=100%"), std::runtime_error);

   // 101% - more events than training set -> error
   EXPECT_THROW(context.GetProcessedValidationOption("ValidationSize=101%"), std::runtime_error);

   // -1% - negative numnber -> error
   EXPECT_THROW(context.GetProcessedValidationOption("ValidationSize=-1%"), std::runtime_error);

   // 0% - no events -> error
   EXPECT_THROW(context.GetProcessedValidationOption("ValidationSize=0%"), std::runtime_error);
}

TEST(MethodDnnValidationSize, RelativeFraction)
{
   const UInt_t nEvents = 100;
   TMVA::TestMethodDNNValidationSize context{nEvents};

   size_t nValidSamples;

   // Values here must be between 0 and 1 exclusive!
   // 0 = no events, 1 = one event

   // 0.2
   nValidSamples = context.GetProcessedValidationOption("ValidationSize=0.2");
   EXPECT_EQ(nValidSamples, (size_t)(nEvents * 0.2));

   // .2
   nValidSamples = context.GetProcessedValidationOption("ValidationSize=.2");
   EXPECT_EQ(nValidSamples, (size_t)(nEvents * 0.2));

   // 2e-1
   nValidSamples = context.GetProcessedValidationOption("ValidationSize=2e-1");
   EXPECT_EQ(nValidSamples, (size_t)(nEvents * 0.2));

   // -0.1 - negative number -> error
   EXPECT_THROW(context.GetProcessedValidationOption("ValidationSize=-0.1"), std::runtime_error);

   // 0.0 - no events -> error
   EXPECT_THROW(context.GetProcessedValidationOption("ValidationSize=0.0"), std::runtime_error);
}

TEST(MethodDnnValidationSize, Absolute)
{
   const UInt_t nEvents = 100;
   TMVA::TestMethodDNNValidationSize context{nEvents};

   size_t nValidSamples;

   // 1
   nValidSamples = context.GetProcessedValidationOption("ValidationSize=1");
   EXPECT_EQ(nValidSamples, (size_t)(1));

   // 1.0
   nValidSamples = context.GetProcessedValidationOption("ValidationSize=1.0");
   EXPECT_EQ(nValidSamples, (size_t)(1.0));

   // 1.5
   nValidSamples = context.GetProcessedValidationOption("ValidationSize=1.5");
   EXPECT_EQ(nValidSamples, (size_t)(1.5));

   // 1e1
   nValidSamples = context.GetProcessedValidationOption("ValidationSize=1e1");
   EXPECT_EQ(nValidSamples, (size_t)(1e1));

   // 100 - uses all events in training set for validation -> error
   EXPECT_THROW(context.GetProcessedValidationOption("ValidationSize=100"), std::runtime_error);

   // -1 - negative number -> error
   EXPECT_THROW(context.GetProcessedValidationOption("ValidationSize=-1"), std::runtime_error);

   // 101 - more events than in training set -> error
   EXPECT_THROW(context.GetProcessedValidationOption("ValidationSize=101"), std::runtime_error);
}

int main(int argc, char *argv[])
{
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
