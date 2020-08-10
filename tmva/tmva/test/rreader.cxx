#include <gtest/gtest.h>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>

#include <TMVA/RReader.hxx>
#include <TMVA/RInferenceUtils.hxx>
#include <TMVA/RTensor.hxx>
#include <TMVA/RTensorUtils.hxx>

using namespace TMVA::Experimental;

// Classification
static const std::string modelClassification = "RReaderClassification/weights/RReaderClassification_BDT.weights.xml";
static const std::string filenameClassification = "http://root.cern/files/tmva_class_example.root";
static const std::vector<std::string> variablesClassification = {"var1", "var2", "var3", "var4"};

void TrainClassificationModel()
{
   // Check for existing training
   if (gSystem->mkdir("RReaderClassification") == -1)
#ifndef _MSC_VER
      return;
#else
      std::cout << "The directory \"RReaderClassification\" exists already...\n";
#endif

   // Create factory
   auto output = TFile::Open("TMVA.root", "RECREATE");
   auto factory = new TMVA::Factory("RReaderClassification",
           output, "Silent:!V:!DrawProgressBar:AnalysisType=Classification");

   // Open trees with signal and background events
   const std::string filename = "http://root.cern/files/tmva_class_example.root";
   auto data = TFile::Open(filename.c_str());
   auto signal = (TTree *)data->Get("TreeS");
   auto background = (TTree *)data->Get("TreeB");

   // Add variables and register the trees with the dataloader
   auto dataloader = new TMVA::DataLoader("RReaderClassification");
   const std::vector<std::string> variables = {"var1", "var2", "var3", "var4"};
   for (const auto &var : variables) {
      dataloader->AddVariable(var);
   }
   dataloader->AddSignalTree(signal, 1.0);
   dataloader->AddBackgroundTree(background, 1.0);
   dataloader->PrepareTrainingAndTestTree("", "");

   // Train a TMVA method
   factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", "!V:!H:NTrees=100:MaxDepth=2");
   factory->TrainAllMethods();
   output->Close();
}

// Regression
static const std::string modelRegression = "RReaderRegression/weights/RReaderRegression_BDTG.weights.xml";
static const std::string filenameRegression = "http://root.cern/files/tmva_reg_example.root";
static const std::vector<std::string> variablesRegression = {"var1", "var2"};

void TrainRegressionModel()
{
   // Check for existing training
   if (gSystem->mkdir("RReaderRegression") == -1)
#ifndef _MSC_VER
      return;
#else
      std::cout << "The directory \"RReaderRegression\" exists already...\n";
#endif

   // Create factory
   auto output = TFile::Open("TMVA.root", "RECREATE");
   auto factory = new TMVA::Factory("RReaderRegression",
           output, "Silent:!V:!DrawProgressBar:AnalysisType=Regression");

   // Open trees with signal and background events
   const std::string filename = "http://root.cern/files/tmva_reg_example.root";
   auto data = TFile::Open(filename.c_str());
   auto tree = (TTree *)data->Get("TreeR");

   // Add variables and register the trees with the dataloader
   auto dataloader = new TMVA::DataLoader("RReaderRegression");
   dataloader->AddVariable("var1");
   dataloader->AddVariable("var2");
   dataloader->AddTarget("fvalue");
   dataloader->AddRegressionTree(tree, 1.0);
   dataloader->PrepareTrainingAndTestTree("", "");

   // Train a TMVA method
   factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG", "!V:!H:NTrees=100:MaxDepth=2");
   factory->TrainAllMethods();
   output->Close();
}

// Multiclass
static const std::string modelMulticlass = "RReaderMulticlass/weights/RReaderMulticlass_BDT.weights.xml";
static const std::string filenameMulticlass = "http://root.cern/files/tmva_multiclass_example.root";
static const std::vector<std::string> variablesMulticlass = variablesClassification;

void TrainMulticlassModel()
{
   // Check for existing training
   if (gSystem->mkdir("RReaderMulticlass") == -1)
#ifndef _MSC_VER
      return;
#else
      std::cout << "The directory \"RReaderMulticlass\" exists already...\n";
#endif

   // Create factory
   auto output = TFile::Open("TMVA.root", "RECREATE");
   auto factory = new TMVA::Factory("RReaderMulticlass",
           output, "Silent:!V:!DrawProgressBar:AnalysisType=Multiclass");

   // Open trees with signal and background events
   const std::string filename = "http://root.cern/files/tmva_multiclass_example.root";
   auto data = TFile::Open(filename.c_str());
   auto signal = (TTree *)data->Get("TreeS");
   auto background0 = (TTree *)data->Get("TreeB0");
   auto background1 = (TTree *)data->Get("TreeB1");
   auto background2 = (TTree *)data->Get("TreeB2");

   // Add variables and register the trees with the dataloader
   auto dataloader = new TMVA::DataLoader("RReaderMulticlass");
   const std::vector<std::string> variables = {"var1", "var2", "var3", "var4"};
   for (const auto &var : variables) {
      dataloader->AddVariable(var);
   }
   dataloader->AddTree(signal, "Signal");
   dataloader->AddTree(background0, "Background_0");
   dataloader->AddTree(background1, "Background_1");
   dataloader->AddTree(background2, "Background_2");
   dataloader->PrepareTrainingAndTestTree("", "");

   // Train a TMVA method
   factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", "!V:!H:NTrees=100:MaxDepth=2:BoostType=Grad");
   factory->TrainAllMethods();
   output->Close();
}

TEST(RReader, ClassificationGetVariables)
{
   TrainClassificationModel();
   RReader model(modelClassification);
   auto vars = model.GetVariableNames();
   EXPECT_EQ(vars.size(), 4ul);
   for (std::size_t i = 0; i < vars.size(); i++) {
      EXPECT_EQ(vars[i], variablesClassification[i]);
   }
}

TEST(RReader, ClassificationComputeVector)
{
   TrainClassificationModel();
   const std::vector<float> x = {1.0, 2.0, 3.0, 4.0};
   RReader model(modelClassification);
   auto y = model.Compute(x);
   EXPECT_EQ(y.size(), 1ul);
}

TEST(RReader, ClassificationComputeTensor)
{
   TrainClassificationModel();
   ROOT::RDataFrame df("TreeS", filenameClassification);
   auto x = AsTensor<float>(df, variablesClassification);

   RReader model(modelClassification);
   auto y = model.Compute(x);

   const auto shapeX = x.GetShape();
   const auto shapeY = y.GetShape();
   EXPECT_EQ(shapeY.size(), 1ul);
   EXPECT_EQ(shapeY[0], shapeX[0]);
}

TEST(RReader, ClassificationComputeDataFrame)
{
   TrainClassificationModel();
   ROOT::RDataFrame df("TreeS", filenameClassification);
   RReader model(modelClassification);
   auto df2 = df.Define("y", Compute<4, float>(model), variablesClassification);
   auto df3 = df2.Filter("y.size() == 1");
   auto c = df3.Count();
   auto y = df2.Take<std::vector<float>>("y");
   EXPECT_EQ(y->size(), *c);
}

TEST(RReader, RegressionGetVariables)
{
   TrainRegressionModel();
   RReader model(modelRegression);
   auto vars = model.GetVariableNames();
   EXPECT_EQ(vars.size(), 2ul);
   for (std::size_t i = 0; i < vars.size(); i++) {
      EXPECT_EQ(vars[i], variablesRegression[i]);
   }
}

TEST(RReader, RegressionComputeVector)
{
   TrainRegressionModel();
   const std::vector<float> x = {1.0, 2.0};
   RReader model(modelRegression);
   auto y = model.Compute(x);
   EXPECT_EQ(y.size(), 1ul);
}

TEST(RReader, RegressionComputeTensor)
{
   TrainRegressionModel();
   ROOT::RDataFrame df("TreeR", filenameRegression);
   auto x = AsTensor<float>(df, variablesRegression);

   RReader model(modelRegression);
   auto y = model.Compute(x);

   const auto shapeX = x.GetShape();
   const auto shapeY = y.GetShape();
   EXPECT_EQ(shapeY.size(), 1ul);
   EXPECT_EQ(shapeY[0], shapeX[0]);
}

TEST(RReader, RegressionComputeDataFrame)
{
   TrainRegressionModel();
   ROOT::RDataFrame df("TreeR", filenameRegression);
   RReader model(modelRegression);
   auto df2 = df.Define("y", Compute<2, float>(model), variablesRegression);
   auto df3 = df2.Filter("y.size() == 1");
   auto c = df3.Count();
   auto y = df2.Take<std::vector<float>>("y");
   EXPECT_EQ(y->size(), *c);
}

TEST(RReader, MulticlassGetVariables)
{
   TrainMulticlassModel();
   RReader model(modelMulticlass);
   auto vars = model.GetVariableNames();
   EXPECT_EQ(vars.size(), 4ul);
   for (std::size_t i = 0; i < vars.size(); i++) {
      EXPECT_EQ(vars[i], variablesMulticlass[i]);
   }
}

TEST(RReader, MulticlassComputeVector)
{
   TrainMulticlassModel();
   const std::vector<float> x = {1.0, 2.0, 3.0, 4.0};
   RReader model(modelMulticlass);
   auto y = model.Compute(x);
   EXPECT_EQ(y.size(), 4ul);
}

TEST(RReader, MulticlassComputeTensor)
{
   TrainMulticlassModel();
   ROOT::RDataFrame df("TreeS", filenameMulticlass);
   auto x = AsTensor<float>(df, variablesMulticlass);

   RReader model(modelMulticlass);
   auto y = model.Compute(x);

   const auto shapeX = x.GetShape();
   const auto shapeY = y.GetShape();
   EXPECT_EQ(shapeY.size(), 2ul);
   EXPECT_EQ(shapeY[0], shapeX[0]);
   EXPECT_EQ(shapeY[1], 4ul);
}

TEST(RReader, MulticlassComputeDataFrame)
{
   TrainMulticlassModel();
   ROOT::RDataFrame df("TreeS", filenameMulticlass);
   RReader model(modelMulticlass);
   auto df2 = df.Define("y", Compute<4, float>(model), variablesMulticlass);
   auto df3 = df2.Filter("y.size() == 4");
   auto c = df3.Count();
   auto y = df2.Take<std::vector<float>>("y");
   EXPECT_EQ(y->size(), *c);
}
