// @(#)root/tmva $Id$
// Author: Omar Zapata

/*************************************************************************
 * Copyright (C) 2017, Omar Zapata                                       *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <gtest/gtest.h>

#include <TMVA/Classification.h>
#include <TMVA/DataLoader.h>
#include <TMVA/ROCCurve.h>
#include <TMVA/DataInputHandler.h>
#include <TFile.h>
#include <TSystem.h>
#include <TString.h>

using namespace TMVA;

class ClassifierInterface : public TMVA::Experimental::Classification {
public:
   explicit ClassifierInterface(DataLoader *loader, TFile *file, TString options)
      : Classification(loader, file, options)
   {
   }
   explicit ClassifierInterface(DataLoader *loader, TString options) : Classification(loader, options) {}

   TString _GetMethodOptions(TString methodname, TString methodtitle)
   {
      return GetMethodOptions(methodname, methodtitle);
   }
   Bool_t _HasMethodObject(TString methodname, TString methodtitle, Int_t &index)
   {
      return HasMethodObject(methodname, methodtitle, index);
   }
   Bool_t _IsCutsMethod(TMVA::MethodBase *method) { return IsCutsMethod(method); }
   TMVA::ROCCurve *
   _GetROC(TMVA::MethodBase *method, UInt_t iClass = 0, TMVA::Types::ETreeType type = TMVA::Types::kTesting)
   {
      return GetROC(method, iClass, type);
   }
   TMVA::ROCCurve *_GetROC(TString methodname, TString methodtitle, UInt_t iClass = 0,
                           TMVA::Types::ETreeType type = TMVA::Types::kTesting)
   {
      return GetROC(methodname, methodtitle, iClass, type);
   }

   Double_t _GetROCIntegral(TString methodname, TString methodtitle, UInt_t iClass = 0)
   {
      return GetROCIntegral(methodname, methodtitle, iClass);
   }

   TMVA::Experimental::ClassificationResult &_GetResults(TString methodname, TString methodtitle)
   {
      return GetResults(methodname, methodtitle);
   }
   void _CopyFrom(TDirectory *src, TFile *file) { CopyFrom(src, file); }
   void _MergeFiles() { MergeFiles(); }
};

// class to load basic configuration for all classification tests,
// like environment variables and data in the DataLoader.
class TestBaseClass {
public:
   ClassifierInterface *cl;
   DataLoader *dl;
   TFile *inputfile;
   TFile *outfile;

   TestBaseClass() : cl(nullptr), dl(nullptr), inputfile(nullptr), outfile(nullptr) { TMVA::Tools::Instance(); }

   void CreateDataLoader()
   {
      TString fname = "./tmva_class_example.root";
      if (!gSystem->AccessPathName(fname)) {
         inputfile = TFile::Open(fname); // check if file in local directory exists
      } else {
         TFile::SetCacheFileDir(".");
         inputfile = TFile::Open("http://root.cern/files/tmva_class_example.root", "CACHEREAD");
      }
      if (!inputfile) {
         std::cout << "ERROR: could not open data file" << std::endl;
         exit(1);
      }

      TTree *signalTree = (TTree *)inputfile->Get("TreeS");
      TTree *background = (TTree *)inputfile->Get("TreeB");

      dl = new TMVA::DataLoader("dataset");
      dl->AddVariable("myvar1 := var1+var2", 'F');
      dl->AddVariable("myvar2 := var1-var2", "Expression 2", "", 'F');
      dl->AddVariable("var3", "Variable 3", "units", 'F');
      dl->AddVariable("var4", "Variable 4", "units", 'F');

      dl->AddSpectator("spec1 := var1*2", "Spectator 1", "units", 'F');
      dl->AddSpectator("spec2 := var1*3", "Spectator 2", "units", 'F');

      Double_t signalWeight = 1.0;
      Double_t backgroundWeight = 1.0;

      dl->AddSignalTree(signalTree, signalWeight);
      dl->AddBackgroundTree(background, backgroundWeight);

      dl->SetBackgroundWeightExpression("weight");
      dl->PrepareTrainingAndTestTree("", "", "nTrain_Signal=1000:nTrain_Background=1000:nTest_Signal=5000:nTest_"
                                             "Background=5000:SplitMode=Random:NormMode=NumEvents:!V");
   }

   ~TestBaseClass()
   {
      if (dl != nullptr)
         delete dl;
      if (inputfile != nullptr)
         inputfile->Close();
      if (outfile != nullptr)
         outfile->Close();
   }
};

// Test calling Train/Test without call evaluation method and without output file
class ClassifierTest1 : public ::testing::Test, public TestBaseClass {
public:
   ClassifierTest1()
   {
      CreateDataLoader();
      cl = new ClassifierInterface(dl, "Jobs=3");
   }
   void Booking()
   {
      cl->BookMethod("BDT", "BDTG", "NTrees=100");
      cl->BookMethod(TMVA::Types::kBDT, "BDTB", "NTrees=100:BoostType=Bagging");
   }
};

// testing methods and Train/Test
TEST_F(ClassifierTest1, BasicTests)
{
   Booking();

   ASSERT_EQ(cl->_GetMethodOptions("BDT", "BDTB"), "NTrees=100:BoostType=Bagging");
   ASSERT_EQ(cl->_GetMethodOptions("BDTB_", "BDTB"), "");
   Int_t index = -1;
   ASSERT_EQ(cl->_HasMethodObject("BDT", "BDTB", index), kFALSE);
   ASSERT_EQ(index, -1);

   auto m = cl->GetMethod("BDT", "BDTB");
   ASSERT_NE(m, nullptr);

   ASSERT_EQ(cl->_HasMethodObject("BDT", "BDTB", index), kTRUE);
   ASSERT_EQ(index, 0);

   ASSERT_EQ(cl->_IsCutsMethod(m), kFALSE);

   cl->TrainMethod("BDT", "BDTB");
   auto m2 = cl->GetMethod("BDT", "BDTB");
   ASSERT_EQ(m, m2);

   cl->TestMethod("BDT", "BDTB");
   EXPECT_NEAR(cl->_GetROCIntegral("BDT", "BDTB"), 0.853, 0.03);

   auto r1 = cl->GetResults();
   ASSERT_EQ(r1.size(), (UInt_t)1);
   ASSERT_EQ(r1[0].GetMethodName(), "BDT");
   ASSERT_EQ(r1[0].GetMethodTitle(), "BDTB");
   auto roc = r1[0].GetROC();
   ASSERT_NE(roc, nullptr);
   EXPECT_NEAR(roc->GetROCIntegral(), 0.853, 0.03);
   ASSERT_EQ(r1[0].GetDataLoaderName(), "dataset");

   ASSERT_EQ(cl->_GetMethodOptions("BDT", "BDTG"), "NTrees=100");
   index = -1;
   ASSERT_EQ(cl->_HasMethodObject("BDT", "BDTG", index), kFALSE);
   ASSERT_EQ(index, -1);

   auto m3 = cl->GetMethod("BDT", "BDTG");
   ASSERT_NE(m3, nullptr);

   ASSERT_EQ(cl->_HasMethodObject("BDT", "BDTG", index), kTRUE);
   ASSERT_EQ(index, 1);

   ASSERT_EQ(cl->_IsCutsMethod(m3), kFALSE);

   cl->TrainMethod("BDT", "BDTG");
   auto m4 = cl->GetMethod("BDT", "BDTG");
   ASSERT_EQ(m3, m4);

   cl->TestMethod("BDT", "BDTG");
   EXPECT_NEAR(cl->_GetROCIntegral("BDT", "BDTG"), 0.897, 0.03);

   auto r2 = cl->GetResults();
   ASSERT_EQ(r2.size(), (UInt_t)2);
   ASSERT_EQ(r2[1].GetMethodName(), "BDT");
   ASSERT_EQ(r2[1].GetMethodTitle(), "BDTG");
   roc = r2[1].GetROC();
   ASSERT_NE(roc, nullptr);
   EXPECT_NEAR(roc->GetROCIntegral(), 0.897, 0.03);
   ASSERT_EQ(r2[0].GetDataLoaderName(), "dataset");
}

// Test calling Train/Test with output file
class ClassifierTest2 : public ::testing::Test, public TestBaseClass {
public:
   ClassifierTest2()
   {
      CreateDataLoader();
      outfile = new TFile("TMVAtmp1.root", "RECREATE");
      cl = new ClassifierInterface(dl, outfile, "Jobs=2");
   }
   void Booking()
   {
      cl->BookMethod("BDT", "BDTG", "NTrees=100");
      cl->BookMethod(TMVA::Types::kBDT, "BDTB", "NTrees=100:BoostType=Bagging");
   }
};

// tests to check that the output in the file is right
TEST_F(ClassifierTest2, TestsOverOutput)
{
   Booking();
   auto ds_ = (TDirectoryFile *)outfile->Get("dataset"); // must not exists yet
   ASSERT_EQ(ds_, nullptr);

   cl->TrainMethod("BDT", "BDTB");
   // looking output for the first trained method
   auto ds = (TDirectoryFile *)outfile->Get("dataset"); // get dataset dir
   ASSERT_NE(ds, nullptr);

   auto Method_BDTB = (TDirectoryFile *)ds->Get("Method_BDT"); // get method 1 training output dir
   ASSERT_NE(Method_BDTB, nullptr);

   auto BDTB = (TDirectoryFile *)Method_BDTB->Get("BDTB");
   ASSERT_NE(BDTB, nullptr);

   ASSERT_NE(BDTB->Get("TrainingPath"), nullptr);
   ASSERT_NE(BDTB->Get("WeightFileName"), nullptr);
   ASSERT_NE(BDTB->Get("MonitorNtuple"), nullptr);

   cl->TestMethod("BDT", "BDTB");
   // looking output for the first tested method

   auto TrainTree = (TTree *)ds->Get("TrainTree");
   auto TestTree = (TTree *)ds->Get("TestTree");
   ASSERT_NE(TrainTree, nullptr);
   ASSERT_NE(TestTree, nullptr);

   ASSERT_EQ(TrainTree->GetEntries(), 2000); // 1000 events for sgn and 1000 for Bkg
   ASSERT_EQ(TestTree->GetEntries(), 10000); // 5000 events for sgn and 5000 for Bkg
}
