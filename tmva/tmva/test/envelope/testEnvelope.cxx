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

#include <TMVA/Envelope.h>
#include <TMVA/DataLoader.h>
#include <TMVA/DataInputHandler.h>
#include <TFile.h>
#include <TSystem.h>
#include <TString.h>

using namespace TMVA;

// Envelope is a pure base class,
// I will to create a basic class to create an object
class BasicEnvelope : public TMVA::Envelope {
public:
   BasicEnvelope(const TString &name, DataLoader *dalaloader, TFile *file, const TString options)
      : Envelope(name, dalaloader, file, options)
   {
   }
   TString GetTransformations() { return fTransformations; }

   Int_t GetJobs() { return fJobs; }

   virtual void Evaluate(){};

   std::vector<OptionMap> &GetMethods() { return fMethods; }

   DataInputHandler &GetDataLoaderDataInput() { return ::Envelope::GetDataLoaderDataInput(); }

   DataSetInfo &GetDataLoaderDataSetInfo() { return ::Envelope::GetDataLoaderDataSetInfo(); }

   DataLoader *GetDataLoader() { return fDataLoader.get(); }
};

// test without DataLoader or output file
class EnvelopeTest1 : public ::testing::Test {
public:
   BasicEnvelope *envelope;
   EnvelopeTest1()
   {
      TMVA::Tools::Instance();
      envelope = new BasicEnvelope("envelope", nullptr, nullptr, "ModelPersistence=1:V:Transformations=I:Jobs=3");
   }
   ~EnvelopeTest1() { delete envelope; }
   void Booking()
   {
      envelope->BookMethod("BDT", "BDTG", "NTrees=100");
      envelope->BookMethod(TMVA::Types::kBDT, "BDTB", "BoostType=Bagging");
   }
};

TEST_F(EnvelopeTest1, ParseOptions)
{
   envelope->ParseOptions();

   EXPECT_EQ(kTRUE, envelope->IsModelPersistence());
   envelope->SetModelPersistence(kFALSE);
   EXPECT_EQ(kFALSE, envelope->IsModelPersistence());

   EXPECT_EQ(kTRUE, envelope->IsVerbose());
   envelope->SetVerbose(kFALSE);
   EXPECT_EQ(kFALSE, envelope->IsVerbose());

   EXPECT_EQ("I", envelope->GetTransformations());
   EXPECT_EQ(3, envelope->GetJobs());

   EXPECT_EQ(kTRUE, envelope->IsSilentFile());
   EXPECT_EQ(nullptr, envelope->GetFile());

   EXPECT_EQ(nullptr, envelope->GetDataLoader());
}

TEST_F(EnvelopeTest1, Booking)
{
   Booking();
   EXPECT_EQ(kTRUE, envelope->HasMethod("BDT", "BDTG"));
   EXPECT_EQ(kFALSE, envelope->HasMethod("SVM", "SVM"));
   EXPECT_EQ(kTRUE, envelope->HasMethod("BDT", "BDTB"));

   envelope->BookMethod(TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm");
   EXPECT_EQ(kTRUE, envelope->HasMethod("SVM", "SVM"));

   auto &meths = envelope->GetMethods();
   for (auto &m : meths) {
      auto mname = m.GetValue<TString>("MethodName");
      auto mtitle = m.GetValue<TString>("MethodTitle");
      auto mopts = m.GetValue<TString>("MethodOptions");
      EXPECT_EQ(kTRUE, envelope->HasMethod(mname.Data(), mtitle.Data()));
      if (mname == "BDT" && mtitle == "BDTG") {
         EXPECT_EQ("NTrees=100", mopts);
      }
      if (mname == "BDT" && mtitle == "BDTB") {
         EXPECT_EQ("BoostType=Bagging", mopts);
      }
      if (mname == "SVM" && mtitle == "SVM") {
         EXPECT_EQ("Gamma=0.25:Tol=0.001:VarTransform=Norm", mopts);
      }
   }
}

// test with DataLoader and output file
class EnvelopeTest2 : public ::testing::Test {
public:
   BasicEnvelope *envelope;
   DataLoader *dl;
   TFile *inputfile;
   TFile *outfile;
   EnvelopeTest2()
   {
      TMVA::Tools::Instance();
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
      dl->PrepareTrainingAndTestTree(
         "", "", "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V");

      outfile = TFile::Open("TMVAClass.root", "RECREATE");

      envelope = new BasicEnvelope("envelope", dl, outfile, "ModelPersistence=1:V:Transformations=I:Jobs=3");
   }
   ~EnvelopeTest2()
   {
      inputfile->Close();
      outfile->Close();
      delete envelope;
   }
};

TEST_F(EnvelopeTest2, FileAndDataLoader)
{
   envelope->ParseOptions();
   EXPECT_EQ(kFALSE, envelope->IsSilentFile());
   ASSERT_NE(nullptr, envelope->GetFile());

   ASSERT_NE(nullptr, envelope->GetDataLoader());
   ASSERT_EQ("dataset", TString(envelope->GetDataLoader()->GetName()));

   EXPECT_EQ((UInt_t)6000, envelope->GetDataLoaderDataInput().GetEntries("Signal"));
   EXPECT_EQ((UInt_t)6000, envelope->GetDataLoaderDataInput().GetEntries("Background"));
   EXPECT_EQ((UInt_t)1, envelope->GetDataLoaderDataInput().GetNTrees("Signal"));
   EXPECT_EQ((UInt_t)1, envelope->GetDataLoaderDataInput().GetNTrees("Background"));

   EXPECT_EQ((UInt_t)4, envelope->GetDataLoaderDataSetInfo().GetNVariables());
   EXPECT_EQ((UInt_t)0, envelope->GetDataLoaderDataSetInfo().GetNTargets());
   EXPECT_EQ((UInt_t)2, envelope->GetDataLoaderDataSetInfo().GetNSpectators());

   envelope->SetDataLoader(nullptr);
   ASSERT_EQ(nullptr, envelope->GetDataLoader());
}
