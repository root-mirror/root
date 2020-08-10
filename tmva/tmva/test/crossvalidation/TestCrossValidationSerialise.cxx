// @(#)root/tmva $Id$
// Author: Kim Albertsson

/*************************************************************************
 * Copyright (C) 2017, Kim Albertsson
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TMVA/CrossValidation.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include "TMVA/Types.h"

#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"

#include <iostream>
#include <memory>
#include <stdexcept>

//
// Sets up a method for cross evaluation and generates output through the
// Factory::EvaluationAllMethods. It then recreates the method using the TMVA
// Reader and verifies that the application phase produces the same output
// using the factory produced TMVA.root as input.
//

class TestContex {
public:
   virtual ~TestContex()
   {
      delete fCrossEvaluate;
      fCrossEvaluate = nullptr;

      delete fTreeClass0;
      fTreeClass0 = nullptr;

      delete fTreeClass1;
      fTreeClass1 = nullptr;
   }

   void setUpCrossValidation(TString jobname, TMVA::Types::EMVA methodType, TString methodName, TString methodOptions);
   void runCrossValidation();
   void setUpApplicationPhase(TString jobname, TString methodName);
   void runApplicationPhase(TString methodName);
   void verifyApplicationPhase();

private:
   const TString fOutputFileName = ("TMVA.root");

   TMVA::CrossValidation *fCrossEvaluate;
   TFile fOutputFile{fOutputFileName, "RECREATE"};
   TMVA::Reader fReader{"!Color:Silent:!V"};

   TTree *fTreeClass0;
   TTree *fTreeClass1;

   Float_t fX;
   Float_t fY;
   Float_t fevNum;
   Float_t fMvaEval;

   std::vector<Float_t> fEvaluationResults;
   std::vector<Float_t> fApplicationResults;
};

// =============================================================================
// === DATAGEN ===
// =============================================================================

std::pair<Float_t, Float_t>
generateCirclePoint(Float_t r, Float_t rsig, Float_t phimin, Float_t phimax, Float_t phisig, TRandom &rng)
{
   Float_t phi = rng.Rndm() * (phimax - phimin) + phimin;
   phi += rng.Gaus() * phisig;

   r += rng.Gaus() * rsig;

   Float_t x = r * cos(TMath::DegToRad() * phi);
   Float_t y = r * sin(TMath::DegToRad() * phi);

   return std::pair<Float_t, Float_t>(x, y);
}

TTree *createCircTree(Int_t nPoints, Float_t radius, Float_t rsig, TString name, UInt_t seed = 100)
{
   TRandom rng(seed);
   std::pair<Float_t, Float_t> p;
   static UInt_t id = 0;
   Float_t x = 0;
   Float_t y = 0;

   TTree *data = new TTree(name, name);
   data->Branch("x", &x, "x/F");
   data->Branch("y", &y, "y/F");
   data->Branch("EventNumber", &id, "EventNumber/I");

   for (Int_t n = 0; n < nPoints; ++n) {
      p = generateCirclePoint(radius, rsig, 0, 160, 0.1, rng);
      x = std::get<0>(p);
      y = std::get<1>(p);
      data->Fill();

      id++;
   }

   data->ResetBranchAddresses();
   return data;
}

// === END DATAGEN ===
// =============================================================================

// =============================================================================
// === RUNNING CROSS EVALUATION ===
// =============================================================================

void TestContex::setUpCrossValidation(TString jobname, TMVA::Types::EMVA methodType, TString methodName, TString methodOptions)
{
   std::cout << "Setting cross validation up" << std::endl;
   fTreeClass0 = createCircTree(500, 3.0, 0.20, "Signal", 100);
   fTreeClass1 = createCircTree(500, 3.5, 0.20, "Background", 101);

   // Dataloader managed by CrossValidation
   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

   dataloader->AddVariable("x", 'F');
   dataloader->AddVariable("y", 'F');
   dataloader->AddSpectator("EventNumber", "EventNumber", "");

   dataloader->AddSignalTree(fTreeClass0);
   dataloader->AddBackgroundTree(fTreeClass1);

   dataloader->PrepareTrainingAndTestTree("", "", "nTest_Signal=0:nTest_Background=0)");

   fCrossEvaluate = new TMVA::CrossValidation(
      jobname, dataloader, &fOutputFile,
      "!V:Silent:!ROC:!FoldFileOutput:AnalysisType=Classification:SplitType=Deterministic:SplitExpr=int([EventNumber])%int([NumFolds])");

   fCrossEvaluate->BookMethod(methodType, methodName, methodOptions);
   std::cout << "Done setting cross validation up" << std::endl;
}

void TestContex::runCrossValidation()
{
   std::cout << "Running cross validation" << std::endl;
   fCrossEvaluate->Evaluate();
   fOutputFile.Flush();
   std::cout << "Done running crossvalidation" << std::endl;
}

// === END CROSS EVALUATION ===
// =============================================================================

// =============================================================================
// === APPLICATION PHASE ===
// =============================================================================

void TestContex::setUpApplicationPhase(TString jobname, TString methodName)
{
   std::cout << "Setting application phase up" << std::endl;

   fReader.AddVariable("x", &fX);
   fReader.AddVariable("y", &fY);
   fReader.AddSpectator("EventNumber", &fevNum);

   TString weightfile = TString("dataset/weights/") + jobname + "_" + methodName + TString(".weights.xml");
   fReader.BookMVA(methodName, weightfile);
   std::cout << "Done setting application phase up" << std::endl;
}

void TestContex::runApplicationPhase(TString methodName)
{
   std::cout << "Running application phase" << std::endl;
   TTree *tree = (TTree *)fOutputFile.Get("dataset/TestTree");
   tree->SetBranchAddress("x", &fX);
   tree->SetBranchAddress("y", &fY);
   tree->SetBranchAddress("EventNumber", &fevNum);
   tree->SetBranchAddress(methodName, &fMvaEval);

   fApplicationResults.reserve(fApplicationResults.size() + tree->GetEntries());
   for (Long64_t ievt = 0; ievt < tree->GetEntries(); ievt++) {
      tree->GetEntry(ievt);

      Float_t val = fReader.EvaluateMVA(methodName);
      fApplicationResults.push_back(val);

      fEvaluationResults.push_back(fMvaEval);
   }

   tree->ResetBranchAddresses();
   delete tree;
   std::cout << "Done running application phase" << std::endl;
}

void TestContex::verifyApplicationPhase()
{
   std::cout << "Running verification phase" << std::endl;
   if (fEvaluationResults.size() != fApplicationResults.size()) {
      throw std::runtime_error("Size not equal!");
   }

   for (UInt_t iEvent = 0; iEvent < fEvaluationResults.size(); ++iEvent) {

      if (!TMath::AreEqualAbs(fEvaluationResults[iEvent], fApplicationResults[iEvent], 1e-5)) {
         std::cout << "eval:appl[" << iEvent << "] -- " << fEvaluationResults[iEvent] << ":"
                   << fApplicationResults[iEvent] << std::endl;
         throw std::runtime_error("Output not equal!");
      }
   }
   std::cout << "Done verification phase" << std::endl;
}

// =============================================================================
// === END APPLICATION PHASE ===
// =============================================================================

void TestCvSerialise(TMVA::Types::EMVA methodType, TString methodName, TString methodOptions)
{
   TString jobname = "test_cv_serialise";

   TestContex tc{};
   tc.setUpCrossValidation(jobname, methodType, methodName, methodOptions);
   tc.runCrossValidation(); // serialises cv method if modelPersitance
   tc.setUpApplicationPhase(jobname, methodName);
   tc.runApplicationPhase(methodName);
   tc.verifyApplicationPhase();
   std::cout << "Done!" << std::endl;
}

void TestCrossValidationSerialise()
{
   {
      TString bdtOptions = "!H:!V:NTrees=10";
      TestCvSerialise(TMVA::Types::kBDT, "BDT", bdtOptions);
   }

#if (defined DNNCPU)
   {
      TString dnnTrainStrat = "TrainingStrategy=LearningRate=1e-1,"
                              "BatchSize=5,ConvergenceSteps=1";
      TString dnnOptions = "!H:!V:" + dnnTrainStrat;
      TestCvSerialise(TMVA::Types::kDNN, "DNN", dnnOptions);
   }
#endif
}

int main()
{
   TestCrossValidationSerialise();
   return 0;
}
