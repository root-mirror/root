#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"


#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodPyRandomForest.h"

void Classification()
{
   TMVA::Tools::Instance();
   TMVA::PyMethodBase::PyInitialize();

   TString outfileName("TMVA.root");
   TFile *outputFile = TFile::Open(outfileName, "RECREATE");

   TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,
         "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

   factory->AddVariable("myvar1 := var1+var2", 'F');
   factory->AddVariable("myvar2 := var1-var2", "Expression 2", "", 'F');
   factory->AddVariable("var3",                "Variable 3", "units", 'F');
   factory->AddVariable("var4",                "Variable 4", "units", 'F');
   factory->AddSpectator("spec1 := var1*2",  "Spectator 1", "units", 'F');
   factory->AddSpectator("spec2 := var1*3",  "Spectator 2", "units", 'F');

   TFile *input(0);
   TString fname = "./tmva_class_example.root";
   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   }
   else {
      TFile::SetCacheFileDir(".");
      input = TFile::Open("http://root.cern/files/tmva_class_example.root", "CACHEREAD");
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }

   std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;

   // --- Register the training and test trees

   TTree *tsignal     = (TTree *)input->Get("TreeS");
   TTree *tbackground = (TTree *)input->Get("TreeB");

   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree(tsignal,     signalWeight);
   factory->AddBackgroundTree(tbackground, backgroundWeight);


   // Set individual event weights (the variables must exist in the original TTree)
   factory->SetBackgroundWeightExpression("weight");


   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   // Tell the factory how to use the training and testing events
   factory->PrepareTrainingAndTestTree(mycuts, mycutb,
                                       "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V");


   ///////////////////
   //Booking         //
   ///////////////////
   // Boosted Decision Trees

   //PyMVA methods
   factory->BookMethod(TMVA::Types::kPyRandomForest, "PyRandomForest",
                       "!V:NEstimators=150:Criterion=gini:MaxFeatures=auto:MaxDepth=3:MinSamplesLeaf=1:MinWeightFractionLeaf=0:Bootstrap=kTRUE");
   factory->BookMethod(TMVA::Types::kPyAdaBoost, "PyAdaBoost",
                       "!V:BaseEstimator=None:NEstimators=100:LearningRate=1:Algorithm=SAMME.R:RandomState=None");
   factory->BookMethod(TMVA::Types::kPyGTB, "PyGTB",
                       "!V:NEstimators=150:Loss=deviance:LearningRate=0.1:Subsample=1:MaxDepth=6:MaxFeatures='auto'");


   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();
   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

}
