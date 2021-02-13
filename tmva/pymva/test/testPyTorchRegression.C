#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/PyMethodBase.h"


int testPyTorchRegression(){
   // Get data file
   std::cout << "Get test data..." << std::endl;
   TString fname = "./tmva_reg_example.root";
   if (gSystem->AccessPathName(fname))  // file does not exist in local directory
      gSystem->Exec("curl -O http://root.cern.ch/files/tmva_reg_example.root");
   TFile *input = TFile::Open(fname);

   // Build model from python file
   std::cout << "Generate PyTorch model..." << std::endl;
   UInt_t ret;
   ret = gSystem->Exec("python generatePyTorchModelRegression.py");
   if(ret!=0){
       std::cout << "[ERROR] Failed to generate model using python" << std::endl;
       return 1;
   }

   // Setup PyMVA and factory
   std::cout << "Setup TMVA..." << std::endl;
   TMVA::PyMethodBase::PyInitialize();
   TFile* outputFile = TFile::Open("ResultsTestPyTorchRegression.root", "RECREATE");
   TMVA::Factory *factory = new TMVA::Factory("testPyTorchRegression", outputFile,
      "!V:Silent:Color:!DrawProgressBar:AnalysisType=Regression");

   // Load data
   TMVA::DataLoader *dataloader = new TMVA::DataLoader("datasetTestPyTorchRegression");

   TTree *tree = (TTree*)input->Get("TreeR");
   dataloader->AddRegressionTree(tree);

   dataloader->AddVariable("var1");
   dataloader->AddVariable("var2");
   dataloader->AddTarget("fvalue");

   dataloader->PrepareTrainingAndTestTree("",
      "SplitMode=Random:NormMode=NumEvents:!V");

   // Book and train method
   factory->BookMethod(dataloader, TMVA::Types::kPyTorch, "PyTorch",
      "!H:!V:VarTransform=D,G:FilenameModel=PyTorchModelRegression.pt:FilenameTrainedModel=trainedPyTorchModelRegression.h5:NumEpochs=10:BatchSize=32:SaveBestOnly=false:UserCode=generatePyTorchModelRegression.py");
   std::cout << "Train model..." << std::endl;
   factory->TrainAllMethods();

   // Clean-up
   delete factory;
   delete dataloader;
   delete outputFile;

   // Setup reader
   UInt_t numEvents = 100;
   std::cout << "Run reader and estimate target of " << numEvents << " events..." << std::endl;
   TMVA::Reader *reader = new TMVA::Reader("!Color:Silent");
   Float_t vars[3];
   reader->AddVariable("var1", vars+0);
   reader->AddVariable("var2", vars+1);
   reader->BookMVA("PyTorch", "datasetTestPyTorchRegression/weights/testPyTorchRegression_PyTorch.weights.xml");

   // Get mean squared error on events
   tree->SetBranchAddress("var1", vars+0);
   tree->SetBranchAddress("var2", vars+1);
   tree->SetBranchAddress("fvalue", vars+2);

   Float_t meanMvaError = 0;
   for(UInt_t i=0; i<numEvents; i++){
      tree->GetEntry(i);
      meanMvaError += std::pow(vars[2]-reader->EvaluateMVA("PyTorch"),2);
   }
   meanMvaError = meanMvaError/float(numEvents);

   // Check whether the response is obviously better than guessing
   std::cout << "Mean squared error: " << meanMvaError << std::endl;
   if(meanMvaError > 30.0){
      std::cout << "[ERROR] Mean squared error is " << meanMvaError << " (>30.0)" << std::endl;
      return 1;
   }

   return 0;
}

int main(){
   int err = testPyTorchRegression();
   return err;
}
