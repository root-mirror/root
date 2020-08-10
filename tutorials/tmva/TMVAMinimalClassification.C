/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// Minimal self-contained example for setting up TMVA with binary
/// classification.
///
/// This is intended as a simple foundation to build on. It assumes you are
/// familiar with TMVA already. As such concepts like the Factory, the DataLoader
/// and others are not explained. For descriptions and tutuorials use the TMVA
/// User's Guide (https://root.cern/root-user-guides-and-manuals under TMVA)
/// or the more detailed examples provided with TMVA e.g. TMVAClassification.C.
///
/// Sets up a minimal binary classification example with two slighly overlapping
/// 2-D gaussian distributions and trains a BDT classifier to discriminate the
/// data.
///
/// - Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TMVAMinimalClassification.C
///
/// \macro_output
/// \macro_code
/// \author Kim Albertsson

#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

//
// Helper function to generate 2-D gaussian data points and fill to a ROOT
// TTree.
//
// Arguments:
//    nPoints Number of points to generate.
//    offset  Mean of the generated numbers
//    scale   Standard deviation of the generated numbers.
//    seed    Seed for random number generator. Use `seed=0` for random
//            seed.
// Returns a TTree ready to be used as input to TMVA.
//
TTree *genTree(Int_t nPoints, Double_t offset, Double_t scale, UInt_t seed = 100)
{
   TRandom rng(seed);
   Double_t x = 0;
   Double_t y = 0;

   TTree *data = new TTree();
   data->Branch("x", &x, "x/D");
   data->Branch("y", &y, "y/D");

   for (Int_t n = 0; n < nPoints; ++n) {
      x = rng.Rndm() * scale;
      y = offset + rng.Rndm() * scale;
      data->Fill();
   }

   // Important: Disconnects the tree from the memory locations of x and y.
   data->ResetBranchAddresses();
   return data;
}

//
// Minimal setup for perfroming binary classification in TMVA.
//
// Modify the setup to your liking and run with
//    `root -l -b -q TMVAMinimalClassification.C`.
// This will generate an output file "out.root" that can be viewed with
//    `root -l -e 'TMVA::TMVAGui("out.root")'`.
//
void TMVAMinimalClassification()
{
   TString outputFilename = "out.root";
   TFile *outFile = new TFile(outputFilename, "RECREATE");

   // Data generatration
   TTree *signalTree = genTree(1000, 0.0, 2.0, 100);
   TTree *backgroundTree = genTree(1000, 1.0, 2.0, 101);

   TString factoryOptions = "AnalysisType=Classification";
   TMVA::Factory factory{"", outFile, factoryOptions};

   TMVA::DataLoader dataloader{"dataset"};

   // Data specification
   dataloader.AddVariable("x", 'D');
   dataloader.AddVariable("y", 'D');

   dataloader.AddSignalTree(signalTree, 1.0);
   dataloader.AddBackgroundTree(backgroundTree, 1.0);

   TCut signalCut = "";
   TCut backgroundCut = "";
   TString datasetOptions = "SplitMode=Random";
   dataloader.PrepareTrainingAndTestTree(signalCut, backgroundCut, datasetOptions);

   // Method specification
   TString methodOptions = "";
   factory.BookMethod(&dataloader, TMVA::Types::kBDT, "BDT", methodOptions);

   // Training and Evaluation
   factory.TrainAllMethods();
   factory.TestAllMethods();
   factory.EvaluateAllMethods();

   // Clean up
   outFile->Close();

   delete outFile;
   delete signalTree;
   delete backgroundTree;
}
