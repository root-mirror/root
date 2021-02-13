/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides an example of how to use TMVA for k-folds cross
/// evaluation in application.
///
/// This requires that CrossValidation was run with a deterministic split, such
/// as `"...:splitExpr=int([eventID])%int([numFolds]):..."`.
///
/// - Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TMVACrossValidationApplication
///
/// \macro_output
/// \macro_code
/// \author Kim Albertsson (adapted from code originally by Andreas Hoecker)

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
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

// Helper function to load data into TTrees.
TTree *fillTree(TTree * tree, Int_t nPoints, Double_t offset, Double_t scale, UInt_t seed = 100)
{
   TRandom3 rng(seed);
   Float_t x = 0;
   Float_t y = 0;
   Int_t eventID = 0;

   tree->SetBranchAddress("x", &x);
   tree->SetBranchAddress("y", &y);
   tree->SetBranchAddress("eventID", &eventID);

   for (Int_t n = 0; n < nPoints; ++n) {
      x = rng.Gaus(offset, scale);
      y = rng.Gaus(offset, scale);

      // For our simple example it is enough that the id's are uniformly
      // distributed and independent of the data.
      ++eventID;

      tree->Fill();
   }

   // Important: Disconnects the tree from the memory locations of x and y.
   tree->ResetBranchAddresses();
   return tree;
}

int TMVACrossValidationApplication()
{
   // This loads the library
   TMVA::Tools::Instance();

   // Set up the TMVA::Reader
   TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent:!V");

   Float_t x;
   Float_t y;
   Int_t eventID;

   reader->AddVariable("x", &x);
   reader->AddVariable("y", &y);
   reader->AddSpectator("eventID", &eventID);

   // Book the serialised methods
   TString jobname("TMVACrossValidation");
   {
      TString methodName = "BDTG";
      TString weightfile = TString("dataset/weights/") + jobname + "_" + methodName + TString(".weights.xml");

      Bool_t weightfileExists = (gSystem->AccessPathName(weightfile) == kFALSE);
      if (weightfileExists) {
         reader->BookMVA(methodName, weightfile);
      } else {
         std::cout << "Weightfile for method " << methodName << " not found."
                      " Did you run TMVACrossValidation with a specified"
                      " splitExpr?" << std::endl;
         exit(0);
      }
      
   }
   {
      TString methodName = "Fisher";
      TString weightfile = TString("dataset/weights/") + jobname + "_" + methodName + TString(".weights.xml");
      
      Bool_t weightfileExists = (gSystem->AccessPathName(weightfile) == kFALSE);
      if (weightfileExists) {
         reader->BookMVA(methodName, weightfile);
      } else {
         std::cout << "Weightfile for method " << methodName << " not found."
                      " Did you run TMVACrossValidation with a specified"
                      " splitExpr?" << std::endl;
         exit(0);
      }
   }

   // Load data
   TTree *tree = new TTree();
   tree->Branch("x", &x, "x/F");
   tree->Branch("y", &y, "y/F");
   tree->Branch("eventID", &eventID, "eventID/I");

   fillTree(tree, 1000, 1.0, 1.0, 100);
   fillTree(tree, 1000, -1.0, 1.0, 101);
   tree->SetBranchAddress("x", &x);
   tree->SetBranchAddress("y", &y);
   tree->SetBranchAddress("eventID", &eventID);

   // Prepare histograms
   Int_t nbin = 100;
   TH1F histBDTG{"BDTG", "BDTG", nbin, -1, 1};
   TH1F histFisher{"Fisher", "Fisher", nbin, -1, 1};

   // Evaluate classifiers
   for (Long64_t ievt = 0; ievt < tree->GetEntries(); ievt++) {
      tree->GetEntry(ievt);

      Double_t valBDTG = reader->EvaluateMVA("BDTG");
      Double_t valFisher = reader->EvaluateMVA("Fisher");

      histBDTG.Fill(valBDTG);
      histFisher.Fill(valFisher);
   }

   tree->ResetBranchAddresses();
   delete tree;

   if (!gROOT->IsBatch()) {
      auto c = new TCanvas();
      c->Divide(2,1);
      c->cd(1);
      histBDTG.DrawClone();
      c->cd(2); 
      histFisher.DrawClone();
   }
   else 
   { // Write histograms to output file
      TFile *target = new TFile("TMVACrossEvaluationApp.root", "RECREATE");
      histBDTG.Write();
      histFisher.Write();
      target->Close();
      delete target;
   }

   delete reader;

   return 0;
}

//
// This is used if the macro is compiled. If run through ROOT with
// `root -l -b -q MACRO.C` or similar it is unused.
//
int main(int argc, char **argv)
{
   TMVACrossValidationApplication();
}
