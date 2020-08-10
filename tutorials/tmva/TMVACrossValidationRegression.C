/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides an example of how to use TMVA for k-folds cross
/// evaluation.
///
/// As input data is used a toy-MC sample consisting of two guassian
/// distributions.
///
/// The output file "TMVA.root" can be analysed with the use of dedicated
/// macros (simply say: root -l <macro.C>), which can be conveniently
/// invoked through a GUI that will appear at the end of the run of this macro.
/// Launch the GUI via the command:
///
/// ```
/// root -l -e 'TMVA::TMVAGui("TMVA.root")'
/// ```
///
/// ## Cross Evaluation
/// Cross evaluation is a special case of k-folds cross validation where the
/// splitting into k folds is computed deterministically. This ensures that the
/// a given event will always end up in the same fold.
///
/// In addition all resulting classifiers are saved and can be applied to new
/// data using `MethodCrossValidation`. One requirement for this to work is a
/// splitting function that is evaluated for each event to determine into what
/// fold it goes (for training/evaluation) or to what classifier (for
/// application).
///
/// ## Split Expression
/// Cross evaluation uses a deterministic split to partition the data into
/// folds called the split expression. The expression can be any valid
/// `TFormula` as long as all parts used are defined.
///
/// For each event the split expression is evaluated to a number and the event
/// is put in the fold corresponding to that number.
///
/// It is recommended to always use `%int([NumFolds])` at the end of the
/// expression.
///
/// The split expression has access to all spectators and variables defined in
/// the dataloader. Additionally, the number of folds in the split can be
/// accessed with `NumFolds` (or `numFolds`).
///
/// ### Example
///  ```
///  "int(fabs([eventID]))%int([NumFolds])"
///  ```
///
/// - Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TMVACrossValidationRegression
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
#include "TMVA/CrossValidation.h"

TFile * getDataFile(TString fname) {
   TFile *input(0);

   if (!gSystem->AccessPathName(fname)) {
      input = TFile::Open(fname); // check if file in local directory exists
   } else {
      // if not: download from ROOT server
      TFile::SetCacheFileDir(".");
      input = TFile::Open("http://root.cern/files/tmva_reg_example.root", "CACHEREAD");
   }

   if (!input) {
      std::cout << "ERROR: could not open data file " << fname << std::endl;
      exit(1);
   }

   return input;
}

int TMVACrossValidationRegression()
{
   // This loads the library
   TMVA::Tools::Instance();

   // --------------------------------------------------------------------------

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName("TMVARegCv.root");
   TFile * outputFile = TFile::Open(outfileName, "RECREATE");

   TString infileName("./files/tmva_reg_example.root");
   TFile * inputFile = getDataFile(infileName);

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   
   dataloader->AddVariable("var1", "Variable 1", "units", 'F');
   dataloader->AddVariable("var2", "Variable 2", "units", 'F');

   // Add the variable carrying the regression target
   dataloader->AddTarget("fvalue");

   TTree * regTree = (TTree*)inputFile->Get("TreeR");
   dataloader->AddRegressionTree(regTree, 1.0);

   // Individual events can be weighted
   // dataloader->SetWeightExpression("weight", "Regression");

   std::cout << "--- TMVACrossValidationRegression: Using input file: " << inputFile->GetName() << std::endl;

   // Bypasses the normal splitting mechanism, CV uses a new system for this.
   // Unfortunately the old system is unhappy if we leave the test set empty so
   // we ensure that there is at least one event by placing the first event in
   // it.
   // You can with the selection cut place a global cut on the defined
   // variables. Only events passing the cut will be using in training/testing.
   // Example: `TCut selectionCut = "var1 < 1";`
   TCut selectionCut = "";
   dataloader->PrepareTrainingAndTestTree(selectionCut, "nTest_Regression=1"
                                                        ":SplitMode=Block"
                                                        ":NormMode=NumEvents"
                                                        ":!V");

   // --------------------------------------------------------------------------

   //
   // This sets up a CrossValidation class (which wraps a TMVA::Factory
   // internally) for 2-fold cross validation. The data will be split into the
   // two folds randomly if `splitExpr` is `""`.
   // 
   // One can also give a deterministic split using spectator variables. An
   // example would be e.g. `"int(fabs([spec1]))%int([NumFolds])"`.
   //
   UInt_t numFolds = 2;
   TString analysisType = "Regression";
   TString splitExpr = "";

   TString cvOptions = Form("!V"
                            ":!Silent"
                            ":ModelPersistence"
                            ":!FoldFileOutput"
                            ":AnalysisType=%s"
                            ":NumFolds=%i"
                            ":SplitExpr=%s",
                            analysisType.Data(), numFolds, splitExpr.Data());

   TMVA::CrossValidation cv{"TMVACrossValidationRegression", dataloader, outputFile, cvOptions};

   // --------------------------------------------------------------------------

   //
   // Books a method to use for evaluation
   //
   cv.BookMethod(TMVA::Types::kBDT, "BDTG",
                 "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.1:"
                 "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3");

   // --------------------------------------------------------------------------

   //
   // Train, test and evaluate the booked methods.
   // Evaluates the booked methods once for each fold and aggregates the result
   // in the specified output file.
   //
   cv.Evaluate();

   // --------------------------------------------------------------------------

   //
   // Save the output
   //
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVACrossValidationRegression is done!" << std::endl;

   // --------------------------------------------------------------------------

   //
   // Launch the GUI for the root macros
   //
   if (!gROOT->IsBatch()) {
      TMVA::TMVAGui(outfileName);
   }

   return 0;
}

//
// This is used if the macro is compiled. If run through ROOT with
// `root -l -b -q MACRO.C` or similar it is unused.
//
int main(int argc, char **argv)
{
   TMVACrossValidationRegression();
}
