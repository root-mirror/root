// @(#)root/tmva $Id$
// Author: Omar Zapata, Thomas James Stevenson and Pourya Vakilipourtakalou
// Modified: Kim Albertsson 2017

/*************************************************************************
 * Copyright (C) 2018, Rene Brun and Fons Rademakers.                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TMVA/CrossValidation.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Config.h"
#include "TMVA/CvSplit.h"
#include "TMVA/DataSet.h"
#include "TMVA/Event.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodCrossValidation.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/ResultsClassification.h"
#include "TMVA/ResultsMulticlass.h"
#include "TMVA/ROCCurve.h"
#include "TMVA/Types.h"

#include "TSystem.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"

#include "ROOT/RMakeUnique.hxx"

#include <iostream>
#include <memory>

//_______________________________________________________________________
TMVA::CrossValidationResult::CrossValidationResult(UInt_t numFolds)
:fROCCurves(new TMultiGraph())
{
   fSigs.resize(numFolds);
   fSeps.resize(numFolds);
   fEff01s.resize(numFolds);
   fEff10s.resize(numFolds);
   fEff30s.resize(numFolds);
   fEffAreas.resize(numFolds);
   fTrainEff01s.resize(numFolds);
   fTrainEff10s.resize(numFolds);
   fTrainEff30s.resize(numFolds);
}

//_______________________________________________________________________
TMVA::CrossValidationResult::CrossValidationResult(const CrossValidationResult &obj)
{
   fROCs=obj.fROCs;
   fROCCurves = obj.fROCCurves;

   fSigs = obj.fSigs;
   fSeps = obj.fSeps;
   fEff01s = obj.fEff01s;
   fEff10s = obj.fEff10s;
   fEff30s = obj.fEff30s;
   fEffAreas = obj.fEffAreas;
   fTrainEff01s = obj.fTrainEff01s;
   fTrainEff10s = obj.fTrainEff10s;
   fTrainEff30s = obj.fTrainEff30s;
}

//_______________________________________________________________________
void TMVA::CrossValidationResult::Fill(CrossValidationFoldResult const & fr)
{
   UInt_t iFold = fr.fFold;

   fROCs[iFold] = fr.fROCIntegral;
   fROCCurves->Add(dynamic_cast<TGraph *>(fr.fROC.Clone()));

   fSigs[iFold] = fr.fSig;
   fSeps[iFold] = fr.fSep;
   fEff01s[iFold] = fr.fEff01;
   fEff10s[iFold] = fr.fEff10;
   fEff30s[iFold] = fr.fEff30;
   fEffAreas[iFold] = fr.fEffArea;
   fTrainEff01s[iFold] = fr.fTrainEff01;
   fTrainEff10s[iFold] = fr.fTrainEff10;
   fTrainEff30s[iFold] = fr.fTrainEff30;
}

//_______________________________________________________________________
TMultiGraph *TMVA::CrossValidationResult::GetROCCurves(Bool_t /*fLegend*/)
{
   return fROCCurves.get();
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Generates a multigraph that contains an average ROC Curve.
///
/// \note You own the returned pointer.
///
/// \param[in] numSamples Number of samples used for generating the average ROC
///                       Curve. Avg. curve will be evaluated only at these
///                       points (using interpolation if necessary).
///

TGraph *TMVA::CrossValidationResult::GetAvgROCCurve(UInt_t numSamples) const
{
   // `numSamples * increment` should equal 1.0!
   Double_t increment = 1.0 / (numSamples-1);
   std::vector<Double_t> x(numSamples), y(numSamples);

   TList *rocCurveList = fROCCurves.get()->GetListOfGraphs();

   for(UInt_t iSample = 0; iSample < numSamples; iSample++) {
      Double_t xPoint = iSample * increment;
      Double_t rocSum = 0;

      for(Int_t iGraph = 0; iGraph < rocCurveList->GetSize(); iGraph++) {
        TGraph *foldROC = static_cast<TGraph *>(rocCurveList->At(iGraph));
        rocSum += foldROC->Eval(xPoint);
      }

      x[iSample] = xPoint;
      y[iSample] = rocSum/rocCurveList->GetSize();
   }

   return new TGraph(numSamples, &x[0], &y[0]);
}

//_______________________________________________________________________
Float_t TMVA::CrossValidationResult::GetROCAverage() const
{
   Float_t avg=0;
   for(auto &roc : fROCs) {
      avg+=roc.second;
   }
   return avg/fROCs.size();
}

//_______________________________________________________________________
Float_t TMVA::CrossValidationResult::GetROCStandardDeviation() const
{
   // NOTE: We are using here the unbiased estimation of the standard deviation.
   Float_t std=0;
   Float_t avg=GetROCAverage();
   for(auto &roc : fROCs) {
      std+=TMath::Power(roc.second-avg, 2);
   }
   return TMath::Sqrt(std/float(fROCs.size()-1.0));
}

//_______________________________________________________________________
void TMVA::CrossValidationResult::Print() const
{
   TMVA::MsgLogger::EnableOutput();
   TMVA::gConfig().SetSilent(kFALSE);

   MsgLogger fLogger("CrossValidation");
   fLogger << kHEADER << " ==== Results ====" << Endl;
   for(auto &item:fROCs) {
      fLogger << kINFO << Form("Fold  %i ROC-Int : %.4f",item.first,item.second) << std::endl;
   }

   fLogger << kINFO << "------------------------" << Endl;
   fLogger << kINFO << Form("Average ROC-Int : %.4f",GetROCAverage()) << Endl;
   fLogger << kINFO << Form("Std-Dev ROC-Int : %.4f",GetROCStandardDeviation()) << Endl;

   TMVA::gConfig().SetSilent(kTRUE);
}

//_______________________________________________________________________
TCanvas* TMVA::CrossValidationResult::Draw(const TString name) const
{
   auto *c = new TCanvas(name.Data());
   fROCCurves->Draw("AL");
   fROCCurves->GetXaxis()->SetTitle(" Signal Efficiency ");
   fROCCurves->GetYaxis()->SetTitle(" Background Rejection ");
   Float_t adjust=1+fROCs.size()*0.01;
   c->BuildLegend(0.15,0.15,0.4*adjust,0.5*adjust);
   c->SetTitle("Cross Validation ROC Curves");
   c->Draw();
   return c;
}

//
TCanvas* TMVA::CrossValidationResult::DrawAvgROCCurve(Bool_t drawFolds, TString title) const
{
   // note this function will create memory leak for the TMultiGraph
   // but it needs to be kept alive in order to display the canvas
   TMultiGraph *rocs = new TMultiGraph();

   // Potentially add the folds
   if (drawFolds) {
      for (auto foldRocObj : *(*fROCCurves).GetListOfGraphs()) {
         TGraph * foldRocGraph = dynamic_cast<TGraph *>(foldRocObj->Clone());
         foldRocGraph->SetLineColor(1);
         foldRocGraph->SetLineWidth(1);
         rocs->Add(foldRocGraph);
      }
   }

   // Add the average roc curve
   TGraph *avgRocGraph = GetAvgROCCurve(100);
   avgRocGraph->SetTitle("Avg ROC Curve");
   avgRocGraph->SetLineColor(2);
   avgRocGraph->SetLineWidth(3);
   rocs->Add(avgRocGraph);

   // Draw
   TCanvas *c = new TCanvas();

   if (title != "") {
      title = "Cross Validation Average ROC Curve";
   }

   rocs->SetName("cv_rocs");
   rocs->SetTitle(title);
   rocs->GetXaxis()->SetTitle("Signal Efficiency");
   rocs->GetYaxis()->SetTitle("Background Rejection");
   rocs->DrawClone("AL");

   // Build legend
   TLegend *leg = new TLegend();
   TList *ROCCurveList = rocs->GetListOfGraphs();

   if (drawFolds) {
      Int_t nCurves = ROCCurveList->GetSize();
      leg->AddEntry(static_cast<TGraph *>(ROCCurveList->At(nCurves-1)),
                    "Avg ROC Curve", "l");
      leg->AddEntry(static_cast<TGraph *>(ROCCurveList->At(0)),
                    "Fold ROC Curves", "l");
      leg->Draw();
   } else {
      c->BuildLegend();
   }

   // Draw Canvas
   c->SetTitle("Cross Validation Average ROC Curve");
   c->Draw();
   return c;
}

/**
* \class TMVA::CrossValidation
* \ingroup TMVA
* \brief

Use html for explicit line breaking<br>
Markdown links? [class reference](#reference)?


~~~{.cpp}
ce->BookMethod(dataloader, options);
ce->Evaluate();
~~~

Cross-evaluation will generate a new training and a test set dynamically from
from `K` folds. These `K` folds are generated by splitting the input training
set. The input test set is currently ignored.

This means that when you specify your DataSet you should include all events
in your training set. One way of doing this would be the following:

~~~{.cpp}
dataloader->AddTree( signalTree, "cls1" );
dataloader->AddTree( background, "cls2" );
dataloader->PrepareTrainingAndTestTree( "", "", "nTest_cls1=1:nTest_cls2=1" );
~~~

## Split Expression
See CVSplit documentation?

*/

////////////////////////////////////////////////////////////////////////////////
///

TMVA::CrossValidation::CrossValidation(TString jobName, TMVA::DataLoader *dataloader, TFile *outputFile,
                                       TString options)
   : TMVA::Envelope(jobName, dataloader, nullptr, options),
     fAnalysisType(Types::kMaxAnalysisType),
     fAnalysisTypeStr("Auto"),
     fSplitTypeStr("Random"),
     fCorrelations(kFALSE),
     fCvFactoryOptions(""),
     fDrawProgressBar(kFALSE),
     fFoldFileOutput(kFALSE),
     fFoldStatus(kFALSE),
     fJobName(jobName),
     fNumFolds(2),
     fNumWorkerProcs(1),
     fOutputFactoryOptions(""),
     fOutputFile(outputFile),
     fSilent(kFALSE),
     fSplitExprString(""),
     fROC(kTRUE),
     fTransformations(""),
     fVerbose(kFALSE),
     fVerboseLevel(kINFO)
{
   InitOptions();
   CrossValidation::ParseOptions();
   CheckForUnusedOptions();
}

////////////////////////////////////////////////////////////////////////////////
///

TMVA::CrossValidation::CrossValidation(TString jobName, TMVA::DataLoader *dataloader, TString options)
   : CrossValidation(jobName, dataloader, nullptr, options)
{
}

////////////////////////////////////////////////////////////////////////////////
///

TMVA::CrossValidation::~CrossValidation() = default;

////////////////////////////////////////////////////////////////////////////////
///

void TMVA::CrossValidation::InitOptions()
{
   // Forwarding of Factory options
   DeclareOptionRef(fSilent, "Silent",
                    "Batch mode: boolean silent flag inhibiting any output from TMVA after the creation of the factory "
                    "class object (default: False)");
   DeclareOptionRef(fVerbose, "V", "Verbose flag");
   DeclareOptionRef(fVerboseLevel = TString("Info"), "VerboseLevel", "VerboseLevel (Debug/Verbose/Info)");
   AddPreDefVal(TString("Debug"));
   AddPreDefVal(TString("Verbose"));
   AddPreDefVal(TString("Info"));

   DeclareOptionRef(fTransformations, "Transformations",
                    "List of transformations to test; formatting example: \"Transformations=I;D;P;U;G,D\", for "
                    "identity, decorrelation, PCA, Uniform and Gaussianisation followed by decorrelation "
                    "transformations");

   DeclareOptionRef(fDrawProgressBar, "DrawProgressBar", "Boolean to show draw progress bar");
   DeclareOptionRef(fCorrelations, "Correlations", "Boolean to show correlation in output");
   DeclareOptionRef(fROC, "ROC", "Boolean to show ROC in output");

   TString analysisType("Auto");
   DeclareOptionRef(fAnalysisTypeStr, "AnalysisType",
                    "Set the analysis type (Classification, Regression, Multiclass, Auto) (default: Auto)");
   AddPreDefVal(TString("Classification"));
   AddPreDefVal(TString("Regression"));
   AddPreDefVal(TString("Multiclass"));
   AddPreDefVal(TString("Auto"));

   // Options specific to CE
   DeclareOptionRef(fSplitTypeStr, "SplitType",
                    "Set the split type (Deterministic, Random, RandomStratified) (default: Random)");
   AddPreDefVal(TString("Deterministic"));
   AddPreDefVal(TString("Random"));
   AddPreDefVal(TString("RandomStratified"));

   DeclareOptionRef(fSplitExprString, "SplitExpr", "The expression used to assign events to folds");
   DeclareOptionRef(fNumFolds, "NumFolds", "Number of folds to generate");
   DeclareOptionRef(fNumWorkerProcs, "NumWorkerProcs",
      "Determines how many processes to use for evaluation. 1 means no"
      " parallelisation. 2 means use 2 processes. 0 means figure out the"
      " number automatically based on the number of cpus available. Default"
      " 1.");

   DeclareOptionRef(fFoldFileOutput, "FoldFileOutput",
                    "If given a TMVA output file will be generated for each fold. Filename will be the same as "
                    "specifed for the combined output with a _foldX suffix. (default: false)");

   DeclareOptionRef(fOutputEnsembling = TString("None"), "OutputEnsembling",
                    "Combines output from contained methods. If None, no combination is performed. (default None)");
   AddPreDefVal(TString("None"));
   AddPreDefVal(TString("Avg"));
}

////////////////////////////////////////////////////////////////////////////////
///

void TMVA::CrossValidation::ParseOptions()
{
   this->Envelope::ParseOptions();

   if (fSplitTypeStr != "Deterministic" && fSplitExprString != "") {
      Log() << kFATAL << "SplitExpr can only be used with Deterministic Splitting" << Endl;
   }

   // Factory options
   fAnalysisTypeStr.ToLower();
   if (fAnalysisTypeStr == "classification") {
      fAnalysisType = Types::kClassification;
   } else if (fAnalysisTypeStr == "regression") {
      fAnalysisType = Types::kRegression;
   } else if (fAnalysisTypeStr == "multiclass") {
      fAnalysisType = Types::kMulticlass;
   } else if (fAnalysisTypeStr == "auto") {
      fAnalysisType = Types::kNoAnalysisType;
   }

   if (fVerbose) {
      fCvFactoryOptions += "V:";
      fOutputFactoryOptions += "V:";
   } else {
      fCvFactoryOptions += "!V:";
      fOutputFactoryOptions += "!V:";
   }

   fCvFactoryOptions += Form("VerboseLevel=%s:", fVerboseLevel.Data());
   fOutputFactoryOptions += Form("VerboseLevel=%s:", fVerboseLevel.Data());

   fCvFactoryOptions += Form("AnalysisType=%s:", fAnalysisTypeStr.Data());
   fOutputFactoryOptions += Form("AnalysisType=%s:", fAnalysisTypeStr.Data());

   if (!fDrawProgressBar) {
      fCvFactoryOptions += "!DrawProgressBar:";
      fOutputFactoryOptions += "!DrawProgressBar:";
   }

   if (fTransformations != "") {
      fCvFactoryOptions += Form("Transformations=%s:", fTransformations.Data());
      fOutputFactoryOptions += Form("Transformations=%s:", fTransformations.Data());
   }

   if (fCorrelations) {
      fCvFactoryOptions += "Correlations:";
      fOutputFactoryOptions += "Correlations:";
   } else {
      fCvFactoryOptions += "!Correlations:";
      fOutputFactoryOptions += "!Correlations:";
   }

   if (fROC) {
      fCvFactoryOptions += "ROC:";
      fOutputFactoryOptions += "ROC:";
   } else {
      fCvFactoryOptions += "!ROC:";
      fOutputFactoryOptions += "!ROC:";
   }

   if (fSilent) {
      fCvFactoryOptions += Form("Silent:");
      fOutputFactoryOptions += Form("Silent:");
   }

   // CE specific options
   if (fFoldFileOutput && fOutputFile == nullptr) {
      Log() << kFATAL << "No output file given, cannot generate per fold output." << Endl;
   }

   // Initialisations

   fFoldFactory = std::make_unique<TMVA::Factory>(fJobName, fCvFactoryOptions);

   // The fOutputFactory should always have !ModelPersistence set since we use a custom code path for this.
   //    In this case we create a special method (MethodCrossValidation) that can only be used by
   //    CrossValidation and the Reader.
   if (fOutputFile == nullptr) {
      fFactory = std::make_unique<TMVA::Factory>(fJobName, fOutputFactoryOptions);
   } else {
      fFactory = std::make_unique<TMVA::Factory>(fJobName, fOutputFile, fOutputFactoryOptions);
   }

   if(fSplitTypeStr == "Random"){
      fSplit = std::unique_ptr<CvSplitKFolds>(new CvSplitKFolds(fNumFolds, fSplitExprString, kFALSE));
   } else if(fSplitTypeStr == "RandomStratified"){
      fSplit = std::unique_ptr<CvSplitKFolds>(new CvSplitKFolds(fNumFolds, fSplitExprString, kTRUE));
   } else {
      fSplit = std::unique_ptr<CvSplitKFolds>(new CvSplitKFolds(fNumFolds, fSplitExprString));
   }

}

//_______________________________________________________________________
void TMVA::CrossValidation::SetNumFolds(UInt_t i)
{
   if (i != fNumFolds) {
      fNumFolds = i;
      fSplit = std::make_unique<CvSplitKFolds>(fNumFolds, fSplitExprString);
      fDataLoader->MakeKFoldDataSet(*fSplit);
      fFoldStatus = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
///

void TMVA::CrossValidation::SetSplitExpr(TString splitExpr)
{
   if (splitExpr != fSplitExprString) {
      fSplitExprString = splitExpr;
      fSplit = std::make_unique<CvSplitKFolds>(fNumFolds, fSplitExprString);
      fDataLoader->MakeKFoldDataSet(*fSplit);
      fFoldStatus = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluates each fold in turn.
///   - Prepares train and test data sets
///   - Trains method
///   - Evalutes on test set
///   - Stores the evaluation internally
///
/// @param iFold fold to evaluate
///

TMVA::CrossValidationFoldResult TMVA::CrossValidation::ProcessFold(UInt_t iFold, const OptionMap & methodInfo)
{
   TString methodTypeName = methodInfo.GetValue<TString>("MethodName");
   TString methodTitle = methodInfo.GetValue<TString>("MethodTitle");
   TString methodOptions = methodInfo.GetValue<TString>("MethodOptions");
   TString foldTitle = methodTitle + TString("_fold") + TString::Format("%i", iFold + 1);

   Log() << kDEBUG << "Processing  " << methodTitle << " fold " << iFold << Endl;

   // Only used if fFoldOutputFile == true
   TFile *foldOutputFile = nullptr;

   if (fFoldFileOutput && fOutputFile != nullptr) {
      TString path = gSystem->GetDirName(fOutputFile->GetName()) + "/" + foldTitle + ".root";
      foldOutputFile = TFile::Open(path, "RECREATE");
      Log() << kINFO << "Creating fold output at:" << path << Endl;
      fFoldFactory = std::make_unique<TMVA::Factory>(fJobName, foldOutputFile, fCvFactoryOptions);
   }

   fDataLoader->PrepareFoldDataSet(*fSplit, iFold, TMVA::Types::kTraining);
   MethodBase *smethod = fFoldFactory->BookMethod(fDataLoader.get(), methodTypeName, foldTitle, methodOptions);

   // Train method (train method and eval train set)
   Event::SetIsTraining(kTRUE);
   smethod->TrainMethod();
   Event::SetIsTraining(kFALSE);

   fFoldFactory->TestAllMethods();
   fFoldFactory->EvaluateAllMethods();

   TMVA::CrossValidationFoldResult result(iFold);

   // Results for aggregation (ROC integral, efficiencies etc.)
   if (fAnalysisType == Types::kClassification || fAnalysisType == Types::kMulticlass) {
      result.fROCIntegral = fFoldFactory->GetROCIntegral(fDataLoader->GetName(), foldTitle);

      TGraph *gr = fFoldFactory->GetROCCurve(fDataLoader->GetName(), foldTitle, true);
      gr->SetLineColor(iFold + 1);
      gr->SetLineWidth(2);
      gr->SetTitle(foldTitle.Data());
      result.fROC = *gr;

      result.fSig = smethod->GetSignificance();
      result.fSep = smethod->GetSeparation();

      if (fAnalysisType == Types::kClassification) {
         Double_t err;
         result.fEff01 = smethod->GetEfficiency("Efficiency:0.01", Types::kTesting, err);
         result.fEff10 = smethod->GetEfficiency("Efficiency:0.10", Types::kTesting, err);
         result.fEff30 = smethod->GetEfficiency("Efficiency:0.30", Types::kTesting, err);
         result.fEffArea = smethod->GetEfficiency("", Types::kTesting, err);
         result.fTrainEff01 = smethod->GetTrainingEfficiency("Efficiency:0.01");
         result.fTrainEff10 = smethod->GetTrainingEfficiency("Efficiency:0.10");
         result.fTrainEff30 = smethod->GetTrainingEfficiency("Efficiency:0.30");
      } else if (fAnalysisType == Types::kMulticlass) {
         // Nothing here for now
      }
   }

   // Per-fold file output
   if (fFoldFileOutput && foldOutputFile != nullptr) {
      foldOutputFile->Close();
   }

   // Clean-up for this fold
   {
      smethod->Data()->DeleteAllResults(Types::kTraining, smethod->GetAnalysisType());
      smethod->Data()->DeleteAllResults(Types::kTesting, smethod->GetAnalysisType());
   }

   fFoldFactory->DeleteAllMethods();
   fFoldFactory->fMethodsMap.clear();

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Does training, test set evaluation and performance evaluation of using
/// cross-evalution.
///

void TMVA::CrossValidation::Evaluate()
{
   // Generate K folds on given dataset
   if (!fFoldStatus) {
      fDataLoader->MakeKFoldDataSet(*fSplit);
      fFoldStatus = kTRUE;
   }

   fResults.reserve(fMethods.size());
   for (auto & methodInfo : fMethods) {
      CrossValidationResult result{fNumFolds};

      TString methodTypeName = methodInfo.GetValue<TString>("MethodName");
      TString methodTitle = methodInfo.GetValue<TString>("MethodTitle");

      if (methodTypeName == "") {
         Log() << kFATAL << "No method booked for cross-validation" << Endl;
      }

      TMVA::MsgLogger::EnableOutput();
      Log() << kINFO << Endl;
      Log() << kINFO << Endl;
      Log() << kINFO << "========================================" << Endl;
      Log() << kINFO << "Processing folds for method " << methodTitle << Endl;
      Log() << kINFO << "========================================" << Endl;
      Log() << kINFO << Endl;

      // Process K folds
      auto nWorkers = fNumWorkerProcs;
      if (nWorkers == 1) {
         // Fall back to global config
         nWorkers = TMVA::gConfig().GetNumWorkers();
      }
      if (nWorkers == 1) {
         for (UInt_t iFold = 0; iFold < fNumFolds; ++iFold) {
            auto fold_result = ProcessFold(iFold, methodInfo);
            result.Fill(fold_result);
         }
      } else {
#ifndef _MSC_VER
         ROOT::TProcessExecutor workers(nWorkers);
         std::vector<CrossValidationFoldResult> result_vector;

         auto workItem = [this, methodInfo](UInt_t iFold) {
            return ProcessFold(iFold, methodInfo);
         };

         result_vector = workers.Map(workItem, ROOT::TSeqI(fNumFolds));

         for (auto && fold_result : result_vector) {
            result.Fill(fold_result);
         }
#endif
      }

      fResults.push_back(result);

      // Serialise the cross evaluated method
      TString options =
         Form("SplitExpr=%s:NumFolds=%i"
              ":EncapsulatedMethodName=%s"
              ":EncapsulatedMethodTypeName=%s"
              ":OutputEnsembling=%s",
              fSplitExprString.Data(), fNumFolds, methodTitle.Data(), methodTypeName.Data(), fOutputEnsembling.Data());

      fFactory->BookMethod(fDataLoader.get(), Types::kCrossValidation, methodTitle, options);

      // Feed EventToFold mapping used when random fold assignments are used
      // (when splitExpr="").
      IMethod *method_interface = fFactory->GetMethod(fDataLoader->GetName(), methodTitle);
      auto *method = dynamic_cast<MethodCrossValidation *>(method_interface);

      method->fEventToFoldMapping = fSplit->fEventToFoldMapping;
   }

   Log() << kINFO << Endl;
   Log() << kINFO << Endl;
   Log() << kINFO << "========================================" << Endl;
   Log() << kINFO << "Folds processed for all methods, evaluating." << Endl;
   Log() << kINFO << "========================================" << Endl;
   Log() << kINFO << Endl;

   // Recombination of data (making sure there is data in training and testing trees).
   fDataLoader->RecombineKFoldDataSet(*fSplit);

   // "Eval" on training set
   for (auto & methodInfo : fMethods) {
      TString methodTypeName = methodInfo.GetValue<TString>("MethodName");
      TString methodTitle = methodInfo.GetValue<TString>("MethodTitle");

      IMethod *method_interface = fFactory->GetMethod(fDataLoader->GetName(), methodTitle);
      auto method = dynamic_cast<MethodCrossValidation *>(method_interface);

      if (fOutputFile != nullptr) {
         fFactory->WriteDataInformation(method->fDataSetInfo);
      }

      Event::SetIsTraining(kTRUE);
      method->TrainMethod();
      Event::SetIsTraining(kFALSE);
   }

   // Eval on Testing set
   fFactory->TestAllMethods();

   // Calc statistics
   fFactory->EvaluateAllMethods();

   Log() << kINFO << "Evaluation done." << Endl;
}

//_______________________________________________________________________
const std::vector<TMVA::CrossValidationResult> &TMVA::CrossValidation::GetResults() const
{
   if (fResults.empty()) {
      Log() << kFATAL << "No cross-validation results available" << Endl;
   }
   return fResults;
}
