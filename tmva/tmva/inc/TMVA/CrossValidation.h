// @(#)root/tmva $Id$
// Author: Omar Zapata, Thomas James Stevenson, Pourya Vakilipourtakalou, Kim Albertsson

/*************************************************************************
 * Copyright (C) 2018, Rene Brun and Fons Rademakers.                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMVA_CROSS_EVALUATION
#define ROOT_TMVA_CROSS_EVALUATION

#include "TGraph.h"
#include "TMultiGraph.h"
#include "TString.h"
#include <vector>
#include <map>

#include "TMVA/IMethod.h"
#include "TMVA/Configurable.h"
#include "TMVA/Types.h"
#include "TMVA/DataSet.h"
#include "TMVA/Event.h"
#include <TMVA/Results.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/OptionMap.h>
#include <TMVA/Envelope.h>

/*! \class TMVA::CrossValidationResult
 * Class to save the results of cross validation,
 * the metric for the classification ins ROC and you can ROC curves
 * ROC integrals, ROC average and ROC standard deviation.
\ingroup TMVA
*/

/*! \class TMVA::CrossValidation
 * Class to perform cross validation, splitting the dataloader into folds.
\ingroup TMVA
*/

namespace TMVA {

class CvSplitKFolds;

using EventCollection_t = std::vector<Event *>;
using EventTypes_t = std::vector<Bool_t>;
using EventOutputs_t = std::vector<Float_t>;
using EventOutputsMulticlass_t = std::vector<std::vector<Float_t>>;

class CrossValidationFoldResult {
public:
   CrossValidationFoldResult() {} // For multi-proc serialisation
   CrossValidationFoldResult(UInt_t iFold)
   : fFold(iFold)
   {}

   UInt_t fFold;

   Float_t fROCIntegral;
   TGraph fROC;

   Double_t fSig;
   Double_t fSep;
   Double_t fEff01;
   Double_t fEff10;
   Double_t fEff30;
   Double_t fEffArea;
   Double_t fTrainEff01;
   Double_t fTrainEff10;
   Double_t fTrainEff30;
};

// Used internally to keep per-fold aggregate statistics
// such as ROC curves, ROC integrals and efficiencies.
class CrossValidationResult {
   friend class CrossValidation;

private:
   std::map<UInt_t, Float_t> fROCs;
   std::shared_ptr<TMultiGraph> fROCCurves;

   std::vector<Double_t> fSigs;
   std::vector<Double_t> fSeps;
   std::vector<Double_t> fEff01s;
   std::vector<Double_t> fEff10s;
   std::vector<Double_t> fEff30s;
   std::vector<Double_t> fEffAreas;
   std::vector<Double_t> fTrainEff01s;
   std::vector<Double_t> fTrainEff10s;
   std::vector<Double_t> fTrainEff30s;

public:
   CrossValidationResult(UInt_t numFolds);
   CrossValidationResult(const CrossValidationResult &);
   ~CrossValidationResult() { fROCCurves = nullptr; }

   std::map<UInt_t, Float_t> GetROCValues() const { return fROCs; }
   Float_t GetROCAverage() const;
   Float_t GetROCStandardDeviation() const;
   TMultiGraph *GetROCCurves(Bool_t fLegend = kTRUE);
   TGraph *GetAvgROCCurve(UInt_t numSamples = 100) const;
   void Print() const;

   TCanvas *Draw(const TString name = "CrossValidation") const;
   TCanvas *DrawAvgROCCurve(Bool_t drawFolds=kFALSE, TString title="") const;

   std::vector<Double_t> GetSigValues() const { return fSigs; }
   std::vector<Double_t> GetSepValues() const { return fSeps; }
   std::vector<Double_t> GetEff01Values() const { return fEff01s; }
   std::vector<Double_t> GetEff10Values() const { return fEff10s; }
   std::vector<Double_t> GetEff30Values() const { return fEff30s; }
   std::vector<Double_t> GetEffAreaValues() const { return fEffAreas; }
   std::vector<Double_t> GetTrainEff01Values() const { return fTrainEff01s; }
   std::vector<Double_t> GetTrainEff10Values() const { return fTrainEff10s; }
   std::vector<Double_t> GetTrainEff30Values() const { return fTrainEff30s; }

private:
   void Fill(CrossValidationFoldResult const & fr);
};

class CrossValidation : public Envelope {

public:
   explicit CrossValidation(TString jobName, TMVA::DataLoader *dataloader, TString options);
   explicit CrossValidation(TString jobName, TMVA::DataLoader *dataloader, TFile *outputFile, TString options);
   ~CrossValidation();

   void InitOptions();
   void ParseOptions();

   void SetNumFolds(UInt_t i);
   void SetSplitExpr(TString splitExpr);

   UInt_t GetNumFolds() { return fNumFolds; }
   TString GetSplitExpr() { return fSplitExprString; }

   Factory &GetFactory() { return *fFactory; }

   const std::vector<CrossValidationResult> &GetResults() const;

   void Evaluate();

private:
   CrossValidationFoldResult ProcessFold(UInt_t iFold, const OptionMap & methodInfo);

   Types::EAnalysisType fAnalysisType;
   TString fAnalysisTypeStr;
   TString fSplitTypeStr;
   Bool_t fCorrelations;
   TString fCvFactoryOptions;
   Bool_t fDrawProgressBar;
   Bool_t fFoldFileOutput; //! If true: generate output file for each fold
   Bool_t fFoldStatus;     //! If true: dataset is prepared
   TString fJobName;
   UInt_t fNumFolds; //! Number of folds to prepare
   UInt_t fNumWorkerProcs; //! Number of processes to use for fold evaluation.
                           //!(Default, no parallel evaluation)
   TString fOutputFactoryOptions;
   TString fOutputEnsembling; //! How to combine output of individual folds
   TFile *fOutputFile;
   Bool_t fSilent;
   TString fSplitExprString;
   std::vector<CrossValidationResult> fResults; //!
   Bool_t fROC;
   TString fTransformations;
   Bool_t fVerbose;
   TString fVerboseLevel;

   std::unique_ptr<Factory> fFoldFactory;
   std::unique_ptr<Factory> fFactory;
   std::unique_ptr<CvSplitKFolds> fSplit;

   ClassDef(CrossValidation, 0);
   };

} // namespace TMVA

#endif // ROOT_TMVA_CROSS_EVALUATION
