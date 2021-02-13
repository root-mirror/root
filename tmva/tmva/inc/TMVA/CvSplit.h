// @(#)root/tmva $Id$
// Author: Kim Albertsson

/*************************************************************************
 * Copyright (C) 2018, Rene Brun and Fons Rademakers.                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMVA_CvSplit
#define ROOT_TMVA_CvSplit

#include "TMVA/Configurable.h"
#include "TMVA/Types.h"

#include <Rtypes.h>
#include <TFormula.h>

#include <memory>
#include <vector>
#include <map>

class TString;

namespace TMVA {

class CrossValidation;
class DataSetInfo;
class Event;

/* =============================================================================
      TMVA::CvSplit
============================================================================= */

class CvSplit : public Configurable {
public:
   CvSplit(UInt_t numFolds);
   virtual ~CvSplit() {}

   virtual void MakeKFoldDataSet(DataSetInfo &dsi) = 0;
   virtual void PrepareFoldDataSet(DataSetInfo &dsi, UInt_t foldNumber, Types::ETreeType tt);
   virtual void RecombineKFoldDataSet(DataSetInfo &dsi, Types::ETreeType tt = Types::kTraining);

   UInt_t GetNumFolds() { return fNumFolds; }
   Bool_t NeedsRebuild() { return fMakeFoldDataSet; }

protected:
   UInt_t fNumFolds;
   Bool_t fMakeFoldDataSet;

   std::vector<std::vector<TMVA::Event *>> fTrainEvents;
   std::vector<std::vector<TMVA::Event *>> fTestEvents;

protected:
   ClassDef(CvSplit, 0);
};

/* =============================================================================
      TMVA::CvSplitKFoldsExpr
============================================================================= */

class CvSplitKFoldsExpr {
public:
   CvSplitKFoldsExpr(DataSetInfo &dsi, TString expr);
   ~CvSplitKFoldsExpr() {}

   UInt_t Eval(UInt_t numFolds, const Event *ev);

   static Bool_t Validate(TString expr);

private:
   UInt_t GetSpectatorIndexForName(DataSetInfo &dsi, TString name);

private:
   DataSetInfo &fDsi;

   std::vector<std::pair<Int_t, Int_t>>
      fFormulaParIdxToDsiSpecIdx; //! Maps parameter indicies in splitExpr to their spectator index in the datasetinfo.
   Int_t fIdxFormulaParNumFolds;  //! Keeps track of the index of reserved par "NumFolds" in splitExpr.
   TString fSplitExpr;     //! Expression used to split data into folds. Should output values between 0 and numFolds.
   TFormula fSplitFormula; //! TFormula for splitExpr.

   std::vector<Double_t> fParValues;
};

/* =============================================================================
      TMVA::CvSplitKFolds
============================================================================= */

class CvSplitKFolds : public CvSplit {

   friend CrossValidation;

public:
   CvSplitKFolds(UInt_t numFolds, TString splitExpr = "", Bool_t stratified = kTRUE, UInt_t seed = 100);
   ~CvSplitKFolds() override {}

   void MakeKFoldDataSet(DataSetInfo &dsi) override;

private:
   std::vector<std::vector<Event *>> SplitSets(std::vector<TMVA::Event *> &oldSet, UInt_t numFolds, UInt_t numClasses);
   std::vector<UInt_t> GetEventIndexToFoldMapping(UInt_t nEntries, UInt_t numFolds, UInt_t seed = 100);

private:
   UInt_t fSeed;
   TString fSplitExprString; //! Expression used to split data into folds. Should output values between 0 and numFolds.
   std::unique_ptr<CvSplitKFoldsExpr> fSplitExpr;
   Bool_t fStratified; // If true, use stratified split. (Balance class presence in each fold).

   // Used for CrossValidation with random splits (not using the
   // CVSplitKFoldsExpr functionality) to communicate Event to fold mapping.
   std::map<const TMVA::Event *, UInt_t> fEventToFoldMapping;

private:
   ClassDefOverride(CvSplitKFolds, 0);
};

} // end namespace TMVA

#endif