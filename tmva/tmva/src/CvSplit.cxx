// @(#)root/tmva $Id$
// Author: Kim Albertsson

/*************************************************************************
 * Copyright (C) 2018, Rene Brun and Fons Rademakers.                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TMVA/CvSplit.h"

#include "TMVA/DataSet.h"
#include "TMVA/DataSetFactory.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Event.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Tools.h"

#include <TString.h>
#include <TFormula.h>

#include <algorithm>
#include <numeric>
#include <stdexcept>

ClassImp(TMVA::CvSplit);
ClassImp(TMVA::CvSplitKFolds);

/* =============================================================================
      TMVA::CvSplit
============================================================================= */

////////////////////////////////////////////////////////////////////////////////
///

TMVA::CvSplit::CvSplit(UInt_t numFolds) : fNumFolds(numFolds), fMakeFoldDataSet(kFALSE) {}

////////////////////////////////////////////////////////////////////////////////
/// \brief Set training and test set vectors of dataset described by `dsi`.
/// \param[in] dsi DataSetInfo for data set to be split
/// \param[in] foldNumber Ordinal of fold to prepare
/// \param[in] tt The set used to prepare fold. If equal to `Types::kTraining`
///               splitting will be based off the original train set. If instead
///               equal to `Types::kTesting` the test set will be used.
///               The original training/test set is the set as defined by
///               `DataLoader::PrepareTrainingAndTestSet`.
///
/// Sets the training and test set vectors of the DataSet described by `dsi` as
/// defined by the split. If `tt` is eqal to `Types::kTraining` the split will
/// be based off of the original training set.
///
/// Note: Requires `MakeKFoldDataSet` to have been called first.
///

void TMVA::CvSplit::PrepareFoldDataSet(DataSetInfo &dsi, UInt_t foldNumber, Types::ETreeType tt)
{
   if (foldNumber >= fNumFolds) {
      Log() << kFATAL << "DataSet prepared for \"" << fNumFolds << "\" folds, requested fold \"" << foldNumber
            << "\" is outside of range." << Endl;
      return;
   }

   auto prepareDataSetInternal = [this, &dsi, foldNumber](std::vector<std::vector<Event *>> vec) {
      UInt_t numFolds = fTrainEvents.size();

      // Events in training set (excludes current fold)
      UInt_t nTotal = std::accumulate(vec.begin(), vec.end(), 0,
                                      [&](UInt_t sum, std::vector<TMVA::Event *> v) { return sum + v.size(); });

      UInt_t nTrain = nTotal - vec.at(foldNumber).size();
      UInt_t nTest = vec.at(foldNumber).size();

      std::vector<Event *> tempTrain;
      std::vector<Event *> tempTest;

      tempTrain.reserve(nTrain);
      tempTest.reserve(nTest);

      // Insert data into training set
      for (UInt_t i = 0; i < numFolds; ++i) {
         if (i == foldNumber) {
            continue;
         }

         tempTrain.insert(tempTrain.end(), vec.at(i).begin(), vec.at(i).end());
      }

      // Insert data into test set
      tempTest.insert(tempTest.end(), vec.at(foldNumber).begin(), vec.at(foldNumber).end());

      Log() << kDEBUG << "Fold prepared, num events in training set: " << tempTrain.size() << Endl;
      Log() << kDEBUG << "Fold prepared, num events in test     set: " << tempTest.size() << Endl;

      // Assign the vectors of the events to rebuild the dataset
      dsi.GetDataSet()->SetEventCollection(&tempTrain, Types::kTraining, false);
      dsi.GetDataSet()->SetEventCollection(&tempTest, Types::kTesting, false);
   };

   if (tt == Types::kTraining) {
      prepareDataSetInternal(fTrainEvents);
   } else if (tt == Types::kTesting) {
      prepareDataSetInternal(fTestEvents);
   } else {
      Log() << kFATAL << "PrepareFoldDataSet can only work with training and testing data sets." << std::endl;
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
///

void TMVA::CvSplit::RecombineKFoldDataSet(DataSetInfo &dsi, Types::ETreeType tt)
{
   if (tt != Types::kTraining) {
      Log() << kFATAL << "Only kTraining is supported for CvSplit::RecombineKFoldDataSet currently." << std::endl;
   }

   std::vector<Event *> *tempVec = new std::vector<Event *>;

   for (UInt_t i = 0; i < fNumFolds; ++i) {
      tempVec->insert(tempVec->end(), fTrainEvents.at(i).begin(), fTrainEvents.at(i).end());
   }

   dsi.GetDataSet()->SetEventCollection(tempVec, Types::kTraining, false);
   dsi.GetDataSet()->SetEventCollection(tempVec, Types::kTesting, false);

   delete tempVec;
}

/* =============================================================================
      TMVA::CvSplitKFoldsExpr
============================================================================= */

////////////////////////////////////////////////////////////////////////////////
///

TMVA::CvSplitKFoldsExpr::CvSplitKFoldsExpr(DataSetInfo &dsi, TString expr)
   : fDsi(dsi), fIdxFormulaParNumFolds(std::numeric_limits<UInt_t>::max()), fSplitFormula("", expr),
     fParValues(fSplitFormula.GetNpar())
{
   if (!fSplitFormula.IsValid()) {
      throw std::runtime_error("Split expression \"" + std::string(fSplitExpr.Data()) + "\" is not a valid TFormula.");
   }

   for (Int_t iFormulaPar = 0; iFormulaPar < fSplitFormula.GetNpar(); ++iFormulaPar) {
      TString name = fSplitFormula.GetParName(iFormulaPar);

      // std::cout << "Found variable with name \"" << name << "\"." << std::endl;

      if (name == "NumFolds" || name == "numFolds") {
         // std::cout << "NumFolds|numFolds is a reserved variable! Adding to context." << std::endl;
         fIdxFormulaParNumFolds = iFormulaPar;
      } else {
         fFormulaParIdxToDsiSpecIdx.push_back(std::make_pair(iFormulaPar, GetSpectatorIndexForName(fDsi, name)));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///

UInt_t TMVA::CvSplitKFoldsExpr::Eval(UInt_t numFolds, const Event *ev)
{
   for (auto &p : fFormulaParIdxToDsiSpecIdx) {
      auto iFormulaPar = p.first;
      auto iSpectator = p.second;

      fParValues.at(iFormulaPar) = ev->GetSpectator(iSpectator);
   }

   if (fIdxFormulaParNumFolds < fSplitFormula.GetNpar()) {
      fParValues[fIdxFormulaParNumFolds] = numFolds;
   }

   // NOTE: We are using a double to represent an integer here. This _will_
   // lead to problems if the norm of the double grows too large. A quick test
   // with python suggests that problems arise at a magnitude of ~1e16.
   Double_t iFold_d = fSplitFormula.EvalPar(nullptr, &fParValues[0]);

   if (iFold_d < 0) {
      throw std::runtime_error("Output of splitExpr must be non-negative.");
   }

   UInt_t iFold = std::lround(iFold_d);
   if (iFold >= numFolds) {
      throw std::runtime_error("Output of splitExpr should be a non-negative"
                               "integer between 0 and numFolds-1 inclusive.");
   }

   return iFold;
}

////////////////////////////////////////////////////////////////////////////////
///

Bool_t TMVA::CvSplitKFoldsExpr::Validate(TString expr)
{
   return TFormula("", expr).IsValid();
}

////////////////////////////////////////////////////////////////////////////////
///

UInt_t TMVA::CvSplitKFoldsExpr::GetSpectatorIndexForName(DataSetInfo &dsi, TString name)
{
   std::vector<VariableInfo> spectatorInfos = dsi.GetSpectatorInfos();

   for (UInt_t iSpectator = 0; iSpectator < spectatorInfos.size(); ++iSpectator) {
      VariableInfo vi = spectatorInfos[iSpectator];
      if (vi.GetName() == name) {
         return iSpectator;
      } else if (vi.GetLabel() == name) {
         return iSpectator;
      } else if (vi.GetExpression() == name) {
         return iSpectator;
      }
   }

   throw std::runtime_error("Spectator \"" + std::string(name.Data()) + "\" not found.");
}

/* =============================================================================
      TMVA::CvSplitKFolds
============================================================================= */

////////////////////////////////////////////////////////////////////////////////
/// \brief Splits a dataset into k folds, ready for use in cross validation.
/// \param[in] numFolds Number of folds to split data into
/// \param[in] stratified If true, use stratified splitting, balancing the
///                       number of events across classes and folds. If false,
///                       no such balancing is done. For
/// \param[in] splitExpr Expression used to split data into folds. If `""` a
///                      random assignment will be done. Otherwise the
///                      expression is fed into a TFormula and evaluated per
///                      event. The resulting value is the the fold assignment.
/// \param[in] seed Used only when using random splitting (i.e. when
///                 `splitExpr` is `""`). Seed is used to initialise the random
///                 number generator when assigning events to folds.
///

TMVA::CvSplitKFolds::CvSplitKFolds(UInt_t numFolds, TString splitExpr, Bool_t stratified, UInt_t seed)
   : CvSplit(numFolds), fSeed(seed), fSplitExprString(splitExpr), fStratified(stratified)
{
   if (!CvSplitKFoldsExpr::Validate(fSplitExprString) && (splitExpr != TString(""))) {
      Log() << kFATAL << "Split expression \"" << fSplitExprString << "\" is not a valid TFormula." << Endl;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// \brief Prepares a DataSet for cross validation

void TMVA::CvSplitKFolds::MakeKFoldDataSet(DataSetInfo &dsi)
{
   // Validate spectator
   // fSpectatorIdx = GetSpectatorIndexForName(dsi, fSpectatorName);

   if (fSplitExprString != TString("")) {
      fSplitExpr = std::unique_ptr<CvSplitKFoldsExpr>(new CvSplitKFoldsExpr(dsi, fSplitExprString));
   }

   // No need to do it again if the sets have already been split.
   if (fMakeFoldDataSet) {
      Log() << kINFO << "Splitting in k-folds has been already done" << Endl;
      return;
   }

   fMakeFoldDataSet = kTRUE;

   UInt_t numClasses = dsi.GetNClasses();

   // Get the original event vectors for testing and training from the dataset.
   std::vector<Event *> trainData = dsi.GetDataSet()->GetEventCollection(Types::kTraining);
   std::vector<Event *> testData = dsi.GetDataSet()->GetEventCollection(Types::kTesting);

   // Split the sets into the number of folds.
   fTrainEvents = SplitSets(trainData, fNumFolds, numClasses);
   fTestEvents = SplitSets(testData, fNumFolds, numClasses);
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Generates a vector of fold assignments
/// \param[in] nEntries Number of events in range
/// \param[in] numFolds Number of folds to split data into
/// \param[in] seed Random seed
///
/// Randomly assigns events to `numFolds` folds. Each fold will hold at most
/// `nEntries / numFolds + 1` events.
///

std::vector<UInt_t> TMVA::CvSplitKFolds::GetEventIndexToFoldMapping(UInt_t nEntries, UInt_t numFolds, UInt_t seed)
{
   // Generate assignment of the pattern `0, 1, 2, 0, 1, 2, 0, 1 ...` for
   // `numFolds = 3`.
   std::vector<UInt_t> fOrigToFoldMapping;
   fOrigToFoldMapping.reserve(nEntries);

   for (UInt_t iEvent = 0; iEvent < nEntries; ++iEvent) {
      fOrigToFoldMapping.push_back(iEvent % numFolds);
   }

   // Shuffle assignment
   TMVA::RandomGenerator<TRandom3> rng(seed);
   std::shuffle(fOrigToFoldMapping.begin(), fOrigToFoldMapping.end(), rng);

   return fOrigToFoldMapping;
}


////////////////////////////////////////////////////////////////////////////////
/// \brief Split sets for into k-folds
/// \param[in] oldSet Original, unsplit, events
/// \param[in] numFolds Number of folds to split data into
///

std::vector<std::vector<TMVA::Event *>>
TMVA::CvSplitKFolds::SplitSets(std::vector<TMVA::Event *> &oldSet, UInt_t numFolds, UInt_t numClasses)
{
   const ULong64_t nEntries = oldSet.size();
   const ULong64_t foldSize = nEntries / numFolds;

   std::vector<std::vector<Event *>> tempSets;
   tempSets.reserve(fNumFolds);
   for (UInt_t iFold = 0; iFold < numFolds; ++iFold) {
      tempSets.emplace_back();
      tempSets.at(iFold).reserve(foldSize);
   }

   Bool_t useSplitExpr = !(fSplitExpr == nullptr || fSplitExprString == "");

   if (useSplitExpr) {
      // Deterministic split
      for (ULong64_t i = 0; i < nEntries; i++) {
         TMVA::Event *ev = oldSet[i];
         UInt_t iFold = fSplitExpr->Eval(numFolds, ev);
         tempSets.at((UInt_t)iFold).push_back(ev);
      }
   } else {
      if(!fStratified){
         // Random split
         std::vector<UInt_t> fOrigToFoldMapping;
         fOrigToFoldMapping = GetEventIndexToFoldMapping(nEntries, numFolds, fSeed);

         for (UInt_t iEvent = 0; iEvent < nEntries; ++iEvent) {
            UInt_t iFold = fOrigToFoldMapping[iEvent];
            TMVA::Event *ev = oldSet[iEvent];
            tempSets.at(iFold).push_back(ev);

            fEventToFoldMapping[ev] = iFold;
         }
      } else {
         // Stratified Split
         std::vector<std::vector<TMVA::Event *>> oldSets;
         oldSets.reserve(numClasses);

         for(UInt_t iClass = 0; iClass < numClasses; iClass++){
            oldSets.emplace_back();
            //find a way to get number of events in each class
            oldSets.reserve(nEntries);
         }

         for(UInt_t iEvent = 0; iEvent < nEntries; ++iEvent){
            // check the class of event and add to its vector of events
            TMVA::Event *ev = oldSet[iEvent];
            UInt_t iClass = ev->GetClass();
            oldSets.at(iClass).push_back(ev);
         }

         for(UInt_t i = 0; i<numClasses; ++i){
            // Shuffle each vector individually
            TMVA::RandomGenerator<TRandom3> rng(fSeed);
            std::shuffle(oldSets.at(i).begin(), oldSets.at(i).end(), rng);
         }

         for(UInt_t i = 0; i<numClasses; ++i) {
            std::vector<UInt_t> fOrigToFoldMapping;
            fOrigToFoldMapping = GetEventIndexToFoldMapping(oldSets.at(i).size(), numFolds, fSeed);

            for (UInt_t iEvent = 0; iEvent < oldSets.at(i).size(); ++iEvent) {
               UInt_t iFold = fOrigToFoldMapping[iEvent];
               TMVA::Event *ev = oldSets.at(i)[iEvent];
               tempSets.at(iFold).push_back(ev);
               fEventToFoldMapping[ev] = iFold;
            }
         }
      }
   }
   return tempSets;
}
