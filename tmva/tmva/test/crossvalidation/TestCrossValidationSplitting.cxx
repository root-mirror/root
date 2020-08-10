/// \file
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Exectuable: TMVACrossValidation
///
///
/// Performs a verification that the cross evaluation splitting was performed as
/// intended. If the input data has ids (EventNumbers)
///
///    0,1,2,3,4,5,6,7,8,9,10
///
/// these will be split into e.g. 3 parts
///
///    Part 0: 0,3,6,9
///    Part 1: 1,4,7,10
///    Part 2: 2,5,8
///
/// This file verifies that each fold is calculated so that the test data for
/// each fold is equal to the equivalent part. That is the test set for fold 0
/// is part 0 and the train set is part 1 + part 2.
///
///    Fold 0:
///       Train: 1,2,4,5,7,8,10,
///       Test : 0,3,6,9
///
/// Et.c.
///

#include "gtest/gtest.h"

#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>

#include "TMVA/CvSplit.h"
#include "TMVA/DataLoader.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/CrossValidation.h"
#include "TMVA/Tools.h"

#include <algorithm>
#include <vector>

using id_vec_t = std::vector<UInt_t>;
using data_t = std::tuple<id_vec_t, TTree *>;
using fold_id_vec_t = std::vector<std::shared_ptr<id_vec_t>>;

namespace TMVA {

/*
 * Creates data for use by the test. Returns a tuple of a TTree and an id
 * vector. The TTree can be used by a dataloader to create a DataSet. The id
 * vector can be used to verify that the splitting was done correctly.
 *
 * \param nPoints Number of data points to generate
 * \param start   Start value for the id, can be used to differentiate trees.
 */
data_t createData(Int_t nPoints, UInt_t start)
{
   std::vector<UInt_t> ids;

   UInt_t id = start;
   Double_t x = 0.;

   TTree *data = new TTree();
   data->Branch("x", &x, "x/D");
   data->Branch("id", &id, "id/I");

   for (Int_t n = 0; n < nPoints; ++n) {
      data->Fill();
      ids.push_back(id);
      x += 0.1;
      id++;
   }
   data->ResetBranchAddresses();

   return std::make_tuple(ids, data);
}

/**
 * Performs the same split as CvSplitKFolds
 *
 * @param ids
 */
fold_id_vec_t getCurrentFoldExternal(id_vec_t ids, UInt_t numFolds, UInt_t iFold)
{
   std::cout << "Entering getCurrentFoldExternal" << std::endl;

   // Generate the individual folds
   fold_id_vec_t fold_vec;
   for (size_t i = 0; i < numFolds; ++i) {
      fold_vec.push_back(std::shared_ptr<id_vec_t>(new id_vec_t()));
   }

   for (auto &val : ids) {
      fold_vec[val % numFolds]->push_back(val);
   }

   for (auto &vec : fold_vec) {
      std::sort(vec->begin(), vec->end());
   }

   for (size_t k = 0; k < numFolds; ++k) {
      std::cout << "Events in fold " << k << ": ";
      for (size_t i = 0; i < fold_vec[k]->size(); ++i) {
         std::cout << fold_vec[k]->at(i) << ", ";
      }
      std::cout << std::endl;
   }

   // Combine folds into a a training and test set
   fold_id_vec_t combined_vec;
   for (size_t i = 0; i < 2; ++i) {
      combined_vec.push_back(std::shared_ptr<id_vec_t>(new id_vec_t()));
   }

   for (size_t i = 0; i < numFolds; ++i) {
      auto fold = fold_vec[i];
      if (i != iFold) {
         combined_vec[0]->insert(combined_vec[0]->end(), fold->begin(), fold->end());
      } else {
         // Fold number iFold is kept as test set
         combined_vec[1]->insert(combined_vec[1]->end(), fold->begin(), fold->end());
      }
   }

   for (auto &vec : combined_vec) {
      std::sort(vec->begin(), vec->end());
   }

   // Print the contents of the sets
   std::cout << "Events in training: ";
   for (size_t i = 0; i < combined_vec[0]->size(); ++i) {
      std::cout << combined_vec[0]->at(i) << ", ";
   }
   std::cout << std::endl;

   std::cout << "Events in testing : ";
   for (size_t i = 0; i < combined_vec[1]->size(); ++i) {
      std::cout << combined_vec[1]->at(i) << ", ";
   }
   std::cout << std::endl;

   return combined_vec;
}

fold_id_vec_t getCurrentFold(DataSet *ds)
{
   std::cout << "Entering getCurrentFold" << std::endl;

   fold_id_vec_t fold_vec;
   for (size_t i = 0; i < 2; ++i) {
      fold_vec.push_back(std::shared_ptr<id_vec_t>(new id_vec_t()));
   }

   for (auto &ev : ds->GetEventCollection(Types::kTraining)) {
      fold_vec[0]->push_back(ev->GetSpectators().at(0));
   }
   for (auto &ev : ds->GetEventCollection(Types::kTesting)) {
      fold_vec[1]->push_back(ev->GetSpectators().at(0));
   }

   for (auto &vec : fold_vec) {
      std::sort(vec->begin(), vec->end());
   }

   // Print the contents of the sets
   std::cout << "Events in training: ";
   for (size_t i = 0; i < fold_vec[0]->size(); ++i) {
      std::cout << fold_vec[0]->at(i) << ", ";
   }
   std::cout << std::endl;

   std::cout << "Events in testing : ";
   for (size_t i = 0; i < fold_vec[1]->size(); ++i) {
      std::cout << fold_vec[1]->at(i) << ", ";
   }
   std::cout << std::endl;

   return fold_vec;
}

void verifySplitExternal(DataSet *ds, id_vec_t ids, UInt_t numFolds, UInt_t iFold)
{
   fold_id_vec_t fold_vec_ext = getCurrentFoldExternal(ids, numFolds, iFold);
   fold_id_vec_t fold_vec = getCurrentFold(ds);

   auto training_ext = fold_vec_ext[0];
   auto training = fold_vec[0];
   auto test_ext = fold_vec_ext[1];
   auto test = fold_vec[1];

   EXPECT_EQ(training_ext->size(), training->size());
   EXPECT_EQ(test_ext->size(), test->size());

   std::cout << "training_ext->size(): " << training_ext->size() << std::endl;
   std::cout << "training->size()    : " << training->size() << std::endl;
   std::cout << "test_ext->size()    : " << test_ext->size() << std::endl;
   std::cout << "test->size()        : " << test->size() << std::endl;

   // Verify Training set
   for (size_t iEvent = 0; iEvent < training->size(); ++iEvent) {
      EXPECT_EQ(training_ext->at(iEvent), training->at(iEvent));
   }

   // Verify Test set
   for (size_t iEvent = 0; iEvent < test->size(); ++iEvent) {
      EXPECT_EQ(test_ext->at(iEvent), test->at(iEvent));
   }
}

/*
 * Checks wether a fold has been prepared successfully. Only does the split on the
 * training set, (d->PrepareFoldDataSet(iFold, Types::kTraining);)
 */
bool testFold(DataLoader *d, id_vec_t ids, CvSplit &split, UInt_t iFold)
{
   DataSet *ds = d->GetDataSetInfo().GetDataSet();
   d->PrepareFoldDataSet(split, iFold, Types::kTraining);

   verifySplitExternal(ds, ids, split.GetNumFolds(), iFold);

   return true;
}

/*
 * Checks that the spread of the number of events of a particular class is at
 * most 1 over all the folds. This is the core of the stratified splitting.
 */
bool testStratified(DataLoader *d, CvSplit &split, UInt_t numFolds)
{
   DataSet *ds = d->GetDataSetInfo().GetDataSet();

   std::vector<UInt_t> nSigFolds;
   std::vector<UInt_t> nBkgFolds;

   for (UInt_t iFold = 0; iFold < numFolds; ++iFold) {
      d->PrepareFoldDataSet(split, iFold, Types::kTraining);

      // Get the number events per class in a fold
      UInt_t nSignal = 0;
      UInt_t nBackground = 0;
      UInt_t nTotal = 0;
      for (auto &ev : ds->GetEventCollection(Types::kTesting)) {
         UInt_t classid = ev->GetClass();
         if (classid == d->GetDataSetInfo().GetSignalClassIndex()) {
            ++nSignal;
         } else {
            ++nBackground;
         }
         ++nTotal;
      }

      nSigFolds.push_back(nSignal);
      nBkgFolds.push_back(nBackground);

      std::cout << "Stats for fold " << iFold << " sig/bkg/tot: " << nSignal
                << "/" << nBackground << "/" << nTotal << std::endl;
   }

   // Check the spread
   Int_t minSig = *std::min_element(nSigFolds.begin(), nSigFolds.end());
   Int_t maxSig = *std::max_element(nSigFolds.begin(), nSigFolds.end());
   Int_t minBkg = *std::min_element(nBkgFolds.begin(), nBkgFolds.end());
   Int_t maxBkg = *std::max_element(nBkgFolds.begin(), nBkgFolds.end());

   EXPECT_LE((maxSig-minSig), 1);
   EXPECT_LE((maxBkg-minBkg), 1);

   return true;
}

} // End namespace TMVA

TEST(CrossValidationSplitting, TrainingSetSplitOnSpectator)
{
   TMVA::Tools::Instance();

   const UInt_t NUM_FOLDS = 3;
   const UInt_t nPointsSig = 11;
   const UInt_t nPointsBkg = 10;

   // Create DataSet
   TMVA::MsgLogger::InhibitOutput();
   data_t data_class0 = TMVA::createData(nPointsSig, 0);
   data_t data_class1 = TMVA::createData(nPointsBkg, 100);

   id_vec_t ids;
   ids.insert(ids.end(), std::get<0>(data_class0).begin(), std::get<0>(data_class0).end());
   ids.insert(ids.end(), std::get<0>(data_class1).begin(), std::get<0>(data_class1).end());

   TMVA::DataLoader *d = new TMVA::DataLoader("dataset");

   d->AddSignalTree(std::get<1>(data_class0));
   d->AddBackgroundTree(std::get<1>(data_class1));

   d->AddVariable("x", 'D');
   d->AddSpectator("id", "id", "");
   d->PrepareTrainingAndTestTree(
      "", Form("SplitMode=Block:nTrain_Signal=%i:nTrain_Background=%i:!V", nPointsSig, nPointsBkg));

   d->GetDataSetInfo().GetDataSet(); // Force creation of dataset.
   TMVA::MsgLogger::EnableOutput();

   TMVA::CvSplitKFolds split{NUM_FOLDS, "int([id])%int([numFolds])", kFALSE, 0};
   d->MakeKFoldDataSet(split);

   // Actual test
   testFold(d, ids, split, 0);
   testFold(d, ids, split, 1);
}

TEST(CrossValidationSplitting, TrainingSetSplitRandomStratified)
{
   TMVA::Tools::Instance();

   // Test for unbalanced classes
   const UInt_t NUM_FOLDS = 3;
   const UInt_t nPointsSig = 110;
   const UInt_t nPointsBkg = 10;

   // Create DataSet
   TMVA::MsgLogger::InhibitOutput();
   data_t data_class0 = TMVA::createData(nPointsSig, 0);
   data_t data_class1 = TMVA::createData(nPointsBkg, 100);

   TMVA::DataLoader *d = new TMVA::DataLoader("dataset");

   d->AddSignalTree(std::get<1>(data_class0));
   d->AddBackgroundTree(std::get<1>(data_class1));

   d->AddVariable("x", 'D');
   d->AddSpectator("id", "id", "");
   d->PrepareTrainingAndTestTree(
      "", Form("SplitMode=Block:nTrain_Signal=%i:nTrain_Background=%i:!V", nPointsSig, nPointsBkg));

   d->GetDataSetInfo().GetDataSet(); // Force creation of dataset.
   TMVA::MsgLogger::EnableOutput();

   TMVA::CvSplitKFolds split{NUM_FOLDS, "", kTRUE, 0};
   d->MakeKFoldDataSet(split);

   // Actual test
   testStratified(d, split, NUM_FOLDS);
}
