// @(#)root/roostats:$Id$
// Authors: Kevin Belasco        17/06/2009
// Authors: Kyle Cranmer         17/06/2009
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////


#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

#ifndef ROOT_TMath
#include "TMath.h"
#endif

#ifndef RooStats_MCMCInterval
#include "RooStats/MCMCInterval.h"
#endif
#ifndef ROOSTATS_MarkovChain
#include "RooStats/MarkovChain.h"
#endif
#ifndef RooStats_Heaviside
#include "RooStats/Heaviside.h"
#endif
#ifndef ROO_DATA_HIST
#include "RooDataHist.h"
#endif
#ifndef ROO_KEYS_PDF
#include "RooNDKeysPdf.h"
#endif
#ifndef ROO_PRODUCT
#include "RooProduct.h"
#endif
#ifndef RooStats_RooStatsUtils
#include "RooStats/RooStatsUtils.h"
#endif
#ifndef ROO_REAL_VAR
#include "RooRealVar.h"
#endif
#ifndef ROO_ARG_LIST
#include "RooArgList.h"
#endif
#ifndef ROOT_TIterator
#include "TIterator.h"
#endif
#ifndef ROOT_TH1
#include "TH1.h"
#endif
#ifndef ROOT_TH1F
#include "TH1F.h"
#endif
#ifndef ROOT_TH2F
#include "TH2F.h"
#endif
#ifndef ROOT_TH3F
#include "TH3F.h"
#endif
#ifndef ROO_MSG_SERVICE
#include "RooMsgService.h"
#endif
#ifndef ROO_GLOBAL_FUNC
#include "RooGlobalFunc.h"
#endif
#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_THnSparse
#include "THnSparse.h"
#endif
#ifndef ROO_NUMBER
#include "RooNumber.h"
#endif
//#ifndef ROOT_TFile
//#include "TFile.h"
//#endif

#include <cstdlib>
#include <string>
#include <algorithm>

ClassImp(RooStats::MCMCInterval);

using namespace RooFit;
using namespace RooStats;
using namespace std;

static const Double_t DEFAULT_EPSILON = 0.01;
static const Double_t DEFAULT_DELTA   = 10e-6;

MCMCInterval::MCMCInterval(const char* name)
   : ConfInterval(name)
{
   fConfidenceLevel = 0.0;
   fHistConfLevel = 0.0;
   fKeysConfLevel = 0.0;
   fTFConfLevel = 0.0;
   fFull = 0.0;
   fChain = NULL;
   fAxes = NULL;
   fDataHist = NULL;
   fSparseHist = NULL;
   fVector.clear();
   fKeysPdf = NULL;
   fProduct = NULL;
   fHeaviside = NULL;
   fKeysDataHist = NULL;
   fCutoffVar = NULL;
   fHist = NULL;
   fNumBurnInSteps = 0;
   fHistCutoff = -1;
   fKeysCutoff = -1;
   fDimension = 1;
   fUseKeys = kFALSE;
   fUseSparseHist = kFALSE;
   fIsHistStrict = kTRUE;
   fEpsilon = DEFAULT_EPSILON;
   fDelta = DEFAULT_DELTA;
   fIntervalType = kShortest;
   fTFLower = -1.0 * RooNumber::infinity();
   fTFUpper = RooNumber::infinity();
   fVecWeight = 0;
   fLeftSideTF = -1;
}

MCMCInterval::MCMCInterval(const char* name,
        const RooArgSet& parameters, MarkovChain& chain) : ConfInterval(name)
{
   fNumBurnInSteps = 0;
   fConfidenceLevel = 0.0;
   fHistConfLevel = 0.0;
   fKeysConfLevel = 0.0;
   fTFConfLevel = 0.0;
   fFull = 0.0;
   fAxes = NULL;
   fChain = &chain;
   fDataHist = NULL;
   fSparseHist = NULL;
   fVector.clear();
   fKeysPdf = NULL;
   fProduct = NULL;
   fHeaviside = NULL;
   fKeysDataHist = NULL;
   fCutoffVar = NULL;
   fHist = NULL;
   fHistCutoff = -1;
   fKeysCutoff = -1;
   fUseKeys = kFALSE;
   fUseSparseHist = kFALSE;
   fIsHistStrict = kTRUE;
   fEpsilon = DEFAULT_EPSILON;
   SetParameters(parameters);
   fDelta = DEFAULT_DELTA;
   fIntervalType = kShortest;
   fTFLower = -1.0 * RooNumber::infinity();
   fTFUpper = RooNumber::infinity();
   fVecWeight = 0;
   fLeftSideTF = -1;
}

MCMCInterval::~MCMCInterval()
{
   // destructor
   delete[] fAxes;
   delete fHist;
   delete fChain;
   // kbelasco: check here for memory management errors
   delete fDataHist;
   delete fSparseHist;
   delete fKeysPdf;
   delete fProduct;
   delete fHeaviside;
   delete fKeysDataHist;
   delete fCutoffVar;
}

struct CompareDataHistBins {
   CompareDataHistBins(RooDataHist* hist) : fDataHist(hist) {}
   bool operator() (Int_t bin1 , Int_t bin2) { 
      fDataHist->get(bin1);
      Double_t n1 = fDataHist->weight();
      fDataHist->get(bin2);
      Double_t n2 = fDataHist->weight();
      return (n1 < n2);
   }
   RooDataHist* fDataHist; 
};

struct CompareSparseHistBins {
   CompareSparseHistBins(THnSparse* hist) : fSparseHist(hist) {}
   bool operator() (Long_t bin1, Long_t bin2) { 
      Double_t n1 = fSparseHist->GetBinContent(bin1);
      Double_t n2 = fSparseHist->GetBinContent(bin2);
      return (n1 < n2);
   }
   THnSparse* fSparseHist; 
};

struct CompareVectorIndices {
   CompareVectorIndices(MarkovChain* chain, RooRealVar* param) :
      fChain(chain), fParam(param) {}
   bool operator() (Int_t i, Int_t j) {
      Double_t xi = fChain->Get(i)->getRealValue(fParam->GetName());
      Double_t xj = fChain->Get(j)->getRealValue(fParam->GetName());
      return (xi < xj);
   }
   MarkovChain* fChain;
   RooRealVar* fParam;
};

// kbelasco: for this method, consider running DetermineInterval() if 
// fKeysPdf==NULL, fSparseHist==NULL, fDataHist==NULL, or fVector.size()==0
// rather than just returning false.  Though this should not be an issue
// because nobody should be able to get an MCMCInterval that has their interval
// or posterior representation NULL/empty since they should only get this
// through the MCMCCalculator
Bool_t MCMCInterval::IsInInterval(const RooArgSet& point) const 
{
   if (fIntervalType == kShortest) {
      if (fUseKeys) {
         if (fKeysPdf == NULL)
            return false;

         // evaluate keyspdf at point and return whether >= cutoff
         RooStats::SetParameters(&point, const_cast<RooArgSet *>(&fParameters));
         return fKeysPdf->getVal(&fParameters) >= fKeysCutoff;
      } else {
         if (fUseSparseHist) {
            if (fSparseHist == NULL)
               return false;

            // evalute sparse hist at bin where point lies and return
            // whether >= cutoff
            RooStats::SetParameters(&point,
                                    const_cast<RooArgSet*>(&fParameters));
            Long_t bin;
            // kbelasco: consider making x static
            Double_t* x = new Double_t[fDimension];
            for (Int_t i = 0; i < fDimension; i++)
               x[i] = fAxes[i]->getVal();
            bin = fSparseHist->GetBin(x, kFALSE);
            Double_t weight = fSparseHist->GetBinContent((Long64_t)bin);
            delete[] x;
            return (weight >= (Double_t)fHistCutoff);
         } else {
            if (fDataHist == NULL)
               return false;

            // evaluate data hist at bin where point lies and return whether
            // >= cutoff
            Int_t bin;
            bin = fDataHist->getIndex(point);
            fDataHist->get(bin);
            return (fDataHist->weight() >= (Double_t)fHistCutoff);
         }
      }
   } else if (fIntervalType == kTailFraction) {
      if (fVector.size() == 0)
         return false;

      // return whether value of point is within the range
      Double_t x = point.getRealValue(fAxes[0]->GetName());
      if (fTFLower <= x && x <= fTFUpper)
         return true;

      return false;
   }

   coutE(InputArguments) << "Error in MCMCInterval::IsInInterval: "
      << "Interval type not set.  Returning false." << endl;
   return false;
}

void MCMCInterval::SetConfidenceLevel(Double_t cl)
{
   fConfidenceLevel = cl;
   DetermineInterval();
}

// kbelasco: update this or just take it out
// kbelasco: consider keeping this around but changing the implementation
// to set the number of bins for each RooRealVar and then reacreating the
// histograms
//void MCMCInterval::SetNumBins(Int_t numBins)
//{
//   if (numBins > 0) {
//      fPreferredNumBins = numBins;
//      for (Int_t d = 0; d < fDimension; d++)
//         fNumBins[d] = numBins;
//   }
//   else {
//      coutE(Eval) << "* Error in MCMCInterval::SetNumBins: " <<
//                     "Negative number of bins given: " << numBins << endl;
//      return;
//   }
//
//   // If the histogram already exists, recreate it with the new bin numbers
//   if (fHist != NULL)
//      CreateHist();
//}

void MCMCInterval::SetAxes(RooArgList& axes)
{
   Int_t size = axes.getSize();
   if (size != fDimension) {
      coutE(InputArguments) << "* Error in MCMCInterval::SetAxes: " <<
                               "number of variables in axes (" << size <<
                               ") doesn't match number of parameters ("
                               << fDimension << ")" << endl;
      return;
   }
   for (Int_t i = 0; i < size; i++)
      fAxes[i] = (RooRealVar*)axes.at(i);
}

void MCMCInterval::CreateKeysPdf()
{
   // kbelasco: check here for memory leak.  does RooNDKeysPdf use
   // the RooArgList passed to it or does it make a clone?
   // also check for memory leak from chain, does RooNDKeysPdf clone that?
   if (fAxes == NULL || fParameters.getSize() == 0) {
      coutE(InputArguments) << "Error in MCMCInterval::CreateKeysPdf: "
         << "parameters have not been set." << endl;
      return;
   }

   if (fNumBurnInSteps >= fChain->Size()) {
      coutE(InputArguments) <<
         "MCMCInterval::CreateKeysPdf: creation of Keys PDF failed: " <<
         "Number of burn-in steps (num steps to ignore) >= number of steps " <<
         "in Markov chain." << endl;
      delete fKeysPdf;
      delete fCutoffVar;
      delete fHeaviside;
      delete fProduct;
      fKeysPdf = NULL;
      fCutoffVar = NULL;
      fHeaviside = NULL;
      fProduct = NULL;
      return;
   }

   RooDataSet* chain = fChain->GetAsDataSet(SelectVars(fParameters),
         EventRange(fNumBurnInSteps, fChain->Size()));
   RooArgList* paramsList = new RooArgList();
   for (Int_t i = 0; i < fDimension; i++)
      paramsList->add(*fAxes[i]);

   // kbelasco: check for memory leaks with chain.  who owns it? does
   // RooNDKeysPdf take ownership?
   fKeysPdf = new RooNDKeysPdf("keysPDF", "Keys PDF", *paramsList, *chain, "a");
   fCutoffVar = new RooRealVar("cutoff", "cutoff", 0);
   fHeaviside = new Heaviside("heaviside", "Heaviside", *fKeysPdf, *fCutoffVar);
   fProduct = new RooProduct("product", "Keys PDF & Heaviside Product",
                                        RooArgSet(*fKeysPdf, *fHeaviside));
}

void MCMCInterval::CreateHist()
{
   if (fAxes == NULL || fChain == NULL) {
      coutE(Eval) << "* Error in MCMCInterval::CreateHist(): " <<
                     "Crucial data member was NULL." << endl;
      coutE(Eval) << "Make sure to fully construct/initialize." << endl;
      return;
   }
   if (fHist != NULL)
      delete fHist;

   if (fNumBurnInSteps >= fChain->Size()) {
      coutE(InputArguments) <<
         "MCMCInterval::CreateHist: creation of histogram failed: " <<
         "Number of burn-in steps (num steps to ignore) >= number of steps " <<
         "in Markov chain." << endl;
      fHist = NULL;
      return;
   }

   if (fDimension == 1)
      fHist = new TH1F("posterior", "MCMC Posterior Histogram",
            fAxes[0]->numBins(), fAxes[0]->getMin(), fAxes[0]->getMax());

   else if (fDimension == 2)
      fHist = new TH2F("posterior", "MCMC Posterior Histogram",
            fAxes[0]->numBins(), fAxes[0]->getMin(), fAxes[0]->getMax(),
            fAxes[1]->numBins(), fAxes[1]->getMin(), fAxes[1]->getMax());

   else if (fDimension == 3) 
      fHist = new TH3F("posterior", "MCMC Posterior Histogram",
            fAxes[0]->numBins(), fAxes[0]->getMin(), fAxes[0]->getMax(),
            fAxes[1]->numBins(), fAxes[1]->getMin(), fAxes[1]->getMax(),
            fAxes[2]->numBins(), fAxes[2]->getMin(), fAxes[2]->getMax());

   else {
      coutE(Eval) << "* Error in MCMCInterval::CreateHist() : " <<
                     "TH1* couldn't handle dimension: " << fDimension << endl;
      return;
   }

   // Fill histogram
   Int_t size = fChain->Size();
   const RooArgSet* entry;
   for (Int_t i = fNumBurnInSteps; i < size; i++) {
      entry = fChain->Get(i);
      if (fDimension == 1)
         ((TH1F*)fHist)->Fill(entry->getRealValue(fAxes[0]->GetName()),
                              fChain->Weight());
      else if (fDimension == 2)
         ((TH2F*)fHist)->Fill(entry->getRealValue(fAxes[0]->GetName()),
                              entry->getRealValue(fAxes[1]->GetName()),
                              fChain->Weight());
      else
         ((TH3F*)fHist)->Fill(entry->getRealValue(fAxes[0]->GetName()),
                              entry->getRealValue(fAxes[1]->GetName()),
                              entry->getRealValue(fAxes[2]->GetName()),
                              fChain->Weight());
   }

   if (fDimension >= 1)
      fHist->GetXaxis()->SetTitle(fAxes[0]->GetName());
   if (fDimension >= 2)
      fHist->GetYaxis()->SetTitle(fAxes[1]->GetName());
   if (fDimension >= 3)
      fHist->GetZaxis()->SetTitle(fAxes[2]->GetName());
}

void MCMCInterval::CreateSparseHist()
{
   if (fAxes == NULL || fChain == NULL) {
      coutE(InputArguments) << "* Error in MCMCInterval::CreateSparseHist(): "
                            << "Crucial data member was NULL." << endl;
      coutE(InputArguments) << "Make sure to fully construct/initialize."
                            << endl;
      return;
   }
   if (fSparseHist != NULL)
      delete fSparseHist;

   Double_t* min = new Double_t[fDimension];
   Double_t* max = new Double_t[fDimension];
   Int_t* bins = new Int_t[fDimension];
   for (Int_t i = 0; i < fDimension; i++) {
      min[i] = fAxes[i]->getMin();
      max[i] = fAxes[i]->getMax();
      bins[i] = fAxes[i]->numBins();
   }
   fSparseHist = new THnSparseF("posterior", "MCMC Posterior Histogram",
         fDimension, bins, min, max);

   delete[] min;
   delete[] max;
   delete[] bins;

   // kbelasco: it appears we need to call Sumw2() just to get the
   // histogram to keep a running total of the weight so that Getsumw doesn't
   // just return 0
   fSparseHist->Sumw2();

   if (fNumBurnInSteps >= fChain->Size()) {
      coutE(InputArguments) <<
         "MCMCInterval::CreateSparseHist: creation of histogram failed: " <<
         "Number of burn-in steps (num steps to ignore) >= number of steps " <<
         "in Markov chain." << endl;
   }

   // Fill histogram
   Int_t size = fChain->Size();
   const RooArgSet* entry;
   Double_t* x = new Double_t[fDimension];
   for (Int_t i = fNumBurnInSteps; i < size; i++) {
      entry = fChain->Get(i);
      for (Int_t ii = 0; ii < fDimension; ii++)
         x[ii] = entry->getRealValue(fAxes[ii]->GetName());
      fSparseHist->Fill(x, fChain->Weight());
   }
   delete[] x;
}

void MCMCInterval::CreateDataHist()
{
   if (fParameters.getSize() == 0 || fChain == NULL) {
      coutE(Eval) << "* Error in MCMCInterval::CreateDataHist(): " <<
                     "Crucial data member was NULL or empty." << endl;
      coutE(Eval) << "Make sure to fully construct/initialize." << endl;
      return;
   }

   if (fNumBurnInSteps >= fChain->Size()) {
      coutE(InputArguments) <<
         "MCMCInterval::CreateDataHist: creation of histogram failed: " <<
         "Number of burn-in steps (num steps to ignore) >= number of steps " <<
         "in Markov chain." << endl;
      fDataHist = NULL;
      return;
   }

   fDataHist = fChain->GetAsDataHist(SelectVars(fParameters),
         EventRange(fNumBurnInSteps, fChain->Size()));
}

void MCMCInterval::CreateVector(RooRealVar* param)
{
   fVector.clear();
   fVecWeight = 0;

   if (fChain == NULL) {
      coutE(InputArguments) << "* Error in MCMCInterval::CreateVector(): " <<
                     "Crucial data member (Markov chain) was NULL." << endl;
      coutE(InputArguments) << "Make sure to fully construct/initialize."
                            << endl;
      return;
   }

   if (fNumBurnInSteps >= fChain->Size()) {
      coutE(InputArguments) <<
         "MCMCInterval::CreateVector: creation of vector failed: " <<
         "Number of burn-in steps (num steps to ignore) >= number of steps " <<
         "in Markov chain." << endl;
   }

   // Fill vector
   Int_t size = fChain->Size() - fNumBurnInSteps;
   fVector.resize(size);
   Int_t i;
   Int_t chainIndex;
   for (i = 0; i < size; i++) {
      chainIndex = i + fNumBurnInSteps;
      fVector[i] = chainIndex;
      fVecWeight += fChain->Weight(chainIndex);
   }

   stable_sort(fVector.begin(), fVector.end(),
               CompareVectorIndices(fChain, param));
}

void MCMCInterval::SetParameters(const RooArgSet& parameters)
{
   fParameters.removeAll();
   fParameters.add(parameters);
   fDimension = fParameters.getSize();
   if (fAxes != NULL)
      delete[] fAxes;
   fAxes = new RooRealVar*[fDimension];
   TIterator* it = fParameters.createIterator();
   Int_t n = 0;
   TObject* obj;
   while ((obj = it->Next()) != NULL) {
      if (dynamic_cast<RooRealVar*>(obj) != NULL)
         fAxes[n] = (RooRealVar*)obj;
      else
         coutE(Eval) << "* Error in MCMCInterval::SetParameters: " <<
                     obj->GetName() << " not a RooRealVar*" << std::endl;
      n++;
   }
   delete it;
}

void MCMCInterval::DetermineInterval()
{
   switch (fIntervalType) {
      case kShortest:
         DetermineShortestInterval();
         break;
      case kTailFraction:
         DetermineTailFractionInterval();
         break;
      default:
         coutE(InputArguments) << "MCMCInterval::DetermineInterval(): " <<
            "Error: Interval type not set" << endl;
         break;
   }
}

void MCMCInterval::DetermineShortestInterval()
{
   if (fUseKeys)
      DetermineByKeys();
   else
      DetermineByHist();
}

void MCMCInterval::DetermineTailFractionInterval()
{
   if (fLeftSideTF < 0 || fLeftSideTF > 1) {
      coutE(InputArguments) << "MCMCInterval::DetermineTailFractionInterval: "
         << "Fraction must be in the range [0, 1].  "
         << fLeftSideTF << "is not allowed." << endl;
      return;
   }

   if (fDimension != 1) {
      coutE(InputArguments) << "MCMCInterval::DetermineTailFractionInterval(): "
         << "Error: Can only find a tail-fraction interval for 1-D intervals"
         << endl;
      return;
   }

   if (fAxes == NULL) {
      coutE(InputArguments) << "MCMCInterval::DetermineTailFractionInterval(): "
                            << "Crucial data member was NULL." << endl;
      coutE(InputArguments) << "Make sure to fully construct/initialize."
                            << endl;
      return;
   }

   // kbelasco: fill in code here to find interval
   //
   // also make changes so that calling GetPosterior...() returns NULL
   // when fIntervalType == kTailFraction, since there really
   // is no posterior for this type of interval determination
   if (fVector.size() == 0)
      CreateVector(fAxes[0]);

   if (fVector.size() == 0 || fVecWeight == 0) {
      // if size is still 0, then creation failed.
      // if fVecWeight == 0, then there are no entries (indicates the same
      // error as fVector.size() == 0 because that only happens when
      // fNumBurnInSteps >= fChain->Size())
      // either way, reset and return
      fVector.clear();
      fTFLower = -1.0 * RooNumber::infinity();
      fTFUpper = RooNumber::infinity();
      fTFConfLevel = 0.0;
      fVecWeight = 0;
      return;
   }

   RooRealVar* param = fAxes[0];

   Double_t c = fConfidenceLevel;
   Double_t leftTailCutoff  = fVecWeight * (1 - c) * fLeftSideTF;
   Double_t rightTailCutoff = fVecWeight * (1 - c) * (1 - fLeftSideTF);
   Double_t leftTailSum  = 0;
   Double_t rightTailSum = 0;

   // kbelasco: consider changing these values to +infinty and -infinity
   Double_t ll = param->getMin();
   Double_t ul = param->getMax();

   Double_t x;
   Double_t w;

   // save a lot of GetName() calls if compiler does not already optimize this
   const char* name = param->GetName();

   // find lower limit
   Int_t i;
   for (i = 0; i < (Int_t)fVector.size(); i++) {
      x = fChain->Get(fVector[i])->getRealValue(name);
      w = fChain->Weight();

      if (TMath::Abs(leftTailSum + w - leftTailCutoff) <
          TMath::Abs(leftTailSum - leftTailCutoff)) {
         // moving the lower limit to x would bring us closer to the desired
         // left tail size
         ll = x;
         leftTailSum += w;
      } else
         break;
   }

   // find upper limit
   for (i = (Int_t)fVector.size() - 1; i >= 0; i--) {
      x = fChain->Get(fVector[i])->getRealValue(name);
      w = fChain->Weight();

      if (TMath::Abs(rightTailSum + w - rightTailCutoff) <
          TMath::Abs(rightTailSum - rightTailCutoff)) {
         // moving the lower limit to x would bring us closer to the desired
         // left tail size
         ul = x;
         rightTailSum += w;
      } else
         break;
   }

   fTFLower = ll;
   fTFUpper = ul;
   fTFConfLevel = 1 - (leftTailSum + rightTailSum) / fVecWeight;
}

void MCMCInterval::DetermineByKeys()
{
   if (fKeysPdf == NULL)
      CreateKeysPdf();

   if (fKeysPdf == NULL) {
      // if fKeysPdf is still NULL, then it means CreateKeysPdf failed
      // so clear all the data members this function would normally determine
      // and return
      fFull = 0.0;
      fKeysCutoff = -1;
      fKeysConfLevel = 0.0;
      return;
   }

   // now we have a keys pdf of the posterior

   Double_t cutoff = 0.0;
   fCutoffVar->setVal(cutoff);
   RooAbsReal* integral = fProduct->createIntegral(fParameters, NormSet(fParameters));
   Double_t full = integral->getVal(fParameters);
   fFull = full;
   delete integral;
   if (full < 0.98) {
      coutW(Eval) << "Warning: Integral of Keys PDF came out to " << full
         << " intead of expected value 1.  Will continue using this "
         << "factor to normalize further integrals of this PDF." << endl;
   }

   // kbelasco: Is there a better way to set the search range?
   // from 0 to max value of Keys
   // kbelasco: how to get max value?
   //Double_t max = product.maxVal(product.getMaxVal(fParameters));

   Double_t volume = 1.0;
   TIterator* it = fParameters.createIterator();
   RooRealVar* var;
   while ((var = (RooRealVar*)it->Next()) != NULL)
      volume *= (var->getMax() - var->getMin());
   delete it;

   Double_t topCutoff = full / volume;
   Double_t bottomCutoff = topCutoff;
   Double_t confLevel = CalcConfLevel(topCutoff, full);
   if (AcceptableConfLevel(confLevel)) {
      fKeysConfLevel = confLevel;
      fKeysCutoff = topCutoff;
      return;
   }
   Bool_t changed = kFALSE;
   // find high end of range
   while (confLevel > fConfidenceLevel) {
      topCutoff *= 2.0;
      confLevel = CalcConfLevel(topCutoff, full);
      if (AcceptableConfLevel(confLevel)) {
         fKeysConfLevel = confLevel;
         fKeysCutoff = topCutoff;
         return;
      }
      changed = kTRUE;
   }
   if (changed) {
      bottomCutoff = topCutoff / 2.0;
   } else {
      changed = kFALSE;
      bottomCutoff /= 2.0;
      confLevel = CalcConfLevel(bottomCutoff, full);
      if (AcceptableConfLevel(confLevel)) {
         fKeysConfLevel = confLevel;
         fKeysCutoff = bottomCutoff;
         return;
      }
      while (confLevel < fConfidenceLevel) {
         bottomCutoff /= 2.0;
         confLevel = CalcConfLevel(bottomCutoff, full);
         if (AcceptableConfLevel(confLevel)) {
            fKeysConfLevel = confLevel;
            fKeysCutoff = bottomCutoff;
            return;
         }
         changed = kTRUE;
      }
      if (changed) {
         topCutoff = bottomCutoff * 2.0;
      }
   }

   coutI(Eval) << "range set: [" << bottomCutoff << ", " << topCutoff << "]"
               << endl;

   cutoff = (topCutoff + bottomCutoff) / 2.0;
   confLevel = CalcConfLevel(cutoff, full);

   // need to use WithinDeltaFraction() because sometimes the integrating the
   // posterior in this binary search seems to not have enough granularity to
   // find an acceptable conf level (small no. of strange cases).
   // WithinDeltaFraction causes the search to terminate when
   // topCutoff is essentially equal to bottomCutoff (compared to the magnitude
   // of their mean).
   while (!AcceptableConfLevel(confLevel) &&
          !WithinDeltaFraction(topCutoff, bottomCutoff)) {
      if (confLevel > fConfidenceLevel)
         bottomCutoff = cutoff;
      else if (confLevel < fConfidenceLevel)
         topCutoff = cutoff;
      cutoff = (topCutoff + bottomCutoff) / 2.0;
      coutI(Eval) << "cutoff range: [" << bottomCutoff << ", "
                  << topCutoff << "]" << endl;
      confLevel = CalcConfLevel(cutoff, full);
   }

   fKeysCutoff = cutoff;
   fKeysConfLevel = confLevel;
}

void MCMCInterval::DetermineByHist()
{
   if (fUseSparseHist)
      DetermineBySparseHist();
   else
      DetermineByDataHist();
}

void MCMCInterval::DetermineBySparseHist()
{
   Long_t numBins;
   if (fSparseHist == NULL)
      CreateSparseHist();

   if (fSparseHist == NULL) {
      // if fSparseHist is still NULL, then CreateSparseHist failed
      fHistCutoff = -1;
      fHistConfLevel = 0.0;
      return;
   }

   numBins = (Long_t)fSparseHist->GetNbins();

   std::vector<Long_t> bins(numBins);
   for (Int_t ibin = 0; ibin < numBins; ibin++)
      bins[ibin] = (Long_t)ibin;
   std::stable_sort(bins.begin(), bins.end(), CompareSparseHistBins(fSparseHist));

   Double_t nEntries = fSparseHist->GetSumw();
   Double_t sum = 0;
   Double_t content;
   Int_t i;
   // see above note on indexing to understand numBins - 3
   for (i = numBins - 1; i >= 0; i--) {
      content = fSparseHist->GetBinContent(bins[i]);
      if ((sum + content) / nEntries >= fConfidenceLevel) {
         fHistCutoff = content;
         if (fIsHistStrict) {
            sum += content;
            i--;
            break;
         } else {
            i++;
            break;
         }
      }
      sum += content;
   }

   if (fIsHistStrict) {
      // keep going to find the sum
      for ( ; i >= 0; i--) {
         content = fSparseHist->GetBinContent(bins[i]);
         if (content == fHistCutoff)
            sum += content;
         else
            break; // content must be < fHistCutoff
      }
   } else {
      // backtrack to find the cutoff and sum
      for ( ; i < numBins; i++) {
         content = fSparseHist->GetBinContent(bins[i]);
         if (content > fHistCutoff) {
            fHistCutoff = content;
            break;
         } else // content == fHistCutoff
            sum -= content;
         if (i == numBins - 1)
            // still haven't set fHistCutoff correctly yet, and we have no bins
            // left, so set fHistCutoff to something higher than the tallest bin
            fHistCutoff = content + 1.0;
      }
   }

   fHistConfLevel = sum / nEntries;
}

void MCMCInterval::DetermineByDataHist()
{
   Int_t numBins;
   if (fDataHist == NULL)
      CreateDataHist();
   if (fDataHist == NULL) {
      // if fDataHist is still NULL, then CreateDataHist failed
      fHistCutoff = -1;
      fHistConfLevel = 0.0;
      return;
   }

   numBins = fDataHist->numEntries();

   std::vector<Int_t> bins(numBins);
   for (Int_t ibin = 0; ibin < numBins; ibin++)
      bins[ibin] = ibin;
   std::stable_sort(bins.begin(), bins.end(), CompareDataHistBins(fDataHist));

   Double_t nEntries = fDataHist->sum(kFALSE);
   Double_t sum = 0;
   Double_t content;
   Int_t i;
   for (i = numBins - 1; i >= 0; i--) {
      fDataHist->get(bins[i]);
      content = fDataHist->weight();
      if ((sum + content) / nEntries >= fConfidenceLevel) {
         fHistCutoff = content;
         if (fIsHistStrict) {
            sum += content;
            i--;
            break;
         } else {
            i++;
            break;
         }
      }
      sum += content;
   }

   if (fIsHistStrict) {
      // keep going to find the sum
      for ( ; i >= 0; i--) {
         fDataHist->get(bins[i]);
         content = fDataHist->weight();
         if (content == fHistCutoff)
            sum += content;
         else
            break; // content must be < fHistCutoff
      }
   } else {
      // backtrack to find the cutoff and sum
      for ( ; i < numBins; i++) {
         fDataHist->get(bins[i]);
         content = fDataHist->weight();
         if (content > fHistCutoff) {
            fHistCutoff = content;
            break;
         } else // content == fHistCutoff
            sum -= content;
         if (i == numBins - 1)
            // still haven't set fHistCutoff correctly yet, and we have no bins
            // left, so set fHistCutoff to something higher than the tallest bin
            fHistCutoff = content + 1.0;
      }
   }

   fHistConfLevel = sum / nEntries;
}

Double_t MCMCInterval::GetActualConfidenceLevel()
{
   if (fIntervalType == kShortest) {
      if (fUseKeys)
         return fKeysConfLevel;
      else
         return fHistConfLevel;
   } else if (fIntervalType == kTailFraction) {
      return fTFConfLevel;
   } else {
      coutE(InputArguments) << "MCMCInterval::GetActualConfidenceLevel: "
         << "not implemented for this type of interval.  Returning 0." << endl;
      return 0;
   }
}

Double_t MCMCInterval::LowerLimit(RooRealVar& param)
{
   switch (fIntervalType) {
      case kShortest:
         return LowerLimitShortest(param);
      case kTailFraction:
         return LowerLimitTailFraction(param);
      default:
         coutE(InputArguments) << "MCMCInterval::LowerLimit(): " <<
            "Error: Interval type not set" << endl;
         return RooNumber::infinity();
   }
}

Double_t MCMCInterval::UpperLimit(RooRealVar& param)
{
   switch (fIntervalType) {
      case kShortest:
         return UpperLimitShortest(param);
      case kTailFraction:
         return UpperLimitTailFraction(param);
      default:
         coutE(InputArguments) << "MCMCInterval::UpperLimit(): " <<
            "Error: Interval type not set" << endl;
         return RooNumber::infinity();
   }
}

Double_t MCMCInterval::LowerLimitTailFraction(RooRealVar& /*param*/)
{
   if (fTFLower == -1.0 * RooNumber::infinity())
      DetermineTailFractionInterval();

   return fTFLower;
}

Double_t MCMCInterval::UpperLimitTailFraction(RooRealVar& /*param*/)
{
   if (fTFUpper == RooNumber::infinity())
      DetermineTailFractionInterval();

   return fTFUpper;
}

Double_t MCMCInterval::LowerLimitShortest(RooRealVar& param)
{
   if (fUseKeys)
      return LowerLimitByKeys(param);
   else
      return LowerLimitByHist(param);
}

Double_t MCMCInterval::UpperLimitShortest(RooRealVar& param)
{
   if (fUseKeys)
      return UpperLimitByKeys(param);
   else
      return UpperLimitByHist(param);
}

// Determine the lower limit for param on this interval
// using the binned data set
Double_t MCMCInterval::LowerLimitByHist(RooRealVar& param)
{
   if (fUseSparseHist)
      return LowerLimitBySparseHist(param);
   else
      return LowerLimitByDataHist(param);
}

// Determine the upper limit for param on this interval
// using the binned data set
Double_t MCMCInterval::UpperLimitByHist(RooRealVar& param)
{
   if (fUseSparseHist)
      return UpperLimitBySparseHist(param);
   else
      return UpperLimitByDataHist(param);
}

// Determine the lower limit for param on this interval
// using the binned data set
Double_t MCMCInterval::LowerLimitBySparseHist(RooRealVar& param)
{
   if (fDimension != 1) {
      coutE(InputArguments) << "In MCMCInterval::LowerLimitBySparseHist: "
         << "Sorry, will not compute lower limit unless dimension == 1" << endl;
      return param.getMin();
   }
   if (fHistCutoff < 0)
      DetermineBySparseHist(); // this initializes fSparseHist

   if (fHistCutoff < 0) {
      // if fHistCutoff < 0 still, then determination of interval failed
      coutE(Eval) << "In MCMCInterval::LowerLimitBySparseHist: "
         << "couldn't determine cutoff.  Check that num burn in steps < num "
         << "steps in the Markov chain.  Returning param.getMin()." << endl;
      return param.getMin();
   }

   std::vector<Int_t> coord(fDimension);
   for (Int_t d = 0; d < fDimension; d++) {
      if (strcmp(fAxes[d]->GetName(), param.GetName()) == 0) {
         Long_t numBins = (Long_t)fSparseHist->GetNbins();
         Double_t lowerLimit = param.getMax();
         Double_t val;
         for (Long_t i = 0; i < numBins; i++) {
            if (fSparseHist->GetBinContent(i, &coord[0]) >= fHistCutoff) {
               val = fSparseHist->GetAxis(d)->GetBinCenter(coord[d]);
               if (val < lowerLimit)
                  lowerLimit = val;
            }
         }
         return lowerLimit;
      }
   }
   return param.getMin();
}

// Determine the lower limit for param on this interval
// using the binned data set
Double_t MCMCInterval::LowerLimitByDataHist(RooRealVar& param)
{
   if (fHistCutoff < 0)
      DetermineByDataHist(); // this initializes fDataHist

   if (fHistCutoff < 0) {
      // if fHistCutoff < 0 still, then determination of interval failed
      coutE(Eval) << "In MCMCInterval::LowerLimitByDataHist: "
         << "couldn't determine cutoff.  Check that num burn in steps < num "
         << "steps in the Markov chain.  Returning param.getMin()." << endl;
      return param.getMin();
   }

   for (Int_t d = 0; d < fDimension; d++) {
      if (strcmp(fAxes[d]->GetName(), param.GetName()) == 0) {
         Int_t numBins = fDataHist->numEntries();
         Double_t lowerLimit = param.getMax();
         Double_t val;
         for (Int_t i = 0; i < numBins; i++) {
            fDataHist->get(i);
            if (fDataHist->weight() >= fHistCutoff) {
               val = fDataHist->get()->getRealValue(param.GetName());
               if (val < lowerLimit)
                  lowerLimit = val;
            }
         }
         return lowerLimit;
      }
   }
   return param.getMin();
}

// Determine the upper limit for param on this interval
// using the binned data set
Double_t MCMCInterval::UpperLimitBySparseHist(RooRealVar& param)
{
   if (fDimension != 1) {
      coutE(InputArguments) << "In MCMCInterval::UpperLimitBySparseHist: "
         << "Sorry, will not compute upper limit unless dimension == 1" << endl;
      return param.getMax();
   }
   if (fHistCutoff < 0)
      DetermineBySparseHist(); // this initializes fSparseHist

   if (fHistCutoff < 0) {
      // if fHistCutoff < 0 still, then determination of interval failed
      coutE(Eval) << "In MCMCInterval::UpperLimitBySparseHist: "
         << "couldn't determine cutoff.  Check that num burn in steps < num "
         << "steps in the Markov chain.  Returning param.getMax()." << endl;
      return param.getMax();
   }

   std::vector<Int_t> coord(fDimension);
   for (Int_t d = 0; d < fDimension; d++) {
      if (strcmp(fAxes[d]->GetName(), param.GetName()) == 0) {
         Long_t numBins = (Long_t)fSparseHist->GetNbins();
         Double_t upperLimit = param.getMin();
         Double_t val;
         for (Long_t i = 0; i < numBins; i++) {
            if (fSparseHist->GetBinContent(i, &coord[0]) >= fHistCutoff) {
               val = fSparseHist->GetAxis(d)->GetBinCenter(coord[d]);
               if (val > upperLimit)
                  upperLimit = val;
            }
         }
         return upperLimit;
      }
   }
   return param.getMax();
}

// Determine the upper limit for param on this interval
// using the binned data set
Double_t MCMCInterval::UpperLimitByDataHist(RooRealVar& param)
{
   if (fHistCutoff < 0)
      DetermineByDataHist(); // this initializes fDataHist

   if (fHistCutoff < 0) {
      // if fHistCutoff < 0 still, then determination of interval failed
      coutE(Eval) << "In MCMCInterval::UpperLimitByDataHist: "
         << "couldn't determine cutoff.  Check that num burn in steps < num "
         << "steps in the Markov chain.  Returning param.getMax()." << endl;
      return param.getMax();
   }

   for (Int_t d = 0; d < fDimension; d++) {
      if (strcmp(fAxes[d]->GetName(), param.GetName()) == 0) {
         Int_t numBins = fDataHist->numEntries();
         Double_t upperLimit = param.getMin();
         Double_t val;
         for (Int_t i = 0; i < numBins; i++) {
            fDataHist->get(i);
            if (fDataHist->weight() >= fHistCutoff) {
               val = fDataHist->get()->getRealValue(param.GetName());
               if (val > upperLimit)
                  upperLimit = val;
            }
         }
         return upperLimit;
      }
   }
   return param.getMax();
}

// Determine the lower limit for param on this interval
// using the keys pdf
Double_t MCMCInterval::LowerLimitByKeys(RooRealVar& param)
{
   if (fKeysCutoff < 0)
      DetermineByKeys();

   if (fKeysDataHist == NULL)
      CreateKeysDataHist();

   if (fKeysCutoff < 0 || fKeysDataHist == NULL) {
      // failure in determination of cutoff and/or creation of histogram
      coutE(Eval) << "in MCMCInterval::LowerLimitByKeys(): "
         << "couldn't find lower limit, check that the number of burn in "
         << "steps < number of total steps in the Markov chain.  Returning "
         << "param.getMin()" << endl;
      return param.getMin();
   }

   for (Int_t d = 0; d < fDimension; d++) {
      if (strcmp(fAxes[d]->GetName(), param.GetName()) == 0) {
         Int_t numBins = fKeysDataHist->numEntries();
         Double_t lowerLimit = param.getMax();
         Double_t val;
         for (Int_t i = 0; i < numBins; i++) {
            fKeysDataHist->get(i);
            if (fKeysDataHist->weight() >= fKeysCutoff) {
               val = fKeysDataHist->get()->getRealValue(param.GetName());
               if (val < lowerLimit)
                  lowerLimit = val;
            }
         }
         return lowerLimit;
      }
   }
   return param.getMin();
}

// Determine the upper limit for param on this interval
// using the keys pdf
Double_t MCMCInterval::UpperLimitByKeys(RooRealVar& param)
{
   if (fKeysCutoff < 0)
      DetermineByKeys();

   if (fKeysDataHist == NULL)
      CreateKeysDataHist();

   if (fKeysCutoff < 0 || fKeysDataHist == NULL) {
      // failure in determination of cutoff and/or creation of histogram
      coutE(Eval) << "in MCMCInterval::UpperLimitByKeys(): "
         << "couldn't find upper limit, check that the number of burn in "
         << "steps < number of total steps in the Markov chain.  Returning "
         << "param.getMax()" << endl;
      return param.getMax();
   }

   for (Int_t d = 0; d < fDimension; d++) {
      if (strcmp(fAxes[d]->GetName(), param.GetName()) == 0) {
         Int_t numBins = fKeysDataHist->numEntries();
         Double_t upperLimit = param.getMin();
         Double_t val;
         for (Int_t i = 0; i < numBins; i++) {
            fKeysDataHist->get(i);
            if (fKeysDataHist->weight() >= fKeysCutoff) {
               val = fKeysDataHist->get()->getRealValue(param.GetName());
               if (val > upperLimit)
                  upperLimit = val;
            }
         }
         return upperLimit;
      }
   }
   return param.getMax();
}

// Determine the approximate maximum value of the Keys PDF
Double_t MCMCInterval::GetKeysMax()
{
   if (fKeysCutoff < 0)
      DetermineByKeys();

   if (fKeysDataHist == NULL)
      CreateKeysDataHist();

   if (fKeysDataHist == NULL) {
      // failure in determination of cutoff and/or creation of histogram
      coutE(Eval) << "in MCMCInterval::KeysMax(): "
         << "couldn't find Keys max value, check that the number of burn in "
         << "steps < number of total steps in the Markov chain.  Returning 0"
         << endl;
      return 0;
   }

   Int_t numBins = fKeysDataHist->numEntries();
   Double_t max = 0;
   Double_t w;
   for (Int_t i = 0; i < numBins; i++) {
      fKeysDataHist->get(i);
      w = fKeysDataHist->weight();
      if (w > max)
         max = w;
   }

   return max;
}

Double_t MCMCInterval::GetHistCutoff()
{
   if (fHistCutoff < 0)
      DetermineByHist();

   return fHistCutoff;
}

Double_t MCMCInterval::GetKeysPdfCutoff()
{
   if (fKeysCutoff < 0)
      DetermineByKeys();

   // kbelasco: if fFull hasn't been set (because Keys creation failed because
   // fNumBurnInSteps >= fChain->Size()) then this will return infinity, which
   // seems ok to me since it will indicate error

   return fKeysCutoff / fFull;
}

Double_t MCMCInterval::CalcConfLevel(Double_t cutoff, Double_t full)
{
   RooAbsReal* integral;
   Double_t confLevel;
   fCutoffVar->setVal(cutoff);
   integral = fProduct->createIntegral(fParameters, NormSet(fParameters));
   confLevel = integral->getVal(fParameters) / full;
   coutI(Eval) << "cutoff = " << cutoff << ", conf = " << confLevel << endl;
   //cout << "tmp: cutoff = " << cutoff << ", conf = " << confLevel << endl;
   delete integral;
   return confLevel;
}

TH1* MCMCInterval::GetPosteriorHist()
{
  if(fConfidenceLevel == 0)
     coutE(InputArguments) << "Error in MCMCInterval::GetPosteriorHist: "
                           << "confidence level not set " << endl;
  if (fHist == NULL)
     CreateHist();

  if (fHist == NULL)
     // if fHist is still NULL, then CreateHist failed
     return NULL;

  return (TH1*) fHist->Clone("MCMCposterior_hist");
}

RooNDKeysPdf* MCMCInterval::GetPosteriorKeysPdf()
{
   if (fConfidenceLevel == 0)
      coutE(InputArguments) << "Error in MCMCInterval::GetPosteriorKeysPdf: "
                            << "confidence level not set " << endl;
   if (fKeysPdf == NULL)
      CreateKeysPdf();

   if (fKeysPdf == NULL)
      // if fKeysPdf is still NULL, then it means CreateKeysPdf failed
      return NULL;

   return (RooNDKeysPdf*) fKeysPdf->Clone("MCMCPosterior_keys");
}

RooProduct* MCMCInterval::GetPosteriorKeysProduct()
{
   if (fConfidenceLevel == 0)
      coutE(InputArguments) << "MCMCInterval::GetPosteriorKeysProduct: "
                            << "confidence level not set " << endl;
   if (fProduct == NULL) {
      CreateKeysPdf();
      DetermineByKeys();
   }

   if (fProduct == NULL)
      // if fProduct is still NULL, then it means CreateKeysPdf failed
      return NULL;

   return (RooProduct*) fProduct->Clone("MCMCPosterior_keysproduct");
}

RooArgSet* MCMCInterval::GetParameters() const
{  
   // returns list of parameters
   return new RooArgSet(fParameters);
}

Bool_t MCMCInterval::AcceptableConfLevel(Double_t confLevel)
{
   return (TMath::Abs(confLevel - fConfidenceLevel) < fEpsilon);
}

Bool_t MCMCInterval::WithinDeltaFraction(Double_t a, Double_t b)
{
   return (TMath::Abs(a - b) < TMath::Abs(fDelta * (a + b)/2));
}

void MCMCInterval::CreateKeysDataHist()
{
   if (fAxes == NULL)
      return;
   if (fProduct == NULL)
      DetermineByKeys();
   if (fProduct == NULL)
      // if fProduct still NULL, then creation failed
      return;

   //RooAbsBinning** savedBinning = new RooAbsBinning*[fDimension];
   Int_t* savedBins = new Int_t[fDimension];
   Int_t i;
   Double_t numBins;
   RooRealVar* var;

   // kbelasco: Note - the accuracy is only increased here if the binning for
   // each RooRealVar is uniform

   // kbelasco: look into why saving the binnings and replacing them doesn't
   // work (replaces with 1 bin always).
   // Note: this code modifies the binning for the parameters (if they are
   // uniform) and sets them back to what they were.  If the binnings are not
   // uniform, this code does nothing.

   // first scan through fAxes to make sure all binnings are uniform, or else
   // we can't change the number of bins because there seems to be an error
   // when setting the binning itself rather than just the number of bins
   Bool_t tempChangeBinning = true;
   for (i = 0; i < fDimension; i++) {
      if (!fAxes[i]->getBinning(NULL, false, false).isUniform()) {
         tempChangeBinning = false;
         break;
      }
   }

   // kbelasco: for 1 dimension this should be fine, but for more dimensions
   // the total nubmer of bins in the histogram increases exponentially with
   // the dimension, so don't do this above 1-D for now.
   if (fDimension >= 2)
      tempChangeBinning = false;

   if (tempChangeBinning) {
      // set high nubmer of bins for high accuracy on lower/upper limit by keys
      for (i = 0; i < fDimension; i++) {
         var = fAxes[i];
         //savedBinning[i] = &var->getBinning("__binning_clone", false, true);
         savedBins[i] = var->getBinning(NULL, false, false).numBins();
         numBins = (var->getMax() - var->getMin()) / fEpsilon;
         var->setBins((Int_t)numBins);
      }
   }

   fKeysDataHist = new RooDataHist("_productDataHist",
         "Keys PDF & Heaviside Product Data Hist", fParameters);
   fKeysDataHist = fProduct->fillDataHist(fKeysDataHist, &fParameters, 1.);

   if (tempChangeBinning) {
      // set the binning back to normal
      for (i = 0; i < fDimension; i++)
         //fAxes[i]->setBinning(*savedBinning[i], NULL);
         //fAxes[i]->setBins(savedBinning[i]->numBins(), NULL);
         fAxes[i]->setBins(savedBins[i], NULL);
   }

   //delete[] savedBinning;
   delete[] savedBins;
}

Bool_t MCMCInterval::CheckParameters(const RooArgSet& parameterPoint) const
{  
   // check that the parameters are correct

   if (parameterPoint.getSize() != fParameters.getSize() ) {
     coutE(Eval) << "MCMCInterval: size is wrong, parameters don't match" << std::endl;
     return kFALSE;
   }
   if ( ! parameterPoint.equals( fParameters ) ) {
     coutE(Eval) << "MCMCInterval: size is ok, but parameters don't match" << std::endl;
     return kFALSE;
   }
   return kTRUE;
}
