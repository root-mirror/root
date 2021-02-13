// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_FeldmanCousins
#define ROOSTATS_FeldmanCousins


#include "Rtypes.h"

#include "RooStats/IntervalCalculator.h"

#include "RooStats/ToyMCSampler.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/PointSetInterval.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"

class RooAbsData;

namespace RooStats {

   class ConfInterval;

   class FeldmanCousins : public IntervalCalculator {

   public:

     //     FeldmanCousins();

     /// Common constructor
     FeldmanCousins(RooAbsData& data, ModelConfig& model);

     virtual ~FeldmanCousins();

      /// Main interface to get a ConfInterval (will be a PointSetInterval)
      virtual PointSetInterval* GetInterval() const;

      /// Get the size of the test (eg. rate of Type I error)
      virtual Double_t Size() const {return fSize;}
      /// Get the Confidence level for the test
      virtual Double_t ConfidenceLevel()  const {return 1.-fSize;}
      /// Set the DataSet
      virtual void SetData(RooAbsData& /*data*/) {
   std::cout << "DEPRECATED, set data in constructor" << std::endl;
      }
      /// Set the Pdf
      virtual void SetPdf(RooAbsPdf& /*pdf*/) {
   std::cout << "DEPRECATED, use ModelConfig" << std::endl;
      }

      /// specify the parameters of interest in the interval
      virtual void SetParameters(const RooArgSet& /*set*/) {
   std::cout << "DEPRECATED, use ModelConfig" << std::endl;
      }

      /// specify the nuisance parameters (eg. the rest of the parameters)
      virtual void SetNuisanceParameters(const RooArgSet& /*set*/) {
   std::cout << "DEPRECATED, use ModelConfig" << std::endl;
      }

      /// User-defined set of points to test
      void SetParameterPointsToTest(RooAbsData& pointsToTest) {
   fPointsToTest = &pointsToTest;
      }

      /// User-defined set of points to test
      void SetPOIPointsToTest(RooAbsData& poiToTest) {
   fPOIToTest = &poiToTest;
      }

      /// set the size of the test (rate of Type I error) ( Eg. 0.05 for a 95% Confidence Interval)
      virtual void SetTestSize(Double_t size) {fSize = size;}
      /// set the confidence level for the interval (eg. 0.95 for a 95% Confidence Interval)
      virtual void SetConfidenceLevel(Double_t cl) {fSize = 1.-cl;}

      virtual void SetModel(const ModelConfig &);

      RooAbsData* GetPointsToScan() {
   if(!fPointsToTest) CreateParameterPoints();
   return fPointsToTest;
      }

      /// Get the confidence belt. This requires that CreateConfBelt() has been set.
      ConfidenceBelt* GetConfidenceBelt() {return fConfBelt;}

      void UseAdaptiveSampling(bool flag=true){fAdaptiveSampling=flag;}

      void AdditionalNToysFactor(double fact){fAdditionalNToysFactor = fact;}

      void SetNBins(Int_t bins) {fNbins = bins;}

      void FluctuateNumDataEntries(bool flag=true){fFluctuateData = flag;}

      void SaveBeltToFile(bool flag=true){
   fSaveBeltToFile = flag;
   if(flag) fCreateBelt = true;
      }
      void CreateConfBelt(bool flag=true){fCreateBelt = flag;}

      /// Returns instance of TestStatSampler. Use to change properties of
      /// TestStatSampler, e.g. GetTestStatSampler.SetTestSize(Double_t size);
      TestStatSampler* GetTestStatSampler() const;


   private:

      /// initializes fPointsToTest data member (mutable)
      void CreateParameterPoints() const;

      /// initializes fTestStatSampler data member (mutable)
      void CreateTestStatSampler() const;

      Double_t fSize; // size of the test (eg. specified rate of Type I error)
      ModelConfig &fModel;
      RooAbsData & fData; // data set

      /*
      RooAbsPdf * fPdf; // common PDF
      RooArgSet fPOI; // RooArgSet specifying  parameters of interest for interval
      RooArgSet fNuisParams;// RooArgSet specifying  nuisance parameters for interval
      RooArgSet fObservables;// RooArgSet specifying  nuisance parameters for interval
      */

      mutable ToyMCSampler* fTestStatSampler; // the test statistic sampler
      mutable RooAbsData* fPointsToTest; // points to perform the construction
      mutable RooAbsData* fPOIToTest; // value of POI points to perform the construction
      mutable ConfidenceBelt* fConfBelt;
      Bool_t fAdaptiveSampling; // controls use of adaptive sampling algorithm
      Double_t fAdditionalNToysFactor; // give user ability to ask for more toys
      Int_t fNbins; // number of samples per variable
      Bool_t fFluctuateData;  // tell ToyMCSampler to fluctuate number of entries in dataset
      Bool_t fDoProfileConstruction; // instead of full construction over nuisance parameters, do profile
      Bool_t fSaveBeltToFile; // controls use if ConfidenceBelt should be saved to a TFile
      Bool_t fCreateBelt; // controls use if ConfidenceBelt should be saved to a TFile

   protected:
      ClassDef(FeldmanCousins,2)   // Interface for tools setting limits (producing confidence intervals)
   };
}


#endif
