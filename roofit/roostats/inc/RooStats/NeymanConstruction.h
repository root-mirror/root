// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_NeymanConstruction
#define ROOSTATS_NeymanConstruction


#include "Rtypes.h"

#include "RooStats/IntervalCalculator.h"

#include "RooStats/TestStatSampler.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/PointSetInterval.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"

#include <map>

class RooAbsData;

namespace RooStats {

   class ConfInterval;

   class NeymanConstruction : public IntervalCalculator{

   public:

     ///     NeymanConstruction();
     NeymanConstruction(RooAbsData& data, ModelConfig& model);

     virtual ~NeymanConstruction();

      /// Main interface to get a ConfInterval (will be a PointSetInterval)
     virtual PointSetInterval* GetInterval() const;

      /// in addition to interface we also need:
      /// Set the TestStatSampler (eg. ToyMC or FFT, includes choice of TestStatistic)
      void SetTestStatSampler(TestStatSampler& sampler) {fTestStatSampler = &sampler;}
      /// fLeftSideTailFraction*fSize defines lower edge of acceptance region.
      /// Unified limits use 0, central limits use 0.5,
      /// for upper/lower limits it is 0/1 depends on sign of test statistic w.r.t. parameter
      void SetLeftSideTailFraction(Double_t leftSideFraction = 0.) {fLeftSideFraction = leftSideFraction;}

      /// User-defined set of points to test
      void SetParameterPointsToTest(RooAbsData& pointsToTest) {
   fPointsToTest = &pointsToTest;
        fConfBelt = new ConfidenceBelt("ConfBelt",pointsToTest);
      }
      /// This class can make regularly spaced scans based on range stored in RooRealVars.
      /// Choose number of steps for a rastor scan (common for each dimension)
      ///      void SetNumSteps(Int_t);
      /// This class can make regularly spaced scans based on range stored in RooRealVars.
      /// Choose number of steps for a rastor scan (specific for each dimension)
      ///      void SetNumSteps(std::map<RooAbsArg, Int_t>)

      /// Get the size of the test (eg. rate of Type I error)
      virtual Double_t Size() const {return fSize;}

      /// Get the Confidence level for the test
      virtual Double_t ConfidenceLevel()  const {return 1.-fSize;}

      /// Set ModelConfig
      virtual void SetModel(const ModelConfig &model) {fModel = model;}

      /// Set the DataSet
      virtual void SetData(RooAbsData& data) { fData = data; }

      /// Set the Pdf, add to the the workspace if not already there
      virtual void SetPdf(RooAbsPdf& /*pdf*/) {
        std::cout << "DEPRECATED, use ModelConfig"<<std::endl;
      }

      /// specify the parameters of interest in the interval
      virtual void SetParameters(const RooArgSet& /*set*/) {
        std::cout << "DEPRECATED, use ModelConfig"<<std::endl;
      }

      /// specify the nuisance parameters (eg. the rest of the parameters)
      virtual void SetNuisanceParameters(const RooArgSet& /*set*/) {
        std::cout << "DEPRECATED, use ModelConfig"<<std::endl;
      }

      /// set the size of the test (rate of Type I error) ( Eg. 0.05 for a 95% Confidence Interval)
      virtual void SetTestSize(Double_t size) {fSize = size;}
      /// set the confidence level for the interval (eg. 0.95 for a 95% Confidence Interval)
      virtual void SetConfidenceLevel(Double_t cl) {fSize = 1.-cl;}

      /// Get confidence belt. This requires that CreateConfBelt() has been called.
      ConfidenceBelt* GetConfidenceBelt() {return fCreateBelt ? fConfBelt : nullptr;}

      /// adaptive sampling algorithm to speed up interval calculation
      void UseAdaptiveSampling(bool flag=true){fAdaptiveSampling=flag;}

      /// give user ability to ask for more toys
      void AdditionalNToysFactor(double fact){fAdditionalNToysFactor = fact;}

      /// save the confidence belt to a file
      void SaveBeltToFile(bool flag=true){
   fSaveBeltToFile = flag;
   if(flag) fCreateBelt = true;
      }
      /// should create confidence belt
      void CreateConfBelt(bool flag=true){fCreateBelt = flag;}

      /// Returns instance of TestStatSampler. Use to change properties of
      /// TestStatSampler, e.g. GetTestStatSampler.SetTestSize(Double_t size);
      TestStatSampler* GetTestStatSampler(void) { return fTestStatSampler; }


   private:

      Double_t fSize; /// size of the test (eg. specified rate of Type I error)
      RooAbsData& fData; /// data set
      ModelConfig &fModel;
      /*
      RooAbsPdf * fPdf; // common PDF
      mutable RooArgSet fPOI; // RooArgSet specifying  parameters of interest for interval
      RooArgSet fNuisParams;// RooArgSet specifying  nuisance parameters for interval
      */

      TestStatSampler* fTestStatSampler;
      RooAbsData* fPointsToTest;
      Double_t fLeftSideFraction;
      ConfidenceBelt* fConfBelt;
      bool fAdaptiveSampling; // controls use of adaptive sampling algorithm
      Double_t fAdditionalNToysFactor; // give user ability to ask for more toys
      bool fSaveBeltToFile; // controls use if ConfidenceBelt should be saved to a TFile
      bool fCreateBelt; // controls use if ConfidenceBelt should be saved to a TFile

   protected:
      ClassDef(NeymanConstruction,1)   // Interface for tools setting limits (producing confidence intervals)
   };
}


#endif
