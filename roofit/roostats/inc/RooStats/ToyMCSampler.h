// @(#)root/roostats:$Id$
// Author: Sven Kreiss and Kyle Cranmer    June 2010
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
// Additions and modifications by Mario Pelliccioni
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_ToyMCSampler
#define ROOSTATS_ToyMCSampler


#include "Rtypes.h"

#include <vector>
#include <list>
#include <string>
#include <sstream>

#include "RooStats/TestStatSampler.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/TestStatistic.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProofConfig.h"

#include "RooWorkspace.h"
#include "RooMsgService.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"

#include "RooDataSet.h"

namespace RooStats {

  class DetailedOutputAggregator;

class NuisanceParametersSampler {

   public:
      NuisanceParametersSampler(RooAbsPdf *prior=NULL, const RooArgSet *parameters=NULL, Int_t nToys=1000, Bool_t asimov=kFALSE) :
         fPrior(prior),
         fParams(parameters),
         fNToys(nToys),
         fExpected(asimov),
         fPoints(NULL),
         fIndex(0)
      {
         if(prior) Refresh();
      }
      virtual ~NuisanceParametersSampler() {
         if(fPoints) { delete fPoints; fPoints = NULL; }
      }

      void NextPoint(RooArgSet& nuisPoint, Double_t& weight);

   protected:
      void Refresh();

   private:
      RooAbsPdf *fPrior;           // prior for nuisance parameters
      const RooArgSet *fParams;    // nuisance parameters
      Int_t fNToys;
      Bool_t fExpected;

      RooAbsData *fPoints;         // generated nuisance parameter points
      Int_t fIndex;                // current index in fPoints array
};

class ToyMCSampler: public TestStatSampler {

   public:

      ToyMCSampler();
      ToyMCSampler(TestStatistic &ts, Int_t ntoys);
      virtual ~ToyMCSampler();

      static void SetAlwaysUseMultiGen(Bool_t flag);

      void SetUseMultiGen(Bool_t flag) { fUseMultiGen = flag ; }

      // main interface
      virtual SamplingDistribution* GetSamplingDistribution(RooArgSet& paramPoint);
      virtual RooDataSet* GetSamplingDistributions(RooArgSet& paramPoint);
      virtual RooDataSet* GetSamplingDistributionsSingleWorker(RooArgSet& paramPoint);

      virtual SamplingDistribution* AppendSamplingDistribution(
         RooArgSet& allParameters,
         SamplingDistribution* last,
         Int_t additionalMC
      );


      // The pdf can be NULL in which case the density from SetPdf()
      // is used. The snapshot and TestStatistic is also optional.
      virtual void AddTestStatistic(TestStatistic* t = NULL) {
         if( t == NULL ) {
            oocoutI((TObject*)0,InputArguments) << "No test statistic given. Doing nothing." << std::endl;
            return;
         }

         //if( t == NULL && fTestStatistics.size() >= 1 ) t = fTestStatistics[0];

         fTestStatistics.push_back( t );
      }

      // generates toy data
      //   without weight
      virtual RooAbsData* GenerateToyData(RooArgSet& paramPoint, RooAbsPdf& pdf) const {
         if(fExpectedNuisancePar) oocoutE((TObject*)NULL,InputArguments) << "ToyMCSampler: using expected nuisance parameters but ignoring weight. Use GetSamplingDistribution(paramPoint, weight) instead." << std::endl;
         double weight;
         return GenerateToyData(paramPoint, weight, pdf);
      }
      virtual RooAbsData* GenerateToyData(RooArgSet& paramPoint) const { return GenerateToyData(paramPoint,*fPdf); }
      //   with weight
      virtual RooAbsData* GenerateToyData(RooArgSet& paramPoint, double& weight, RooAbsPdf& pdf) const;
      virtual RooAbsData* GenerateToyData(RooArgSet& paramPoint, double& weight) const { return GenerateToyData(paramPoint,weight,*fPdf); }

      // generate global observables
      virtual void GenerateGlobalObservables(RooAbsPdf& pdf) const;


      // Main interface to evaluate the test statistic on a dataset
      virtual Double_t EvaluateTestStatistic(RooAbsData& data, RooArgSet& nullPOI, int i ) {
         return fTestStatistics[i]->Evaluate(data, nullPOI);
      }
      virtual Double_t EvaluateTestStatistic(RooAbsData& data, RooArgSet& nullPOI) { return EvaluateTestStatistic( data,nullPOI, 0 ); }
      virtual RooArgList* EvaluateAllTestStatistics(RooAbsData& data, const RooArgSet& poi);


      virtual TestStatistic* GetTestStatistic(unsigned int i) const {
         if( fTestStatistics.size() <= i ) return NULL;
         return fTestStatistics[i];
      }
      virtual TestStatistic* GetTestStatistic(void) const { return GetTestStatistic(0); }

      virtual Double_t ConfidenceLevel() const { return 1. - fSize; }
      virtual void Initialize(
         RooAbsArg& /*testStatistic*/,
         RooArgSet& /*paramsOfInterest*/,
         RooArgSet& /*nuisanceParameters*/
      ) {}

      virtual Int_t GetNToys(void) { return fNToys; }
      virtual void SetNToys(const Int_t ntoy) { fNToys = ntoy; }
      /// Forces the generation of exactly `n` events even for extended PDFs. Set to 0 to
      /// use the Poisson-distributed events from the extended PDF.
      virtual void SetNEventsPerToy(const Int_t nevents) {
         fNEvents = nevents;
      }


      // Set the Pdf, add to the the workspace if not already there
      virtual void SetParametersForTestStat(const RooArgSet& nullpoi) {
         if( fParametersForTestStat ) delete fParametersForTestStat;
         fParametersForTestStat = (const RooArgSet*)nullpoi.snapshot();
      }

      virtual void SetPdf(RooAbsPdf& pdf) { fPdf = &pdf; ClearCache(); }

      // How to randomize the prior. Set to NULL to deactivate randomization.
      virtual void SetPriorNuisance(RooAbsPdf* pdf) {
         fPriorNuisance = pdf;
         if (fNuisanceParametersSampler) {
            delete fNuisanceParametersSampler;
            fNuisanceParametersSampler = NULL;
         }
      }
      // specify the nuisance parameters (eg. the rest of the parameters)
      virtual void SetNuisanceParameters(const RooArgSet& np) { fNuisancePars = &np; }
      // specify the observables in the dataset (needed to evaluate the test statistic)
      virtual void SetObservables(const RooArgSet& o) { fObservables = &o; }
      // specify the conditional observables
      virtual void SetGlobalObservables(const RooArgSet& o) { fGlobalObservables = &o; }


      // set the size of the test (rate of Type I error) ( Eg. 0.05 for a 95% Confidence Interval)
      virtual void SetTestSize(Double_t size) { fSize = size; }
      // set the confidence level for the interval (eg. 0.95 for a 95% Confidence Interval)
      virtual void SetConfidenceLevel(Double_t cl) { fSize = 1. - cl; }

      // Set the TestStatistic (want the argument to be a function of the data & parameter points
      virtual void SetTestStatistic(TestStatistic *testStatistic, unsigned int i) {
         if( fTestStatistics.size() < i ) {
            oocoutE((TObject*)NULL,InputArguments) << "Cannot set test statistic for this index." << std::endl;
            return;
         }
    if( fTestStatistics.size() == i)
       fTestStatistics.push_back(testStatistic);
    else
       fTestStatistics[i] = testStatistic;
      }
      virtual void SetTestStatistic(TestStatistic *t) { return SetTestStatistic(t,0); }

      virtual void SetExpectedNuisancePar(Bool_t i = kTRUE) { fExpectedNuisancePar = i; }
      virtual void SetAsimovNuisancePar(Bool_t i = kTRUE) { fExpectedNuisancePar = i; }

      // Checks for sufficient information to do a GetSamplingDistribution(...).
      Bool_t CheckConfig(void);

      // control to use bin data generation (=> see RooFit::AllBinned() option)
      void SetGenerateBinned(bool binned = true) { fGenerateBinned = binned; }
      // name of the tag for individual components to be generated binned (=> see RooFit::GenBinned() option)
      void SetGenerateBinnedTag( const char* binnedTag = "" ) { fGenerateBinnedTag = binnedTag; }
      // set auto binned generation (=> see RooFit::AutoBinned() option)
      void SetGenerateAutoBinned( Bool_t autoBinned = kTRUE ) { fGenerateAutoBinned = autoBinned; }

      // Set the name of the sampling distribution used for plotting
      void SetSamplingDistName(const char* name) { if(name) fSamplingDistName = name; }
      std::string GetSamplingDistName(void) { return fSamplingDistName; }

      // This option forces a maximum number of total toys.
      void SetMaxToys(Double_t t) { fMaxToys = t; }

      void SetToysLeftTail(Double_t toys, Double_t threshold) {
         fToysInTails = toys;
         fAdaptiveLowLimit = threshold;
         fAdaptiveHighLimit = RooNumber::infinity();
      }
      void SetToysRightTail(Double_t toys, Double_t threshold) {
         fToysInTails = toys;
         fAdaptiveHighLimit = threshold;
         fAdaptiveLowLimit = -RooNumber::infinity();
      }
      void SetToysBothTails(Double_t toys, Double_t low_threshold, Double_t high_threshold) {
         fToysInTails = toys;
         fAdaptiveHighLimit = high_threshold;
         fAdaptiveLowLimit = low_threshold;
      }

      // calling with argument or NULL deactivates proof
      void SetProofConfig(ProofConfig *pc = NULL) { fProofConfig = pc; }

      void SetProtoData(const RooDataSet* d) { fProtoData = d; }

   protected:

      const RooArgList* EvaluateAllTestStatistics(RooAbsData& data, const RooArgSet& poi, DetailedOutputAggregator& detOutAgg);

      // helper for GenerateToyData
      RooAbsData* Generate(RooAbsPdf &pdf, RooArgSet &observables, const RooDataSet *protoData=NULL, int forceEvents=0) const;

      // helper method for clearing  the cache
      virtual void ClearCache();


      // densities, snapshots, and test statistics to reweight to
      RooAbsPdf *fPdf; // model (can be alt or null)
      const RooArgSet* fParametersForTestStat;
      std::vector<TestStatistic*> fTestStatistics;

      std::string fSamplingDistName; // name of the model
      RooAbsPdf *fPriorNuisance; // prior pdf for nuisance parameters
      const RooArgSet *fNuisancePars;
      const RooArgSet *fObservables;
      const RooArgSet *fGlobalObservables;
      Int_t fNToys; // number of toys to generate
      Int_t fNEvents; // number of events per toy (may be ignored depending on settings)
      Double_t fSize;
      Bool_t fExpectedNuisancePar; // whether to use expectation values for nuisance parameters (ie Asimov data set)
      Bool_t fGenerateBinned;
      TString fGenerateBinnedTag;
      Bool_t fGenerateAutoBinned;

      // minimum no of toys in tails for adaptive sampling
      // (taking weights into account, therefore double)
      // Default: 0.0 which means no adaptive sampling
      Double_t fToysInTails;
      // maximum no of toys
      // (taking weights into account, therefore double)
      Double_t fMaxToys;
      // tails
      Double_t fAdaptiveLowLimit;
      Double_t fAdaptiveHighLimit;

      const RooDataSet *fProtoData; // in dev

      ProofConfig *fProofConfig;   //!

      mutable NuisanceParametersSampler *fNuisanceParametersSampler; //!

      // objects below cache information and are mutable and non-persistent
      mutable RooArgSet* _allVars ; //!
      mutable std::list<RooAbsPdf*> _pdfList ; //!
      mutable std::list<RooArgSet*> _obsList ; //!
      mutable std::list<RooAbsPdf::GenSpec*> _gsList ; //!
      mutable RooAbsPdf::GenSpec* _gs1 ; //! GenSpec #1
      mutable RooAbsPdf::GenSpec* _gs2 ; //! GenSpec #2
      mutable RooAbsPdf::GenSpec* _gs3 ; //! GenSpec #3
      mutable RooAbsPdf::GenSpec* _gs4 ; //! GenSpec #4

      static Bool_t fgAlwaysUseMultiGen ;  // Use PrepareMultiGen always
      Bool_t fUseMultiGen ; // Use PrepareMultiGen?

   protected:
   ClassDef(ToyMCSampler,3) // A simple implementation of the TestStatSampler interface
};
}


#endif
