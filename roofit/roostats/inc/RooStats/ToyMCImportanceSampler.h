// @(#)root/roostats:$Id$
// Author: Sven Kreiss and Kyle Cranmer    January 2012
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_ToyMCImportanceSampler
#define ROOSTATS_ToyMCImportanceSampler

#include "RooStats/ToyMCSampler.h"
#include <vector>

namespace RooStats {

enum toysStrategies { EQUALTOYSPERDENSITY, EXPONENTIALTOYDISTRIBUTION };

class ToyMCImportanceSampler: public ToyMCSampler {

   public:
      ToyMCImportanceSampler() :
         ToyMCSampler()
      {
         // Proof constructor. Do not use.

         fIndexGenDensity = 0;
         fGenerateFromNull = true;
         fApplyVeto = true;
         fReuseNLL = true;
         fToysStrategy = EQUALTOYSPERDENSITY;
      }
      ToyMCImportanceSampler(TestStatistic &ts, Int_t ntoys) :
         ToyMCSampler(ts, ntoys)
      {
         fIndexGenDensity = 0;
         fGenerateFromNull = true;
         fApplyVeto = true;
         fReuseNLL = true;
         fToysStrategy = EQUALTOYSPERDENSITY;
      }

      virtual ~ToyMCImportanceSampler();

      // overwrite GetSamplingDistributionsSingleWorker(paramPoint) with a version that loops
      // over nulls and importance densities, but calls the parent
      // ToyMCSampler::GetSamplingDistributionsSingleWorker(paramPoint).
      virtual RooDataSet* GetSamplingDistributionsSingleWorker(RooArgSet& paramPoint);

      using ToyMCSampler::GenerateToyData;
      virtual RooAbsData* GenerateToyData(RooArgSet& paramPoint, double& weight) const;
      virtual RooAbsData* GenerateToyData(RooArgSet& paramPoint, double& weight, std::vector<double>& impNLLs, double& nullNLL) const;
      virtual RooAbsData* GenerateToyData(std::vector<double>& weights) const;
      virtual RooAbsData* GenerateToyData(std::vector<double>& weights, std::vector<double>& nullNLLs, std::vector<double>& impNLLs) const;


      /// specifies the pdf to sample from
      void SetDensityToGenerateFromByIndex(unsigned int i, bool fromNull = false) {
         if( (fromNull  &&  i >= fNullDensities.size())  ||
             (!fromNull &&  i >= fImportanceDensities.size())
         ) {
            oocoutE((TObject*)0,InputArguments) << "Index out of range. Requested index: "<<i<<
               " , but null densities: "<<fNullDensities.size()<<
               " and importance densities: "<<fImportanceDensities.size() << std::endl;
         }

         fIndexGenDensity = i;
         fGenerateFromNull = fromNull;

         ClearCache();
      }

      // For importance sampling with multiple densities/snapshots:
      // This is used to check the current Likelihood against Likelihoods from
      // other importance densities apart from the one given as importance snapshot.
      // The pdf can be NULL in which case the density from SetImportanceDensity()
      // is used. The snapshot is also optional.
      void AddImportanceDensity(RooAbsPdf* p, const RooArgSet* s) {
         if( p == NULL && s == NULL ) {
            oocoutI((TObject*)0,InputArguments) << "Neither density nor snapshot given. Doing nothing." << std::endl;
            return;
         }
         if( p == NULL && fPdf == NULL ) {
            oocoutE((TObject*)0,InputArguments) << "No density given, but snapshot is there. Aborting." << std::endl;
            return;
         }

         if( p == NULL ) p = fPdf;

         if( s ) s = (const RooArgSet*)s->snapshot();

         fImportanceDensities.push_back( p );
         fImportanceSnapshots.push_back( s );
         fImpNLLs.push_back( NULL );
      }

      // The pdf can be NULL in which case the density from SetPdf()
      // is used. The snapshot and TestStatistic is also optional.
      void AddNullDensity(RooAbsPdf* p, const RooArgSet* s = NULL) {
         if( p == NULL && s == NULL ) {
            oocoutI((TObject*)0,InputArguments) << "Neither density nor snapshot nor test statistic given. Doing nothing." << std::endl;
            return;
         }

         if( p == NULL && fNullDensities.size() >= 1 ) p = fNullDensities[0];
         if( s == NULL ) s = fParametersForTestStat;
         if( s ) s = (const RooArgSet*)s->snapshot();

         fNullDensities.push_back( p );
         fNullSnapshots.push_back( s );
         fNullNLLs.push_back( NULL );
         ClearCache();
      }
      // overwrite from ToyMCSampler
      virtual void SetPdf(RooAbsPdf& pdf) {
         ToyMCSampler::SetPdf(pdf);

         if( fNullDensities.size() == 1 ) { fNullDensities[0] = &pdf; }
         else if( fNullDensities.size() == 0) AddNullDensity( &pdf );
         else{
            oocoutE((TObject*)0,InputArguments) << "Cannot use SetPdf() when already multiple null densities are specified. Please use AddNullDensity()." << std::endl;
         }
      }
      // overwrite from ToyMCSampler
      void SetParametersForTestStat(const RooArgSet& nullpoi) {
         ToyMCSampler::SetParametersForTestStat(nullpoi);
         if( fNullSnapshots.size() == 0 ) AddNullDensity( NULL, &nullpoi );
         else if( fNullSnapshots.size() == 1 ) {
            oocoutI((TObject*)0,InputArguments) << "Overwriting snapshot for the only defined null density." << std::endl;
            if( fNullSnapshots[0] ) delete fNullSnapshots[0];
            fNullSnapshots[0] = (const RooArgSet*)nullpoi.snapshot();
         }else{
            oocoutE((TObject*)0,InputArguments) << "Cannot use SetParametersForTestStat() when already multiple null densities are specified. Please use AddNullDensity()." << std::endl;
         }
      }

      // When set to true, this sets the weight of all toys to zero that
      // do not have the largest likelihood under the density it was generated
      // compared to the other densities.
      void SetApplyVeto(bool b = true) { fApplyVeto = b; }

      void SetReuseNLL(bool r = true) { fReuseNLL = r; }

     // set the conditional observables which will be used when creating the NLL
     // so the pdf's will not be normalized on the conditional observables when computing the NLL
     // Since the class use a NLL we need to set the conditional observables if they exist in the model
     virtual void SetConditionalObservables(const RooArgSet& set) {fConditionalObs.removeAll(); fConditionalObs.add(set);}

      int CreateNImpDensitiesForOnePOI(
         RooAbsPdf& pdf,
         const RooArgSet& allPOI,
         RooRealVar& poi,
         int n,
         double poiValueForBackground = 0.0
      );
      int CreateImpDensitiesForOnePOIAdaptively(
         RooAbsPdf& pdf,
         const RooArgSet& allPOI,
         RooRealVar& poi,
         double nStdDevOverlap = 0.5,
         double poiValueForBackground = 0.0
      );

      void SetEqualNumToysPerDensity( void ) { fToysStrategy = EQUALTOYSPERDENSITY; }
      void SetExpIncreasingNumToysPerDensity( void ) { fToysStrategy = EXPONENTIALTOYDISTRIBUTION; }

   protected:

      // helper method for clearing  the cache
      virtual void ClearCache();

      unsigned int fIndexGenDensity;
      bool fGenerateFromNull;
      bool fApplyVeto;

   RooArgSet fConditionalObs;  // set of conditional observables

      // support multiple null densities
      std::vector<RooAbsPdf*> fNullDensities;
      mutable std::vector<const RooArgSet*> fNullSnapshots;

      // densities and snapshots to generate from
      std::vector<RooAbsPdf*> fImportanceDensities;
      std::vector<const RooArgSet*> fImportanceSnapshots;

      bool fReuseNLL;

      toysStrategies fToysStrategy;

      mutable std::vector<RooAbsReal*> fNullNLLs;    //!
      mutable std::vector<RooAbsReal*> fImpNLLs;     //!

   protected:
   ClassDef(ToyMCImportanceSampler,2) // An implementation of importance sampling
};
}

#endif
