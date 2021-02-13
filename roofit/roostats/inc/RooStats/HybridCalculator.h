// @(#)root/roostats:$Id$
// Author: Sven Kreiss, Kyle Cranmer   Nov 2010
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_HybridCalculator
#define ROOSTATS_HybridCalculator

#include "RooStats/HypoTestCalculatorGeneric.h"

#include "RooStats/ToyMCSampler.h"



namespace RooStats {

   class HybridCalculator : public HypoTestCalculatorGeneric {

   public:
      HybridCalculator(
                        const RooAbsData &data,
                        const ModelConfig &altModel,
                        const ModelConfig &nullModel,
                        TestStatSampler* sampler=0
      ) :
         HypoTestCalculatorGeneric(data, altModel, nullModel, sampler),
         fPriorNuisanceNull(MakeNuisancePdf(nullModel, "PriorNuisanceNull")),
         fPriorNuisanceAlt(MakeNuisancePdf(altModel, "PriorNuisanceAlt")),
         fPriorNuisanceNullExternal(false),
         fPriorNuisanceAltExternal(false),
         fNToysNull(-1),
         fNToysAlt(-1),
         fNToysNullTail(0),
         fNToysAltTail(0)
      {
      }

      ~HybridCalculator() {
         if(!fPriorNuisanceNullExternal) delete fPriorNuisanceNull;
         if(!fPriorNuisanceAltExternal) delete fPriorNuisanceAlt;
      }


      /// Override the distribution used for marginalizing nuisance parameters that is inferred from ModelConfig
      virtual void ForcePriorNuisanceNull(RooAbsPdf& priorNuisance) {
         if(!fPriorNuisanceNullExternal) delete fPriorNuisanceNull;
         fPriorNuisanceNull = &priorNuisance; 
         fPriorNuisanceNullExternal = true;
      }
      virtual void ForcePriorNuisanceAlt(RooAbsPdf& priorNuisance) {
         if(!fPriorNuisanceAltExternal) delete fPriorNuisanceAlt;
         fPriorNuisanceAlt = &priorNuisance; 
         fPriorNuisanceAltExternal = true;
      }

      virtual void SetNullModel(const ModelConfig &nullModel) {
         fNullModel = &nullModel;
         if(!fPriorNuisanceNullExternal) delete fPriorNuisanceNull;
         fPriorNuisanceNull = MakeNuisancePdf(nullModel, "PriorNuisanceNull");
         fPriorNuisanceAltExternal = false;
      }

      virtual void SetAlternateModel(const ModelConfig &altModel) {
         fAltModel = &altModel;
         if(!fPriorNuisanceAltExternal) delete fPriorNuisanceAlt;
         fPriorNuisanceAlt = MakeNuisancePdf(altModel, "PriorNuisanceAlt");
         fPriorNuisanceAltExternal = false; 
      }

      /// set number of toys
      void SetToys(int toysNull, int toysAlt) { fNToysNull = toysNull; fNToysAlt = toysAlt; }

      /// set least number of toys in tails
      void SetNToysInTails(int toysNull, int toysAlt) { fNToysNullTail = toysNull; fNToysAltTail = toysAlt; }

   protected:
      /// check whether all input is consistent
      int CheckHook(void) const;

      /// configure TestStatSampler for the Null run
      int PreNullHook(RooArgSet* /*parameterPoint*/, double obsTestStat) const;

      /// configure TestStatSampler for the Alt run
      int PreAltHook(RooArgSet* /*parameterPoint*/, double obsTestStat) const;

   protected:
      RooAbsPdf *fPriorNuisanceNull;
      RooAbsPdf *fPriorNuisanceAlt;

      // these flags tell us if the nuisance pdfs came from an external resource (via ForcePriorNuisance)
      // or were created internally and should be deleted
      Bool_t fPriorNuisanceNullExternal;
      Bool_t fPriorNuisanceAltExternal;

      // different number of toys for null and alt
      int fNToysNull;
      int fNToysAlt;

      // adaptive sampling
      int fNToysNullTail;
      int fNToysAltTail;

   protected:
      ClassDef(HybridCalculator,2)
   };
}

#endif
