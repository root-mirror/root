// @(#)root/roostats:$Id$
// Authors: Kyle Cranmer, Sven Kreiss    June 2010
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_RatioOfProfiledLikelihoodsTestStat
#define ROOSTATS_RatioOfProfiledLikelihoodsTestStat


#include "Rtypes.h"

#include "RooNLLVar.h"

#include "RooStats/TestStatistic.h"

#include "RooStats/ProfileLikelihoodTestStat.h"


namespace RooStats {

class RatioOfProfiledLikelihoodsTestStat: public TestStatistic {

  public:

   RatioOfProfiledLikelihoodsTestStat() :
      fNullProfile(),
      fAltProfile(),
      fAltPOI(NULL),
      fSubtractMLE(true),
      fDetailedOutputEnabled(false),
      fDetailedOutput(NULL)
   {
      // Proof constructor. Don't use.
   }

  RatioOfProfiledLikelihoodsTestStat(RooAbsPdf& nullPdf, RooAbsPdf& altPdf,
                 const RooArgSet* altPOI=0) :
    fNullProfile(nullPdf),
    fAltProfile(altPdf),
    fSubtractMLE(true),
    fDetailedOutputEnabled(false),
    fDetailedOutput(NULL)
      {
         //  Calculates the ratio of profiled likelihoods.

   if(altPOI)
     fAltPOI = (RooArgSet*) altPOI->snapshot();
   else
     fAltPOI = new RooArgSet(); // empty set

      }

    //__________________________________________
    ~RatioOfProfiledLikelihoodsTestStat(void) {
      if(fAltPOI) delete fAltPOI;
      if(fDetailedOutput) delete fDetailedOutput;
    }


   // returns -logL(poi, conditional MLE of nuisance params)
   // it does not subtract off the global MLE
   // because  nuisance parameters of null and alternate may not
   // be the same.
   Double_t ProfiledLikelihood(RooAbsData& data, RooArgSet& poi, RooAbsPdf& pdf);

    // evaluate the ratio of profile likelihood
   virtual Double_t Evaluate(RooAbsData& data, RooArgSet& nullParamsOfInterest);

   virtual void EnableDetailedOutput( bool e=true ) {
      fDetailedOutputEnabled = e;
      fNullProfile.EnableDetailedOutput(fDetailedOutputEnabled);
      fAltProfile.EnableDetailedOutput(fDetailedOutputEnabled);
   }

   static void SetAlwaysReuseNLL(Bool_t flag);

   void SetReuseNLL(Bool_t flag) {
      fNullProfile.SetReuseNLL(flag);
      fAltProfile.SetReuseNLL(flag);
   }

   void SetMinimizer(const char* minimizer){
      fNullProfile.SetMinimizer(minimizer);
      fAltProfile.SetMinimizer(minimizer);
   }
   void SetStrategy(Int_t strategy){
      fNullProfile.SetStrategy(strategy);
      fAltProfile.SetStrategy(strategy);
   }
   void SetTolerance(Double_t tol){
      fNullProfile.SetTolerance(tol);
      fAltProfile.SetTolerance(tol);
   }
   void SetPrintLevel(Int_t printLevel){
      fNullProfile.SetPrintLevel(printLevel);
      fAltProfile.SetPrintLevel(printLevel);
   }

     // set the conditional observables which will be used when creating the NLL
     // so the pdf's will not be normalized on the conditional observables when computing the NLL
     virtual void SetConditionalObservables(const RooArgSet& set) {
        fNullProfile.SetConditionalObservables(set);
        fAltProfile.SetConditionalObservables(set);
    }

     virtual const RooArgSet* GetDetailedOutput(void) const {
        // Returns detailed output. The value returned by this function is updated after each call to Evaluate().
        // The returned RooArgSet contains the following for the alternative and null hypotheses:
        //  - the minimum nll, fitstatus and convergence quality for each fit
        //  - for each fit and for each non-constant parameter, the value, error and pull of the parameter are stored
        return fDetailedOutput;
     }




   virtual const TString GetVarName() const { return "log(L(#mu_{1},#hat{#nu}_{1}) / L(#mu_{0},#hat{#nu}_{0}))"; }

    //    const bool PValueIsRightTail(void) { return false; } // overwrites default

    void SetSubtractMLE(bool subtract){fSubtractMLE = subtract;}

  private:

    ProfileLikelihoodTestStat fNullProfile;
    ProfileLikelihoodTestStat fAltProfile;

    RooArgSet* fAltPOI;
    Bool_t fSubtractMLE;
   static Bool_t fgAlwaysReuseNll ;

    bool fDetailedOutputEnabled;
    RooArgSet* fDetailedOutput;


  protected:
    ClassDef(RatioOfProfiledLikelihoodsTestStat,3)  // implements the ratio of profiled likelihood as test statistic
};

}


#endif
