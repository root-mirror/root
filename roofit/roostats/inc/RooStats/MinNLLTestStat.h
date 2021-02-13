// @(#)root/roostats:$Id: MinNLLTestStat.h 43035 2012-02-16 16:48:57Z sven $
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
// Additional Contributions: Giovanni Petrucciani
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_MinNLLTestStat
#define ROOSTATS_MinNLLTestStat


#include "Rtypes.h"

#include "RooStats/RooStatsUtils.h"

//#include "RooStats/DistributionCreator.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/TestStatistic.h"

#include "RooStats/RooStatsUtils.h"

#include "RooRealVar.h"
#include "RooProfileLL.h"
#include "RooNLLVar.h"
#include "RooMsgService.h"

#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"
#include "ProfileLikelihoodTestStat.h"

namespace RooStats {

/** \class MinNLLTestStat
   \ingroup Roostats

MinNLLTestStat is an implementation of the TestStatistic interface that
calculates the minimum value of the negative log likelihood
function and returns it as a test statistic.
Internally it operates by delegating to a MinNLLTestStat object.

*/

  class MinNLLTestStat : public TestStatistic{

   public:
     MinNLLTestStat() {
        // Proof constructor. Do not use.
   fProflts = 0;
     }
     MinNLLTestStat(RooAbsPdf& pdf) {
   fProflts = new ProfileLikelihoodTestStat(pdf);
     }

     MinNLLTestStat(const MinNLLTestStat& rhs) : TestStatistic(rhs), fProflts(0) {
        RooAbsPdf * pdf = rhs.fProflts->GetPdf();
        if (pdf)  fProflts = new ProfileLikelihoodTestStat(*pdf);
     }

     MinNLLTestStat & operator=(const MinNLLTestStat& rhs)  {
        if (this == &rhs) return *this;
        RooAbsPdf * pdf = rhs.fProflts->GetPdf();
        if (fProflts) delete fProflts;
        fProflts = NULL;
        if (pdf)  fProflts = new ProfileLikelihoodTestStat(*pdf);
        return *this;
     }

     virtual ~MinNLLTestStat() {
   delete fProflts;
     }

     void SetOneSided(Bool_t flag=true) {fProflts->SetOneSided(flag);}
     void SetOneSidedDiscovery(Bool_t flag=true) {fProflts->SetOneSidedDiscovery(flag);}
     void SetReuseNLL(Bool_t flag) { fProflts->SetReuseNLL(flag); }
     void SetMinimizer(const char* minimizer){ fProflts->SetMinimizer(minimizer); }
     void SetStrategy(Int_t strategy){ fProflts->SetStrategy(strategy); }
     void SetTolerance(double tol){ fProflts->SetTolerance(tol); }
     void SetPrintLevel(Int_t printlevel){ fProflts->SetPrintLevel(printlevel); }
     void SetLOffset(Bool_t flag=kTRUE) { fProflts->SetLOffset(flag) ; }

     // Main interface to evaluate the test statistic on a dataset
     virtual Double_t Evaluate(RooAbsData& data, RooArgSet& paramsOfInterest) {
       return fProflts->EvaluateProfileLikelihood(1, data, paramsOfInterest); //find unconditional NLL minimum
     }

     virtual void EnableDetailedOutput( bool e=true ) { fProflts->EnableDetailedOutput(e); }

     virtual const RooArgSet* GetDetailedOutput(void) const {
        // Returns detailed output. The value returned by this function is updated after each call to Evaluate().
        // The returned RooArgSet contains the following:
        //
        //  - the minimum nll, fitstatus and convergence quality for each fit </li>
        //  - for all non-constant parameters their value, error and pull </li>
        return fProflts->GetDetailedOutput();
     }

     virtual void SetVarName(const char* name) { fProflts->SetVarName(name); }

     virtual const TString GetVarName() const { return fProflts->GetVarName(); }

   private:
     ProfileLikelihoodTestStat* fProflts;

   protected:
      ClassDef(MinNLLTestStat,1)   // implements the minimum NLL as a test statistic to be used with several tools
   };
}


#endif
