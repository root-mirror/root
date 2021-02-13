// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_DebuggingTestStat
#define ROOSTATS_DebuggingTestStat

/** \class DebuggingTestStat
    \ingroup Roostats

DebuggingTestStat is a simple implementation of the DistributionCreator interface used for debugging.
The sampling distribution is uniformly random between [0,1] and is INDEPENDENT of the data.  So it is not useful
for true statistical tests, but it is useful for debugging.

*/

#include "Rtypes.h"

//#include "RooStats/DistributionCreator.h"
#include "RooStats/TestStatistic.h"
#include "RooStats/ToyMCSampler.h"

#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "SamplingDistribution.h"
#include "TRandom.h"

namespace RooStats {

  class DebuggingTestStat : public TestStatistic {

   public:
     DebuggingTestStat() {
       fTestStatistic = new RooRealVar("UniformTestStatistic","UniformTestStatistic",0,0,1);
       fRand = new TRandom();
     }
     virtual ~DebuggingTestStat() {
       //       delete fRand;
       //       delete fTestStatistic;
     }

     // Main interface to evaluate the test statistic on a dataset
     virtual Double_t Evaluate(RooAbsData& /*data*/, RooArgSet& /*paramsOfInterest*/)  {
       //data = data; // avoid warning
       //paramsOfInterest = paramsOfInterest; //avoid warning
       return fRand->Uniform();
     }




   private:

      RooRealVar* fTestStatistic;
      TRandom* fRand;

   protected:
      ClassDef(DebuggingTestStat,1)   // A concrete implementation of the TestStatistic interface, useful for debugging.
   };

}


#endif
