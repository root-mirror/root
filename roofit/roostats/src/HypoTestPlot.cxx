// @(#)root/roostats:$Id$
// Author: Sven Kreiss   June 2010
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class RooStats::HypoTestPlot
    \ingroup Roostats

This class provides the plots for the result of a study performed with any of the
HypoTestCalculatorGeneric (e.g. HybridCalculator or FrequentistCalculator)  class.

   */

#include "RooStats/HypoTestPlot.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/SamplingDistribution.h"

ClassImp(RooStats::HypoTestPlot);

using namespace RooStats;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

HypoTestPlot::HypoTestPlot(HypoTestResult& result, Int_t bins, Option_t* opt) :
   SamplingDistPlot(bins),
   fHypoTestResult(&result)
{
   ApplyResult(result, opt);
}
HypoTestPlot::HypoTestPlot(HypoTestResult& result, Int_t bins, Double_t min, Double_t max, Option_t* opt) :
   SamplingDistPlot(bins,min,max),
   fHypoTestResult(&result)
{
   ApplyResult(result, opt);
}

////////////////////////////////////////////////////////////////////////////////

void HypoTestPlot::ApplyResult(HypoTestResult& result, Option_t* opt) {
   fLegend = new TLegend(0.55,0.95-0.3*0.66,0.95,0.95);

   const SamplingDistribution *alt = result.GetAltDistribution();
   const SamplingDistribution *null = result.GetNullDistribution();
   if(!result.HasTestStatisticData()) {
      if(alt) AddSamplingDistribution(alt, opt);
      if(null) AddSamplingDistribution(null, opt);
   }else{
      if(result.GetPValueIsRightTail()) {
         if(alt) AddSamplingDistributionShaded(alt, result.GetTestStatisticData(), RooNumber::infinity(), opt);
         if(null) AddSamplingDistributionShaded(null, result.GetTestStatisticData(), RooNumber::infinity() , opt);
      }else{
         if(alt) AddSamplingDistributionShaded(alt, -RooNumber::infinity(), result.GetTestStatisticData(), opt);
         if(null) AddSamplingDistributionShaded(null, - RooNumber::infinity(), result.GetTestStatisticData() , opt);
      }
   }

   if(result.HasTestStatisticData()) {
      Double_t theMin(0.), theMax(0.), theYMax(0.);
      GetAbsoluteInterval(theMin, theMax, theYMax);

      AddLine(result.GetTestStatisticData(), 0, result.GetTestStatisticData(), theYMax*0.66, "test statistic data");
   }

   ApplyDefaultStyle();
}

////////////////////////////////////////////////////////////////////////////////

void HypoTestPlot::ApplyDefaultStyle(void) {
   if(!fHypoTestResult) return;

   const SamplingDistribution *alt = fHypoTestResult->GetAltDistribution();
   const SamplingDistribution *null = fHypoTestResult->GetNullDistribution();

   if(alt) {
      SetLineWidth(2, alt);
      SetLineColor(kBlue, alt);
   }
   if(null) {
      SetLineWidth(2, null);
      SetLineColor(kRed, null);
   }
}
