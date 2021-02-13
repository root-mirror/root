// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke, Sven Kreiss
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*****************************************************************************
 * Project: RooStats
 * Package: RooFit/RooStats
 * @(#)root/roofit/roostats:$Id$
 * Authors:
 *   Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke, Sven Kreiss
 *
 *****************************************************************************/


/** \class RooStats::HypoTestResult
    \ingroup Roostats

HypoTestResult is a base class for results from hypothesis tests.
Any tool inheriting from HypoTestCalculator can return a HypoTestResult.
As such, it stores a p-value for the null-hypothesis (eg. background-only)
and an alternate hypothesis (eg. signal+background).
The p-values can also be transformed into confidence levels
(\f$CL_{b}\f$, \f$CL_{s+b}\f$) in a trivial way.
The ratio of the \f$CL_{s+b}\f$ to \f$CL_{b}\f$ is often called
\f$CL_{s}\f$, and is considered useful, though it is not a probability.
Finally, the p-value of the null can be transformed into a number of
equivalent Gaussian sigma using the Significance method.

The p-value of the null for a given test statistic is rigorously defined and
this is the starting point for the following conventions.

### Conventions used in this class

The p-value for the null and alternate are on the **same side** of the
observed value of the test statistic. This is the more standard
convention and avoids confusion when doing inverted tests.

For exclusion, we also want the formula \f$CL_{s} = CL_{s+b} / CL_{b}\f$
to hold which therefore defines our conventions for \f$CL_{s+b}\f$ and
\f$CL_{b}\f$. \f$CL_{s}\f$ was specifically invented for exclusion
and therefore all quantities need be related through the assignments
as they are for exclusion: \f$CL_{s+b} = p_{s+b}\f$; \f$CL_{b} = p_{b}\f$. This
is derived by considering the scenarios of a powerful and not powerful
inverted test, where for the not so powerful test, \f$CL_{s}\f$ must be
close to one.

For results of Hypothesis tests,
\f$CL_{s}\f$ has no similar direct interpretation as for exclusion and can
be larger than one.

*/

#include "RooStats/HypoTestResult.h"
#include "RooStats/SamplingDistribution.h"
#include "RooAbsReal.h"

#include "RooStats/RooStatsUtils.h"

#include <limits>
#define NaN numeric_limits<float>::quiet_NaN()
#define IsNaN(a) TMath::IsNaN(a)

ClassImp(RooStats::HypoTestResult); ;

using namespace RooStats;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

HypoTestResult::HypoTestResult(const char* name) :
   TNamed(name,name),
   fNullPValue(NaN), fAlternatePValue(NaN),
   fNullPValueError(0), fAlternatePValueError(0),
   fTestStatisticData(NaN),
   fAllTestStatisticsData(NULL),
   fNullDistr(NULL), fAltDistr(NULL),
   fNullDetailedOutput(NULL), fAltDetailedOutput(NULL),
   fPValueIsRightTail(kTRUE),
   fBackgroundIsAlt(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Alternate constructor

HypoTestResult::HypoTestResult(const char* name, Double_t nullp, Double_t altp) :
   TNamed(name,name),
   fNullPValue(nullp), fAlternatePValue(altp),
   fNullPValueError(0), fAlternatePValueError(0),
   fTestStatisticData(NaN),
   fAllTestStatisticsData(NULL),
   fNullDistr(NULL), fAltDistr(NULL),
   fNullDetailedOutput(NULL), fAltDetailedOutput(NULL),
   fPValueIsRightTail(kTRUE),
   fBackgroundIsAlt(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// copy constructor

HypoTestResult::HypoTestResult(const HypoTestResult& other) :
   TNamed(other),
   fNullPValue(NaN), fAlternatePValue(NaN),
   fNullPValueError(0), fAlternatePValueError(0),
   fTestStatisticData(NaN),
   fAllTestStatisticsData(NULL),
   fNullDistr(NULL), fAltDistr(NULL),
   fNullDetailedOutput(NULL), fAltDetailedOutput(NULL),
   fPValueIsRightTail( other.GetPValueIsRightTail() ),
   fBackgroundIsAlt( other.GetBackGroundIsAlt() )
{
   this->Append( &other );
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

HypoTestResult::~HypoTestResult()
{
   if( fNullDistr ) delete fNullDistr;
   if( fAltDistr ) delete fAltDistr;

   if( fNullDetailedOutput ) delete fNullDetailedOutput;
   if( fAltDetailedOutput ) delete fAltDetailedOutput;

   if( fAllTestStatisticsData ) delete fAllTestStatisticsData;
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

HypoTestResult & HypoTestResult::operator=(const HypoTestResult& other) {
   if (this == &other) return *this;
   SetName(other.GetName());
   SetTitle(other.GetTitle());
   fNullPValue = other.fNullPValue;
   fAlternatePValue = other.fAlternatePValue;
   fNullPValueError = other.fNullPValueError;
   fAlternatePValueError = other.fAlternatePValueError;
   fTestStatisticData = other.fTestStatisticData;

   if( fAllTestStatisticsData ) delete fAllTestStatisticsData;
   fAllTestStatisticsData = NULL;
   if( fNullDistr ) { delete fNullDistr; fNullDistr = NULL; }
   if( fAltDistr ) { delete fAltDistr; fAltDistr = NULL; }
   if( fNullDetailedOutput ) { delete fNullDetailedOutput; fNullDetailedOutput = NULL; }
   if( fAltDetailedOutput ) { delete fAltDetailedOutput;  fAltDetailedOutput = NULL; }
   fFitInfo = nullptr;

   fPValueIsRightTail =  other.GetPValueIsRightTail();
   fBackgroundIsAlt = other.GetBackGroundIsAlt();

   this->Append( &other );

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Add additional toy-MC experiments to the current results.
/// Use the data test statistics of the added object if it is not already
/// set (otherwise, ignore the new one).

void HypoTestResult::Append(const HypoTestResult* other) {
   if(fNullDistr)
      fNullDistr->Add(other->GetNullDistribution());
   else
      if(other->GetNullDistribution()) fNullDistr = new SamplingDistribution( *other->GetNullDistribution() );

   if(fAltDistr)
      fAltDistr->Add(other->GetAltDistribution());
   else
      if(other->GetAltDistribution()) fAltDistr = new SamplingDistribution( *other->GetAltDistribution() );


   if( fNullDetailedOutput ) {
      if( other->GetNullDetailedOutput() ) fNullDetailedOutput->append( *other->GetNullDetailedOutput() );
   }else{
      if( other->GetNullDetailedOutput() ) fNullDetailedOutput = new RooDataSet( *other->GetNullDetailedOutput() );
   }

   if( fAltDetailedOutput ) {
      if( other->GetAltDetailedOutput() ) fAltDetailedOutput->append( *other->GetAltDetailedOutput() );
   }else{
      if( other->GetAltDetailedOutput() ) fAltDetailedOutput = new RooDataSet( *other->GetAltDetailedOutput() );
   }

   if( fFitInfo ) {
      if( other->GetFitInfo() ) fFitInfo->append( *other->GetFitInfo() );
   }else{
      if( other->GetFitInfo() ) fFitInfo.reset(new RooDataSet( *other->GetFitInfo() ));
   }

   // if no data is present use the other HypoTestResult's data
   if(IsNaN(fTestStatisticData)) fTestStatisticData = other->GetTestStatisticData();

   UpdatePValue(fNullDistr, fNullPValue, fNullPValueError, kTRUE);
   UpdatePValue(fAltDistr, fAlternatePValue, fAlternatePValueError, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////

void HypoTestResult::SetAltDistribution(SamplingDistribution *alt) {
   fAltDistr = alt;
   UpdatePValue(fAltDistr, fAlternatePValue, fAlternatePValueError, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////

void HypoTestResult::SetNullDistribution(SamplingDistribution *null) {
   fNullDistr = null;
   UpdatePValue(fNullDistr, fNullPValue, fNullPValueError, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////

void HypoTestResult::SetTestStatisticData(const Double_t tsd) {
   fTestStatisticData = tsd;

   UpdatePValue(fNullDistr, fNullPValue, fNullPValueError, kTRUE);
   UpdatePValue(fAltDistr, fAlternatePValue, fAlternatePValueError, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////

void HypoTestResult::SetAllTestStatisticsData(const RooArgList* tsd) {
   if (fAllTestStatisticsData) {
      delete fAllTestStatisticsData;
      fAllTestStatisticsData = 0;
   }
   if (tsd) fAllTestStatisticsData = (const RooArgList*)tsd->snapshot();

   if( fAllTestStatisticsData  &&  fAllTestStatisticsData->getSize() > 0 ) {
      RooRealVar* firstTS = (RooRealVar*)fAllTestStatisticsData->at(0);
      if( firstTS ) SetTestStatisticData( firstTS->getVal() );
   }
}

////////////////////////////////////////////////////////////////////////////////

void HypoTestResult::SetPValueIsRightTail(Bool_t pr) {
   fPValueIsRightTail = pr;

   UpdatePValue(fNullDistr, fNullPValue, fNullPValueError, kTRUE);
   UpdatePValue(fAltDistr, fAlternatePValue, fAlternatePValueError, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////

Bool_t HypoTestResult::HasTestStatisticData(void) const {
   return !IsNaN(fTestStatisticData);
}

////////////////////////////////////////////////////////////////////////////////

Double_t HypoTestResult::NullPValueError() const {
   // compute error on Null pvalue
   return fNullPValueError;
}

////////////////////////////////////////////////////////////////////////////////
/// compute \f$CL_{b}\f$ error
/// \f$CL_{b}\f$ = 1 - NullPValue()
/// must use opposite condition that routine above

Double_t HypoTestResult::CLbError() const {
   return fBackgroundIsAlt ? fAlternatePValueError : fNullPValueError;
}

////////////////////////////////////////////////////////////////////////////////

Double_t HypoTestResult::CLsplusbError() const {
   return fBackgroundIsAlt ? fNullPValueError : fAlternatePValueError;
}

////////////////////////////////////////////////////////////////////////////////
/// Taylor expansion series approximation for standard deviation (error propagation)

Double_t HypoTestResult::SignificanceError() const {
   return NullPValueError() / ROOT::Math::normal_pdf(Significance());
}

////////////////////////////////////////////////////////////////////////////////
/// Returns an estimate of the error on \f$CL_{s}\f$ through combination of the
/// errors on \f$CL_{b}\f$ and \f$CL_{s+b}\f$:
/// \f[
/// \sigma_{CL_s} = CL_s
/// \sqrt{\left( \frac{\sigma_{CL_{s+b}}}{CL_{s+b}} \right)^2 + \left( \frac{\sigma_{CL_{b}}}{CL_{b}} \right)^2}
/// \f]

Double_t HypoTestResult::CLsError() const {
   if(!fAltDistr || !fNullDistr) return 0.0;

   // unsigned const int n_b = fNullDistr->GetSamplingDistribution().size();
   // unsigned const int n_sb = fAltDistr->GetSamplingDistribution().size();

   // if CLb() == 0 CLs = -1 so return a -1 error
   if (CLb() == 0 ) return -1;

   double cl_b_err2 = pow(CLbError(),2);
   double cl_sb_err2 = pow(CLsplusbError(),2);

   return TMath::Sqrt(cl_sb_err2 + cl_b_err2 * pow(CLs(),2))/CLb();
}

////////////////////////////////////////////////////////////////////////////////
/// updates the pvalue if sufficient data is available

void HypoTestResult::UpdatePValue(const SamplingDistribution* distr, Double_t &pvalue, Double_t &perror, Bool_t /*isNull*/) {
   if(IsNaN(fTestStatisticData)) return;
   if(!distr) return;

   /* Got to be careful for discrete distributions:
    * To get the right behaviour for limits, the p-value must
    * include the value of fTestStatistic both for Alt and Null cases
    */
   if(fPValueIsRightTail) {
      pvalue = distr->IntegralAndError(perror, fTestStatisticData, RooNumber::infinity(), kTRUE,
                                       kTRUE , kTRUE );   // always closed interval [ fTestStatistic, inf ]

   }else{
      pvalue = distr->IntegralAndError(perror, -RooNumber::infinity(), fTestStatisticData, kTRUE,
                                       kTRUE,  kTRUE  ); // // always closed  [ -inf, fTestStatistic ]
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print out some information about the results
/// Note: use Alt/Null labels for the hypotheses here as the Null
/// might be the s+b hypothesis.

void HypoTestResult::Print(Option_t * ) const
{
   bool fromToys = (fAltDistr || fNullDistr);

   std::cout << std::endl << "Results " << GetName() << ": " << endl;
   std::cout << " - Null p-value = " << NullPValue();
   if (fromToys) std::cout << " +/- " << NullPValueError();
   std::cout << std::endl;
   std::cout << " - Significance = " << Significance();
   if (fromToys) std::cout << " +/- " << SignificanceError() << " sigma";
   std::cout << std::endl;
   if(fAltDistr)
      std::cout << " - Number of Alt toys: " << fAltDistr->GetSize() << std::endl;
   if(fNullDistr)
      std::cout << " - Number of Null toys: " << fNullDistr->GetSize() << std::endl;

   if (HasTestStatisticData() ) std::cout << " - Test statistic evaluated on data: " << fTestStatisticData << std::endl;
   std::cout << " - CL_b: " << CLb();
   if (fromToys) std::cout << " +/- " << CLbError();
   std::cout << std::endl;
   std::cout << " - CL_s+b: " << CLsplusb();
   if (fromToys) std::cout << " +/- " << CLsplusbError();
   std::cout << std::endl;
   std::cout << " - CL_s: " << CLs();
   if (fromToys) std::cout << " +/- " << CLsError();
   std::cout << std::endl;

   return;
}
