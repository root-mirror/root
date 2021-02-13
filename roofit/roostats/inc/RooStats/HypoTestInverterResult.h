// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_HypoTestInverterResult
#define ROOSTATS_HypoTestInverterResult

#include "RooStats/SimpleInterval.h"

#include "RooStats/HypoTestResult.h"

#include <vector>

class RooRealVar;

namespace RooStats {

class SamplingDistribution;

class HypoTestInverterResult : public SimpleInterval {

public:

   /// default constructor
   explicit HypoTestInverterResult(const char* name = 0);

   /// constructor
   HypoTestInverterResult( const char* name,
                           const RooRealVar& scannedVariable,
                           double cl ) ;

   HypoTestInverterResult( const HypoTestInverterResult& other, const char* name );

   /// destructor
   virtual ~HypoTestInverterResult();

   /// operator =
   HypoTestInverterResult& operator = (const HypoTestInverterResult& other);

   /// remove points that appear to have failed.
   int ExclusionCleanup();

   /// merge with the content of another HypoTestInverterResult object
   bool Add( const HypoTestInverterResult& otherResult );

   ///add the result of a single point (an HypoTestRsult)
   bool Add( Double_t x, const HypoTestResult & result );

   /// function to return the value of the parameter of interest for the i^th entry in the results
   double GetXValue( int index ) const ;

   /// function to return the value of the confidence level for the i^th entry in the results
   double GetYValue( int index ) const ;

   /// function to return the estimated error on the value of the confidence level for the i^th entry in the results
   double GetYError( int index ) const ;

   /// return the observed CLsplusb value  for the i-th entry
   double CLsplusb( int index) const;

   /// return the observed CLb value  for the i-th entry
   double CLb( int index) const;

   /// return the observed CLb value  for the i-th entry
   double CLs( int index) const;

   /// return the observed CLsplusb value  for the i-th entry
   double CLsplusbError( int index) const;

   /// return the observed CLb value  for the i-th entry
   double CLbError( int index) const;

   /// return the observed CLb value  for the i-th entry
   double CLsError( int index) const;

   /// return a pointer to the i^th result object
   HypoTestResult* GetResult( int index ) const ;

   double GetLastYValue( ) const  { return GetYValue(  fXValues.size()-1); }

   double GetLastXValue( ) const  { return GetXValue(  fXValues.size()-1); }

   double GetLastYError( ) const  { return GetYError(  fXValues.size()-1); }

   HypoTestResult * GetLastResult( ) const  { return GetResult(  fXValues.size()-1); }

   /// number of entries in the results array
   int ArraySize() const { return fXValues.size(); };

   int FindIndex(double xvalue) const;

   /// set the size of the test (rate of Type I error) (eg. 0.05 for a 95% Confidence Interval)
   virtual void SetTestSize( Double_t size ) { fConfidenceLevel = 1.-size; }

   /// set the confidence level for the interval (eg. 0.95 for a 95% Confidence Interval)
   virtual void SetConfidenceLevel( Double_t cl ) { fConfidenceLevel = cl; }

   /// set CLs threshold for exclusion cleanup function
   inline void SetCLsCleanupThreshold( Double_t th ) { fCLsCleanupThreshold = th; }

   /// flag to switch between using CLsb (default) or CLs as confidence level
   void UseCLs( bool on = true ) { fUseCLs = on; }

   /// query if one sided result
   bool IsOneSided() const { return !fIsTwoSided; }
   /// query if two sided result
   bool IsTwoSided() const { return fIsTwoSided; }

   /// lower and upper bound of the confidence interval (to get upper/lower limits, multiply the size( = 1-confidence level ) by 2
   Double_t LowerLimit();
   Double_t UpperLimit();

   /// rough estimation of the error on the computed bound of the confidence interval
   /// Estimate of lower limit error
   ///function evaluates only a rough error on the lower limit. Be careful when using this estimation
   Double_t LowerLimitEstimatedError();

   /// Estimate of lower limit error
   ///function evaluates only a rough error on the lower limit. Be careful when using this estimation
   Double_t UpperLimitEstimatedError();

   /// return expected distribution of p-values (Cls or Clsplusb)

   SamplingDistribution * GetExpectedPValueDist(int index) const;

   SamplingDistribution * GetBackgroundTestStatDist(int index ) const;

   SamplingDistribution * GetSignalAndBackgroundTestStatDist(int index) const;

   /// same in terms of alt and null
   SamplingDistribution * GetNullTestStatDist(int index) const {
      return  GetSignalAndBackgroundTestStatDist(index);
   }
   SamplingDistribution * GetAltTestStatDist(int index) const {
      return  GetBackgroundTestStatDist(index);
   }

   /// get expected lower limit distributions
   /// implemented using interpolation
   ///  The size for the sampling distribution is given (by default is given by the average number of toy/point)
   SamplingDistribution* GetLowerLimitDistribution() const { return GetLimitDistribution(true); }

   /// get expected upper limit distributions
   /// implemented using interpolation
   SamplingDistribution* GetUpperLimitDistribution() const { return GetLimitDistribution(false); }

   /// get Limit value corresponding at the desired nsigma level (0) is median -1 sigma is 1 sigma
   double GetExpectedLowerLimit(double nsig = 0, const char * opt = "" ) const ;

   /// get Limit value corresponding at the desired nsigma level (0) is median -1 sigma is 1 sigma
   double GetExpectedUpperLimit(double nsig = 0, const char * opt = "") const ;


   double FindInterpolatedLimit(double target, bool lowSearch = false, double xmin=1, double xmax=0);

   enum InterpolOption_t { kLinear, kSpline };

   /// set the interpolation option, linear (kLinear ) or spline  (kSpline)
   void SetInterpolationOption( InterpolOption_t opt) { fInterpolOption = opt; }

   InterpolOption_t GetInterpolationOption() const { return fInterpolOption; }

private:


   double CalculateEstimatedError(double target, bool lower = true, double xmin = 1, double xmax = 0);

   int FindClosestPointIndex(double target, int mode = 0, double xtarget = 0);

   SamplingDistribution* GetLimitDistribution(bool lower ) const;

   double GetExpectedLimit(double nsig, bool lower, const char * opt = "" ) const ;

   double GetGraphX(const TGraph & g, double y0, bool lowSearch, double &xmin, double &xmax) const;
   double GetGraphX(const TGraph & g, double y0, bool lowSearch = true) const {
      double xmin=1; double xmax = 0;
      return GetGraphX(g,y0,lowSearch,xmin,xmax);
   }


protected:

   bool fUseCLs;
   bool fIsTwoSided;                  /// two sided scan (look for lower/upper limit)
   bool fInterpolateLowerLimit;
   bool fInterpolateUpperLimit;
   bool fFittedLowerLimit;
   bool fFittedUpperLimit;
   InterpolOption_t fInterpolOption;  /// interpolation option (linear or spline)

   double fLowerLimitError;
   double fUpperLimitError;

   double fCLsCleanupThreshold;

   static double fgAsymptoticMaxSigma;  /// max sigma value used to scan asymptotic expected p values
   static int fgAsymptoticNumPoints;    /// number of points used to build expected p-values

   std::vector<double> fXValues;

   TList fYObjects;       /// list of HypoTestResult for each point
   TList fExpPValues;     /// list of expected sampling distribution for each point

   friend class HypoTestInverter;
   friend class HypoTestInverterPlot;
   friend class HypoTestInverterOriginal;

   ClassDef(HypoTestInverterResult,5)  /// HypoTestInverterResult class
};
}

#endif
