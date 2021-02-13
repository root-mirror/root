// @(#)root/mathcore:$Id$
// Authors: Bartolomeu Rabacal    05/2010
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
// Header file for GoFTest

#ifndef ROOT_Math_GoFTest
#define ROOT_Math_GoFTest

#include "Math/WrappedFunction.h"
#include "TMath.h"

#include <memory>
#include <vector>

/*
*/

namespace ROOT {

   namespace Fit { 
      class BinData; 
   }
namespace Math {

/////  @defgroup GoFClasses Goodness of Fit Statistical Tests Tools 
   
/*
  Class for Goodness of Fit tests implementing the Anderson-Darling and Kolmogorov-Smirnov 1- and 2-Samples Goodness of Fit Tests.
  @ingroup MathCore

 */

   
class GoFTest {
public:

   enum EDistribution { // H0 distributions for using only with 1-sample tests
      kUndefined,       // Default value for non templated 1-sample test. Set with SetDistribution
      kUserDefined,     // For internal use only within the class's template constructor
      kGaussian,
      kLogNormal,
      kExponential
   };

   enum EUserDistribution { // User input distribution option
      kCDF,
      kPDF                  // Default value
   };

   enum ETestType { // Goodness of Fit test types for using with the class's unary funtions as a shorthand for the in-built methods
      kAD,   // Anderson-Darling Test. Default value
      kAD2s, // Anderson-Darling 2-Samples Test
      kKS,   // Kolmogorov-Smirnov Test
      kKS2s  // Kolmogorov-Smirnov 2-Samples Test
   };

   /* Constructor for using only with 2-samples tests */
   GoFTest(UInt_t sample1Size, const Double_t* sample1, UInt_t sample2Size, const Double_t* sample2);

   /* Constructor for using only with 1-sample tests with a specified distribution */
   GoFTest(UInt_t sampleSize, const Double_t* sample, EDistribution dist = kUndefined);

   /* Templated constructor for using only with 1-sample tests with a user specified distribution */
   template<class Dist>
   GoFTest(UInt_t sampleSize, const Double_t* sample, Dist& dist, EUserDistribution userDist = kPDF,
           Double_t xmin = 1, Double_t xmax = 0)
   {
      Instantiate(sample, sampleSize);
      SetUserDistribution<Dist>(dist, userDist, xmin, xmax);
   }

   /* Specialization using IGenFunction interface */
   GoFTest(UInt_t sampleSize, const Double_t* sample, const IGenFunction& dist, EUserDistribution userDist = kPDF,
           Double_t xmin = 1, Double_t xmax = 0)
   {
      Instantiate(sample, sampleSize);
      SetUserDistribution(dist, userDist, xmin, xmax);
   }

   /* Sets the user input distribution function for 1-sample tests. */
   template<class Dist>
   void SetUserDistribution(Dist& dist, EUserDistribution userDist = kPDF, Double_t xmin = 1, Double_t xmax = 0) {
      WrappedFunction<Dist&> wdist(dist);
      SetDistributionFunction(wdist, userDist, xmin, xmax);
   }

   /* Template specialization to set the user input distribution for 1-sample tests */
   void SetUserDistribution(const IGenFunction& dist, GoFTest::EUserDistribution userDist = kPDF, Double_t xmin = 1, Double_t xmax = 0) {
      SetDistributionFunction(dist, userDist, xmin, xmax);
   }

   /* Sets the user input distribution as a probability density function for 1-sample tests */
   template<class Dist>
   void SetUserPDF(Dist& pdf, Double_t xmin = 1, Double_t xmax = 0) {
      SetUserDistribution<Dist>(pdf, kPDF, xmin, xmax);
   }

   /* Template specialization to set the user input distribution as a probability density function for 1-sample tests */
   void SetUserPDF(const IGenFunction& pdf, Double_t xmin = 1, Double_t xmax = 0) {
      SetUserDistribution(pdf, kPDF, xmin, xmax);
   }

   /* Sets the user input distribution as a cumulative distribution function for 1-sample tests
      The CDF must return zero
    */
   template<class Dist>
   void SetUserCDF(Dist& cdf, Double_t xmin = 1, Double_t xmax = 0) {
      SetUserDistribution<Dist>(cdf, kCDF, xmin, xmax);
   }

   /* Template specialization to set the user input distribution as a cumulative distribution function for 1-sample tests */
   void SetUserCDF(const IGenFunction& cdf, Double_t xmin = 1, Double_t xmax = 0)  {
      SetUserDistribution(cdf, kCDF, xmin, xmax);
   }


   /* Sets the distribution for the predefined distribution types  */
   void SetDistribution(EDistribution dist);


   virtual ~GoFTest();

/*
  The Anderson-Darling K-Sample Test algorithm is described and taken from
  http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/andeksam.htm
  and described and taken from
   (1) Scholz F.W., Stephens M.A. (1987), K-sample Anderson-Darling Tests, Journal of the American Statistical Association, 82, 918–924. (2-samples variant implemented)
*/ void AndersonDarling2SamplesTest(Double_t& pvalue, Double_t& testStat) const;
   Double_t AndersonDarling2SamplesTest(const Char_t* option = "p") const; // Returns default p-value; option "t" returns the test statistic value "A2"

/*
  The Anderson-Darling 1-Sample Test algorithm for a specific distribution is described at
  http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/andedarl.htm
  and described and taken from (2)
  Marsaglia J.C.W., Marsaglia G. (2004), Evaluating the Anderson-Darling Distribution, Journal of Statistical Software, Volume 09, Issue i02.
  and described and taken from (3)
  Lewis P.A.W. (1961), The Annals of Mathematical Statistics, Distribution of the Anderson-Darling Statistic, Volume 32, Number 4, 1118-1124.
*/ void AndersonDarlingTest(Double_t& pvalue, Double_t& testStat) const;
   Double_t AndersonDarlingTest(const Char_t* option = "p") const; // Returns default p-value; option "t" returns the test statistic value "A2"

/*
  The Kolmogorov-Smirnov 2-Samples Test algorithm is described at
  http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/ks2samp.htm
  and described and taken from
  http://root.cern.ch/root/html/TMath.html#TMath:KolmogorovTest
*/ void KolmogorovSmirnov2SamplesTest(Double_t& pvalue, Double_t& testStat) const;
   Double_t KolmogorovSmirnov2SamplesTest(const Char_t* option = "p") const; // Returns default p-value; option "t" returns the test statistic value "Dn"

/*
  The Kolmogorov-Smirnov 1-Sample Test algorithm for a specific distribution is described at
  http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/kstest.htm
  and described and taken from (4)
  Press W. H., Teukolsky S.A., Vetterling W.T., Flannery B.P. (2007), Numerical Recipes - The Art of Scientific Computing (Third Edition), Cambridge Univerdity Press
*/ void KolmogorovSmirnovTest(Double_t& pvalue, Double_t& testStat) const;
   Double_t KolmogorovSmirnovTest(const Char_t* option = "p") const; // Returns default p-value; option "t" returns the test statistic value "Dn"

   // The class's unary functions
   void operator()(ETestType test, Double_t& pvalue, Double_t& testStat) const;

   // Returns default Anderson Darling 1-Sample Test and default p-value; option "t" returns the test statistic value
   // specific to the test type
   Double_t operator()(ETestType test = kAD, const Char_t* option = "p") const; 
 
   // Computation of the K-Sample Anderson-Darling Test's p-value as described in (1) 
   // given a normalized test statistic. The first variant described in the paper is used 
   static Double_t PValueADKSamples(UInt_t nsamples, Double_t A2 ); 

   // Compute The 2-Sample Anderson Darling test for binned data
   static void  AndersonDarling2SamplesTest(const ROOT::Fit::BinData & data1, const ROOT::Fit::BinData & data2, Double_t& pvalue, Double_t& testStat);

private:

   GoFTest();                       // Disallowed default constructor
   GoFTest(GoFTest& gof);           // Disallowed copy constructor
   GoFTest operator=(GoFTest& gof); // Disallowed assign operator

   std::unique_ptr<IGenFunction> fCDF;


   EDistribution fDist;

   Double_t fMean;
   Double_t fSigma;

   std::vector<Double_t> fCombinedSamples;

   std::vector<std::vector<Double_t> > fSamples;

   Bool_t fTestSampleFromH0;

   void SetCDF();
   void SetDistributionFunction(const IGenFunction& cdf, Bool_t isPDF, Double_t xmin, Double_t xmax);

   void Instantiate(const Double_t* sample, UInt_t sampleSize);


   Double_t LogNormalCDF(Double_t x) const;
   Double_t GaussianCDF(Double_t x) const;
   Double_t ExponentialCDF(Double_t x) const;

   static Double_t GetSigmaN(const std::vector<UInt_t> & ns, UInt_t N); // Computation of sigma_N as described in (1)

   static Double_t InterpolatePValues(int nsamples,Double_t A2); // Linear interpolation used in GoFTest::PValueAD2Samples


   Double_t PValueAD1Sample(Double_t A2) const; // Computation of the 1-Sample Anderson-Darling Test's p-value

   void LogSample(); // Applies the logarithm to the sample when the specified distribution to test is LogNormal

   void SetSamples(std::vector<const Double_t*> samples, const std::vector<UInt_t> samplesSizes);

   void SetParameters(); // Sets the estimated mean and standard-deviation from the samples
}; // end GoFTest class


} // ROOT namespace
} // Math namespace
#endif
