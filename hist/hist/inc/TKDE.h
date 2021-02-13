// @(#)root/hist:$Id$
// Authors: Bartolomeu Rabacal    07/2010
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
// Header file for TKDE

#ifndef ROOT_TKDE
#define ROOT_TKDE

#include "Math/WrappedFunction.h"

#include "TNamed.h"

#include "Math/Math.h"

#include <string>
#include <vector>
#include <memory>

class TGraphErrors;
class TF1;

/*
   Kernel Density Estimation class. The three main references are (1) "Scott DW, Multivariate Density Estimation.
Theory, Practice and Visualization. New York: Wiley", (2) "Jann Ben - ETH Zurich, Switzerland -, Univariate kernel density estimation document for KDENS: Stata module for univariate kernel density estimation." and (3) "Hardle W, Muller M, Sperlich S, Werwatz A, Nonparametric and Semiparametric Models. Springer."
   The algorithm is briefly described in (4) "Cranmer KS, Kernel Estimation in High-Energy
Physics. Computer Physics Communications 136:198-207,2001" - e-Print Archive: hep ex/0011057.
   A binned version is also implemented to address the performance issue due to its data size dependence.
*/
class TKDE : public TNamed  {
public:

   /// Types of Kernel functions
   /// They can be set using the function SetKernelType()
   // or as a string in the constructor
   enum EKernelType { // Kernel function type option
      kGaussian,
      kEpanechnikov,
      kBiweight,
      kCosineArch,
      kUserDefined, // Internal use only for the class's template constructor
      kTotalKernels // Internal use only for member initialization
   };

   /// Iteration types. They can be set using SetIteration()
   enum EIteration {
      kAdaptive,
      kFixed
   };

   /// Data "mirroring" option to address the probability "spill out" boundary effect
   /// They can be set using SetMirror()
   enum EMirror {
      kNoMirror,
      kMirrorLeft,
      kMirrorRight,
      kMirrorBoth,
      kMirrorAsymLeft,
      kMirrorRightAsymLeft,
      kMirrorAsymRight,
      kMirrorLeftAsymRight,
      kMirrorAsymBoth
   };

   /// Data binning option.
   /// They can be set using SetBinning()
   enum EBinning{
      kUnbinned,
      kRelaxedBinning, // The algorithm is allowed to use binning if the data is large enough
      kForcedBinning
   };

   ///  default constructor used only by I/O
   TKDE();

   /// Constructor for unweighted data
   /// Varius option for TKDE can be passed in the option string as below.
   /// Note that min and max will define the plotting range but will not restrict the data in the unbinned case
   /// Instead when use binning, only the data in the range will be considered.
   /// Note also, that when some data exists outside the range, one should not use the mirror option with unbinned.
   /// Adaptive will be soon very slow especially for Nevents > 10000.
   /// For this reason, by default for Nevents >=10000, the data are automatically binned  in
   /// nbins=Min(10000,Nevents/10)
   /// In case of ForceBinning option the default number of bins is 1000
   TKDE(UInt_t events, const Double_t* data, Double_t xMin = 0.0, Double_t xMax = 0.0, const Option_t* option =
                 "KernelType:Gaussian;Iteration:Adaptive;Mirror:noMirror;Binning:RelaxedBinning", Double_t rho = 1.0) {
      Instantiate( nullptr,  events, data, nullptr, xMin, xMax, option, rho);
   }

   /// Constructor for weighted data
   TKDE(UInt_t events, const Double_t* data, const Double_t* dataWeight, Double_t xMin = 0.0, Double_t xMax = 0.0, const Option_t* option =
        "KernelType:Gaussian;Iteration:Adaptive;Mirror:noMirror;Binning:RelaxedBinning", Double_t rho = 1.0) {
      Instantiate( nullptr,  events, data, dataWeight, xMin, xMax, option, rho);
   }

   /// Constructor for unwweighted data and a user defined kernel function
   template<class KernelFunction>
   TKDE(const Char_t* /*name*/, const KernelFunction& kernfunc, UInt_t events, const Double_t* data, Double_t xMin = 0.0, Double_t xMax = 0.0, const Option_t* option = "KernelType:UserDefined;Iteration:Adaptive;Mirror:noMirror;Binning:RelaxedBinning", Double_t rho = 1.0)  {
      Instantiate(new ROOT::Math::WrappedFunction<const KernelFunction&>(kernfunc), events, data, nullptr, xMin, xMax, option, rho);
   }
   /// Constructor for weighted data and a user defined kernel function
   template<class KernelFunction>
   TKDE(const Char_t* /*name*/, const KernelFunction& kernfunc, UInt_t events, const Double_t* data, const Double_t * dataWeight, Double_t xMin = 0.0, Double_t xMax = 0.0, const Option_t* option = "KernelType:UserDefined;Iteration:Adaptive;Mirror:noMirror;Binning:RelaxedBinning", Double_t rho = 1.0)  {
      Instantiate(new ROOT::Math::WrappedFunction<const KernelFunction&>(kernfunc), events, data, dataWeight, xMin, xMax, option, rho);
   }

   virtual ~TKDE();

   void Fill(Double_t data);
   void Fill(Double_t data, Double_t weight);
   void SetKernelType(EKernelType kern);
   void SetIteration(EIteration iter);
   void SetMirror(EMirror mir);
   void SetBinning(EBinning);
   void SetNBins(UInt_t nbins);
   void SetUseBinsNEvents(UInt_t nEvents);
   void SetTuneFactor(Double_t rho);
   void SetRange(Double_t xMin, Double_t xMax); // By default computed from the data

   virtual void Draw(const Option_t* option = "");

   Double_t operator()(Double_t x) const;
   Double_t operator()(const Double_t* x, const Double_t* p=0) const;  // Needed for creating TF1

   Double_t GetValue(Double_t x) const { return (*this)(x); }
   Double_t GetError(Double_t x) const;

   Double_t GetBias(Double_t x) const;
   Double_t GetMean() const;
   Double_t GetSigma() const;
   Double_t GetRAMISE() const;

   Double_t GetFixedWeight() const;

   TF1* GetFunction(UInt_t npx = 100, Double_t xMin = 1.0, Double_t xMax = 0.0);
   TF1* GetUpperFunction(Double_t confidenceLevel = 0.95, UInt_t npx = 100, Double_t xMin = 1.0, Double_t xMax = 0.0);
   TF1* GetLowerFunction(Double_t confidenceLevel = 0.95, UInt_t npx = 100, Double_t xMin = 1.0, Double_t xMax = 0.0);
   TF1* GetApproximateBias(UInt_t npx = 100, Double_t xMin = 1.0, Double_t xMax = 0.0);
   TGraphErrors * GetGraphWithErrors(UInt_t npx = 100, Double_t xMin = 1.0, Double_t xMax = 0.0);

   // get the drawn object to chanage settings
   // These objects are managed by TKDE and should not be deleted by the user
   TF1 * GetDrawnFunction() { return fPDF;}
   TF1 * GetDrawnUpperFunction() { return fUpperPDF;}
   TF1 * GetDrawnLowerFunction() { return fLowerPDF;}
   TGraphErrors * GetDrawnGraph() { return fGraph;}

   const Double_t * GetAdaptiveWeights() const;


public:

   class TKernel {
      TKDE *fKDE;
      UInt_t fNWeights;               // Number of kernel weights (bandwidth as vectorized for binning)
      std::vector<Double_t> fWeights; // Kernel weights (bandwidth)
   public:
      TKernel(Double_t weight, TKDE *kde);
      void ComputeAdaptiveWeights();
      Double_t operator()(Double_t x) const;
      Double_t GetWeight(Double_t x) const;
      Double_t GetFixedWeight() const;
      const std::vector<Double_t> &GetAdaptiveWeights() const;
   };

   friend class TKernel;

private:

   TKDE(TKDE& kde);           // Disallowed copy constructor
   TKDE operator=(TKDE& kde); // Disallowed assign operator

   // Kernel funciton pointer. It is managed by class for internal kernels or exernally for user defined kernels
   typedef ROOT::Math::IBaseFunctionOneDim* KernelFunction_Ptr;
   KernelFunction_Ptr fKernelFunction;  ///<! pointer to kernel function

   std::unique_ptr<TKernel> fKernel;             ///<! internal kernel class. Transient because it is recreated after reading from a file

   std::vector<Double_t> fData;   // Data events
   std::vector<Double_t> fEvents; // Original data storage
   std::vector<Double_t> fEventWeights; // Original data weights

   TF1* fPDF;             //! Output Kernel Density Estimation PDF function
   TF1* fUpperPDF;        //! Output Kernel Density Estimation upper confidence interval PDF function
   TF1* fLowerPDF;        //! Output Kernel Density Estimation lower confidence interval PDF function
   TF1* fApproximateBias; //! Output Kernel Density Estimation approximate bias
   TGraphErrors* fGraph;  //! Graph with the errors

   EKernelType fKernelType;
   EIteration fIteration;
   EMirror fMirror;
   EBinning fBinning;


   Bool_t fUseMirroring, fMirrorLeft, fMirrorRight, fAsymLeft, fAsymRight;
   Bool_t fUseBins;
   Bool_t fNewData;        // flag to control when new data are given
   Bool_t fUseMinMaxFromData; // flag top control if min and max must be used from data

   UInt_t fNBins;          // Number of bins for binned data option
   UInt_t fNEvents;        // Data's number of events
   Double_t fSumOfCounts; // Data sum of weights
   UInt_t fUseBinsNEvents; // If the algorithm is allowed to use automatic (relaxed) binning this is the minimum number of events to do so

   Double_t fMean;  // Data mean
   Double_t fSigma; // Data std deviation
   Double_t fSigmaRob; // Data std deviation (robust estimation)
   Double_t fXMin;  // Data minimum value
   Double_t fXMax;  // Data maximum value
   Double_t fRho;   // Adjustment factor for sigma
   Double_t fAdaptiveBandwidthFactor; // Geometric mean of the kernel density estimation from the data for adaptive iteration

   Double_t fWeightSize; // Caches the weight size

   std::vector<Double_t> fCanonicalBandwidths;
   std::vector<Double_t> fKernelSigmas2;

   std::vector<Double_t> fBinCount; // Number of events per bin for binned data option

   std::vector<Bool_t> fSettedOptions; // User input options flag

   struct KernelIntegrand;
   friend struct KernelIntegrand;

   void Instantiate(KernelFunction_Ptr kernfunc, UInt_t events, const Double_t* data, const Double_t* weight,
                    Double_t xMin, Double_t xMax, const Option_t* option, Double_t rho);

   inline Double_t GaussianKernel(Double_t x) const {
      // Returns the kernel evaluation at x
      Double_t k2_PI_ROOT_INV = 0.398942280401432703; // (2 * M_PI)**-0.5
      return (x > -9. && x < 9.) ? k2_PI_ROOT_INV * std::exp(-.5 * x * x) : 0.0;
   }
   inline Double_t EpanechnikovKernel(Double_t x) const {
      return (x > -1. &&  x < 1.) ? 3. / 4. * (1. - x * x) : 0.0;
   }
   inline Double_t BiweightKernel(Double_t x) const {
      // Returns the kernel evaluation at x
      return (x > -1. &&  x < 1.) ? 15. / 16. * (1. - x * x) * (1. - x * x) : 0.0;
   }
   inline Double_t CosineArchKernel(Double_t x) const {
      // Returns the kernel evaluation at x
      return (x > -1. &&  x < 1.) ? M_PI_4 * std::cos(M_PI_2 * x) : 0.0;
   }
   Double_t UpperConfidenceInterval(const Double_t* x, const Double_t* p) const; // Valid if the bandwidth is small compared to nEvents**1/5
   Double_t LowerConfidenceInterval(const Double_t* x, const Double_t* p) const; // Valid if the bandwidth is small compared to nEvents**1/5
   Double_t ApproximateBias(const Double_t* x, const Double_t* ) const { return GetBias(*x); }
   Double_t ComputeKernelL2Norm() const;
   Double_t ComputeKernelSigma2() const;
   Double_t ComputeKernelMu() const;
   Double_t ComputeKernelIntegral() const;
   Double_t ComputeMidspread() ;
   void ComputeDataStats() ;

   UInt_t Index(Double_t x) const;

   void SetBinCentreData(Double_t xmin, Double_t xmax);
   void SetBinCountData();
   void CheckKernelValidity();
   void SetUserCanonicalBandwidth();
   void SetUserKernelSigma2();
   void SetCanonicalBandwidths();
   void SetKernelSigmas2();
   void SetHistogram();
   void SetUseBins();
   void SetMirror();
   void SetMean();
   void SetSigma(Double_t R);
   void SetKernel();
   void SetKernelFunction(KernelFunction_Ptr kernfunc = 0);
   void SetOptions(const Option_t* option, Double_t rho);
   void CheckOptions(Bool_t isUserDefinedKernel = kFALSE);
   void GetOptions(std::string optionType, std::string option);
   void AssureOptions();
   void SetData(const Double_t* data, const Double_t * weights);
   void ReInit();
   void InitFromNewData();
   void SetMirroredEvents();
   void SetDrawOptions(const Option_t* option, TString& plotOpt, TString& drawOpt);
   void DrawErrors(TString& drawOpt);
   void DrawConfidenceInterval(TString& drawOpt, double cl=0.95);

   TF1* GetKDEFunction(UInt_t npx = 100, Double_t xMin = 1.0, Double_t xMax = 0.0);
   TF1* GetKDEApproximateBias(UInt_t npx = 100, Double_t xMin = 1.0, Double_t xMax = 0.0);
   // The density to estimate should be at least twice differentiable.
   TF1* GetPDFUpperConfidenceInterval(Double_t confidenceLevel = 0.95, UInt_t npx = 100, Double_t xMin = 1.0, Double_t xMax = 0.0);
   TF1* GetPDFLowerConfidenceInterval(Double_t confidenceLevel = 0.95, UInt_t npx = 100, Double_t xMin = 1.0, Double_t xMax = 0.0);

   ClassDef(TKDE, 3) // One dimensional semi-parametric Kernel Density Estimation

};

#endif
