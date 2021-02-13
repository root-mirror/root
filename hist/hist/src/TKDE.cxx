// // @(#)root/hist:$Id$
// Authors: Bartolomeu Rabacal    07/2010
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006 , ROOT MathLib Team                             *
 *                                                                    *
 * For the licensing terms see $ROOTSYS/LICENSE.                      *
 * For the list of contributors see $ROOTSYS/README/CREDITS.          *
 *                                                                    *
 **********************************************************************/
//////////////////////////////////////////////////////////////////////////////
/** \class TKDE
    \ingroup Hist
 Kernel Density Estimation class.
 The three main references are:
 1. "Scott DW, Multivariate Density Estimation. Theory, Practice and Visualization. New York: Wiley",
 2. "Jann Ben - ETH Zurich, Switzerland -, Univariate kernel density estimation document for KDENS:
    Stata module for univariate kernel density estimation."
 3. "Hardle W, Muller M, Sperlich S, Werwatz A, Nonparametric and Semiparametric Models. Springer."
 4. "Cranmer KS, Kernel Estimation in High-Energy
 Physics. Computer Physics Communications 136:198-207,2001" - e-Print Archive: hep ex/0011057.

 The algorithm is briefly described in (4). A binned version is also implemented to address the
 performance issue due to its data size dependance.
 */


#include <functional>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cassert>

#include "Math/Error.h"
#include "TMath.h"
#include "Math/Functor.h"
#include "Math/Integrator.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/RichardsonDerivator.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TVirtualPad.h"
#include "TKDE.h"

// for make_unique
#include "ROOT/RMakeUnique.hxx"

ClassImp(TKDE);


struct TKDE::KernelIntegrand {
   enum EIntegralResult{kNorm, kMu, kSigma2, kUnitIntegration};
   KernelIntegrand(const TKDE* kde, EIntegralResult intRes);
   Double_t operator()(Double_t x) const;
private:
   const TKDE* fKDE;
   EIntegralResult fIntegralResult;
};

///////////////////////////////////////////////////////////////////////
/// default constructor used by I/O
TKDE::TKDE() :
   fKernelFunction(nullptr),
   fKernel(nullptr),
   fPDF(nullptr),
   fUpperPDF(nullptr),
   fLowerPDF(nullptr),
   fApproximateBias(nullptr),
   fGraph(nullptr),
   fUseMirroring(false), fMirrorLeft(false), fMirrorRight(false), fAsymLeft(false), fAsymRight(false),
   fUseBins(false), fNewData(false), fUseMinMaxFromData(false),
   fNBins(0), fNEvents(0), fSumOfCounts(0), fUseBinsNEvents(0),
   fMean(0.),fSigma(0.), fSigmaRob(0.), fXMin(0.), fXMax(0.),
   fRho(0.), fAdaptiveBandwidthFactor(0.), fWeightSize(0)
{
}

TKDE::~TKDE() {
   //Class destructor
   if (fPDF)              delete fPDF;
   if (fUpperPDF)         delete fUpperPDF;
   if (fLowerPDF)         delete fLowerPDF;
   if (fGraph)         delete fGraph;
   if (fApproximateBias)  delete fApproximateBias;
   // delete fKernelFunction if it is not user defined
   if (fKernelFunction && fKernelType != kUserDefined)
      delete fKernelFunction;
}

////////////////////////////////////////////////////////////////////
//// Internal function to instantiate a TKDE object
void TKDE::Instantiate(KernelFunction_Ptr kernfunc, UInt_t events, const Double_t* data, const Double_t* dataWeights, Double_t xMin, Double_t xMax, const Option_t* option, Double_t rho) {

   fPDF = 0;
   fKernelFunction = 0;
   fUpperPDF = 0;
   fLowerPDF = 0;
   fApproximateBias = 0;
   fGraph = 0;
   fNewData = false;
   fUseMirroring = false; fMirrorLeft = false; fMirrorRight = false;
   fAsymLeft = false; fAsymRight = false;
   fNBins = events < 10000 ? 1000 : std::min(10000, int(events / 100)*10);
   fNEvents = events;
   fUseBinsNEvents = 10000;
   fMean = 0.0;
   fSigma = 0.0;
   fXMin = xMin;
   fXMax = xMax;
   fUseMinMaxFromData = (fXMin >= fXMax);
   fSumOfCounts = 0;
   fAdaptiveBandwidthFactor = 1.;
   fRho = rho;
   fWeightSize = 0;
   fCanonicalBandwidths = std::vector<Double_t>(kTotalKernels, 0.0);
   fKernelSigmas2 = std::vector<Double_t>(kTotalKernels, -1.0);
   fSettedOptions = std::vector<Bool_t>(4, kFALSE);
   SetOptions(option, rho);
   CheckOptions(kTRUE);
   SetMirror();
   SetUseBins();
   SetData(data, dataWeights);
   SetKernelFunction(kernfunc);

}

void TKDE::SetOptions(const Option_t* option, Double_t rho) {
   //Sets User defined construction options
   TString opt = option;
   opt.ToLower();
   std::string options = opt.Data();
   size_t numOpt = 4;
   std::vector<std::string> voption(numOpt, "");
   for (std::vector<std::string>::iterator it = voption.begin(); it != voption.end() && !options.empty(); ++it) {
      size_t pos = options.find_last_of(';');
      if (pos == std::string::npos) {
         *it = options;
         break;
      }
      *it = options.substr(pos + 1);
      options = options.substr(0, pos);
   }
   for (std::vector<std::string>::iterator it = voption.begin(); it != voption.end(); ++it) {
      size_t pos = (*it).find(':');
      if (pos != std::string::npos) {
         GetOptions((*it).substr(0, pos) , (*it).substr(pos + 1));
      }
   }
   AssureOptions();
   fRho = rho;
}

void TKDE::SetDrawOptions(const Option_t* option, TString& plotOpt, TString& drawOpt) {
   // Sets User defined drawing options
   size_t numOpt = 2;
   std::string options = TString(option).Data();
   std::vector<std::string> voption(numOpt, "");
   for (std::vector<std::string>::iterator it = voption.begin(); it != voption.end() && !options.empty(); ++it) {
      size_t pos = options.find_last_of(';');
      if (pos == std::string::npos) {
         *it = options;
         break;
      }
      *it = options.substr(pos + 1);
      options = options.substr(0, pos);
   }
   Bool_t foundPlotOPt = kFALSE;
   Bool_t foundDrawOPt = kFALSE;
   for (std::vector<std::string>::iterator it = voption.begin(); it != voption.end() && !options.empty(); ++it) {
      size_t pos = (*it).find(':');
      if (pos == std::string::npos) break;
      TString optionType = (*it).substr(0, pos);
      TString optionInstance = (*it).substr(pos + 1);
      optionType.ToLower();
      optionInstance.ToLower();
      if (optionType.Contains("plot")) {
         foundPlotOPt = kTRUE;
         if (optionInstance.Contains("estimate") || optionInstance.Contains("errors") || optionInstance.Contains("confidenceinterval"))
            plotOpt = optionInstance;
         else {
            this->Warning("SetDrawOptions", "Unknown plotting option %s: setting to KDE estimate plot.",optionInstance.Data());
            this->Info("SetDrawOptions", "Possible plotting options are: Estimate,Errors,ConfidenceInterval");
         }
      } else if (optionType.Contains("drawoptions")) {
         foundDrawOPt = kTRUE;
         drawOpt = optionInstance;
      }
   }
   if (!foundPlotOPt) {
      this->Warning("SetDrawOptions", "No plotting option: setting to KDE estimate plot.");
      plotOpt = "estimate";
   }
   if (!foundDrawOPt) {
      this->Warning("SetDrawOptions", "No drawing options: setting to default ones.");
      drawOpt = "apl4";
   }
}

void TKDE::GetOptions(std::string optionType, std::string option) {
   // Gets User defined KDE construction options
   if (optionType.compare("kerneltype") == 0) {
      fSettedOptions[0] = kTRUE;
      if (option.compare("gaussian") == 0) {
         fKernelType = kGaussian;
      } else if (option.compare("epanechnikov") == 0) {
         fKernelType = kEpanechnikov;
      } else if (option.compare("biweight") == 0) {
         fKernelType = kBiweight;
      } else if (option.compare("cosinearch") == 0) {
         fKernelType = kCosineArch;
      } else if (option.compare("userdefined") == 0) {
         fKernelType = kUserDefined;
      } else {
         this->Warning("GetOptions", "Unknown kernel type option %s: setting to Gaussian",option.c_str());
         this->Info("GetOptions", "Possible kernel type options are: Gaussian, Epanechnikov, Biweight, Cosinearch, Userdefined");
         fKernelType = kGaussian;
      }
   } else if (optionType.compare("iteration") == 0) {
      fSettedOptions[1] = kTRUE;
      if (option.compare("adaptive") == 0) {
         fIteration = kAdaptive;
      } else if (option.compare("fixed") == 0) {
         fIteration = kFixed;
      } else {
         this->Warning("GetOptions", "Unknown iteration option %s: setting to Adaptive",option.c_str());
         this->Info("GetOptions","Possible iteration type options are: Adaptive, Fixed");
         fIteration = kAdaptive;
      }
   } else if (optionType.compare("mirror") == 0) {
      fSettedOptions[2] = kTRUE;
      if (option.compare("nomirror") == 0) {
         fMirror = kNoMirror;
      } else if (option.compare("mirrorleft") == 0) {
         fMirror = kMirrorLeft;
      } else if (option.compare("mirrorright") == 0) {
         fMirror = kMirrorRight;
      } else if (option.compare("mirrorboth") == 0) {
         fMirror = kMirrorBoth;
      } else if (option.compare("mirrorasymleft") == 0) {
         fMirror = kMirrorAsymLeft;
      } else if (option.compare("mirrorrightasymleft") == 0) {
         fMirror = kMirrorRightAsymLeft;
      } else if (option.compare("mirrorasymright") == 0) {
         fMirror = kMirrorAsymRight;
      } else if (option.compare("mirrorleftasymright") == 0) {
         fMirror = kMirrorLeftAsymRight;
      } else if (option.compare("mirrorasymboth") == 0) {
         fMirror = kMirrorAsymBoth;
      } else {
         this->Warning("GetOptions", "Unknown mirror option %s: setting to NoMirror",option.c_str());
         this->Info("GetOptions", "Possible mirror type options are: NoMirror, MirrorLeft, MirrorRight, MirrorAsymLeft,"
                                  "MirrorAsymRight, MirrorRightAsymLeft, MirrorLeftAsymRight, MirrorAsymBoth");
         fMirror = kNoMirror;
      }
   } else if (optionType.compare("binning") == 0) {
      fSettedOptions[3] = kTRUE;
      if (option.compare("unbinned") == 0) {
         fBinning = kUnbinned;
      } else if (option.compare("relaxedbinning") == 0) {
         fBinning = kRelaxedBinning;
      } else if (option.compare("forcedbinning") == 0) {
         fBinning = kForcedBinning;
      } else {
         this->Warning("GetOptions", "Unknown binning option %s: setting to RelaxedBinning", option.c_str());
         this->Info("GetOptions", "Possible binning type options are: Unbinned, ForcedBinning, RelaxedBinning");
         fBinning = kRelaxedBinning;
      }
   }
}

void TKDE::AssureOptions() {
   // Sets missing construction options to default ones
   if (!fSettedOptions[0]) {
      fKernelType = kGaussian;
   }
   if (!fSettedOptions[1]) {
      fIteration = kAdaptive;
   }
   if (!fSettedOptions[2]) {
      fMirror = kNoMirror;
   }
   if (!fSettedOptions[3]) {
      fBinning = kRelaxedBinning;
   }
}

void TKDE::CheckOptions(Bool_t isUserDefinedKernel) {
   // Sets User global options
   if (!(isUserDefinedKernel) && !(fKernelType >= kGaussian && fKernelType < kUserDefined)) {
      Error("CheckOptions", "Illegal user kernel type input! Use template constructor for user defined kernel.");
   }
   if (fIteration != kAdaptive && fIteration != kFixed) {
      Warning("CheckOptions", "Illegal user iteration type input - use default value !");
      fIteration = kAdaptive;
   }
   if (!(fMirror >= kNoMirror && fMirror <= kMirrorAsymBoth)) {
      Warning("CheckOptions", "Illegal user mirroring type input - use default value !");
      fMirror = kNoMirror;
   }
   if (!(fBinning >= kUnbinned && fBinning <= kForcedBinning)) {
      Warning("CheckOptions", "Illegal user binning type input - use default value !");
      fBinning = kRelaxedBinning;
   }
   if (fRho <= 0.0) {
      Warning("CheckOptions", "Tuning factor rho cannot be non-positive - use default value !");
      fRho = 1.0;
   }
}

void TKDE::SetKernelType(EKernelType kern) {
   // Sets User option for the choice of kernel estimator
   if (fKernelFunction && fKernelType != kUserDefined) {
      delete fKernelFunction;
      fKernelFunction = nullptr;
   }
   fKernelType = kern;
   CheckOptions();
   // resetting fKernel automatically calls ReInit when needed
   fKernel.reset();
}

void TKDE::SetIteration(EIteration iter) {
   // Sets User option for fixed or adaptive iteration
   fIteration = iter;
   CheckOptions();
   fKernel.reset();
}

void TKDE::SetMirror(EMirror mir) {
   // Sets User option for mirroring the data
   fMirror = mir;
   CheckOptions();
   SetMirror();
   if (fUseMirroring) {
      SetMirroredEvents();
   }
   fKernel.reset();
}

void TKDE::SetBinning(EBinning bin) {
   // Sets User option for binning the weights
   fBinning = bin;
   CheckOptions();
   SetUseBins();
}

void TKDE::SetNBins(UInt_t nbins) {
   // Sets User option for number of bins
   if (!nbins) {
      Error("SetNBins", "Number of bins must be greater than zero.");
      return;
   }

   fNBins = nbins;

   SetUseBins();
   if (!fUseBins) {
      if (fBinning == kUnbinned)
         Warning("SetNBins", "Bin type using SetBinning must be set for using a binned evaluation");
      else
         Warning("SetNBins", "Bin type using SetBinning or with SetUseBinsNEvents must be set for using a binned evaluation");
   }
}

void TKDE::SetUseBinsNEvents(UInt_t nEvents) {
   // Sets User option for the minimum number of events for allowing automatic binning
   fUseBinsNEvents = nEvents;
   SetUseBins();
}

void TKDE::SetTuneFactor(Double_t rho) {
   // Factor which can be used to tune the smoothing.
   // It is used as multiplicative factor for the fixed and adaptive bandwidth.
   // A value < 1 will reproduce better the tails but oversmooth the peak
   // while a factor > 1 will overestimate the tail
   fRho = rho;
   CheckOptions();
   fKernel.reset();
}

void TKDE::SetRange(Double_t xMin, Double_t xMax) {
   // Sets minimum range value and maximum range value
   if (xMin >= xMax) {
      Error("SetRange", "Minimum range cannot be bigger or equal than the maximum range! Present range values remain the same.");
      return;
   }
   fXMin = xMin;
   fXMax = xMax;
   fUseMinMaxFromData = false;
   fKernel.reset();
}

// private methods

void TKDE::SetUseBins() {
   // Sets User option for using binned weights
   switch (fBinning) {
      default:
      case kRelaxedBinning:
         if (fNEvents >= fUseBinsNEvents) {
            fUseBins = kTRUE;
         } else {
            fUseBins = kFALSE;
         }
         break;
      case kForcedBinning:
         fUseBins = kTRUE;
         break;
      case kUnbinned:
         fUseBins = kFALSE;
   }

   // during initialization we don't need to recompute the bins
   // it is done within TKDE::SetData
   // in this case fEvents is empty
   if (fEvents.empty()) return;

   // in case we need to recompute bins structure
   // needed when called from the TKDE setters method
   if (fUseBins) {
      fWeightSize = fNBins / (fXMax - fXMin);
      SetBinCentreData(fXMin, fXMax);
      SetBinCountData();
   }
   else {
      // unbinned case
      fWeightSize = 0.;
      fBinCount.clear();
      fData = fEvents;
   }
   fKernel.reset();
}

void TKDE::SetMirror() {
   // Sets the mirroring
   fMirrorLeft   = fMirror == kMirrorLeft      || fMirror == kMirrorBoth          || fMirror == kMirrorLeftAsymRight;
   fMirrorRight  = fMirror == kMirrorRight     || fMirror == kMirrorBoth          || fMirror == kMirrorRightAsymLeft;
   fAsymLeft     = fMirror == kMirrorAsymLeft  || fMirror == kMirrorRightAsymLeft || fMirror == kMirrorAsymBoth;
   fAsymRight    = fMirror == kMirrorAsymRight || fMirror == kMirrorLeftAsymRight || fMirror == kMirrorAsymBoth;
   fUseMirroring = fMirrorLeft                 || fMirrorRight ;
}

void TKDE::SetData(const Double_t* data, const Double_t* wgts) {
   // Sets the data events input sample or bin centres for binned option and computes basic estimators
   if (!data) {
      if (fNEvents) fData.reserve(fNEvents);
      return;
   }
   fEvents.assign(data, data + fNEvents);
   if (wgts) fEventWeights.assign(wgts, wgts + fNEvents);

   if (fUseMinMaxFromData) {
      fXMin = *std::min_element(fEvents.begin(), fEvents.end());
      fXMax = *std::max_element(fEvents.begin(), fEvents.end());
   }

   if (fUseBins) {
      if (fNBins >= fNEvents) {
         this->Warning("SetData", "Default number of bins is greater or equal to number of events. Use SetNBins(UInt_t) to set the appropriate number of bins");
      }
      fWeightSize = fNBins / (fXMax - fXMin);
      SetBinCentreData(fXMin, fXMax);
   } else {
      fWeightSize = 0;  //fNEvents / (fXMax - fXMin);
      fData = fEvents;
   }

   ComputeDataStats();
   if (fUseMirroring) {
      // set mirror events called automatically SetBinCount
      SetMirroredEvents();
   }
   else {
      // to set fBinCount and fSumOfCounts
      SetBinCountData();
   }
}

void TKDE::ReInit() {
   // re-initialize the KDE in case of new data or in case of reading from a file

   // in case of new data re-compute the statistics (works only for unbinned case)
   if (fNewData)  {
      InitFromNewData();
      return;
   }

   if (fEvents.size() == 0) {
      Error("ReInit","TKDE does not contain any data !");
      return;
   }

   //  when of reading from a file fKernelFunction is a nullptr
   // we need to recreate Kernel class (with adaptive weights if needed) and
   // recreate kernel function pointer
   if (!fKernelFunction) SetKernelFunction(nullptr);

   SetKernel();
}

void TKDE::InitFromNewData() {
   // re-initialize when new data have been filled in TKDE
   // re-compute kernel quantities and statistics
   if (fUseBins) {
      Error("InitFromNewData","Re-felling is not supported with binning");
      return;
   }

   fNewData = false;
   // in case of unbinned data
   if (!fUseBins) fEvents = fData;
   if (fUseMinMaxFromData) {
      fXMin = *std::min_element(fEvents.begin(), fEvents.end());
      fXMax = *std::max_element(fEvents.begin(), fEvents.end());
   }
   ComputeDataStats();
   //
   fWeightSize = fNEvents / (fXMax - fXMin);
   if (fUseMirroring) {
      SetMirroredEvents();
   }
   // in case of I/O reset the kernel
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// Intgernal function to mirror the data
void TKDE::SetMirroredEvents() {
   assert(fUseMirroring);
   if (fUseBins) {
      // binned case
      // first fill the bin data in [fXmin,fXmax]
      // bin center should have been already filled before
      SetBinCentreData(fXMin, fXMax);
      SetBinCountData();

      // now mirror the bins
      assert(fNBins == fData.size());
      if (fMirrorLeft) {
         fData.insert(fData.begin(), fNBins, 0.);
         fBinCount.insert(fBinCount.begin(), fNBins, 0.);
         for (UInt_t i = 0; i < fNBins; ++i) {
            fData[i] = fData[i + fNBins] - (fXMax - fXMin);
            fBinCount[fNBins - i - 1] = fBinCount[i + fNBins];
         }
      }
      if (fMirrorRight) {
         fData.insert(fData.end(), fNBins, 0.);
         fBinCount.insert(fBinCount.end(), fNBins, 0.);
         int i0 = (fMirrorLeft) ? fNBins : 0;
         int iLast = i0 + 2 * fNBins - 1;
         for (UInt_t i = 0; i < fNBins; ++i) {
            fData[i0 + fNBins + i] = fData[i0 + i] + (fXMax - fXMin);
            fBinCount[iLast - i] = fBinCount[i0 + i];
         }
      }
   }
   else {
      // unbinned case (transform events)
      std::vector<Double_t> originalEvents = fEvents;
      std::vector<Double_t> originalWeights = fEventWeights;
      if (fMirrorLeft) {
         fEvents.resize(2 * fNEvents, 0.0);
         transform(fEvents.begin(), fEvents.begin() + fNEvents, fEvents.begin() + fNEvents,
                   std::bind(std::minus<Double_t>(), 2 * fXMin, std::placeholders::_1));
      }
      if (fMirrorRight) {
         fEvents.resize((fMirrorLeft + 2) * fNEvents, 0.0);
         transform(fEvents.begin(), fEvents.begin() + fNEvents, fEvents.begin() + (fMirrorLeft + 1) * fNEvents,
                   std::bind(std::minus<Double_t>(), 2 * fXMax, std::placeholders::_1));
      }
      if (!fEventWeights.empty() && (fMirrorLeft || fMirrorRight)) {
         // copy weights too
         fEventWeights.insert(fEventWeights.end(), fEventWeights.begin(), fEventWeights.begin() + fNEvents);
         if (fMirrorLeft && fMirrorRight)
            fEventWeights.insert(fEventWeights.end(), fEventWeights.begin(), fEventWeights.begin() + fNEvents);
      }

      fData = fEvents;
      SetBinCountData();
      fEvents = originalEvents;
      fEventWeights = originalWeights;
   }

}

void TKDE::SetMean() {
   // Computes input data's mean
   fMean = std::accumulate(fEvents.begin(), fEvents.end(), 0.0) / fEvents.size();
}

void TKDE::SetSigma(Double_t R) {
   // Computes input data's sigma
   fSigma = std::sqrt(1. / (fEvents.size() - 1.) * (std::inner_product(fEvents.begin(), fEvents.end(), fEvents.begin(), 0.0) - fEvents.size() * std::pow(fMean, 2.)));
   fSigmaRob = std::min(fSigma, R / 1.349); // Sigma's robust estimator
}

void TKDE::SetKernel() {
   // Sets the kernel density estimator
   // n here should be fNEvents or the total size including the mirrored events ???
   // it should not be fData.size() that in binned case is number of bins
   UInt_t n =  (fUseBins) ? fNBins : fNEvents;
   if (n == 0) return;
   // Optimal bandwidth (Silverman's rule of thumb with assumed Gaussian density)
   Double_t weight = fCanonicalBandwidths[kGaussian] * fSigmaRob * std::pow(3. / (8. * std::sqrt(M_PI)) * n, -0.2);
   weight *= fRho * fCanonicalBandwidths[fKernelType] / fCanonicalBandwidths[kGaussian];

   fKernel = std::make_unique<TKernel>(weight, this);

   if (fIteration == kAdaptive) {
      fKernel->ComputeAdaptiveWeights();
   }
   if (gDebug) {
      if (fIteration != kAdaptive)
         Info("SetKernel",
              "Using a fix kernel - bandwidth = %f - using n = %d, rho = %f , sigmaRob = %f , mean = %f , sigma = %f "
              ", canonicalBandwidth= %f",
              weight, n, fRho, fSigmaRob, fMean, fSigma, fCanonicalBandwidths[fKernelType]);
      else
         Info("SetKernel",
              "Using an adaptive kernel - weight = %f - using n = %d, rho = %f , sigmaRob = %f , mean = %f , sigma = %f "
              ", canonicalBandwidth= %f",
              weight, n, fRho, fSigmaRob, fMean, fSigma, fCanonicalBandwidths[fKernelType]);
   }
}

void TKDE::SetKernelFunction(KernelFunction_Ptr kernfunc) {
   // sets the kernel function pointer. It creates and manages the pointer for internal classes
   if (fKernelFunction) {
      // kernel function pointer must be set to null before calling SetKernelFunction to avoid memory leaks
      Fatal("SetKernelFunction", "Kernel function pointer is not null");
   }
   switch (fKernelType) {
      case kGaussian :
         fKernelFunction = new ROOT::Math::WrappedMemFunction<TKDE, Double_t (TKDE::*)(Double_t) const>(*this, &TKDE::GaussianKernel);
         break;
      case kEpanechnikov :
         fKernelFunction = new ROOT::Math::WrappedMemFunction<TKDE, Double_t (TKDE::*)(Double_t) const>(*this, &TKDE::EpanechnikovKernel);
         break;
      case kBiweight :
         fKernelFunction = new ROOT::Math::WrappedMemFunction<TKDE, Double_t (TKDE::*)(Double_t) const>(*this, &TKDE::BiweightKernel);
         break;
      case kCosineArch :
         fKernelFunction = new ROOT::Math::WrappedMemFunction<TKDE, Double_t (TKDE::*)(Double_t) const>(*this, &TKDE::CosineArchKernel);
         break;
      case kUserDefined :
         fKernelFunction = kernfunc;
         if (fKernelFunction)  CheckKernelValidity();
         break;
      case kTotalKernels :
      default:
         /// for user defined kernels
         fKernelFunction = kernfunc;
         fKernelType = kUserDefined;
   }

   if (fKernelType == kUserDefined) {
      if (fKernelFunction) {
         CheckKernelValidity();
         SetUserCanonicalBandwidth();
         SetUserKernelSigma2();
      }
      else {
         Error("SetKernelFunction", "User kernel function is not defined !");
         return;
      }
   }
   assert(fKernelFunction);
   SetKernelSigmas2();
   SetCanonicalBandwidths();
   SetKernel();
}

void TKDE::SetCanonicalBandwidths() {
   // Sets the canonical bandwidths according to the kernel type
   fCanonicalBandwidths[kGaussian] = 0.7764;     // Checked in Mathematica
   fCanonicalBandwidths[kEpanechnikov] = 1.7188; // Checked in Mathematica
   fCanonicalBandwidths[kBiweight] = 2.03617;    // Checked in Mathematica
   fCanonicalBandwidths[kCosineArch] = 1.7663;   // Checked in Mathematica
   fCanonicalBandwidths[kUserDefined] = 1.0;     // To be Checked
}

void TKDE::SetKernelSigmas2() {
   // Sets the kernel sigmas2 according to the kernel type
   fKernelSigmas2[kGaussian] = 1.0;
   fKernelSigmas2[kEpanechnikov] = 1.0 / 5.0;
   fKernelSigmas2[kBiweight] = 1.0 / 7.0;
   fKernelSigmas2[kCosineArch] = 1.0 - 8.0 / std::pow(M_PI, 2);
}

TF1* TKDE::GetFunction(UInt_t npx, Double_t xMin, Double_t xMax) {
   // Returns the PDF estimate as a function sampled in npx between xMin and xMax
   // the KDE is not re-normalized to the xMin/xMax range.
   // The user manages the returned function
   // For getting a non-sampled TF1, one can create a TF1 directly from the TKDE by doing
   // TF1 * f1  = new TF1("f1",kde,xMin,xMax,0);
   return GetKDEFunction(npx,xMin,xMax);
}

TF1* TKDE::GetUpperFunction(Double_t confidenceLevel, UInt_t npx, Double_t xMin, Double_t xMax) {
   // Returns the PDF upper estimate (upper confidence interval limit)
   return GetPDFUpperConfidenceInterval(confidenceLevel,npx,xMin,xMax);
}

TF1* TKDE::GetLowerFunction(Double_t confidenceLevel, UInt_t npx, Double_t xMin, Double_t xMax) {
   // Returns the PDF lower estimate (lower confidence interval limit)
   return GetPDFLowerConfidenceInterval(confidenceLevel,npx,xMin,xMax);
}

TF1* TKDE::GetApproximateBias(UInt_t npx, Double_t xMin, Double_t xMax) {
   // Returns the PDF estimate bias
   return GetKDEApproximateBias(npx,xMin,xMax);
}

void TKDE::Fill(Double_t data) {
   // Fills data member with User input data event for the unbinned option
   if (fUseBins) {
      this->Warning("Fill", "Cannot fill data with data binned option. Data input ignored.");
      return;
   }
   fData.push_back(data);
   fNEvents++;
   fNewData = kTRUE;
}

void TKDE::Fill(Double_t data, Double_t weight) {
   // Fills data member with User input data event for the unbinned option
   if (fUseBins) {
      this->Warning("Fill", "Cannot fill data with data binned option. Data input ignored.");
      return;
   }
   fData.push_back(data);  // should not be here fEvent ??
   fEventWeights.push_back(weight);
   fNEvents++;
   fNewData = kTRUE;
}

Double_t TKDE::operator()(const Double_t* x, const Double_t*) const {
   // The class's unary function: returns the kernel density estimate
   return (*this)(*x);
}

Double_t TKDE::operator()(Double_t x) const {
   // The class's unary function: returns the kernel density estimate
   if (!fKernel) {
      (const_cast<TKDE*>(this))->ReInit();
      // in case of failed re-initialization
      if (!fKernel) return TMath::QuietNaN();
   }
   return (*fKernel)(x);
}

Double_t TKDE::GetMean() const {
   // return the mean of the data
   if (fNewData) (const_cast<TKDE*>(this))->InitFromNewData();
   return fMean;
}

Double_t TKDE::GetSigma() const {
   // return the standard deviation  of the data
   if (fNewData) (const_cast<TKDE*>(this))->InitFromNewData();
   return fSigma;
}

Double_t TKDE::GetRAMISE() const {
   // Returns the Root Asymptotic Mean Integrated Squared Error according to Silverman's rule of thumb with assumed Gaussian density
   Double_t result = 5. / 4. * fKernelSigmas2[fKernelType] * std::pow(fCanonicalBandwidths[fKernelType], 4) * std::pow(3. / (8. * std::sqrt(M_PI)) , -0.2) * fSigmaRob * std::pow(fNEvents, -0.8);
   return std::sqrt(result);
}

TKDE::TKernel::TKernel(Double_t weight, TKDE* kde) :
// Internal class constructor
fKDE(kde),
fNWeights(kde->fData.size()),
fWeights(1, weight)
{}

void TKDE::TKernel::ComputeAdaptiveWeights() {
   // Gets the adaptive weights (bandwidths) for TKernel internal computation
   unsigned int n = fKDE->fData.size();
   Double_t minWeight = fWeights[0] * 0.05;
   // we will store computed adaptive weights in weights
   std::vector<Double_t> weights(n, fWeights[0]);
   bool useDataWeights = (fKDE->fBinCount.size() == n);
   Double_t f = 0.0;
   for (unsigned int i = 0; i < n; ++i) {
      // for negative or null bin contents use the fixed weight value (fWeights[0])
      if (useDataWeights && fKDE->fBinCount[i] <= 0) {
         weights[i] = fWeights[0];
         continue; // skip negative or null weights
      }
      f = (*fKDE->fKernel)(fKDE->fData[i]);
      if (f <= 0) {
         // this can happen when data are outside range and fAsymLeft or fAsymRight is on
         fKDE->Warning("ComputeAdativeWeights","function value is zero or negative for x = %f w = %f - set their bandwidth to zero",
                       fKDE->fData[i],(useDataWeights) ? fKDE->fBinCount[i] : 1.);
         // set bandwidth for these points to zero.
         weights[i] = 0;
         continue;
      }
      // use weight if it is larger
      weights[i] = std::max(weights[i] /= std::sqrt(f), minWeight);
      fKDE->fAdaptiveBandwidthFactor += std::log(f);
      // printf("(f = %f w = %f af = %f ),",f,*weight,fKDE->fAdaptiveBandwidthFactor);
   }
   Double_t kAPPROX_GEO_MEAN = 0.241970724519143365; // 1 / TMath::Power(2 * TMath::Pi(), .5) * TMath::Exp(-.5). Approximated geometric mean over pointwise data (the KDE function is substituted by the "real Gaussian" pdf) and proportional to sigma. Used directly when the mirroring is enabled, otherwise computed from the data
   // not sure for this special case for mirror. This results in a much smaller bandwidth for mirror case
   fKDE->fAdaptiveBandwidthFactor = fKDE->fUseMirroring ? kAPPROX_GEO_MEAN / fKDE->fSigmaRob : std::sqrt(std::exp(fKDE->fAdaptiveBandwidthFactor / fKDE->fData.size()));
   // set adaptive weights in fWeights matrix
   fWeights.resize(n);
   transform(weights.begin(), weights.end(), fWeights.begin(),
             std::bind(std::multiplies<Double_t>(), std::placeholders::_1, fKDE->fAdaptiveBandwidthFactor));
   //printf("adaptive bandwidth factor % f weight 0 %f , %f \n",fKDE->fAdaptiveBandwidthFactor, weights[0],fWeights[0] );
}

Double_t TKDE::TKernel::GetWeight(Double_t x) const {
   // Returns the bandwidth
   return fWeights[fKDE->Index(x)];
}

void TKDE::SetBinCentreData(Double_t xmin, Double_t xmax) {
   // Returns the bins' centres from the data for using with the binned option
   fData.assign(fNBins, 0.0);
   Double_t binWidth((xmax - xmin) / fNBins);
   for (UInt_t i = 0; i < fNBins; ++i) {
      fData[i] = xmin + (i + 0.5) * binWidth;
   }
}

void TKDE::SetBinCountData() {
   // Returns the bins' count from the data for using with the binned option
   // or set the bin count to the weights in case of weighted data
   if (fUseBins) {
      fBinCount.assign(fNBins, 0.0);
      fSumOfCounts = 0;
      // note this function should be called before the data have been mirrored
      UInt_t nevents = fNEvents;
      assert(fEvents.size() == nevents);
      // case of weighted events
      if (!fEventWeights.empty() ) {
         assert(nevents == fEventWeights.size());
         for (UInt_t i = 0; i < nevents; ++i) {
            if (fEvents[i] >= fXMin && fEvents[i] < fXMax) {
               fBinCount[Index(fEvents[i])] += fEventWeights[i];
               fSumOfCounts += fEventWeights[i];
               //printf("sum of counts %f - bin count %d - %f \n",fSumOfCounts, Index(fEvents[i]), fBinCount[Index(fEvents[i])] );
            }
         }
      }
      // case of unweighted data
      else {
         for (UInt_t i = 0; i < nevents; ++i) {
            if (fEvents[i] >= fXMin && fEvents[i] < fXMax) {
               fBinCount[Index(fEvents[i])] += 1;
               fSumOfCounts += 1;
            }
         }
      }
   }
   else if (!fEventWeights.empty() ) {
      //unbinned weighted data
      fBinCount = fEventWeights;
      fSumOfCounts = 0;
      for (UInt_t i = 0; i < fNEvents; ++i) {
         if (fEvents[i] >= fXMin && fEvents[i] < fXMax) {
            fSumOfCounts += fEventWeights[i];
         }
      }
   }
   else {
      // unbinned unweighted data
      //fSumOfCounts = fNEvents;
      fSumOfCounts = 0;
      for (UInt_t i = 0; i < fNEvents; ++i) {
         if (fEvents[i] >= fXMin && fEvents[i] < fXMax) {
            fSumOfCounts += 1;
         }
      }
      fBinCount.clear();
   }
}

////////////////////////////////////////////////////////////////////////////////
///  Draws either the KDE functions or its errors
//    @param opt  : Drawing options:
//                  -  ""  (default) - draw just the kde
//                  - "same" draw on top of existing pad
//                  - "Errors" draw a TGraphErrors with the point and errors
//                  -"confidenceinterval" draw KDE + conf interval functions (default is 95%)
//                  -"confidenceinterval@0.90" draw KDE + conf interval functions at 90%
//                  - Extra options can be passed in opt for the drawing of the corresponding TF1 or TGraph
//
// NOTE:  The functions GetDrawnFunction(), GetDrawnUpperFunction(), GetDrawnLowerFunction()
//  and GetGraphWithErrors() return the corresponding drawn objects (which are maneged by the TKDE)
// They can be used to changes style, color, etc..
////
void TKDE::Draw(const Option_t* opt) {

   TString plotOpt = opt;
   plotOpt.ToLower();
   TString drawOpt = plotOpt;
   if(gPad && !plotOpt.Contains("same")) {
      gPad->Clear();
   }
   if (plotOpt.Contains("errors"))  {
      drawOpt.ReplaceAll("errors","");
      DrawErrors(drawOpt);
   }
   else if (plotOpt.Contains("confidenceinterval") ||
            plotOpt.Contains("confinterval")) {
      // parse level option
      drawOpt.ReplaceAll("confidenceinterval","");
      drawOpt.ReplaceAll("confinterval","");
      Double_t level = 0.95;
      const char * s = strstr(plotOpt.Data(),"interval@");
      // coverity [secure_coding : FALSE]
      if (s != 0) sscanf(s,"interval@%lf",&level);
      if((level <= 0) || (level >= 1)) {
         Warning("Draw","given confidence level %.3lf is invalid - use default 0.95",level);
         level = 0.95;
      }
      DrawConfidenceInterval(drawOpt,level);
   }
   else {
      if (fPDF) delete fPDF;
      fPDF = GetKDEFunction();
      fPDF->Draw(drawOpt);
   }
}

///////////////////////////////////////////////////////////////////////
///  Draws a TGraphErrors wih KDE values and errors
void TKDE::DrawErrors(TString& drawOpt)
{
   if (fGraph) delete fGraph;
   fGraph = GetGraphWithErrors();
   fGraph->Draw(drawOpt.Data());
}
///////////////////////////////////////////////////////////////////////
/// return a TGraphErrors with the KDE values and errors
/// The return object is managed by the user
TGraphErrors* TKDE::GetGraphWithErrors(UInt_t npx, double xmin, double xmax) {
   if (xmin>= xmax) { xmin = fXMin; xmax = fXMax; }
   UInt_t n = npx;
   Double_t* x = new Double_t[n + 1];
   Double_t* ex = new Double_t[n + 1];
   Double_t* y = new Double_t[n + 1];
   Double_t* ey = new Double_t[n + 1];
   for (UInt_t i = 0; i <= n; ++i) {
      x[i] = xmin + i * (xmax - xmin) / n;
      y[i] = (*this)(x[i]);
      ex[i] = 0;
      ey[i] = this->GetError(x[i]);
   }
   TGraphErrors* ge = new TGraphErrors(n, &x[0], &y[0], &ex[0], &ey[0]);
   ge->SetName("kde_graph_error");
   ge->SetTitle("Errors");
   delete [] x;
   delete [] ex;
   delete [] y;
   delete [] ey;
   return ge;
}

//////////////////////////////////////////////////////////////
/// // Draws the KDE and its confidence interval
void TKDE::DrawConfidenceInterval(TString& drawOpt,double cl) {
   GetKDEFunction()->Draw(drawOpt.Data());
   TF1* upper = GetPDFUpperConfidenceInterval(cl);
   upper->SetLineColor(kBlue);
   upper->Draw(("same" + drawOpt).Data());
   TF1* lower = GetPDFLowerConfidenceInterval(cl);
   lower->SetLineColor(kRed);
   lower->Draw(("same" + drawOpt).Data());
   if (fUpperPDF) delete fUpperPDF;
   if (fLowerPDF) delete fLowerPDF;
   fUpperPDF = upper;
   fLowerPDF = lower;
}

Double_t TKDE::GetFixedWeight() const {
   // Returns the bandwidth for the non adaptive KDE
   Double_t result = -1.;
   if (fIteration == TKDE::kAdaptive) {
      this->Warning("GetFixedWeight()", "Fixed iteration option not enabled. Returning %f.", result);
   } else {
      result = fKernel->GetFixedWeight();
   }
   return result;
}

const Double_t *  TKDE::GetAdaptiveWeights() const {
   // Returns the bandwidths for the adaptive KDE
   if (fIteration != TKDE::kAdaptive) {
      this->Warning("GetFixedWeight()", "Adaptive iteration option not enabled. Returning a NULL pointer<");
      return 0;
   }
   if (fNewData) (const_cast<TKDE*>(this))->InitFromNewData();
   return &(fKernel->GetAdaptiveWeights()).front();
}

Double_t TKDE::TKernel::GetFixedWeight() const {
   // Returns the bandwidth for the non adaptive KDE
   return fWeights[0];
}

const std::vector<Double_t> & TKDE::TKernel::GetAdaptiveWeights() const {
   // Returns the bandwidth for the adaptive KDE
   return fWeights;
}

Double_t TKDE::TKernel::operator()(Double_t x) const {
   // The internal class's unary function: returns the kernel density estimate
   Double_t result(0.0);
   UInt_t n = fKDE->fData.size();
   // case of bins or weighted data
   Bool_t useCount = (fKDE->fBinCount.size() == n);
   // also in case of unbinned unweighted data fSumOfCounts is sum of events in range
   // events outside range should be used to normalize the TKDE ??
   Double_t nSum = fKDE->fSumOfCounts; //(useBins) ? fKDE->fSumOfCounts : fKDE->fNEvents;
   //if (!useCount) nSum = fKDE->fNEvents;
   // in case of non-adaptive fWeights is a vector of size 1
   Bool_t hasAdaptiveWeights = (fWeights.size() == n);
   Double_t invWeight = (!hasAdaptiveWeights) ? 1. / fWeights[0] : 0;
   for (UInt_t i = 0; i < n; ++i) {
      Double_t binCount = (useCount) ? fKDE->fBinCount[i] : 1.0;
      // uncommenting following line slows down so keep computation for
      // zero bincounts
      //if (binCount <= 0) continue;
      if (hasAdaptiveWeights) {
         // skip data points that have 0 bandwidth (this can happen, see TKernel::ComputeAdaptiveWeight)
         if (fWeights[i] == 0) continue;
         invWeight = 1. / fWeights[i];
      }
      result += binCount * invWeight * (*fKDE->fKernelFunction)((x - fKDE->fData[i]) * invWeight );
      if (fKDE->fAsymLeft) {
         result += binCount * invWeight * (*fKDE->fKernelFunction)((x - (2. * fKDE->fXMin - fKDE->fData[i])) * invWeight);
      }
      if (fKDE->fAsymRight) {
         result += binCount * invWeight * (*fKDE->fKernelFunction)((x - (2. * fKDE->fXMax - fKDE->fData[i])) * invWeight);
      }
      // printf("data point %i  %f  %f  count %f weight % f result % f\n",i,fKDE->fData[i],fKDE->fEvents[i],binCount,fWeights[i], result);
   }
   if ( TMath::IsNaN(result) ) {
      fKDE->Warning("operator()","Result is NaN for  x %f \n",x);
   }
   return result / nSum;
}

////////////////////////////////////////////////////
/// compute the bin index given a data point x
UInt_t TKDE::Index(Double_t x) const {
   // Returns the indices (bins) for the data point
   // fWeightSize is the inverse of the bin width
   Int_t bin = Int_t((x - fXMin) * fWeightSize);
   if (bin == (Int_t)fData.size()) return --bin;
   if (bin > (Int_t)fData.size()) {
      return (Int_t)(fData.size()) - 1;
   } else if (bin <= 0) {
      return 0;
   }
   return bin;
}

Double_t TKDE::UpperConfidenceInterval(const Double_t* x, const Double_t* p) const {
   // Returns the pointwise upper estimated density
   Double_t f = (*this)(x);
   Double_t sigma = GetError(*x);
   Double_t prob = 1. - (1.-*p)/2;   // this is 1.-alpha/2
   Double_t z = ROOT::Math::normal_quantile(prob, 1.0);
   return f + z * sigma;
}

Double_t TKDE::LowerConfidenceInterval(const Double_t* x, const Double_t* p) const {
   // Returns the pointwise lower estimated density
   Double_t f = (*this)(x);
   Double_t sigma = GetError(*x);
   Double_t prob = (1.-*p)/2;    // this is alpha/2
   Double_t z = ROOT::Math::normal_quantile(prob, 1.0);
   return f + z * sigma;
}


Double_t TKDE::GetBias(Double_t x) const {
   // Returns the pointwise approximate estimated density bias
   ROOT::Math::WrappedFunction<const TKDE&> kern(*this);
   ROOT::Math::RichardsonDerivator rd;
   rd.SetFunction(kern);
   Double_t df2 = rd.Derivative2(x);
   Double_t weight = fKernel->GetWeight(x); // Bandwidth
   return  0.5 * fKernelSigmas2[fKernelType] * std::pow(weight, 2) * df2;
}
Double_t TKDE::GetError(Double_t x) const {
   // Returns the pointwise sigma of estimated density
   Double_t kernelL2Norm = ComputeKernelL2Norm();
   Double_t f = (*this)(x);
   Double_t weight = fKernel->GetWeight(x); // Bandwidth
   Double_t result = f * kernelL2Norm / (fNEvents * weight);
   return std::sqrt(result);
}

void TKDE::CheckKernelValidity() {
   // Checks if kernel has unit integral, mu = 0 and positive finite sigma conditions
   Double_t valid = kTRUE;
   Double_t unity = ComputeKernelIntegral();
   valid = valid && unity == 1.;
   if (!valid) {
      Error("CheckKernelValidity", "Kernel's integral is %f",unity);
   }
   Double_t mu = ComputeKernelMu();
   valid = valid && mu == 0.;
   if (!valid) {
      Error("CheckKernelValidity", "Kernel's mu is %f" ,mu);
   }
   Double_t sigma2 = ComputeKernelSigma2();
   valid = valid && sigma2 > 0 && sigma2 != std::numeric_limits<Double_t>::infinity();
   if (!valid) {
      Error("CheckKernelValidity", "Kernel's sigma2 is %f",sigma2);
   }
   if (!valid) {
      Error("CheckKernelValidity", "Validation conditions: the kernel's integral must be 1, the kernel's mu must be zero and the kernel's sigma2 must be finite positive to be a suitable kernel.");
      //exit(EXIT_FAILURE);
   }
}

Double_t TKDE::ComputeKernelL2Norm() const {
   // Computes the kernel's L2 norm
   ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kGAUSS);
   KernelIntegrand kernel(this, TKDE::KernelIntegrand::kNorm);
   ig.SetFunction(kernel);
   Double_t result = ig.Integral();
   return result;
}

Double_t TKDE::ComputeKernelSigma2() const {
   // Computes the kernel's sigma squared
   ROOT::Math::IntegratorOneDim ig( ROOT::Math::IntegrationOneDim::kGAUSS);
   KernelIntegrand kernel(this, TKDE::KernelIntegrand::kSigma2);
   ig.SetFunction(kernel);
   Double_t result = ig.Integral();
   return result;
}

Double_t TKDE::ComputeKernelMu() const {
   // Computes the kernel's mu
   ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kGAUSS);
   KernelIntegrand kernel(this, TKDE::KernelIntegrand::kMu);
   ig.SetFunction(kernel);
   Double_t result = ig.Integral();
   return result;
}

Double_t TKDE::ComputeKernelIntegral() const {
   // Computes the kernel's integral which ought to be unity
   ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kGAUSS);
   KernelIntegrand kernel(this, TKDE::KernelIntegrand::kUnitIntegration);
   ig.SetFunction(kernel);
   Double_t result = ig.Integral();
   return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Internal funciton to compute statistics (mean,stddev) using always all the provided data (i.e. no binning)
void TKDE::ComputeDataStats() {
   /// in case of weights use
   if (!fEventWeights.empty() ) {
      // weighted data
      double x1 = fXMin - 0.001*(fXMax-fXMin);
      double x2 = fXMax + 0.001*(fXMax-fXMin);
      TH1D h1("temphist","", 500, x1, x2);
      h1.FillN(fEvents.size(), fEvents.data(), fEventWeights.data() );
      assert (h1.GetSumOfWeights() > 0) ;
      fMean = h1.GetMean();
      fSigma = h1.GetRMS();
      // compute robust sigma using midspread
      Double_t quantiles[2] = {0.0, 0.0};
      Double_t prob[2] = {0.25, 0.75};
      h1.GetQuantiles(2, quantiles, prob);
      Double_t midspread = quantiles[1] - quantiles[0];
      fSigmaRob = std::min(fSigma, midspread / 1.349); // Sigma's robust estimator
      //printf("weight case - stat: m = %f, s= %f, sr = %f \n",fMean, fSigma, midspread);
      return;
   }
   else {
      // compute statistics using the data
      SetMean();
      Double_t midspread = ComputeMidspread();
      SetSigma(midspread);
      //printf("un-weight case - stat: m = %f, s= %f, sr = %f \n",fMean, fSigma, fSigmaRob);
   }
}

Double_t TKDE::ComputeMidspread () {
   // Computes the inter-quartile range from the data
   std::sort(fEvents.begin(), fEvents.end());
   Double_t quantiles[2] = {0.0, 0.0};
   Double_t prob[2] = {0.25, 0.75};
   TMath::Quantiles(fEvents.size(), 2, &fEvents[0], quantiles, prob);
   Double_t lowerquartile = quantiles[0];
   Double_t upperquartile = quantiles[1];
   return upperquartile - lowerquartile;
}

void TKDE::SetUserCanonicalBandwidth() {
   // Computes the user's input kernel function canonical bandwidth
   fCanonicalBandwidths[kUserDefined] = std::pow(ComputeKernelL2Norm() / std::pow(ComputeKernelSigma2(), 2), 1. / 5.);
}

void TKDE::SetUserKernelSigma2() {
   // Computes the user's input kernel function sigma2
   fKernelSigmas2[kUserDefined] = ComputeKernelSigma2();
}

TKDE::KernelIntegrand::KernelIntegrand(const TKDE* kde, EIntegralResult intRes) : fKDE(kde), fIntegralResult(intRes) {
   // Internal class constructor
}

Double_t TKDE::KernelIntegrand::operator()(Double_t x) const {
   // Internal class unary function
   if (fIntegralResult == kNorm) {
      return std::pow((*fKDE->fKernelFunction)(x), 2);
   } else if (fIntegralResult == kMu) {
      return x * (*fKDE->fKernelFunction)(x);
   } else if (fIntegralResult == kSigma2) {
      return std::pow(x, 2) * (*fKDE->fKernelFunction)(x);
   } else if (fIntegralResult == kUnitIntegration) {
      return (*fKDE->fKernelFunction)(x);
   } else {
      return -1;
   }
}

TF1* TKDE::GetKDEFunction(UInt_t npx, Double_t xMin , Double_t xMax) {
   //Returns the estimated density
   TString name = "KDEFunc_"; name+= GetName();
   TString title = "KDE "; title+= GetTitle();
   if (xMin >= xMax) { xMin = fXMin; xMax = fXMax; }
   //Do not register the TF1 to global list
   bool previous = TF1::DefaultAddToGlobalList(kFALSE);
   TF1 * pdf = new TF1(name.Data(), this, xMin, xMax, 0);
   TF1::DefaultAddToGlobalList(previous);
   if (npx > 0) pdf->SetNpx(npx);
   pdf->SetTitle(title);
   return pdf;
}

TF1* TKDE::GetPDFUpperConfidenceInterval(Double_t confidenceLevel, UInt_t npx, Double_t xMin , Double_t xMax) {
   // Returns the upper estimated density
   TString name;
   name.Form("KDE_UpperCL%f5.3_%s",confidenceLevel,GetName());
   if (xMin >= xMax) { xMin = fXMin; xMax = fXMax; }
   TF1 * upperPDF = new TF1(name, this, &TKDE::UpperConfidenceInterval, xMin, xMax, 1);
   upperPDF->SetParameter(0, confidenceLevel);
   if (npx > 0) upperPDF->SetNpx(npx);
   TF1 * f =  (TF1*)upperPDF->Clone();
   delete upperPDF;
   return f;
}

TF1* TKDE::GetPDFLowerConfidenceInterval(Double_t confidenceLevel, UInt_t npx, Double_t xMin , Double_t xMax) {
   // Returns the upper estimated density
   TString name;
   name.Form("KDE_LowerCL%f5.3_%s",confidenceLevel,GetName());
   if (xMin >= xMax) { xMin = fXMin; xMax = fXMax; }
   TF1 * lowerPDF = new TF1(name, this, &TKDE::LowerConfidenceInterval, xMin, xMax, 1);
   lowerPDF->SetParameter(0, confidenceLevel);
   if (npx > 0) lowerPDF->SetNpx(npx);
   TF1 * f = (TF1*)lowerPDF->Clone();
   delete lowerPDF;
   return f;
}

TF1* TKDE::GetKDEApproximateBias(UInt_t npx, Double_t xMin , Double_t xMax) {
   // Returns the approximate bias
   TString name = "KDE_Bias_"; name += GetName();
   if (xMin >= xMax) { xMin = fXMin; xMax = fXMax; }
   TF1 * approximateBias = new TF1(name, this, &TKDE::ApproximateBias, xMin, xMax, 0);
   if (npx > 0) approximateBias->SetNpx(npx);
   TF1 * f =  (TF1*)approximateBias->Clone();
   delete approximateBias;
   return f;
}
