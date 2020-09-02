/// \file ROOT/RFitImpl.hxx
/// \ingroup Hist ROOT7
/// \author Claire Guyot
/// \date 2020-07
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT7_RFitImpl
#define ROOT7_RFitImpl

template <int DIMENSIONS, class PRECISION, template <int D_, class P_> class... STAT>
class RHist;

#include "HFitInterface.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TF1.h"

#include "Fit/DataOptions.h"
#include "Fit/DataRange.h"
#include "Fit/FitConfig.h"
#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "Fit/FitExecutionPolicy.h"

#include "TMath.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/IParamFunction.h"

#include "TError.h"

namespace ROOT {
namespace Experimental {
namespace RFit {

   bool AdjustError(const ROOT::Fit::DataOptions & option, double & error, double value = 1);
   int CheckFitFunction(const TF1 * f1, int dim);
   void GetFunctionRange(const TF1 & f1, ROOT::Fit::DataRange & range);

   template <int DIMENSIONS, class PRECISION, template <int D_, class P_> class... STAT>
   void BinContentToBinData(RHist<DIMENSIONS, PRECISION, STAT...> & hist, ROOT::Fit::BinData & fitData,
                                    TF1 *f1, const ROOT::Fit::DataOptions & fitOption, const ROOT::Fit::DataRange & range)
   {
      int ndim = hist.GetNDim();
      int nPoints = hist.GetImpl()->GetNBins();
      bool useRange = (range.Size(0) > 0);

      // Get the error on bin data to initialize
      ROOT::Fit::BinData::ErrorType errorType = ROOT::Fit::BinData::kValueError;
      if (fitOption.fErrors1) {
         errorType =  ROOT::Fit::BinData::kNoError;
      }
      fitData.Initialize(nPoints, ndim, errorType);
      
      double xmin = 0, ymin = 0, zmin = 0;
      double xmax = 0, ymax = 0, zmax = 0;
      switch (ndim) {
         case 1:
            range.GetRange(xmin,xmax);
            break;
         case 2:
            range.GetRange(xmin,xmax,ymin,ymax);
            break;
         case 3:
            range.GetRange(xmin,xmax,ymin,ymax,zmin,zmax);
         default:
            break;
      }
      
      // Convert RHist bin data to BinData
      for (auto bin: hist) {
         // Don't add bin if out of range
         switch (ndim) {
            case 1:
               if (useRange && (bin.GetFrom()[0] < xmin || bin.GetTo()[0] > xmax))
                  continue;
               break;
            case 2:
               if (useRange && (bin.GetFrom()[0] < xmin || bin.GetTo()[0] > xmax || 
                                 bin.GetFrom()[1] < ymin || bin.GetTo()[1] > ymax ))
                  continue;
               break;
            case 3:
               if (useRange && (bin.GetFrom()[0] < xmin || bin.GetTo()[0] > xmax || 
                                 bin.GetFrom()[1] < ymin || bin.GetTo()[1] > ymax ||
                                 bin.GetFrom()[2] < zmin || bin.GetTo()[2] > zmax ))
                  continue;
            default:
               break;
         }

         // Don't add bin if rejected by TF1
         if (f1) {
            TF1::RejectPoint(false);
            (*f1)(bin);
            if (TF1::RejectedPoint())
               continue;
         }

         // Add bin data depending on error requested
         if (fitOption.fErrors1) {
            fitData.Add(bin.GetCenter(), bin.GetContent());
         }
         else if (errorType == ROOT::Fit::BinData::kValueError) {
            if (!RFit::AdjustError(fitOption, bin.GetUncertainty(), bin.GetContent()))
               continue;
            fitData.Add(bin.GetCenter(), bin.GetContent(), bin.GetUncertainty());
         }
      }
   }

   /// For now, only for dim <= 3 (due to the use of functions restricting dimension like in DataRange) 
   template <int DIMENSIONS, class PRECISION,
            template <int D_, class P_> class... STAT>
   TFitResultPtr Fit(RHist<DIMENSIONS, PRECISION, STAT...> & hist, TF1 *f1, ROOT::Fit::DataOptions & fitOption, ROOT::Fit::FitConfig & fitConfig)
   {
      // Make sure function and histogram are compatible for fitting
      int ndim = hist.GetNDim();
      if (ndim == 0 || ndim > 3)
         // Arbitrary error return
         return -1;
      int checkResult = RFit::CheckFitFunction(f1, ndim);
      if (checkResult != 0)
         return checkResult;

      // If function dimension is less than hist dimension, then integral option is not possible
      if (f1->GetNdim() < ndim && fitOption.fIntegral)
         throw std::runtime_error("Attempted to fit a model function with a lesser dimension than that of the data object");

      // If specified, use range of function when fitting
      ROOT::Fit::DataRange range(ndim);
      if (fitOption.fUseRange) {
         RFit::GetFunctionRange(*f1,range);
      }

      Int_t special = f1->GetNumber();
      Bool_t linear = f1->IsLinear();
      Int_t npar = f1->GetNpar();
      // If polynomial function, make linear
      if (special == 299 + npar)  linear = kTRUE;
      // If option "integral" true, make non linear
      if (fitOption.fIntegral)
         linear = kFALSE;

      // Create an empty TFitResult, result of the fitting
      std::shared_ptr<TFitResult> tFitResult(new TFitResult());
      // Create the fitter from an empty fit result
      std::shared_ptr<ROOT::Fit::Fitter> fitter(new ROOT::Fit::Fitter(std::static_pointer_cast<ROOT::Fit::FitResult>(tFitResult)));
      // Set config options for fitter
      ROOT::Fit::FitConfig & fitterConfig = fitter->Config();
      fitterConfig = fitConfig;

      // Option special cases
      if (fitOption.fExpErrors) fitOption.fUseEmpty = true;  // use empty bins in log-likelihood fits
      if (special == 300) fitOption.fCoordErrors = false; // no need to use coordinate errors in a pol0 fit
      if (!fitOption.fErrors1) fitOption.fUseEmpty = true; // use empty bins with weight=1

      // Fill data for fitting
      std::shared_ptr<ROOT::Fit::BinData> fitData(new ROOT::Fit::BinData(fitOption,range));
      // TODO: check the different errors/uncertainty wanted for implementation

      BinContentToBinData(hist, fitData, f1, fitOption, range);
      
      if (fitData->Size() == 0) {
         Warning("Fit","Fit data is empty ");
         return -1;
      }

      // Switch off linear fitting in case data has coordinate errors and the option is set
      if (fitData->GetErrorType() == ROOT::Fit::BinData::kCoordError && fitData->Opt().fCoordErrors ) linear = false;
      // Linear fit cannot be done also in case of asymmetric errors
      if (fitData->GetErrorType() == ROOT::Fit::BinData::kAsymError && fitData->Opt().fAsymErrors ) linear = false;

      // TShis functions use the TVirtualFitter
      if (special != 0 && !linear) {
         if      (special == 100)      ROOT::Fit::InitGaus  (*fitData,f1); // gaussian
         else if (special == 110 || special == 112)   ROOT::Fit::Init2DGaus(*fitData,f1); // 2D gaussians ( xygaus or bigaus)
         else if (special == 400)      ROOT::Fit::InitGaus  (*fitData,f1); // landau (use the same)
         else if (special == 410)      ROOT::Fit::Init2DGaus(*fitData,f1); // 2D landau (use the same)
         else if (special == 200)      ROOT::Fit::InitExpo  (*fitData, f1); // exponential

      }

      // Set the fit function
      if (linear)
         fitter->SetFunction(ROOT::Math::WrappedMultiTF1(*f1));
      else
         fitter->SetFunction(static_cast<const ROOT::Math::IParamMultiFunction &>(ROOT::Math::WrappedMultiTF1(*f1)));
      // Create the wrapped TF1 function to transform into IParamFunction
      //ROOT::Math::WrappedTF1 wrappedTF1(*f1);
      // Create the IParamFunction
      //ROOT::Math::IParamFunction & paramFunction = wrappedTF1;

      // Error normalization in case of zero error in the data
      if (fitData->GetErrorType() == ROOT::Fit::BinData::kNoError) fitterConfig.SetNormErrors(true);
      // Error normalization also in case of weights = 1
      if (fitData->Opt().fErrors1)  fitterConfig.SetNormErrors(true);
      // Normalize errors also in case you are fitting a Ndim histo with a N-1 function
      if (int(fitData->NDim())  == ndim -1 ) fitterConfig.SetNormErrors(true);

      // Parameter settings and transfer of the parameters values, names and limits from the functions
      // are done automatically in Fitter.cxx
      for (int i = 0; i < npar; ++i) {
         ROOT::Fit::ParameterSettings & parSettings = fitterConfig.ParSettings(i);

         // check limits
         double plow,pup;
         f1->GetParLimits(i,plow,pup);
         if (plow*pup != 0 && plow >= pup) { // this is a limitation - cannot fix a parameter to zero value
            parSettings.Fix();
         }
         else if (plow < pup ) {
            if (!TMath::Finite(pup) && TMath::Finite(plow) )
               parSettings.SetLowerLimit(plow);
            else if (!TMath::Finite(plow) && TMath::Finite(pup) )
               parSettings.SetUpperLimit(pup);
            else
               parSettings.SetLimits(plow,pup);
         }

         // set the parameter step size (by default are set to 0.3 of value)
         // if function provides meaningful error values
         double err = f1->GetParError(i);
         if ( err > 0)
            parSettings.SetStepSize(err);
         else if (plow < pup && TMath::Finite(plow) && TMath::Finite(pup) ) { // in case of limits improve step sizes
            double step = 0.1 * (pup - plow);
            // check if value is not too close to limit otherwise trim value
            if (  parSettings.Value() < pup && pup - parSettings.Value() < 2 * step  )
               step = (pup - parSettings.Value() ) / 2;
            else if ( parSettings.Value() > plow && parSettings.Value() - plow < 2 * step )
               step = (parSettings.Value() - plow ) / 2;

            parSettings.SetStepSize(step);
         }
      }

      // Set all default minimizer options (tolerance, max iterations, etc..)
      ROOT::Math::MinimizerOptions minimizerOpts = fitConfig.MinimizerOptions();
      fitterConfig.SetMinimizerOptions(minimizerOpts);

      // Specific minimizer options depending on minimizer
      if (linear) {
         fitterConfig.SetMinimizer("Linear","");
      }

      // Run fitting
      bool fitDone = false;

      fitDone = fitter->Fit(fitData, ROOT::Fit::ExecutionPolicy::kSerial);

      checkResult |= !fitDone;

      const ROOT::Fit::FitResult & fitResult = fitter->Result();

      // Set directly the fit result in TF1
      checkResult = fitResult.Status();
      if (!fitResult.IsEmpty()) {
         f1->SetChisquare(fitResult.Chi2());
         f1->SetNDF(fitResult.Ndf());
         f1->SetNumberFitPoints(fitData->Size());

         assert((Int_t)fitResult.Parameters().size() >= f1->GetNpar());
         f1->SetParameters( const_cast<double*>(&(fitResult.Parameters().front())));
         if (int(fitResult.Errors().size()) >= f1->GetNpar())
            f1->SetParErrors( &(fitResult.Errors().front()));
      }

      // Print the result of the fitting
      if (fitter->GetMinimizer() && fitConfig.MinimizerType() == "Minuit" &&
               !fitConfig.NormalizeErrors()) {
            fitter->GetMinimizer()->PrintResults();
      }
      else {
         fitResult.PrintCovMatrix(std::cout);
         fitResult.Print(std::cout);
      }

      return TFitResultPtr(checkResult);
   }

}// namespace RFit
}// namespace Experimental
}// namespace ROOT

#endif
