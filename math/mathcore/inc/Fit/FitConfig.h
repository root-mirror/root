// @(#)root/mathcore:$Id$
// Author: L. Moneta Thu Sep 21 16:21:29 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class FitConfig

#ifndef ROOT_Fit_FitConfig
#define ROOT_Fit_FitConfig


#include "Fit/ParameterSettings.h"

#include "Math/MinimizerOptions.h"

#include "Math/IParamFunctionfwd.h"

#include "TMath.h"

#include <vector>
#include <string>

namespace ROOT {

   namespace Math {

      class Minimizer;
      class MinimizerOptions;
   }

   namespace Fit {

      class FitResult;

//___________________________________________________________________________________
/**
   Class describing the configuration of the fit, options and parameter settings
   using the ROOT::Fit::ParameterSettings class

   @ingroup FitMain
*/
class FitConfig {

public:

   /**
      Default constructor
   */
   FitConfig (unsigned int npar = 0);


   /*
     Copy constructor
    */
   FitConfig(const FitConfig & rhs);

   /**
      Destructor
   */
   ~FitConfig ();

   /*
     Assignment operator
   */
   FitConfig & operator= (const FitConfig & rhs);


   /**
      get the parameter settings for the i-th parameter (const method)
   */
   const ParameterSettings & ParSettings(unsigned int i) const { return fSettings.at(i); }

   /**
      get the parameter settings for the i-th parameter (non-const method)
   */
   ParameterSettings & ParSettings(unsigned int i) { return fSettings.at(i); }

   /**
      get the vector of parameter settings  (const method)
   */
   const std::vector<ROOT::Fit::ParameterSettings> & ParamsSettings() const { return fSettings; }

   /**
      get the vector of parameter settings  (non-const method)
   */
   std::vector<ROOT::Fit::ParameterSettings> & ParamsSettings() { return fSettings; }

   /**
      number of parameters settings
    */
   unsigned int NPar() const { return fSettings.size(); }

   /**
      return a vector of stored parameter values (i.e initial fit parameters)
    */
   std::vector<double> ParamsValues() const;


   /**
      set the parameter settings from a model function.
      Create always new parameter setting list from a given model function
   */
   template <class T>
   void CreateParamsSettings(const ROOT::Math::IParamMultiFunctionTempl<T> &func) {
      // initialize from model function
      // set the parameters values from the function
      unsigned int npar = func.NPar();
      const double *begin = func.Parameters();
      if (begin == 0) {
         fSettings = std::vector<ParameterSettings>(npar);
         return;
      }

      fSettings.clear();
      fSettings.reserve(npar);
      const double *end = begin + npar;
      unsigned int i = 0;
      for (const double *ipar = begin; ipar != end; ++ipar) {
         double val = *ipar;
         double step = 0.3 * fabs(val); // step size is 30% of par value
         // double step = 2.0*fabs(val);   // step size is 30% of par value
         if (val == 0) step = 0.3;

         fSettings.push_back(ParameterSettings(func.ParameterName(i), val, step));
#ifdef DEBUG
         std::cout << "FitConfig: add parameter " << func.ParameterName(i) << " val = " << val << std::endl;
#endif
         i++;
      }
   }

   /**
      set the parameter settings from number of parameters and a vector of values and optionally step values. If there are not existing or number of parameters does not match existing one, create a new parameter setting list.
   */
   void SetParamsSettings(unsigned int npar, const double * params, const double * vstep = 0);

   /*
     Set the parameter settings from a vector of parameter settings
   */
   void SetParamsSettings (const std::vector<ROOT::Fit::ParameterSettings>& pars ) {
      fSettings = pars;
   }


   /*
     Set the parameter settings from a fit Result
   */
   void SetFromFitResult (const FitResult & rhs);



   /**
      create a new minimizer according to chosen configuration
   */
   ROOT::Math::Minimizer * CreateMinimizer();



   /**
      access to the minimizer  control parameter (non const method)
   */
   ROOT::Math::MinimizerOptions & MinimizerOptions()  { return fMinimizerOpts; }


#ifndef __CINT__   // this method fails on Windows
   /**
      set all the minimizer options using class MinimizerOptions
    */
   void SetMinimizerOptions(const ROOT::Math::MinimizerOptions & minopt);
#endif


   /**
      set minimizer type
   */
   void SetMinimizer(const char * type, const char * algo = 0) {
      if (type) fMinimizerOpts.SetMinimizerType(type);
      if (algo) fMinimizerOpts.SetMinimizerAlgorithm(algo);
   }

   /**
      return type of minimizer package
   */
   const std::string & MinimizerType() const { return fMinimizerOpts.MinimizerType(); }

   /**
      return type of minimizer algorithms
   */
   const std::string & MinimizerAlgoType() const { return fMinimizerOpts.MinimizerAlgorithm(); }

   /**
    * return Minimizer full name (type / algorithm)
    */
   std::string MinimizerName() const;

   /**
      flag to check if resulting errors are be normalized according to chi2/ndf
   */
   bool NormalizeErrors() const { return fNormErrors; }

   ///do analysis for parabolic errors
   bool ParabErrors() const { return fParabErrors; }

   ///do minos errros analysis on the  parameters
   bool MinosErrors() const { return fMinosErrors; }

   ///Update configuration after a fit using the FitResult
   bool UpdateAfterFit() const { return fUpdateAfterFit; }

   ///Apply Weight correction for error matrix computation
   bool UseWeightCorrection() const { return fWeightCorr; }


   /// return vector of parameter indeces for which the Minos Error will be computed
   const std::vector<unsigned int> & MinosParams() const { return fMinosParams; }

   /**
      set the option to normalize the error on the result  according to chi2/ndf
   */
   void SetNormErrors(bool on = true) { fNormErrors= on; }

   ///set parabolic erros
   void SetParabErrors(bool on = true) { fParabErrors = on; }

   ///set Minos erros computation to be performed after fitting
   void SetMinosErrors(bool on = true) { fMinosErrors = on; }

   ///apply the weight correction for error matric computation
   void SetWeightCorrection(bool on = true) { fWeightCorr = on; }

   /// set parameter indeces for running Minos
   /// this can be used for running Minos on a subset of parameters - otherwise is run on all of them
   /// if MinosErrors() is set
   void SetMinosErrors(const std::vector<unsigned int> & paramInd ) {
      fMinosErrors = true;
      fMinosParams = paramInd;
   }

   ///Update configuration after a fit using the FitResult
   void SetUpdateAfterFit(bool on = true) { fUpdateAfterFit = on; }


   /**
      static function to control default minimizer type and algorithm
   */
   static void SetDefaultMinimizer(const char * type, const char * algo = 0);




protected:


private:

   bool fNormErrors;       // flag for error normalization
   bool fParabErrors;      // get correct parabolic errors estimate (call Hesse after minimizing)
   bool fMinosErrors;      // do full error analysis using Minos
   bool fUpdateAfterFit;   // update the configuration after a fit using the result
   bool fWeightCorr;       // apply correction to errors for weights fits

   std::vector<ROOT::Fit::ParameterSettings> fSettings;  // vector with the parameter settings
   std::vector<unsigned int> fMinosParams;               // vector with the parameter indeces for running Minos

   ROOT::Math::MinimizerOptions fMinimizerOpts;   //minimizer control parameters including name and algo type

};

   } // end namespace Fit

} // end namespace ROOT


#endif /* ROOT_Fit_FitConfig */
