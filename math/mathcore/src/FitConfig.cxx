// @(#)root/mathcore:$Id$
// Author: L. Moneta Thu Sep 21 16:21:29 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for class FitConfig

#include "Fit/FitConfig.h"

#include "Fit/FitResult.h"

#include "Math/IParamFunction.h"
#include "Math/Util.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include <cmath>

#include <string>
#include <sstream>

#include "Math/Error.h"

//#define DEBUG
#ifdef DEBUG
#endif

namespace ROOT {

namespace Fit {



FitConfig::FitConfig(unsigned int npar) :
   fNormErrors(false),
   fParabErrors(false), // ensure that in any case correct parabolic errors are estimated
   fMinosErrors(false),    // do full Minos error analysis for all parameters
   fUpdateAfterFit(true),    // update after fit
   fWeightCorr(false),
   fSettings(std::vector<ParameterSettings>(npar) )
{
   // constructor implementation
}


FitConfig::~FitConfig()
{
   // destructor implementation. No Operations
}

FitConfig::FitConfig(const FitConfig &rhs) {
   // Implementation of copy constructor
   (*this) = rhs;
}

FitConfig & FitConfig::operator = (const FitConfig &rhs) {
   // Implementation of assignment operator.
   if (this == &rhs) return *this;  // time saving self-test

   fNormErrors = rhs.fNormErrors;
   fParabErrors = rhs.fParabErrors;
   fMinosErrors = rhs.fMinosErrors;
   fUpdateAfterFit = rhs.fUpdateAfterFit;
   fWeightCorr     = rhs.fWeightCorr;

   fSettings = rhs.fSettings;
   fMinosParams = rhs.fMinosParams;

   fMinimizerOpts = rhs.fMinimizerOpts;

   return *this;
}

void FitConfig::SetFromFitResult(const FitResult &result) {
   // Implementation of setting of parameters from the result of the fit
   // all the other options will stay the same.
   // If the size of parameters do not match they will be re-created
   // but in that case the bound on the parameter will be lost

   unsigned int npar = result.NPar();
   if (fSettings.size() !=  npar) {
      fSettings.clear();
      fSettings.resize(npar);
   }
   // fill the parameter settings
   for (unsigned int i = 0; i < npar; ++i) {
      if (result.IsParameterFixed(i) )
         fSettings[i].Set(result.ParName(i), result.Value(i) );
      else {
         fSettings[i].Set( result.ParName(i), result.Value(i), result.Error(i) );
         // check if parameter is bound
         double lower = 0;
         double upper = 0;
         if (result.ParameterBounds(i,lower,upper) ) {
            if (lower == -std::numeric_limits<double>::infinity()) fSettings[i].SetUpperLimit(upper);
            else if (upper ==  std::numeric_limits<double>::infinity()) fSettings[i].SetLowerLimit(lower);
            else fSettings[i].SetLimits(lower,upper);
         }

         // query if parameter needs to run Minos
         if (result.HasMinosError(i) ) {
            if (fMinosParams.size() == 0) {
               fMinosErrors = true;
               fMinosParams.reserve(npar-i);
            }
            fMinosParams.push_back(i);
         }
      }
   }

   // set information about errors
   SetNormErrors( result.NormalizedErrors() );

   // set also minimizer type
   // algorithm is after " / "
   const std::string & minname = result.MinimizerType();
   size_t pos = minname.find(" / ");
   if (pos != std::string::npos) {
      std::string minimType = minname.substr(0,pos);
      std::string algoType = minname.substr(pos+3,minname.length() );
      SetMinimizer(minimType.c_str(), algoType.c_str() );
   }
   else {
      SetMinimizer(minname.c_str());
   }
}


void FitConfig::SetParamsSettings(unsigned int npar, const double *params, const double * vstep ) {
   // initialize FitConfig from given parameter values and step sizes
   // if npar different than existing one - clear old one and create new ones
   if (params == 0) {
      fSettings =  std::vector<ParameterSettings>(npar);
      return;
   }
   // if a vector of parameters is given and parameters are not existing or are of different size
   bool createNew = false;
   if (npar != fSettings.size() ) {
      fSettings.clear();
      fSettings.reserve(npar);
      createNew = true;
   }
   unsigned int i = 0;
   const double * end = params+npar;
   for (const double * ipar = params; ipar !=  end; ++ipar) {
      double val = *ipar;
      double step = 0;
      if (vstep == 0) {
         step = 0.3*std::fabs(val);   // step size is 30% of par value
         //double step = 2.0*std::fabs(val);   // step size is 30% of par value
         if (val ==  0) step  =  0.3;
      }
      else
         step = vstep[i];

      if (createNew)
         fSettings.push_back( ParameterSettings("Par_" + ROOT::Math::Util::ToString(i), val, step ) );
      else {
         fSettings[i].SetValue(val);
         fSettings[i].SetStepSize(step);
      }

      i++;
   }
}

ROOT::Math::Minimizer * FitConfig::CreateMinimizer() {
   // create minimizer according to the chosen configuration using the
   // plug-in manager

   // in case of empty string usesd default values
   if (fMinimizerOpts.MinimizerType().empty())
      fMinimizerOpts.SetMinimizerType(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str());
   if (fMinimizerOpts.MinimizerAlgorithm().empty())
      fMinimizerOpts.SetMinimizerAlgorithm(ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());

   const std::string &minimType = fMinimizerOpts.MinimizerType();
   const std::string & algoType  = fMinimizerOpts.MinimizerAlgorithm();

   std::string  defaultMinim = ROOT::Math::MinimizerOptions::DefaultMinimizerType();

   ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minimType, algoType);
   // check if a different minimizer is used (in case a default value is passed, then set correctly in FitConfig)
   const std::string & minim_newDefault = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
   if (defaultMinim != minim_newDefault )  fMinimizerOpts.SetMinimizerType(minim_newDefault.c_str());

   if (min == 0) {
      // if creation of minimizer failed force the use by default of Minuit
      std::string minim2 = "Minuit";
      if (minimType == "Minuit") minim2 = "Minuit2";
      if (minimType != minim2 ) {
         std::string msg = "Could not create the " + minimType + " minimizer. Try using the minimizer " + minim2;
         MATH_WARN_MSG("FitConfig::CreateMinimizer",msg.c_str());
         min = ROOT::Math::Factory::CreateMinimizer(minim2,"Migrad");
         if (min == 0) {
            MATH_ERROR_MSG("FitConfig::CreateMinimizer","Could not create the Minuit2 minimizer");
            return 0;
         }
         SetMinimizer( minim2.c_str(),"Migrad");
      }
      else {
         std::string msg = "Could not create the Minimizer " + minimType;
         MATH_ERROR_MSG("FitConfig::CreateMinimizer",msg.c_str());
         return 0;
      }
   }

   // set default max of function calls according to the number of parameters
   // formula from Minuit2 (adapted)
   if (fMinimizerOpts.MaxFunctionCalls() == 0) {
      unsigned int npar =  fSettings.size();
      int maxfcn = 1000 + 100*npar + 5*npar*npar;
      fMinimizerOpts.SetMaxFunctionCalls(maxfcn);
   }


   // set default minimizer control parameters
   min->SetPrintLevel( fMinimizerOpts.PrintLevel() );
   min->SetMaxFunctionCalls( fMinimizerOpts.MaxFunctionCalls() );
   min->SetMaxIterations( fMinimizerOpts.MaxIterations() );
   min->SetTolerance( fMinimizerOpts.Tolerance() );
   min->SetPrecision( fMinimizerOpts.Precision() );
   min->SetValidError( fParabErrors );
   min->SetStrategy( fMinimizerOpts.Strategy() );
   min->SetErrorDef( fMinimizerOpts.ErrorDef() );


   return min;
}

std::string FitConfig::MinimizerName() const
{
   // set minimizer type
   std::string name = MinimizerType();

   // append algorithm name for minimizer that support it
   if ((name.find("Fumili") == std::string::npos) && (name.find("GSLMultiFit") == std::string::npos)) {
      if (MinimizerAlgoType() != "")
         name += " / " + MinimizerAlgoType();
   }
   return name;
}

void FitConfig::SetDefaultMinimizer(const char * type, const char *algo ) {
   // set the default minimizer type and algorithms
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer(type, algo);
}

void FitConfig::SetMinimizerOptions(const ROOT::Math::MinimizerOptions & minopt) {
   // set all the minimizer options
   fMinimizerOpts = minopt;
}

std::vector<double> FitConfig::ParamsValues() const {

   std::vector<double> params(NPar() );
   for (unsigned int i = 0; i < params.size();  ++i) {
      params[i] = fSettings[i].Value();
   }
   return params;
}
   } // end namespace Fit

} // end namespace ROOT
