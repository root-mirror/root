// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
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
 *   Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
 *
 *****************************************************************************/


/** \class RooStats::LikelihoodInterval
    \ingroup Roostats

   LikelihoodInterval is a concrete implementation of the RooStats::ConfInterval interface.
   It implements a connected N-dimensional intervals based on the contour of a likelihood ratio.
   The boundary of the interval is equivalent to a MINUIT/MINOS contour about the maximum likelihood estimator

   The interval does not need to be an ellipse (eg. it is not the HESSE error matrix).
   The level used to make the contour is the same as that used in MINOS, eg. it uses Wilks' theorem,
   which states that under certain regularity conditions the function -2* log (profile likelihood ratio) is asymptotically distributed as a chi^2 with N-dof, where
   N is the number of parameters of interest.


   Note, a boundary on the parameter space (eg. s>= 0) or a degeneracy (eg. mass of signal if Nsig = 0) can lead to violations of the conditions necessary for Wilks'
   theorem to be true.

   Also note, one can use any RooAbsReal as the function that will be used in the contour; however, the level of the contour
   is based on Wilks' theorem as stated above.


#### References

*  1. F. James., Minuit.Long writeup D506, CERN, 1998.

*/


#include "RooStats/LikelihoodInterval.h"
#include "RooStats/RooStatsUtils.h"

#include "RooAbsReal.h"
#include "RooMsgService.h"

#include "Math/WrappedFunction.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/MinimizerOptions.h"
#include "RooFunctor.h"
#include "RooProfileLL.h"

#include "TMinuitMinimizer.h"

#include <string>
#include <algorithm>
#include <functional>
#include <ctype.h>   // need to use c version of toupper defined here

/*
// for debugging
#include "RooNLLVar.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
*/

ClassImp(RooStats::LikelihoodInterval); ;

using namespace RooStats;
using namespace std;


////////////////////////////////////////////////////////////////////////////////
/// Default constructor with name and title

LikelihoodInterval::LikelihoodInterval(const char* name) :
   ConfInterval(name), fBestFitParams(0), fLikelihoodRatio(0), fConfidenceLevel(0.95)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Alternate constructor taking a pointer to the profile likelihood ratio, parameter of interest and
/// optionally a snapshot of best parameter of interest for interval

LikelihoodInterval::LikelihoodInterval(const char* name, RooAbsReal* lr, const RooArgSet* params,  RooArgSet * bestParams) :
   ConfInterval(name),
   fParameters(*params),
   fBestFitParams(bestParams),
   fLikelihoodRatio(lr),
   fConfidenceLevel(0.95)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

LikelihoodInterval::~LikelihoodInterval()
{
   if (fBestFitParams) delete fBestFitParams;
   if (fLikelihoodRatio) delete fLikelihoodRatio;
}


////////////////////////////////////////////////////////////////////////////////
/// This is the main method to satisfy the RooStats::ConfInterval interface.
/// It returns true if the parameter point is in the interval.

Bool_t LikelihoodInterval::IsInInterval(const RooArgSet &parameterPoint) const
{
   RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  // Method to determine if a parameter point is in the interval
  if( !this->CheckParameters(parameterPoint) ) {
    std::cout << "parameters don't match" << std::endl;
    RooMsgService::instance().setGlobalKillBelow(msglevel);
    return false;
  }

  // make sure likelihood ratio is set
  if(!fLikelihoodRatio) {
    std::cout << "likelihood ratio not set" << std::endl;
    RooMsgService::instance().setGlobalKillBelow(msglevel);
    return false;
  }



  // set parameters
  SetParameters(&parameterPoint, fLikelihoodRatio->getVariables() );


  // evaluate likelihood ratio, see if it's bigger than threshold
  if (fLikelihoodRatio->getVal()<0){
    std::cout << "The likelihood ratio is < 0, indicates a bad minimum or numerical precision problems.  Will return true" << std::endl;
    RooMsgService::instance().setGlobalKillBelow(msglevel);
    return true;
  }


  // here we use Wilks' theorem.
  if ( TMath::Prob( 2* fLikelihoodRatio->getVal(), parameterPoint.getSize()) < (1.-fConfidenceLevel) ){
    RooMsgService::instance().setGlobalKillBelow(msglevel);
    return false;
  }


  RooMsgService::instance().setGlobalKillBelow(msglevel);

  return true;

}

////////////////////////////////////////////////////////////////////////////////
/// returns list of parameters

RooArgSet* LikelihoodInterval::GetParameters() const
{
   return new RooArgSet(fParameters);
}

////////////////////////////////////////////////////////////////////////////////
/// check that the parameters are correct

Bool_t LikelihoodInterval::CheckParameters(const RooArgSet &parameterPoint) const
{
  if (parameterPoint.getSize() != fParameters.getSize() ) {
    std::cout << "size is wrong, parameters don't match" << std::endl;
    return false;
  }
  if ( ! parameterPoint.equals( fParameters ) ) {
    std::cout << "size is ok, but parameters don't match" << std::endl;
    return false;
  }
  return true;
}



////////////////////////////////////////////////////////////////////////////////
/// Compute lower limit, check first if limit has been computed
/// status is a boolean flag which will b set to false in case of error
/// and is true if calculation is successful
/// in case of error return also a lower limit value of zero

Double_t LikelihoodInterval::LowerLimit(const RooRealVar& param, bool & status)
{
   double lower = 0;
   double upper = 0;
   status = FindLimits(param, lower, upper);
   return lower;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute upper limit, check first if limit has been computed
/// status is a boolean flag which will b set to false in case of error
/// and is true if calculation is successful
/// in case of error return also a lower limit value of zero

Double_t LikelihoodInterval::UpperLimit(const RooRealVar& param, bool & status)
{
   double lower = 0;
   double upper = 0;
   status = FindLimits(param, lower, upper);
   return upper;
}


void LikelihoodInterval::ResetLimits() {
   // reset map with cached limits - called every time the test size or CL has been changed
   fLowerLimits.clear();
   fUpperLimits.clear();
}


bool LikelihoodInterval::CreateMinimizer() {
   // internal function to create minimizer object needed to find contours or interval limits
   // (running MINOS).
   // Minimizer must be Minuit or Minuit2

   RooProfileLL * profilell = dynamic_cast<RooProfileLL*>(fLikelihoodRatio);
   if (!profilell) return false;

   RooAbsReal & nll  = profilell->nll();
   // bind the nll function in the right interface for the Minimizer class
   // as a function of only the parameters (poi + nuisance parameters)

   RooArgSet * partmp = profilell->getVariables();
   // need to remove constant parameters
   RemoveConstantParameters(partmp);

   RooArgList params(*partmp);
   delete partmp;

   // need to restore values and errors for POI
   if (fBestFitParams) {
      for (int i = 0; i < params.getSize(); ++i) {
         RooRealVar & par =  (RooRealVar &) params[i];
         RooRealVar * fitPar =  (RooRealVar *) (fBestFitParams->find(par.GetName() ) );
         if (fitPar) {
            par.setVal( fitPar->getVal() );
            par.setError( fitPar->getError() );
         }
      }
   }

   const auto& config = GetGlobalRooStatsConfig();

   // now do binding of NLL with a functor for Minimizer
   if (config.useLikelihoodOffset) {
      ccoutI(InputArguments) << "LikelihoodInterval: using nll offset - set all RooAbsReal to hide the offset  " << std::endl;
      RooAbsReal::setHideOffset(kFALSE); // need to keep this false
   }
   fFunctor = std::make_shared<RooFunctor>(nll, RooArgSet(), params);

   std::string minimType =  ROOT::Math::MinimizerOptions::DefaultMinimizerType();
   std::transform(minimType.begin(), minimType.end(), minimType.begin(), (int(*)(int)) tolower );
   *minimType.begin() = toupper( *minimType.begin());

   if (minimType != "Minuit" && minimType != "Minuit2") {
      ccoutE(InputArguments) << minimType << " is wrong type of minimizer for getting interval limits or contours - must use Minuit or Minuit2" << std::endl;
      return false;
   }
   // do not use static instance of TMInuit which could interfere with RooFit
   if (minimType == "Minuit")  TMinuitMinimizer::UseStaticMinuit(false);
   // create minimizer class
   fMinimizer = std::shared_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer(minimType, "Migrad"));

   if (!fMinimizer.get()) return false;

   fMinFunc = std::static_pointer_cast<ROOT::Math::IMultiGenFunction>(
      std::make_shared<ROOT::Math::WrappedMultiFunction<RooFunctor &>>(*fFunctor, fFunctor->nPar()) );
   fMinimizer->SetFunction(*fMinFunc);

   // set minimizer parameters
   assert( params.getSize() == int(fMinFunc->NDim()) );

   for (unsigned int i = 0; i < fMinFunc->NDim(); ++i) {
      RooRealVar & v = (RooRealVar &) params[i];
      fMinimizer->SetLimitedVariable( i, v.GetName(), v.getVal(), v.getError(), v.getMin(), v.getMax() );
   }
   // for finding the contour need to find first global minimum
   bool iret = fMinimizer->Minimize();
   if (!iret || fMinimizer->X() == 0) {
      ccoutE(Minimization) << "Error: Minimization failed  " << std::endl;
      return false;
   }

   //std::cout << "print minimizer result..........." << std::endl;
   //fMinimizer->PrintResults();

   return true;
}

bool LikelihoodInterval::FindLimits(const RooRealVar & param, double &lower, double & upper)
{
   // Method to find both lower and upper limits using MINOS
   // If cached values exist (limits have been already found) return them in that case
   // check first if limit has been computed
   // otherwise compute limit using MINOS
   // in case of failure lower and upper will maintain previous value (will not be modified)

   std::map<std::string, double>::const_iterator itrl = fLowerLimits.find(param.GetName());
   std::map<std::string, double>::const_iterator itru = fUpperLimits.find(param.GetName());
   if ( itrl != fLowerLimits.end() && itru != fUpperLimits.end() ) {
      lower = itrl->second;
      upper = itru->second;
      return true;
   }


   RooArgSet * partmp = fLikelihoodRatio->getVariables();
   RemoveConstantParameters(partmp);
   RooArgList params(*partmp);
   delete partmp;
   int ix = params.index(&param);
   if (ix < 0 ) {
      ccoutE(InputArguments) << "Error - invalid parameter " << param.GetName() << " specified for finding the interval limits " << std::endl;
      return false;
   }

   bool ret = true;
   if (!fMinimizer.get()) ret = CreateMinimizer();
   if (!ret) {
      ccoutE(Eval) << "Error returned from minimization of likelihood function - cannot find interval limits " << std::endl;
      return false;
   }

   assert(fMinimizer.get());

   // getting a 1D interval so ndf = 1
   double err_level = TMath::ChisquareQuantile(ConfidenceLevel(),1); // level for -2log LR
   err_level = err_level/2; // since we are using -log LR
   fMinimizer->SetErrorDef(err_level);

   unsigned int ivarX = ix;

   double elow = 0;
   double eup = 0;
   ret = fMinimizer->GetMinosError(ivarX, elow, eup );
   if (!ret)  {
      ccoutE(Minimization) << "Error  running Minos for parameter " << param.GetName() << std::endl;
      return false;
   }

   // WHEN error is zero normally is at limit
   if (elow == 0) {
      lower = param.getMin();
      ccoutW(Minimization) << "Warning: lower value for " << param.GetName() << " is at limit " << lower << std::endl;
   }
   else
      lower = fMinimizer->X()[ivarX] + elow;  // elow is negative

   if (eup == 0) {
      ccoutW(Minimization) << "Warning: upper value for " << param.GetName() << " is at limit " << upper << std::endl;
      upper = param.getMax();
   }
   else
      upper = fMinimizer->X()[ivarX] + eup;

   // store limits in the map
   // minos return error limit = minValue +/- error
   fLowerLimits[param.GetName()] = lower;
   fUpperLimits[param.GetName()] = upper;

   return true;
}


Int_t LikelihoodInterval::GetContourPoints(const RooRealVar & paramX, const RooRealVar & paramY, Double_t * x, Double_t *y, Int_t npoints ) {
   // use Minuit to find the contour of the likelihood function at the desired CL

   // check the parameters
   // variable index in minimizer
   // is index in the RooArgList obtained from the profileLL variables
   RooArgSet * partmp = fLikelihoodRatio->getVariables();
   RemoveConstantParameters(partmp);
   RooArgList params(*partmp);
   delete partmp;
   int ix = params.index(&paramX);
   int iy = params.index(&paramY);
   if (ix < 0 || iy < 0) {
      coutE(InputArguments) << "LikelihoodInterval - Error - invalid parameters specified for finding the contours; parX = " << paramX.GetName()
             << " parY = " << paramY.GetName() << std::endl;
         return 0;
   }

   bool ret = true;
   if (!fMinimizer.get()) ret = CreateMinimizer();
   if (!ret) {
      coutE(Eval) << "LikelihoodInterval - Error returned creating minimizer for likelihood function - cannot find contour points " << std::endl;
      return 0;
   }

   assert(fMinimizer.get());

   // getting a 2D contour so ndf = 2
   double cont_level = TMath::ChisquareQuantile(ConfidenceLevel(),2); // level for -2log LR
   cont_level = cont_level/2; // since we are using -log LR
   fMinimizer->SetErrorDef(cont_level);

   unsigned int ncp = npoints;
   unsigned int ivarX = ix;
   unsigned int ivarY = iy;
   coutI(Minimization)  << "LikelihoodInterval - Finding the contour of " << paramX.GetName() << " ( " << ivarX << " ) and " << paramY.GetName() << " ( " << ivarY << " ) " << std::endl;
   ret = fMinimizer->Contour(ivarX, ivarY, ncp, x, y );
   if (!ret) {
      coutE(Minimization) << "LikelihoodInterval - Error finding contour for parameters " << paramX.GetName() << " and " << paramY.GetName()  << std::endl;
      return 0;
   }
   if (int(ncp) < npoints) {
      coutW(Minimization) << "LikelihoodInterval -Warning - Less points calculated in contours np = " << ncp << " / " << npoints << std::endl;
   }

   return ncp;
 }
