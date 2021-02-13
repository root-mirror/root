// @(#)root/mathmore:$Id$
// Author: L. Moneta Oct 2012

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for class BasicMinimizer

#include "Math/BasicMinimizer.h"

#include "Math/IFunction.h"

#include "Math/IFunctionfwd.h"

#include "Math/IParamFunctionfwd.h"

#include "Math/FitMethodFunction.h"

#include "Math/MinimTransformFunction.h"

#include "Math/Error.h"

#include "Fit/ParameterSettings.h"

#include <cassert>

#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <cctype>   // need to use c version of tolower defined here
#include <limits>

namespace ROOT {

   namespace Math {


BasicMinimizer::BasicMinimizer( ) :
   fDim(0),
   fObjFunc(0),
   fMinVal(0)
{
   fValues.reserve(10);
   fNames.reserve(10);
   fSteps.reserve(10);

   int niter = ROOT::Math::MinimizerOptions::DefaultMaxIterations();
   if (niter <=0 ) niter = 1000;
   SetMaxIterations(niter);
   SetPrintLevel(ROOT::Math::MinimizerOptions::DefaultPrintLevel());
}


BasicMinimizer::~BasicMinimizer () {
   if (fObjFunc) delete fObjFunc;
}

bool BasicMinimizer::SetVariable(unsigned int ivar, const std::string & name, double val, double step) {
   // set variable in minimizer - support only free variables
   // no transformation implemented - so far
   if (ivar > fValues.size() ) return false;
   if (ivar == fValues.size() ) {
      fValues.push_back(val);
      fNames.push_back(name);
      fSteps.push_back(step);
      fVarTypes.push_back(kDefault);
   }
   else {
      fValues[ivar] = val;
      fNames[ivar] = name;
      fSteps[ivar] = step;
      fVarTypes[ivar] = kDefault;

      // remove bounds if needed
      std::map<unsigned  int, std::pair<double, double> >::iterator iter = fBounds.find(ivar);
      if ( iter !=  fBounds.end() ) fBounds.erase (iter);

   }

   return true;
}

bool BasicMinimizer::SetLowerLimitedVariable(unsigned int ivar, const std::string & name, double val, double step, double lower) {
   // set lower limited variable
   bool ret = SetVariable(ivar, name, val, step);
   if (!ret) return false;
   const double upper = std::numeric_limits<double>::infinity();
   fBounds[ivar] = std::make_pair( lower, upper);
   fVarTypes[ivar] = kLowBound;
   return true;
}
bool BasicMinimizer::SetUpperLimitedVariable(unsigned int ivar, const std::string & name, double val, double step, double upper ) {
   // set upper limited variable
   bool ret = SetVariable(ivar, name, val, step);
   if (!ret) return false;
   const double lower = -std::numeric_limits<double>::infinity();
   fBounds[ivar] = std::make_pair( lower, upper);
   fVarTypes[ivar] = kUpBound;
   return true;
}

bool BasicMinimizer::SetLimitedVariable(unsigned int ivar, const std::string & name, double val, double step, double lower, double upper) {
   // set double bounded variable
   bool ret = SetVariable(ivar, name, val, step);
   if (!ret) return false;
   fBounds[ivar] = std::make_pair( lower, upper);
   fVarTypes[ivar] = kBounds;
   return true;
}

bool BasicMinimizer::SetFixedVariable(unsigned int ivar , const std::string & name , double val ) {
   /// set fixed variable
   bool ret = SetVariable(ivar, name, val, 0.);
   if (!ret) return false;
   fVarTypes[ivar] = kFix;
   return true;
}


bool BasicMinimizer::SetVariableValue(unsigned int ivar, double val) {
   // set variable value in minimizer
   // no change to transformation or variable status
   if (ivar >= fValues.size() ) return false;
   fValues[ivar] = val;
   return true;
}

bool BasicMinimizer::SetVariableValues( const double * x) {
   // set all variable values in minimizer
   if (x == 0) return false;
   std::copy(x,x+fValues.size(), fValues.begin() );
   return true;
}

bool BasicMinimizer::SetVariableStepSize(unsigned int ivar, double step) {
   // set step size
   if (ivar > fValues.size() ) return false;
   fSteps[ivar] = step;
   return true;
}

bool BasicMinimizer::SetVariableLowerLimit(unsigned int ivar, double lower) {
   // set variable lower limit
   double upper =  (fBounds.count(ivar)) ? fBounds[ivar].second : std::numeric_limits<double>::infinity();
   return SetVariableLimits(ivar, lower, upper);
}

bool BasicMinimizer::SetVariableUpperLimit(unsigned int ivar, double upper) {
   // set variable upper limit
   double lower =  (fBounds.count(ivar)) ? fBounds[ivar].first : - std::numeric_limits<double>::infinity();
   return SetVariableLimits(ivar, lower, upper);
}

bool BasicMinimizer::SetVariableLimits(unsigned int ivar, double lower, double upper) {
   // set variable limits (remove limits if lower >= upper)
   if (ivar > fVarTypes.size() ) return false;
   // if limits do not exists add them or update
   fBounds[ivar] = std::make_pair( lower, upper);
   if (lower > upper || (lower == - std::numeric_limits<double>::infinity() &&
                         upper ==   std::numeric_limits<double>::infinity() ) ) {
      fBounds.erase(ivar);
      fVarTypes[ivar] = kDefault;
   }
   else if (lower == upper)
      FixVariable(ivar);
   else {
      if (lower == - std::numeric_limits<double>::infinity() )
         fVarTypes[ivar] = kLowBound;
      else if (upper == std::numeric_limits<double>::infinity() )
         fVarTypes[ivar] = kUpBound;
      else
         fVarTypes[ivar] = kBounds;
   }
   return true;
}

bool BasicMinimizer::FixVariable(unsigned int ivar) {
   // fix variable
   if (ivar >= fVarTypes.size() ) return false;
   fVarTypes[ivar] = kFix;
   return true;
}

bool BasicMinimizer::ReleaseVariable(unsigned int ivar) {
   // fix variable
   if (ivar >= fVarTypes.size() ) return false;
   if (fBounds.count(ivar) == 0)  {
      fVarTypes[ivar] = kDefault;
      return true;
   }
   if (fBounds[ivar].first == - std::numeric_limits<double>::infinity() )
      fVarTypes[ivar] = kLowBound;
   else if (fBounds[ivar].second == std::numeric_limits<double>::infinity() )
      fVarTypes[ivar] = kUpBound;
   else
      fVarTypes[ivar] = kBounds;

   return true;
}

bool BasicMinimizer::IsFixedVariable(unsigned int ivar) const {
   if (ivar >= fVarTypes.size() ) return false;
   return (fVarTypes[ivar] == kFix ) ;
}

bool BasicMinimizer::GetVariableSettings(unsigned int ivar, ROOT::Fit::ParameterSettings & varObj) const {
   if (ivar >= fValues.size() ) return false;
   assert(fValues.size() == fNames.size() && fValues.size() == fVarTypes.size() );
   varObj.Set(fNames[ivar],fValues[ivar],fSteps[ivar]);
   std::map< unsigned int , std::pair< double, double> >::const_iterator itr = fBounds.find(ivar);
   if (itr != fBounds.end() )  {
      double lower = (itr->second).first;
      double upper = (itr->second).second;
      if (fVarTypes[ivar] == kLowBound) varObj.SetLowerLimit( lower );
      if (fVarTypes[ivar] == kUpBound) varObj.SetUpperLimit( upper );
      else varObj.SetLimits( lower,upper);
   }
   if (fVarTypes[ivar] == kFix ) varObj.Fix();
   return true;
}

std::string BasicMinimizer::VariableName(unsigned int ivar) const {
   if (ivar >= fNames.size() ) return "";
   return fNames[ivar];
}

int BasicMinimizer::VariableIndex(const std::string & name) const {
   std::vector<std::string>::const_iterator itr = std::find( fNames.begin(), fNames.end(), name);
   if (itr == fNames.end() ) return -1;
   return itr - fNames.begin();
}



void BasicMinimizer::SetFunction(const ROOT::Math::IMultiGenFunction & func) {
   // set the function to minimizer
   fObjFunc = func.Clone();
   fDim = fObjFunc->NDim();
}

void BasicMinimizer::SetFunction(const ROOT::Math::IMultiGradFunction & func) {
   // set the function to minimize
   fObjFunc = dynamic_cast<const ROOT::Math::IMultiGradFunction *>( func.Clone());
   assert(fObjFunc != 0);
   fDim = fObjFunc->NDim();
}


bool BasicMinimizer::CheckDimension() const {
   unsigned int npar = fValues.size();
   if (npar == 0 || npar < fDim  ) {
      MATH_ERROR_MSGVAL("BasicMinimizer::CheckDimension","Wrong number of parameters",npar);
      return false;
   }
   return true;
}

bool BasicMinimizer::CheckObjFunction() const {
   if (fObjFunc == 0) {
      MATH_ERROR_MSG("BasicMinimizer::CheckFunction","Function has not been set");
      return false;
   }
   return true;
}


MinimTransformFunction * BasicMinimizer::CreateTransformation(std::vector<double> & startValues, const ROOT::Math::IMultiGradFunction * func) {

   bool doTransform = (fBounds.size() > 0);
   unsigned int ivar = 0;
   while (!doTransform && ivar < fVarTypes.size() ) {
      doTransform = (fVarTypes[ivar++] != kDefault );
   }

   startValues = std::vector<double>(fValues.begin(), fValues.end() );

   MinimTransformFunction * trFunc  = 0;

   // in case of transformation wrap objective function in a new transformation function
   // and transform from external variables  to internals one
   // Transformations are supported only for gradient function
   const IMultiGradFunction * gradObjFunc = (func) ? func : dynamic_cast<const IMultiGradFunction *>(fObjFunc);
   doTransform &= (gradObjFunc != 0);

   if (doTransform)   {
      // minim transform function manages the passed function pointer (gradObjFunc)
      trFunc =  new MinimTransformFunction ( gradObjFunc, fVarTypes, fValues, fBounds );
      // transform from external to internal
      trFunc->InvTransformation(&fValues.front(), &startValues[0]);
      // size can be different since internal parameter can have smaller size
      // if there are fixed parameters
      startValues.resize( trFunc->NDim() );
      // no need to save fObjFunc since trFunc will manage it
      fObjFunc = trFunc;
   }
   else {
      if (func) fObjFunc = func;  // to manege the passed function object
   }

//    std::cout << " f has transform " << doTransform << "  " << fBounds.size() << "   " << startValues.size() <<  " ndim " << fObjFunc->NDim() << std::endl;   std::cout << "InitialValues external : ";
//    for (int i = 0; i < fValues.size(); ++i) std::cout << fValues[i] << "  ";
//    std::cout << "\n";
//    std::cout << "InitialValues internal : ";
//    for (int i = 0; i < startValues.size(); ++i) std::cout << startValues[i] << "  ";
//    std::cout << "\n";


   return trFunc;
}

bool BasicMinimizer::Minimize() {

   // do nothing
   return false;
}

void BasicMinimizer::SetFinalValues(const double * x) {
   // check to see if a transformation need to be applied
   const MinimTransformFunction * trFunc = TransformFunction();
   if (trFunc) {
      assert(fValues.size() >= trFunc->NTot() );
      trFunc->Transformation(x, &fValues[0]);
   }
   else {
      // case of no transformation applied
      assert( fValues.size() >= NDim() );
      std::copy(x, x + NDim(),  fValues.begin() );
   }
}

void BasicMinimizer::PrintResult() const {
   int pr = std::cout.precision(18);
   std::cout << "FVAL         = " << fMinVal << std::endl;
   std::cout.precision(pr);
//      std::cout << "Edm   = " << fState.Edm() << std::endl;
   std::cout << "Niterations  = " << NIterations() << std::endl;
   unsigned int ncalls = NCalls();
   if (ncalls) std::cout << "NCalls     = " << ncalls << std::endl;
   for (unsigned int i = 0; i < fDim; ++i)
      std::cout << fNames[i] << "\t  = " << fValues[i] << std::endl;
}

const ROOT::Math::IMultiGradFunction * BasicMinimizer::GradObjFunction() const {
      return  dynamic_cast<const ROOT::Math::IMultiGradFunction *>(fObjFunc);
}

const MinimTransformFunction * BasicMinimizer::TransformFunction() const {
   return dynamic_cast<const MinimTransformFunction *>(fObjFunc);
}

unsigned int BasicMinimizer::NFree() const {
   // number of free variables
   unsigned int nfree = fValues.size();
   for (unsigned int i = 0; i < fVarTypes.size(); ++i)
      if (fVarTypes[i] == kFix) nfree--;
   return nfree;
}


   } // end namespace Math

} // end namespace ROOT

