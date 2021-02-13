// @(#)root/mathcore:$Id$
// Author: L. Moneta Fri Aug 17 14:29:24 2007

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2007  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class LogLikelihoodFCN

#ifndef ROOT_Fit_LogLikelihoodFCN
#define ROOT_Fit_LogLikelihoodFCN

#include "ROOT/EExecutionPolicy.hxx"
#include "Fit/BasicFCN.h"
#include "Fit/FitUtil.h"
#include "Fit/UnBinData.h"
#include "Math/IParamFunction.h"

#include <memory>
#include <vector>

namespace ROOT {

   namespace Fit {


//___________________________________________________________________________________
/**
   LogLikelihoodFCN class
   for likelihood fits

   it is template to distinguish gradient and non-gradient case

   @ingroup  FitMethodFunc
*/
template<class DerivFunType,class ModelFunType = ROOT::Math::IParamMultiFunction>
class LogLikelihoodFCN : public BasicFCN<DerivFunType,ModelFunType,UnBinData>  {

public:

   typedef typename ModelFunType::BackendType T;
   typedef  BasicFCN<DerivFunType,ModelFunType,UnBinData> BaseFCN;

   typedef  ::ROOT::Math::BasicFitMethodFunction<DerivFunType> BaseObjFunction;
   typedef typename  BaseObjFunction::BaseFunction BaseFunction;

   typedef  ::ROOT::Math::IParamMultiFunctionTempl<T> IModelFunction;
   typedef typename BaseObjFunction::Type_t Type_t;


   /**
      Constructor from unbin data set and model function (pdf)
   */
   LogLikelihoodFCN (const std::shared_ptr<UnBinData> & data, const std::shared_ptr<IModelFunction> & func, int weight = 0, bool extended = false, const ::ROOT::EExecutionPolicy &executionPolicy = ::ROOT::EExecutionPolicy::kSequential) :
      BaseFCN( data, func),
      fIsExtended(extended),
      fWeight(weight),
      fNEffPoints(0),
      fGrad ( std::vector<double> ( func->NPar() ) ),
      fExecutionPolicy(executionPolicy)
   {}

      /**
      Constructor from unbin data set and model function (pdf) for object managed by users
   */
   LogLikelihoodFCN (const UnBinData & data, const IModelFunction & func, int weight = 0, bool extended = false, const ::ROOT::EExecutionPolicy &executionPolicy = ::ROOT::EExecutionPolicy::kSequential) :
      BaseFCN(std::shared_ptr<UnBinData>(const_cast<UnBinData*>(&data), DummyDeleter<UnBinData>()), std::shared_ptr<IModelFunction>(dynamic_cast<IModelFunction*>(func.Clone() ) ) ),
      fIsExtended(extended),
      fWeight(weight),
      fNEffPoints(0),
      fGrad ( std::vector<double> ( func.NPar() ) ),
      fExecutionPolicy(executionPolicy)
   {}

   /**
      Destructor (no operations)
   */
   virtual ~LogLikelihoodFCN () {}

   /**
      Copy constructor
   */
   LogLikelihoodFCN(const LogLikelihoodFCN & f) :
      BaseFCN(f.DataPtr(), f.ModelFunctionPtr() ),
      fIsExtended(f.fIsExtended ),
      fWeight( f.fWeight ),
      fNEffPoints( f.fNEffPoints ),
      fGrad( f.fGrad),
      fExecutionPolicy(f.fExecutionPolicy)
   {  }


   /**
      Assignment operator
   */
   LogLikelihoodFCN & operator = (const LogLikelihoodFCN & rhs) {
      SetData(rhs.DataPtr() );
      SetModelFunction(rhs.ModelFunctionPtr() );
      fNEffPoints = rhs.fNEffPoints;
      fGrad = rhs.fGrad;
      fIsExtended = rhs.fIsExtended;
      fWeight = rhs.fWeight;
      fExecutionPolicy = rhs.fExecutionPolicy;
      return *this;
   }


   /// clone the function (need to return Base for Windows)
   virtual BaseFunction * Clone() const { return  new LogLikelihoodFCN(*this); }


   //using BaseObjFunction::operator();

   // effective points used in the fit
   virtual unsigned int NFitPoints() const { return fNEffPoints; }

   /// i-th likelihood contribution and its gradient
   virtual double DataElement(const double * x, unsigned int i, double * g) const {
      if (i==0) this->UpdateNCalls();
      return FitUtil::EvaluatePdf(BaseFCN::ModelFunction(), BaseFCN::Data(), x, i, g);
   }

   // need to be virtual to be instantited
   virtual void Gradient(const double *x, double *g) const {
      // evaluate the chi2 gradient
      FitUtil::Evaluate<typename BaseFCN::T>::EvalLogLGradient(BaseFCN::ModelFunction(), BaseFCN::Data(), x, g,
                                                               fNEffPoints, fExecutionPolicy);
   }

   /// get type of fit method function
   virtual  typename BaseObjFunction::Type_t Type() const { return BaseObjFunction::kLogLikelihood; }


   // Use sum of the weight squared in evaluating the likelihood
   // (this is needed for calculating the errors)
   void UseSumOfWeightSquare(bool on = true) {
      if (fWeight == 0) return; // do nothing if it was not weighted
      if (on) fWeight = 2;
      else fWeight = 1;
   }



protected:


private:

   /**
      Evaluation of the  function (required by interface)
    */
   virtual double DoEval (const double * x) const {
      this->UpdateNCalls();
      return FitUtil::Evaluate<T>::EvalLogL(BaseFCN::ModelFunction(), BaseFCN::Data(), x, fWeight, fIsExtended, fNEffPoints, fExecutionPolicy);
   }

   // for derivatives
   virtual double  DoDerivative(const double * x, unsigned int icoord ) const {
      Gradient(x, &fGrad[0]);
      return fGrad[icoord];
   }


      //data member
   bool fIsExtended;  // flag for indicating if likelihood is extended
   int  fWeight;  // flag to indicate if needs to evaluate using weight or weight squared (default weight = 0)


   mutable unsigned int fNEffPoints;  // number of effective points used in the fit

   mutable std::vector<double> fGrad; // for derivatives

   ::ROOT::EExecutionPolicy fExecutionPolicy; // Execution policy
};
      // define useful typedef's
      // using LogLikelihoodFunction_v = LogLikelihoodFCN<ROOT::Math::IMultiGenFunction, ROOT::Math::IParametricFunctionMultiDimTempl<T>>;
      typedef LogLikelihoodFCN<ROOT::Math::IMultiGenFunction, ROOT::Math::IParamMultiFunction>  LogLikelihoodFunction;
      typedef LogLikelihoodFCN<ROOT::Math::IMultiGradFunction, ROOT::Math::IParamMultiFunction> LogLikelihoodGradFunction;

   } // end namespace Fit

} // end namespace ROOT


#endif /* ROOT_Fit_LogLikelihoodFCN */
