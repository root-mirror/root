// @(#)root/mathcore:$Id$
// Author: L. Moneta Tue Sep  5 09:13:32 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class Chi2FCN

#ifndef ROOT_Fit_Chi2FCN
#define ROOT_Fit_Chi2FCN

#include "ROOT/EExecutionPolicy.hxx"
#include "Fit/BasicFCN.h"
#include "Fit/BinData.h"
#include "Fit/FitUtil.h"
#include "Math/IFunction.h"
#include "Math/IFunctionfwd.h"
#include "Math/IParamFunction.h"

#include <memory>
#include <vector>

/**
@defgroup FitMethodFunc Fit Method Classes

Classes describing Fit Method functions
@ingroup Fit
*/

namespace ROOT {


   namespace Fit {

//___________________________________________________________________________________
/**
   Chi2FCN class for binnned fits using the least square methods

   @ingroup  FitMethodFunc
*/
template<class DerivFunType, class ModelFunType = ROOT::Math::IParamMultiFunction>
class Chi2FCN : public BasicFCN<DerivFunType, ModelFunType, BinData> {

public:

   typedef typename ModelFunType::BackendType T;
   typedef  BasicFCN<DerivFunType, ModelFunType, BinData> BaseFCN;

   typedef  ::ROOT::Math::BasicFitMethodFunction<DerivFunType> BaseObjFunction;
   typedef typename  BaseObjFunction::BaseFunction BaseFunction;

   //typedef  typename ::ROOT::Math::ParamFunctionTrait<FunType>::PFType IModelFunction;
   typedef  ::ROOT::Math::IParamMultiFunctionTempl<T> IModelFunction;
   typedef typename BaseObjFunction::Type_t Type_t;

   /**
      Constructor from data set (binned ) and model function
   */
   Chi2FCN (const std::shared_ptr<BinData> & data, const std::shared_ptr<IModelFunction> & func, const ::ROOT::EExecutionPolicy &executionPolicy = ::ROOT::EExecutionPolicy::kSequential) :
      BaseFCN( data, func),
      fNEffPoints(0),
      fGrad ( std::vector<double> ( func->NPar() ) ),
      fExecutionPolicy(executionPolicy)
   { }

   /**
      Same Constructor from data set (binned ) and model function but now managed by the user
      we clone the function but not the data
   */
   Chi2FCN ( const BinData & data, const IModelFunction & func, const ::ROOT::EExecutionPolicy &executionPolicy = ::ROOT::EExecutionPolicy::kSequential) :
      BaseFCN(std::shared_ptr<BinData>(const_cast<BinData*>(&data), DummyDeleter<BinData>()), std::shared_ptr<IModelFunction>(dynamic_cast<IModelFunction*>(func.Clone() ) ) ),
      fNEffPoints(0),
      fGrad ( std::vector<double> ( func.NPar() ) ),
      fExecutionPolicy(executionPolicy)
   { }

   /**
      Destructor (no operations)
   */
   virtual ~Chi2FCN ()  {}
   /**
      Copy constructor
   */
   Chi2FCN(const Chi2FCN & f) :
      BaseFCN(f.DataPtr(), f.ModelFunctionPtr() ),
      fNEffPoints( f.fNEffPoints ),
      fGrad( f.fGrad),
      fExecutionPolicy(f.fExecutionPolicy)
   {  }

   /**
      Assignment operator
   */
   Chi2FCN & operator = (const Chi2FCN & rhs) {
      SetData(rhs.DataPtr() );
      SetModelFunction(rhs.ModelFunctionPtr() );
      fNEffPoints = rhs.fNEffPoints;
      fGrad = rhs.fGrad;
   }

   /*
      clone the function
    */
   virtual BaseFunction * Clone() const {
      return new Chi2FCN(*this);
   }



   using BaseObjFunction::operator();


   /// i-th chi-square residual
   virtual double DataElement(const double *x, unsigned int i, double *g) const {
      if (i==0) this->UpdateNCalls();
      return FitUtil::Evaluate<T>::EvalChi2Residual(BaseFCN::ModelFunction(), BaseFCN::Data(), x, i, g);
   }

   // need to be virtual to be instantiated
   virtual void Gradient(const double *x, double *g) const {
      // evaluate the chi2 gradient
      FitUtil::Evaluate<T>::EvalChi2Gradient(BaseFCN::ModelFunction(), BaseFCN::Data(), x, g, fNEffPoints,
                                             fExecutionPolicy);
   }

   /// get type of fit method function
   virtual  typename BaseObjFunction::Type_t Type() const { return BaseObjFunction::kLeastSquare; }


protected:

   /// set number of fit points (need to be called in const methods, make it const)
   virtual void SetNFitPoints(unsigned int n) const { fNEffPoints = n; }

private:

   /**
      Evaluation of the  function (required by interface)
    */
   virtual double DoEval (const double * x) const {
      this->UpdateNCalls();
      if (BaseFCN::Data().HaveCoordErrors() || BaseFCN::Data().HaveAsymErrors())
         return FitUtil::Evaluate<T>::EvalChi2Effective(BaseFCN::ModelFunction(), BaseFCN::Data(), x, fNEffPoints);
      else
         return FitUtil::Evaluate<T>::EvalChi2(BaseFCN::ModelFunction(), BaseFCN::Data(), x, fNEffPoints, fExecutionPolicy);
   }

   // for derivatives
   virtual double  DoDerivative(const double * x, unsigned int icoord ) const {
      Gradient(x, fGrad.data());
      return fGrad[icoord];
   }


   mutable unsigned int fNEffPoints;  // number of effective points used in the fit

   mutable std::vector<double> fGrad; // for derivatives
   ::ROOT::EExecutionPolicy fExecutionPolicy;

};

      // define useful typedef's
      typedef Chi2FCN<ROOT::Math::IMultiGenFunction,ROOT::Math::IParamMultiFunction> Chi2Function;
      typedef Chi2FCN<ROOT::Math::IMultiGradFunction, ROOT::Math::IParamMultiFunction> Chi2GradFunction;


   } // end namespace Fit

} // end namespace ROOT


#endif /* ROOT_Fit_Chi2FCN */
