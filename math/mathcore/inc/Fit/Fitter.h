// @(#)root/mathcore:$Id$
// Author: L. Moneta Wed Aug 30 11:05:19 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class Fitter

#ifndef ROOT_Fit_Fitter
#define ROOT_Fit_Fitter

/**
@defgroup Fit Fitting and Parameter Estimation

Classes used for fitting (regression analysis) and estimation of parameter values given a data sample.

@ingroup MathCore

*/

#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "Fit/FitConfig.h"
#include "ROOT/EExecutionPolicy.hxx"
#include "Fit/FitResult.h"
#include "Math/IParamFunction.h"
#include <memory>

namespace ROOT {


   namespace Math {
      class Minimizer;

      // should maybe put this in a FitMethodFunctionfwd file
      template<class FunctionType> class BasicFitMethodFunction;

      // define the normal and gradient function
      typedef BasicFitMethodFunction<ROOT::Math::IMultiGenFunction>  FitMethodFunction;
      typedef BasicFitMethodFunction<ROOT::Math::IMultiGradFunction> FitMethodGradFunction;

   }

   /**
      Namespace for the fitting classes
      @ingroup Fit
    */

   namespace Fit {

/**
   @defgroup FitMain User Fitting classes

   Main Classes used for fitting a given data set
   @ingroup Fit
*/


//___________________________________________________________________________________
/**
   Fitter class, entry point for performing all type of fits.
   Fits are performed using the generic ROOT::Fit::Fitter::Fit method.
   The inputs are the data points and a model function (using a ROOT::Math::IParamFunction)
   The result of the fit is returned and kept internally in the  ROOT::Fit::FitResult class.
   The configuration of the fit (parameters, options, etc...) are specified in the
   ROOT::Math::FitConfig class.
   After fitting the config of the fit will be modified to have the new values the resulting
   parameter of the fit with step sizes equal to the errors. FitConfig can be preserved with
   initial parameters by calling FitConfig.SetUpdateAfterFit(false);

   @ingroup FitMain
*/
class Fitter {

public:

   typedef ROOT::Math::IParamMultiFunction                 IModelFunction;
   template <class T>
   using IModelFunctionTempl =                             ROOT::Math::IParamMultiFunctionTempl<T>;
#ifdef R__HAS_VECCORE
   typedef ROOT::Math::IParametricFunctionMultiDimTempl<ROOT::Double_v>  IModelFunction_v;
   typedef ROOT::Math::IParamMultiGradFunctionTempl<ROOT::Double_v> IGradModelFunction_v;
#else
   typedef ROOT::Math::IParamMultiFunction                 IModelFunction_v;
   typedef ROOT::Math::IParamMultiGradFunction IGradModelFunction_v;
#endif
   typedef ROOT::Math::IParamMultiGradFunction             IGradModelFunction;
   typedef ROOT::Math::IParamFunction                      IModel1DFunction;
   typedef ROOT::Math::IParamGradFunction                  IGradModel1DFunction;

   typedef ROOT::Math::IMultiGenFunction BaseFunc;
   typedef ROOT::Math::IMultiGradFunction BaseGradFunc;


   /**
      Default constructor
   */
   Fitter ();

   /**
      Constructor from a result
   */
   Fitter (const std::shared_ptr<FitResult> & result);


   /**
      Destructor
   */
   ~Fitter ();

private:

   /**
      Copy constructor (disabled, class is not copyable)
   */
   Fitter(const Fitter &);

   /**
      Assignment operator (disabled, class is not copyable)
   */
   Fitter & operator = (const Fitter & rhs);


public:

   /**
       fit a data set using any  generic model  function
       If data set is binned a least square fit is performed
       If data set is unbinned a maximum likelihood fit (not extended) is done
       Pre-requisite on the function:
       it must implement the 1D or multidimensional parametric function interface
   */
   template <class Data, class Function,
             class cond = typename std::enable_if<!(std::is_same<Function, ROOT::EExecutionPolicy>::value ||
                                                    std::is_same<Function, int>::value),
                                                  Function>::type>
   bool Fit(const Data &data, const Function &func,
            const ROOT::EExecutionPolicy &executionPolicy = ROOT::EExecutionPolicy::kSequential)
   {
      SetFunction(func);
      return Fit(data, executionPolicy);
   }

   /**
       Fit a binned data set using a least square fit (default method)
   */
   bool Fit(const BinData & data, const ROOT::EExecutionPolicy &executionPolicy = ROOT::EExecutionPolicy::kSequential) {
      SetData(data);
      return DoLeastSquareFit(executionPolicy);
   }
   bool Fit(const std::shared_ptr<BinData> & data, const ROOT::EExecutionPolicy &executionPolicy = ROOT::EExecutionPolicy::kSequential) {
      SetData(data);
      return DoLeastSquareFit(executionPolicy);
   }

   /**
       Fit a binned data set using a least square fit
   */
   bool LeastSquareFit(const BinData & data) {
      return Fit(data);
   }

   /**
       fit an unbinned data set using loglikelihood method
   */
   bool Fit(const UnBinData & data, bool extended = false, const ROOT::EExecutionPolicy &executionPolicy = ROOT::EExecutionPolicy::kSequential) {
      SetData(data);
      return DoUnbinnedLikelihoodFit(extended, executionPolicy);
   }

   /**
      Binned Likelihood fit. Default is extended
    */
   bool LikelihoodFit(const BinData &data, bool extended = true,
                      const ROOT::EExecutionPolicy &executionPolicy = ROOT::EExecutionPolicy::kSequential) {
      SetData(data);
      return DoBinnedLikelihoodFit(extended, executionPolicy);
   }

   bool LikelihoodFit(const std::shared_ptr<BinData> &data, bool extended = true,
                      const ROOT::EExecutionPolicy &executionPolicy = ROOT::EExecutionPolicy::kSequential) {
      SetData(data);
      return DoBinnedLikelihoodFit(extended, executionPolicy);
   }
   /**
      Unbinned Likelihood fit. Default is not extended
    */
   bool LikelihoodFit(const UnBinData & data, bool extended = false, const ROOT::EExecutionPolicy &executionPolicy = ROOT::EExecutionPolicy::kSequential) {
      SetData(data);
      return DoUnbinnedLikelihoodFit(extended, executionPolicy);
   }
   bool LikelihoodFit(const std::shared_ptr<UnBinData> & data, bool extended = false, const ROOT::EExecutionPolicy &executionPolicy = ROOT::EExecutionPolicy::kSequential) {
      SetData(data);
      return DoUnbinnedLikelihoodFit(extended, executionPolicy);
   }


   /**
       fit a data set using any  generic model  function
       Pre-requisite on the function:
   */
   template < class Data , class Function>
   bool LikelihoodFit( const Data & data, const Function & func, bool extended) {
      SetFunction(func);
      return LikelihoodFit(data, extended);
   }

   /**
      do a linear fit on a set of bin-data
    */
   bool LinearFit(const BinData & data) {
      SetData(data);
      return DoLinearFit();
   }
   bool LinearFit(const std::shared_ptr<BinData> & data) {
      SetData(data);
      return DoLinearFit();
   }


   /**
      Fit using the a generic FCN function as a C++ callable object implementing
      double () (const double *)
      Note that the function dimension (i.e. the number of parameter) is needed in this case
      For the options see documentation for following methods FitFCN(IMultiGenFunction & fcn,..)
    */
   template <class Function>
   bool FitFCN(unsigned int npar, Function  & fcn, const double * params = 0, unsigned int dataSize = 0, bool chi2fit = false);

   /**
      Set a generic FCN function as a C++ callable object implementing
      double () (const double *)
      Note that the function dimension (i.e. the number of parameter) is needed in this case
      For the options see documentation for following methods FitFCN(IMultiGenFunction & fcn,..)
    */
   template <class Function>
   bool SetFCN(unsigned int npar, Function  & fcn, const double * params = 0, unsigned int dataSize = 0, bool chi2fit = false);

   /**
      Fit using the given FCN function represented by a multi-dimensional function interface
      (ROOT::Math::IMultiGenFunction).
      Give optionally the initial arameter values, data size to have the fit Ndf correctly
      set in the FitResult and flag specifying if it is a chi2 fit.
      Note that if the parameters values are not given (params=0) the
      current parameter settings are used. The parameter settings can be created before
      by using the FitConfig::SetParamsSetting. If they have not been created they are created
      automatically when the params pointer is not zero.
      Note that passing a params != 0 will set the parameter settings to the new value AND also the
      step sizes to some pre-defined value (stepsize = 0.3 * abs(parameter_value) )
    */
   bool FitFCN(const ROOT::Math::IMultiGenFunction &fcn, const double *params = 0, unsigned int dataSize = 0, bool chi2fit = false);

   /**
       Fit using a FitMethodFunction interface. Same as method above, but now extra information
       can be taken from the function class
   */
   bool FitFCN(const ROOT::Math::FitMethodFunction & fcn, const double * params = 0);

   /**
      Set the FCN function represented by a multi-dimensional function interface
      (ROOT::Math::IMultiGenFunction) and optionally the initial parameters
      See also note above for the initial parameters for FitFCN
    */
   bool SetFCN(const ROOT::Math::IMultiGenFunction &fcn, const double *params = 0, unsigned int dataSize = 0, bool chi2fit = false);

   /**
      Set the FCN function represented by a multi-dimensional function interface
     (ROOT::Math::IMultiGenFunction) and optionally the initial parameters
      See also note above for the initial parameters for FitFCN
      With this interface we pass in addition a ModelFunction that will be attached to the FitResult and
      used to compute confidence interval of the fit
   */
   bool SetFCN(const ROOT::Math::IMultiGenFunction &fcn, const IModelFunction & func, const double *params = 0,
               unsigned int dataSize = 0, bool chi2fit = false);

   /**
       Set the objective function (FCN)  using a FitMethodFunction interface.
       Same as method above, but now extra information can be taken from the function class
   */
   bool SetFCN(const ROOT::Math::FitMethodFunction & fcn, const double * params = 0);

   /**
      Fit using the given FCN function representing a multi-dimensional gradient function
      interface (ROOT::Math::IMultiGradFunction). In this case the minimizer will use the
      gradient information provided by the function.
      For the options same consideration as in the previous method
    */
   bool FitFCN(const ROOT::Math::IMultiGradFunction &fcn, const double *params = 0, unsigned int dataSize = 0, bool chi2fit = false);

   /**
       Fit using a FitMethodGradFunction interface. Same as method above, but now extra information
       can be taken from the function class
   */
   bool FitFCN(const ROOT::Math::FitMethodGradFunction & fcn, const double * params = 0);

   /**
      Set the FCN function represented by a multi-dimensional gradient function interface
      (ROOT::Math::IMultiGradFunction) and optionally the initial parameters
      See also note above for the initial parameters for FitFCN
    */
   bool SetFCN(const ROOT::Math::IMultiGradFunction &fcn, const double *params = 0, unsigned int dataSize = 0, bool chi2fit = false);

   /**
      Set the FCN function represented by a multi-dimensional gradient function interface
     (ROOT::Math::IMultiGradFunction) and optionally the initial parameters
      See also note above for the initial parameters for FitFCN
      With this interface we pass in addition a ModelFunction that will be attached to the FitResult and
      used to compute confidence interval of the fit
   */
   bool SetFCN(const ROOT::Math::IMultiGradFunction &fcn, const IModelFunction &func, const double *params = 0,
               unsigned int dataSize = 0, bool chi2fit = false);

   /**
       Set the objective function (FCN)  using a FitMethodGradFunction interface.
       Same as method above, but now extra information can be taken from the function class
   */
   bool SetFCN(const ROOT::Math::FitMethodGradFunction & fcn, const double * params = 0);


   /**
      fit using user provided FCN with Minuit-like interface
      If npar = 0 it is assumed that the parameters are specified in the parameter settings created before
      For the options same consideration as in the previous method
    */
   typedef  void (* MinuitFCN_t )(int &npar, double *gin, double &f, double *u, int flag);
   bool FitFCN( MinuitFCN_t fcn, int npar = 0, const double * params = 0, unsigned int dataSize = 0, bool chi2fit = false);

   /**
      set objective function using user provided FCN with Minuit-like interface
      If npar = 0 it is assumed that the parameters are specified in the parameter settings created before
      For the options same consideration as in the previous method
    */
   bool SetFCN( MinuitFCN_t fcn, int npar = 0, const double * params = 0, unsigned int dataSize = 0, bool chi2fit = false);

   /**
      Perform a fit with the previously set FCN function. Require SetFCN before
    */
   bool FitFCN();

   /**
      Perform a simple FCN evaluation. FitResult will be modified and contain  the value of the FCN
    */
   bool EvalFCN();



   /**
       Set the fitted function (model function) from a parametric function interface
   */
   void  SetFunction(const IModelFunction & func, bool useGradient = false);

   /**
       Set the fitted function (model function) from a vectorized parametric function interface
   */
#ifdef R__HAS_VECCORE
   template <class NotCompileIfScalarBackend = std::enable_if<!(std::is_same<double, ROOT::Double_v>::value)>>
   void SetFunction(const IModelFunction_v &func, bool useGradient = false);

   template <class NotCompileIfScalarBackend = std::enable_if<!(std::is_same<double, ROOT::Double_v>::value)>>
   void SetFunction(const IGradModelFunction_v &func, bool useGradient = true);
#endif
   /**
      Set the fitted function from a parametric 1D function interface
    */
   void  SetFunction(const IModel1DFunction & func, bool useGradient = false);

   /**
       Set the fitted function (model function) from a parametric gradient function interface
   */
   void  SetFunction(const IGradModelFunction & func, bool useGradient = true);
   /**
      Set the fitted function from 1D gradient parametric function interface
    */
   void  SetFunction(const IGradModel1DFunction & func, bool useGradient = true);


   /**
      get fit result
   */
   const FitResult & Result() const {
      assert( fResult.get() );
      return *fResult;
   }


   /**
      perform an error analysis on the result using the Hessian
      Errors are obtaied from the inverse of the Hessian matrix
      To be called only after fitting and when a minimizer supporting the Hessian calculations is used
      otherwise an error (false) is returned.
      A new  FitResult with the Hessian result will be produced
    */
   bool CalculateHessErrors();

   /**
      perform an error analysis on the result using MINOS
      To be called only after fitting and when a minimizer supporting MINOS is used
      otherwise an error (false) is returned.
      The result will be appended in the fit result class
      Optionally a vector of parameter indeces can be passed for selecting
      the parameters to analyse using FitConfig::SetMinosErrors
    */
   bool CalculateMinosErrors();

   /**
      access to the fit configuration (const method)
   */
   const FitConfig & Config() const { return fConfig; }

   /**
      access to the configuration (non const method)
   */
   FitConfig & Config() { return fConfig; }

   /**
      query if fit is binned. In cse of false teh fit can be unbinned
      or is not defined (like in case of fitting through a ::FitFCN)
    */
   bool IsBinFit() const { return fBinFit; }

   /**
      return pointer to last used minimizer
      (is NULL in case fit is not yet done)
      This pointer is guranteed to be valid as far as the fitter class is valid and a new fit is not redone.
      To be used only after fitting.
      The pointer should not be stored and will be invalided after performing a new fitting.
      In this case a new instance of ROOT::Math::Minimizer will be re-created and can be
      obtained calling again GetMinimizer()
    */
   ROOT::Math::Minimizer * GetMinimizer() const { return fMinimizer.get(); }

   /**
      return pointer to last used objective function
      (is NULL in case fit is not yet done)
      This pointer will be valid as far as the fitter class
      has not been deleted. To be used after the fitting.
      The pointer should not be stored and will be invalided after performing a new fitting.
      In this case a new instance of the function pointer will be re-created and can be
      obtained calling again GetFCN()
    */
   ROOT::Math::IMultiGenFunction * GetFCN() const { return fObjFunction.get(); }


   /**
      apply correction in the error matrix for the weights for likelihood fits
      This method can be called only after a fit. The
      passed function (loglw2) is a log-likelihood function impelemented using the
      sum of weight squared
      When using FitConfig.SetWeightCorrection() this correction is applied
      automatically when doing a likelihood fit (binned or unbinned)
   */
   bool ApplyWeightCorrection(const ROOT::Math::IMultiGenFunction & loglw2, bool minimizeW2L=false);


protected:


   /// least square fit
   bool DoLeastSquareFit(const ROOT::EExecutionPolicy &executionPolicy = ROOT::EExecutionPolicy::kSequential);
   /// binned likelihood fit
   bool DoBinnedLikelihoodFit(bool extended = true, const ROOT::EExecutionPolicy &executionPolicy = ROOT::EExecutionPolicy::kSequential);
   /// un-binned likelihood fit
   bool DoUnbinnedLikelihoodFit( bool extended = false, const ROOT::EExecutionPolicy &executionPolicy = ROOT::EExecutionPolicy::kSequential);
   /// linear least square fit
   bool DoLinearFit();

   // initialize the minimizer
   bool DoInitMinimizer();
   /// do minimization
   bool DoMinimization(const BaseFunc & f, const ROOT::Math::IMultiGenFunction * chifunc = 0);
   // do minimization after having set obj function
   bool DoMinimization(const ROOT::Math::IMultiGenFunction * chifunc = 0);
   // update config after fit
   void DoUpdateFitConfig();
   // update minimizer options for re-fitting
   bool DoUpdateMinimizerOptions(bool canDifferentMinim = true);
   // get function calls from the FCN
   int GetNCallsFromFCN();

   //set data for the fit
   void SetData(const FitData & data) {
      fData = std::shared_ptr<FitData>(const_cast<FitData*>(&data),DummyDeleter<FitData>());
   }
   // set data and function without cloning them
   template <class T>
   void SetFunctionAndData(const IModelFunctionTempl<T> & func, const FitData & data) {
      SetData(data);
      fFunc = std::shared_ptr<IModelFunctionTempl<T>>(const_cast<IModelFunctionTempl<T>*>(&func),DummyDeleter<IModelFunctionTempl<T>>());
   }

   //set data for the fit using a shared ptr
   template <class Data>
   void SetData(const std::shared_ptr<Data> & data) {
      fData = std::static_pointer_cast<Data>(data);
   }

   /// look at the user provided FCN and get data and model function is
   /// they derive from ROOT::Fit FCN classes
   void ExamineFCN();


   /// internal functions to get data set and model function from FCN
   /// useful for fits done with customized FCN classes
   template <class ObjFuncType>
   bool GetDataFromFCN();


private:

   bool fUseGradient;       // flag to indicate if using gradient or not

   bool fBinFit;            // flag to indicate if fit is binned
                            // in case of false the fit is unbinned or undefined)
                            // flag it is used to compute chi2 for binned likelihood fit

   int fFitType;   // type of fit   (0 undefined, 1 least square, 2 likelihood)

   int fDataSize;  // size of data sets (need for Fumili or LM fitters)

   FitConfig fConfig;       // fitter configuration (options and parameter settings)

   std::shared_ptr<IModelFunction_v> fFunc_v;  //! copy of the fitted  function containing on output the fit result

   std::shared_ptr<IModelFunction> fFunc;  //! copy of the fitted  function containing on output the fit result

   std::shared_ptr<ROOT::Fit::FitResult>  fResult;  //! pointer to the object containing the result of the fit

   std::shared_ptr<ROOT::Math::Minimizer>  fMinimizer;  //! pointer to used minimizer

   std::shared_ptr<ROOT::Fit::FitData>  fData;  //! pointer to the fit data (binned or unbinned data)

   std::shared_ptr<ROOT::Math::IMultiGenFunction>  fObjFunction;  //! pointer to used objective function

};


// internal functions to get data set and model function from FCN
// useful for fits done with customized FCN classes
template <class ObjFuncType>
bool Fitter::GetDataFromFCN()  {
   ObjFuncType * objfunc = dynamic_cast<ObjFuncType*>(fObjFunction.get() );
   if (objfunc) {
      fFunc = objfunc->ModelFunctionPtr();
      fData = objfunc->DataPtr();
      return true;
   }
   else {
      return false;
   }
}

#ifdef R__HAS_VECCORE
template <class NotCompileIfScalarBackend>
void Fitter::SetFunction(const IModelFunction_v &func, bool useGradient)
{
   fUseGradient = useGradient;
   if (fUseGradient) {
      const IGradModelFunction_v *gradFunc = dynamic_cast<const IGradModelFunction_v *>(&func);
      if (gradFunc) {
         SetFunction(*gradFunc, true);
         return;
      } else {
         MATH_WARN_MSG("Fitter::SetFunction",
                       "Requested function does not provide gradient - use it as non-gradient function ");
      }
   }

   //  set the fit model function (clone the given one and keep a copy )
   //  std::cout << "set a non-grad function" << std::endl;
   fUseGradient = false;
   fFunc_v = std::shared_ptr<IModelFunction_v>(dynamic_cast<IModelFunction_v *>(func.Clone()));
   assert(fFunc_v);

   // creates the parameter  settings
   fConfig.CreateParamsSettings(*fFunc_v);
   fFunc.reset();
}

template <class NotCompileIfScalarBackend>
void Fitter::SetFunction(const IGradModelFunction_v &func, bool useGradient)
{
   fUseGradient = useGradient;

   //  set the fit model function (clone the given one and keep a copy )
   fFunc_v = std::shared_ptr<IModelFunction_v>(dynamic_cast<IGradModelFunction_v *>(func.Clone()));
   assert(fFunc_v);

   // creates the parameter  settings
   fConfig.CreateParamsSettings(*fFunc_v);
   fFunc.reset();
}
#endif

   } // end namespace Fit

} // end namespace ROOT

// implementation of inline methods


#ifndef __CINT__

#include "Math/WrappedFunction.h"

template<class Function>
bool ROOT::Fit::Fitter::FitFCN(unsigned int npar, Function & f, const double * par, unsigned int datasize,bool chi2fit) {
   ROOT::Math::WrappedMultiFunction<Function &> wf(f,npar);
   return FitFCN(wf,par,datasize,chi2fit);
}
template<class Function>
bool ROOT::Fit::Fitter::SetFCN(unsigned int npar, Function & f, const double * par, unsigned int datasize,bool chi2fit) {
   ROOT::Math::WrappedMultiFunction<Function &> wf(f,npar);
   return SetFCN(wf,par,datasize,chi2fit);
}




#endif  // endif __CINT__

#endif /* ROOT_Fit_Fitter */
